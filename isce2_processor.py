"""
isce2_processor.py  v2.0
========================
Processes Sentinel-1 SLC pairs using ISCE2 topsApp.py + snaphu to generate
interferograms ready for MintPy.

Output per pair — everything flat in the pair folder:
    interferograms/
      20250106_20250118/
        20250106_20250118_fine.unw          ← unwrapped phase (float32, 2-band BIL)
        20250106_20250118_fine.unw.xml      ← ISCE2 XML metadata
        20250106_20250118_fine.unw.rsc      ← ROI_PAC metadata (DATE12 etc.)
        20250106_20250118_fine.unw.conncomp ← connected components (uint8)
        20250106_20250118_fine.cor          ← coherence (float32)
        20250106_20250118_fine.cor.xml
        20250106_20250118_fine.cor.rsc
        hgt.rdr + .xml                      ← DEM heights
        los.rdr + .xml                      ← LOS geometry (2-band)
        lat.rdr + .xml                      ← latitude lookup
        lon.rdr + .xml                      ← longitude lookup
        incidenceAngle.rdr                  ← extracted from los.rdr band1
        azimuthAngle.rdr                    ← extracted from los.rdr band2
        topsApp.xml / topsApp.log / snaphu.log
      mintpy_config.cfg                     ← auto-generated MintPy template

Known issues / fixes applied:
    - WSL2 Bus error: use sequential processing (num_workers=1)
    - polarization must be 'vv' lowercase (Linux case-sensitive)
    - snaphu_mcf bug in ISCE2 2.6.3: call snaphu binary directly
    - snaphu -g flag for conncomp (not --conncomp)
    - runTopo.py boxes ndim bug: patch separately (see SETUP section)
    - lat/lon/hgt not merged by topsApp --end=filter: merge manually
    - MintPy needs DATE12 in .rsc files (not set by prep_isce.py for topsApp)
    - MintPy mixes ascending/descending (different widths): run separately

SETUP (one-time):
    conda activate isce2
    sudo apt install snaphu
    python3 -c "
    f='/home/pc/miniconda3/envs/isce2/lib/python3.9/site-packages/isce/components/isceobj/TopsProc/runTopo.py'
    c=open(f).read()
    old='    boxes = np.array(boxes)\\n    bbox'
    new='    boxes = np.array(boxes)\\n    if boxes.ndim == 1:\\n        boxes = boxes.reshape(1, -1)\\n    bbox'
    open(f,'w').write(c.replace(old,new)); print('Patched')
    "
"""

import glob
import json
import logging
import math
import os
import re
import shutil
import subprocess
import time

import pandas as pd
from tqdm import tqdm

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────
# DEFAULTS
# ─────────────────────────────────────────────

PAIRS_CSV  = os.path.join("outputs", "sbas_pairs.csv")
SLC_DIR    = os.path.join("outputs", "slc_products")
OUTPUT_DIR = os.path.join("outputs", "interferograms")
DEM_PATH   = ""
SWATHS     = [2]

# Sentinel-1 C-band constants
S1_WAVELENGTH         = 0.05546576  # m
S1_RANGE_PIXEL_SIZE   = 2.329562    # m  (IW SLC)
S1_AZIMUTH_PIXEL_SIZE = 14.1        # m  (7 looks azimuth)
S1_EARTH_RADIUS       = 6371000.0   # m
S1_HEIGHT             = 693000.0    # m  (orbit altitude)
S1_STARTING_RANGE     = 845960.0    # m
S1_RLOOKS             = 19
S1_ALOOKS             = 7


# ─────────────────────────────────────────────
# 1. SAFE SYMLINKS
# ─────────────────────────────────────────────

def ensure_safe_symlinks(slc_dir):
    """Create .SAFE symlinks for downloaded SLC scenes (ISCE2 requires .SAFE suffix).
    
    On NTFS mounts (WSL2 → Windows), symlinks may fail with EPERM.
    In that case we skip silently — find_slc_path() handles both
    bare directories and .SAFE directories.
    """
    for sub in ["ascending", "descending"]:
        subdir = os.path.join(slc_dir, sub)
        if not os.path.isdir(subdir):
            continue
        for name in os.listdir(subdir):
            if name.startswith("S1") and not name.endswith(".SAFE"):
                src  = os.path.join(subdir, name)
                dest = src + ".SAFE"
                if os.path.isdir(src) and not os.path.exists(dest):
                    try:
                        os.symlink(src, dest)
                        logger.info("Created symlink: %s.SAFE", name)
                    except OSError as e:
                        # NTFS doesn't support symlinks from WSL2 without
                        # Windows Developer Mode enabled — skip gracefully
                        logger.debug("Symlink skipped (NTFS): %s — %s", name, e)


# ─────────────────────────────────────────────
# 2. FIND SLC ON DISK
# ─────────────────────────────────────────────

def find_slc_path(scene_name, slc_dir):
    """Find SLC directory. Returns path WITHOUT .SAFE suffix.
    
    Handles three cases:
    1. Directory exists as-is (bare name)
    2. Directory exists with .SAFE suffix already
    3. Symlink exists with .SAFE suffix
    """
    clean = scene_name.replace(".SAFE", "")
    for sub in ["ascending", "descending", "unknown"]:
        # Case 1: bare directory
        candidate = os.path.join(slc_dir, sub, clean)
        if os.path.isdir(candidate):
            return candidate
        # Case 2: already has .SAFE suffix
        if os.path.isdir(candidate + ".SAFE"):
            return candidate   # return without .SAFE — build_topsapp_xml adds it
    return None


# ─────────────────────────────────────────────
# 3. DEM
# ─────────────────────────────────────────────

def get_dem(dem_path, aoi, output_dir):
    """Return DEM path. Uses provided file or downloads via ISCE2 dem.py."""
    if dem_path and os.path.isfile(dem_path):
        logger.info("Using provided DEM: %s", dem_path)
        return dem_path

    dem_output = os.path.join(output_dir, "dem", "dem.dem")
    os.makedirs(os.path.dirname(dem_output), exist_ok=True)

    if os.path.isfile(dem_output):
        logger.info("DEM already exists: %s", dem_output)
        return dem_output

    logger.info("Downloading DEM via ISCE2 dem.py ...")
    cmd = [
        "dem.py", "-a", "stitch",
        "-b",
        str(int(math.floor(aoi["min_lat"]))),
        str(int(math.ceil(aoi["max_lat"]))),
        str(int(math.floor(aoi["min_lon"]))),
        str(int(math.ceil(aoi["max_lon"]))),
        "-r", "-s", "1", "-c", "-o", dem_output,
    ]
    subprocess.run(cmd, check=True)
    logger.info("DEM downloaded: %s", dem_output)
    return dem_output


# ─────────────────────────────────────────────
# 4. BUILD topsApp.xml
# ─────────────────────────────────────────────

def build_topsapp_xml(ref_path, sec_path, dem_path, work_dir, swaths):
    """
    Write topsApp.xml for one pair.

    Key lessons:
        - reference/secondary MUST be inside topsinsar component
        - polarization must be 'vv' lowercase
        - all paths must be absolute
        - do unwrap = False (we run snaphu directly)
        - swaths=[2] reduces memory vs [1,2,3]
    """
    # DEM must be referenced by filename only (no path) because ISCE2 looks
    # for the .vrt sidecar in the pair work_dir. We copy the DEM files there
    # in process_single_pair() before calling this function.
    dem_basename = os.path.basename(dem_path)

    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<topsApp>
  <component name="topsinsar">
    <property name="sensor name">SENTINEL1</property>
    <property name="demFilename">{dem_basename}</property>
    <property name="do unwrap">False</property>
    <property name="do dense offsets">False</property>
    <property name="do ESD">False</property>
    <property name="swaths">{swaths}</property>

    <component name="reference">
      <property name="safe">{os.path.abspath(ref_path)}.SAFE</property>
      <property name="output directory">reference</property>
    </component>

    <component name="secondary">
      <property name="safe">{os.path.abspath(sec_path)}.SAFE</property>
      <property name="output directory">secondary</property>
    </component>
  </component>
</topsApp>
"""
    path = os.path.join(work_dir, "topsApp.xml")
    with open(path, "w") as f:
        f.write(xml)
    return os.path.abspath(path)


# ─────────────────────────────────────────────
# 5. MERGE lat/lon/hgt FROM geom_reference
# ─────────────────────────────────────────────

def merge_lat_lon_hgt(work_dir):
    """
    Merge per-burst lat/lon/hgt files from geom_reference/IW{n}/ into merged/.

    ISCE2 topo step writes per-burst files to geom_reference/ but only merges
    them in the geocode step. Since we stop at --end=filter, we merge manually.

    Critical dtype rules (must match ISCE2 native output):
        lat.rdr → DOUBLE (float64) in geom_reference, but MintPy reads as FLOAT
                  → convert to float32 for MintPy compatibility
        lon.rdr → same as lat
        hgt.rdr → FLOAT (float32) — must match radar geometry size (~10MB not 50MB)

    Strategy:
        1. gdal_merge.py with explicit -ot Float32 — correct dtype + all bursts
        2. numpy manual merge — fallback
        3. Copy first burst — last resort approximation
    """
    merged_dir = os.path.join(work_dir, "merged")
    geom_dir   = os.path.join(work_dir, "geom_reference")

    if not os.path.isdir(geom_dir):
        logger.warning("geom_reference/ not found — lat/lon/hgt will be missing")
        return

    # All three must be float32 for MintPy compatibility
    for geo_type in ["lat", "lon", "hgt"]:
        out_file = os.path.join(merged_dir, f"{geo_type}.rdr")
        if os.path.isfile(out_file):
            continue

        burst_files = sorted(glob.glob(
            os.path.join(geom_dir, "IW*", f"{geo_type}_*.rdr")
        ))
        if not burst_files:
            logger.warning("No %s burst files in geom_reference/", geo_type)
            continue

        # Try gdal_merge.py with explicit float32 output
        try:
            r = subprocess.run(
                ["gdal_merge.py", "-ot", "Float32", "-o", out_file] + burst_files,
                capture_output=True, timeout=120
            )
            if r.returncode == 0 and os.path.isfile(out_file):
                size_mb = os.path.getsize(out_file) / 1024 / 1024
                logger.info("Merged %s.rdr (%d bursts, %.1f MB)", geo_type, len(burst_files), size_mb)
                continue
            else:
                logger.debug("gdal_merge.py stderr: %s", r.stderr.decode()[:200])
        except Exception as e:
            logger.debug("gdal_merge.py failed for %s: %s", geo_type, e)

        # Fallback: numpy merge — read each burst XML for dims, concatenate rows
        try:
            import numpy as np
            arrays = []
            for bf in burst_files:
                bxml = bf + ".xml"
                bw, bl = _read_image_dims(bxml)
                if bw and bl:
                    arr = np.fromfile(bf, dtype=np.float64).reshape(bl, bw)
                    arrays.append(arr.astype(np.float32))
            if arrays:
                merged = np.vstack(arrays)
                merged.tofile(out_file)
                size_mb = os.path.getsize(out_file) / 1024 / 1024
                logger.info("Merged %s.rdr via numpy (%d bursts, %.1f MB)", geo_type, len(arrays), size_mb)
                continue
        except Exception as e:
            logger.debug("numpy merge failed for %s: %s", geo_type, e)

        # Last resort: copy first burst, convert to float32
        try:
            import numpy as np
            bxml = burst_files[0] + ".xml"
            bw, bl = _read_image_dims(bxml)
            if bw and bl:
                arr = np.fromfile(burst_files[0], dtype=np.float64).reshape(bl, bw)
                arr.astype(np.float32).tofile(out_file)
                logger.warning("%s.rdr — first burst only (approximate, float32)", geo_type)
            else:
                shutil.copy2(burst_files[0], out_file)
                logger.warning("%s.rdr — raw copy of first burst", geo_type)
        except Exception:
            shutil.copy2(burst_files[0], out_file)
            logger.warning("%s.rdr — raw copy of first burst (fallback)", geo_type)


# ─────────────────────────────────────────────
# 6. EXTRACT incidenceAngle / azimuthAngle
# ─────────────────────────────────────────────

def extract_geometry_angles(merged_dir):
    """
    Split los.rdr (2-band: incidence + heading) into separate single-band files.
    MintPy prefers separate files for incidenceAngle and azimuthAngle.
    """
    import numpy as np

    los_file = os.path.join(merged_dir, "los.rdr")
    if not os.path.isfile(los_file):
        logger.warning("los.rdr not found — skipping angle extraction")
        return

    try:
        # Use direct XML parsing — avoids ISCE2 2.6.3 "filename attribute" exception
        los_width, length = _read_image_dims(los_file + ".xml")
        if los_width is None or length is None:
            # Fallback: los.rdr is FLOAT 2-band (8 bytes/pixel)
            los_width, length = _read_image_dims_from_file(los_file, bands=2, dtype_bytes=4)
        if los_width is None or length is None:
            logger.warning("Cannot determine los.rdr dimensions — skipping angle extraction")
            return

        los_data  = np.fromfile(los_file, dtype=np.float32)
        los_data  = los_data.reshape(length, los_width, 2)
        inc_angle = los_data[:, :, 0]
        az_angle  = los_data[:, :, 1]

        inc_angle.tofile(os.path.join(merged_dir, "incidenceAngle.rdr"))
        az_angle.tofile(os.path.join(merged_dir, "azimuthAngle.rdr"))

        logger.info(
            "Angles extracted — incidence mean=%.1f°  azimuth mean=%.1f°",
            inc_angle[inc_angle > 0].mean(),
            az_angle[az_angle != 0].mean()
        )
    except Exception as e:
        logger.warning("Angle extraction failed: %s", e)


# ─────────────────────────────────────────────
# 7. XML HELPERS
# ─────────────────────────────────────────────

def _write_imageFile_xml(file_path, width, length, dtype, bands, scheme="BIL"):
    """
    Write an ISCE2 imageFile-format XML sidecar.

    MintPy's read_isce_xml() requires root.tag.startswith('image') and
    properties as direct children of root — NOT inside <component>.

    dtype values: FLOAT, DOUBLE, CFLOAT, INT, BYTE
    """
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<imageFile>
    <property name="width"><value>{width}</value></property>
    <property name="length"><value>{length}</value></property>
    <property name="number_bands"><value>{bands}</value></property>
    <property name="data_type"><value>{dtype}</value></property>
    <property name="scheme"><value>{scheme}</value></property>
    <property name="byte_order"><value>l</value></property>
    <property name="filename"><value>{os.path.abspath(file_path)}</value></property>
</imageFile>
"""
    with open(file_path + ".xml", "w") as f:
        f.write(xml)


def _read_image_dims(xml_path):
    """
    Read (width, length) from an ISCE2 XML sidecar without using the ISCE2
    Image API.

    The ISCE2 Image.load() method logs a non-fatal error:
        "The attribute corresponding to the key 'filename' is not present ..."
    when the XML contains a <property name="filename"> tag that no longer
    maps to an attribute in the Image class (ISCE2 2.6.3 regression).
    Although the error is non-fatal for topsApp itself, it is raised as an
    exception inside our try/except block, causing process_single_pair() to
    return early before snaphu runs.

    This function parses the XML directly using stdlib xml.etree and is
    therefore immune to ISCE2 version quirks.

    Returns (width, length) as ints, or (None, None) if the file is absent
    or unparseable.
    """
    import xml.etree.ElementTree as ET

    if not os.path.isfile(xml_path):
        return None, None

    try:
        root  = ET.parse(xml_path).getroot()
        props = {}
        # Handle both root-level properties and component-wrapped properties
        for prop in root.iter("property"):
            name = prop.get("name", "")
            val  = prop.findtext("value")
            if val is not None:
                props[name] = val.strip()

        width  = int(props["width"])   if "width"  in props else None
        length = int(props["length"])  if "length" in props else None
        return width, length
    except Exception as e:
        logger.debug("_read_image_dims(%s): %s", xml_path, e)
        return None, None


def _read_image_dims_from_file(data_path, bands=2, dtype_bytes=4):
    """
    Infer (width, length) from raw file size when the XML is missing or broken.

    Assumes BIL interleave: total_bytes = width * length * bands * dtype_bytes
    Tries known Sentinel-1 IW widths for this AOI.

    Returns (width, length) or (None, None).
    """
    if not os.path.isfile(data_path):
        return None, None

    size = os.path.getsize(data_path)
    # Common lengths for Sentinel-1 IW SLC bursts processed with 7 az-looks
    for length in [1941, 1960, 2000, 1800, 2100]:
        for width in [1359, 1372, 1380, 1400, 1450]:
            if size == width * length * bands * dtype_bytes:
                logger.debug(
                    "_read_image_dims_from_file: matched %s w=%d l=%d",
                    os.path.basename(data_path), width, length
                )
                return width, length
    return None, None


def _fix_xml_filename(xml_path, new_data_path):
    """Update the filename property inside an existing XML sidecar."""
    if not os.path.isfile(xml_path):
        return
    try:
        with open(xml_path) as f:
            content = f.read()
        content = re.sub(
            r'<property name="filename">.*?</property>',
            f'<property name="filename">{os.path.abspath(new_data_path)}</property>',
            content,
            flags=re.DOTALL
        )
        with open(xml_path, "w") as f:
            f.write(content)
    except Exception as e:
        logger.debug("Could not fix XML filename in %s: %s", xml_path, e)


# ─────────────────────────────────────────────
# 8. CREATE snaphu OUTPUT XMLs
# ─────────────────────────────────────────────

def create_snaphu_xmls(merged_dir, width, length):
    """
    Create XML sidecars for snaphu outputs (snaphu doesn't generate them).

    filt_topophase.unw      → FLOAT, 2 bands (amplitude + phase), BIL
    filt_topophase.unw.conncomp → BYTE, 1 band, BIL
    """
    # Unwrapped phase
    unw_file = os.path.join(merged_dir, "filt_topophase.unw")
    if os.path.isfile(unw_file) and not os.path.isfile(unw_file + ".xml"):
        # Write directly — UnwImage.dump() triggers the same "filename attribute"
        # ISCE2 2.6.3 exception as Image.load(), so we bypass the API entirely.
        _write_imageFile_xml(unw_file, width, length, "FLOAT", 2)
        logger.info("Created: filt_topophase.unw.xml")

    # Connected components
    cc_file = os.path.join(merged_dir, "filt_topophase.unw.conncomp")
    if os.path.isfile(cc_file) and not os.path.isfile(cc_file + ".xml"):
        _write_imageFile_xml(cc_file, width, length, "BYTE", 1)
        logger.info("Created: filt_topophase.unw.conncomp.xml")


# ─────────────────────────────────────────────
# 9. WRITE .rsc FILES
# ─────────────────────────────────────────────

def write_rsc_files(pair_dir, pair_id, width, length, orbit_direction="DESCENDING"):
    """
    Write ROI_PAC-style .rsc metadata files for all interferogram products.

    MintPy's load_data reads DATE12 from these files when processing
    ISCE2 topsApp output (prep_isce.py cannot handle topsApp.xml as metadata).

    DATE12 format: YYYYMMDD_YYYYMMDD  (full 8-digit dates with underscore)
    """
    ref8   = pair_id[:8]
    sec8   = pair_id[9:]
    date12 = f"{ref8}_{sec8}"

    rsc = (
        f"WIDTH                    {width}\n"
        f"FILE_LENGTH              {length}\n"
        f"XMIN                     0\n"
        f"XMAX                     {width - 1}\n"
        f"YMIN                     0\n"
        f"YMAX                     {length - 1}\n"
        f"DATE12                   {date12}\n"
        f"DATE                     {ref8}\n"
        f"DATE2                    {sec8}\n"
        f"WAVELENGTH               {S1_WAVELENGTH}\n"
        f"PLATFORM                 SENTINEL1\n"
        f"SENSOR                   SENTINEL1\n"
        f"ORBIT_DIRECTION          {orbit_direction}\n"
        f"RANGE_PIXEL_SIZE         {S1_RANGE_PIXEL_SIZE}\n"
        f"AZIMUTH_PIXEL_SIZE       {S1_AZIMUTH_PIXEL_SIZE}\n"
        f"EARTH_RADIUS             {S1_EARTH_RADIUS}\n"
        f"HEIGHT                   {S1_HEIGHT}\n"
        f"STARTING_RANGE           {S1_STARTING_RANGE}\n"
        f"RLOOKS                   {S1_RLOOKS}\n"
        f"ALOOKS                   {S1_ALOOKS}\n"
    )

    for fname in [
        f"{pair_id}_fine.unw",
        f"{pair_id}_fine.cor",
        f"{pair_id}_fine.unw.conncomp",
    ]:
        rsc_path = os.path.join(pair_dir, fname + ".rsc")
        with open(rsc_path, "w") as f:
            f.write(rsc)

    logger.info("Written .rsc files for %s  DATE12=%s", pair_id, date12)


# ─────────────────────────────────────────────
# 10. FINALIZE PAIR FOLDER
# ─────────────────────────────────────────────

def finalize_pair(pair_id, work_dir, orbit_direction="DESCENDING"):
    """
    Move all outputs from merged/ to pair root, rename to MintPy names,
    write XML/RSC sidecars, remove intermediate directories.

    Final pair folder contents:
        {pair_id}_fine.unw + .xml + .rsc
        {pair_id}_fine.unw.conncomp + .xml
        {pair_id}_fine.cor + .xml + .rsc
        hgt.rdr + .xml
        los.rdr + .xml
        lat.rdr + .xml
        lon.rdr + .xml
        incidenceAngle.rdr
        azimuthAngle.rdr
        topsApp.xml / topsApp.log / snaphu.log
    """
    merged  = os.path.join(work_dir, "merged")
    unw_xml = os.path.join(merged, "filt_topophase.unw.xml")

    # Read dimensions from unw XML (before moving) — use direct XML parsing to
    # avoid the ISCE2 2.6.3 "filename attribute not present" exception.
    width, length = _read_image_dims(unw_xml)
    if width is None or length is None:
        unw_flat = os.path.join(merged, "filt_topophase.flat")
        width, length = _read_image_dims_from_file(unw_flat, bands=1, dtype_bytes=8)

    # -- Move + rename interferogram files --
    renames = {
        "filt_topophase.unw":          f"{pair_id}_fine.unw",
        "filt_topophase.unw.xml":      f"{pair_id}_fine.unw.xml",
        "filt_topophase.unw.conncomp": f"{pair_id}_fine.unw.conncomp",
        "filt_topophase.unw.conncomp.xml": f"{pair_id}_fine.unw.conncomp.xml",
        "topophase.cor":               f"{pair_id}_fine.cor",
        "topophase.cor.xml":           f"{pair_id}_fine.cor.xml",
        "filt_topophase.flat":         f"{pair_id}_fine.flat",
        "filt_topophase.flat.xml":     f"{pair_id}_fine.flat.xml",
    }
    for src_name, dst_name in renames.items():
        src = os.path.join(merged, src_name)
        dst = os.path.join(work_dir, dst_name)
        if os.path.isfile(src) and not os.path.exists(dst):
            shutil.move(src, dst)

    # Fix filename references inside moved XMLs
    for fname in [f"{pair_id}_fine.unw", f"{pair_id}_fine.cor",
                  f"{pair_id}_fine.unw.conncomp"]:
        _fix_xml_filename(
            os.path.join(work_dir, fname + ".xml"),
            os.path.join(work_dir, fname)
        )

    # -- Move geometry files --
    geom_files = [
        "hgt.rdr", "los.rdr", "lat.rdr", "lon.rdr",
        "incidenceAngle.rdr", "azimuthAngle.rdr",
        "hgt.rdr.xml", "los.rdr.xml", "lat.rdr.xml", "lon.rdr.xml",
        "hgt.rdr.vrt", "los.rdr.vrt",
    ]
    for fname in geom_files:
        src = os.path.join(merged, fname)
        dst = os.path.join(work_dir, fname)
        if os.path.isfile(src) and not os.path.exists(dst):
            shutil.move(src, dst)

    # Fix geometry XML filename references
    for fname in ["hgt.rdr", "los.rdr", "lat.rdr", "lon.rdr"]:
        _fix_xml_filename(
            os.path.join(work_dir, fname + ".xml"),
            os.path.join(work_dir, fname)
        )

    # -- Write clean geometry XMLs (always rewrite, never just patch) --
    # Native topsApp XMLs contain doc-only <property> elements without <value>
    # children that crash MintPy's read_isce_xml(). Always emit our minimal
    # 7-property imageFile format so MintPy never sees the raw topsApp XML.
    if width and length:
        geom_specs = {
            "hgt.rdr":            ("FLOAT", 1),
            "los.rdr":            ("FLOAT", 2),
            "lat.rdr":            ("FLOAT", 1),   # float32 after merge_lat_lon_hgt
            "lon.rdr":            ("FLOAT", 1),   # float32 after merge_lat_lon_hgt
            "incidenceAngle.rdr": ("FLOAT", 1),
            "azimuthAngle.rdr":   ("FLOAT", 1),
        }
        for fname, (dtype, bands) in geom_specs.items():
            fpath = os.path.join(work_dir, fname)
            if os.path.isfile(fpath):
                _write_imageFile_xml(fpath, width, length, dtype, bands)
                logger.debug("Wrote clean XML: %s.xml", fname)

    # -- Write .rsc files with DATE12 (required by MintPy) --
    if width and length:
        write_rsc_files(work_dir, pair_id, width, length, orbit_direction)

    # -- Remove intermediate directories --
    for subdir in [
        "merged", "reference", "secondary",
        "fine_coreg", "fine_offsets", "geom_reference",
        "PICKLE", "coarse_coreg", "coarse_offsets",
        "fine_interferogram", "overlap",
    ]:
        path = os.path.join(work_dir, subdir)
        if os.path.isdir(path):
            shutil.rmtree(path)

    # -- Verify required files --
    required = {
        f"{pair_id}_fine.unw":              "unwrapped phase",
        f"{pair_id}_fine.unw.xml":          "unw XML",
        f"{pair_id}_fine.unw.rsc":          "unw RSC (DATE12)",
        f"{pair_id}_fine.unw.conncomp":     "connected components",
        f"{pair_id}_fine.cor":              "coherence",
        "hgt.rdr":                          "DEM heights",
        "los.rdr":                          "LOS geometry",
        "lat.rdr":                          "latitude lookup",
        "lon.rdr":                          "longitude lookup",
    }
    missing = [f"{f} ({d})" for f, d in required.items()
               if not os.path.isfile(os.path.join(work_dir, f))]

    if missing:
        logger.warning("Missing in %s: %s", pair_id, missing)
    else:
        logger.info("Pair %s complete ✅", pair_id)


# ─────────────────────────────────────────────
# 11. PROCESS ONE PAIR
# ─────────────────────────────────────────────

def process_single_pair(pair, slc_dir, output_dir, dem_path, swaths):
    """
    Full processing pipeline for one interferogram pair:
        topsApp.py --end=filter → snaphu → finalize

    Returns a result dict with status: success / already_exists / failed
    """
    ref_name = pair["reference_scene"]
    sec_name = pair["secondary_scene"]
    ref_date = pair["reference_date"].replace("-", "")
    sec_date = pair["secondary_date"].replace("-", "")
    pair_id  = f"{ref_date}_{sec_date}"
    orbit    = str(pair.get("path", "")).strip()

    result = {
        "pair_id": pair_id, "status": "failed",
        "output_dir": "", "reason": "",
    }

    work_dir  = os.path.join(output_dir, pair_id)
    done_file = os.path.join(work_dir, f"{pair_id}_fine.unw")

    if os.path.isfile(done_file):
        result["status"]     = "already_exists"
        result["output_dir"] = work_dir
        logger.info("Already exists — skipping: %s", pair_id)
        return result

    os.makedirs(work_dir, exist_ok=True)

    # Find SLC files
    ref_path = find_slc_path(ref_name, slc_dir)
    sec_path = find_slc_path(sec_name, slc_dir)
    if not ref_path:
        result["reason"] = f"Reference SLC not found: {ref_name}"
        return result
    if not sec_path:
        result["reason"] = f"Secondary SLC not found: {sec_name}"
        return result

    # Determine orbit direction from path number or CSV
    flight_dir = str(pair.get("flight_dir", "")).upper()
    if not flight_dir:
        # Infer from path: Sentinel-1 descending paths typically < 100
        try:
            path_num = int(orbit)
            flight_dir = "ASCENDING" if path_num > 100 else "DESCENDING"
        except (ValueError, TypeError):
            flight_dir = "DESCENDING"

    # Copy DEM files to pair work_dir — ISCE2 looks for the .vrt sidecar
    # in the current working directory (cwd=work_dir), not by absolute path.
    # Use copyfile (not copy2) to avoid NTFS permission metadata errors on WSL2.
    import shutil as _shutil
    dem_basename = os.path.basename(dem_path)
    for ext in ["", ".xml", ".vrt", ".hdr"]:
        src = dem_path + ext if ext else dem_path
        dst = os.path.join(work_dir, dem_basename + (ext if ext else ""))
        if os.path.isfile(src) and not os.path.isfile(dst):
            try:
                _shutil.copyfile(src, dst)
            except Exception as cp_e:
                logger.warning("DEM copy %s: %s", ext or ".hgt", cp_e)

    # Generate ISCE2-native DEM XML in the pair work_dir
    dem_local = os.path.join(work_dir, dem_basename)
    dem_xml_dst = dem_local + ".xml"
    try:
        import isce
        from isceobj.Image import createDemImage
        img = createDemImage()
        img.initImage(dem_local, "read", 3601)
        img.setLength(7201)
        img.setWidth(3601)
        img.firstLatitude   = 1.00013888888889
        img.firstLongitude  = -92.0001388888889
        img.deltaLatitude   = -0.000277777777777778
        img.deltaLongitude  =  0.000277777777777778
        img.renderHdr()
        img.dump(dem_xml_dst)
        logger.debug("DEM XML written: %s", dem_xml_dst)
    except Exception as dem_e:
        logger.warning("DEM XML via ISCE2 API failed (%s) — using existing", dem_e)

    # Build XML
    xml_path = build_topsapp_xml(ref_path, sec_path, dem_path, work_dir, swaths)

    import sys as _sys
    PYTHON   = _sys.executable
    TOPSAPP  = os.path.join(os.path.dirname(os.path.dirname(PYTHON)),
                            "lib", "python3.9", "site-packages", "isce",
                            "applications", "topsApp.py")
    if not os.path.isfile(TOPSAPP):
        # fallback: search on PATH
        import shutil as _sh
        _found = _sh.which("topsApp.py")
        TOPSAPP = _found if _found else TOPSAPP

    log_file = os.path.join(work_dir, "topsApp.log")
    logger.info("Running topsApp.py for %s ...", pair_id)
    try:
        subprocess.run(
            [PYTHON, TOPSAPP, os.path.basename(xml_path), "--end=filter"],
            cwd=work_dir,
            stdout=open(log_file, "w"),
            stderr=subprocess.STDOUT,
            timeout=10800,
        )
    except subprocess.TimeoutExpired:
        result["reason"] = "topsApp.py timeout (3h)"
        return result
    except FileNotFoundError:
        result["reason"] = "topsApp.py not found — is ISCE2 on PATH?"
        return result

    # Check output file — ignore exit code (ISCE2 prints warnings to stderr
    # that cause non-zero exit even on success, e.g. "filename" attribute warning)
    merged_dir = os.path.join(work_dir, "merged")
    ifg_flat   = os.path.join(merged_dir, "filt_topophase.flat")
    if not os.path.isfile(ifg_flat):
        result["reason"] = "filt_topophase.flat not found — topsApp.py failed"
        return result
    logger.info("topsApp.py completed — filt_topophase.flat found ✅")

    # Read interferogram dimensions — parse XML directly to avoid the ISCE2
    # 2.6.3 "filename attribute not present" exception from Image.load().
    width, length = _read_image_dims(ifg_flat + ".xml")
    if width is None or length is None:
        # Fallback: infer from file size (filt_topophase.flat is CFLOAT = 8 bytes/pixel, 1 band)
        width, length = _read_image_dims_from_file(ifg_flat, bands=1, dtype_bytes=8)
    if width is None or length is None:
        result["reason"] = (
            f"Cannot read interferogram dimensions from "
            f"{ifg_flat}.xml — check topsApp.log"
        )
        return result
    logger.info("Interferogram dimensions: width=%d  length=%d", width, length)

    # Run snaphu
    unw_file = os.path.join(merged_dir, "filt_topophase.unw")
    cor_file = os.path.join(merged_dir, "topophase.cor")
    cc_file  = os.path.join(merged_dir, "filt_topophase.unw.conncomp")
    snap_log = os.path.join(work_dir, "snaphu.log")

    logger.info("Running snaphu for %s  width=%d  length=%d ...", pair_id, width, length)
    try:
        with open(snap_log, "w") as log:
            subprocess.run(
                [
                    "snaphu", ifg_flat, str(width),
                    "-o", unw_file,
                    "-c", cor_file,
                    "-g", cc_file,    # conncomp (-g, not --conncomp)
                    "--mcf", "-v",
                ],
                check=True, cwd=work_dir,
                stdout=log, stderr=subprocess.STDOUT,
                timeout=7200,
            )
    except subprocess.TimeoutExpired:
        result["reason"] = "snaphu timeout (2h)"
        return result
    except subprocess.CalledProcessError as e:
        result["reason"] = f"snaphu exit code {e.returncode}"
        return result
    except FileNotFoundError:
        result["reason"] = "snaphu not found — sudo apt install snaphu"
        return result

    # Post-processing
    create_snaphu_xmls(merged_dir, width, length)
    merge_lat_lon_hgt(work_dir)
    extract_geometry_angles(merged_dir)
    finalize_pair(pair_id, work_dir, orbit_direction=flight_dir)

    result["status"]     = "success"
    result["output_dir"] = work_dir
    return result


# ─────────────────────────────────────────────
# 12. SEQUENTIAL PROCESSING
# ─────────────────────────────────────────────

def run_sequential(pairs, slc_dir, output_dir, dem_path, swaths):
    """
    Process all pairs one at a time.

    WHY SEQUENTIAL: WSL2 crashes (SIGBUS) when multiple subprocesses read
    large SLC files from /mnt/d/ simultaneously. On native Linux, parallel
    processing works fine.
    """
    results = []
    bar = tqdm(total=len(pairs), desc="Interferograms", unit="pair")

    for idx, pair in enumerate(pairs, 1):
        pid = (pair.get("reference_date", "?").replace("-", "") + "_" +
               pair.get("secondary_date", "?").replace("-", ""))
        logger.info("Pair %d/%d — %s", idx, len(pairs), pid)

        result = process_single_pair(pair, slc_dir, output_dir, dem_path, swaths)
        results.append(result)
        bar.update(1)
        bar.set_postfix_str(f"{result['pair_id']}  [{result['status']}]")

        if result["status"] == "failed":
            logger.warning("FAILED %s — %s", result["pair_id"], result["reason"])

    bar.close()
    return results


# ─────────────────────────────────────────────
# 13. SUMMARY + LOG
# ─────────────────────────────────────────────

def save_log(results, output_dir):
    path = os.path.join(output_dir, "isce2_processing_log.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Log saved → %s", path)


def print_summary(results):
    statuses = [r["status"] for r in results]
    print("\n" + "=" * 50)
    print("  ISCE2 PROCESSING SUMMARY")
    print("=" * 50)
    print(f"  Total pairs    : {len(results)}")
    print(f"  Success        : {statuses.count('success')}")
    print(f"  Already exists : {statuses.count('already_exists')}")
    print(f"  Failed         : {statuses.count('failed')}")
    failed = [r for r in results if r["status"] == "failed"]
    if failed:
        print("\n  Failed pairs:")
        for r in failed:
            print(f"    {r['pair_id']} — {r['reason']}")
    print("=" * 50 + "\n")


# ─────────────────────────────────────────────
# 14. MAIN ENTRY POINT
# ─────────────────────────────────────────────

def run(
    pairs_csv  = PAIRS_CSV,
    slc_dir    = SLC_DIR,
    output_dir = OUTPUT_DIR,
    dem_path   = DEM_PATH,
    num_workers = 1,   # ignored — always sequential on WSL2
    batch_size  = 1,   # ignored
    aoi        = None,
    swaths     = None,
):
    """Called by main_v3.py Step 3."""
    if swaths is None:
        swaths = SWATHS

    os.makedirs(output_dir, exist_ok=True)
    ensure_safe_symlinks(slc_dir)

    logger.info("Loading pairs from %s ...", pairs_csv)
    df = pd.read_csv(pairs_csv)
    if df.empty or "reference_scene" not in df.columns:
        raise ValueError("sbas_pairs.csv is empty or invalid.")

    pairs = df.to_dict("records")
    logger.info("Loaded %d pairs.", len(pairs))

    dem = get_dem(dem_path, aoi or {}, output_dir)

    start   = time.time()
    results = run_sequential(pairs, slc_dir, output_dir, dem, swaths)
    elapsed = time.time() - start

    save_log(results, output_dir)
    print_summary(results)
    logger.info("ISCE2 finished in %.1f minutes.", elapsed / 60)

    return results


# ─────────────────────────────────────────────
# 15. CLI
# ─────────────────────────────────────────────

def main():
    import argparse
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    p = argparse.ArgumentParser(
        description="ISCE2 Interferogram Processor v2.0 — Sentinel-1 SBAS",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--pairs-csv",  default=PAIRS_CSV,  help="SBAS pairs CSV")
    p.add_argument("--slc-dir",    default=SLC_DIR,    help="SLC root directory")
    p.add_argument("--output-dir", default=OUTPUT_DIR, help="Output directory (Linux fs on WSL2)")
    p.add_argument("--dem-path",   default=DEM_PATH,   help="DEM .hgt file")
    p.add_argument("--swaths",     default="2",        help="Swaths: 2 or 1,2,3")
    args = p.parse_args()

    run(
        pairs_csv  = args.pairs_csv,
        slc_dir    = args.slc_dir,
        output_dir = args.output_dir,
        dem_path   = args.dem_path,
        swaths     = [int(s.strip()) for s in args.swaths.split(",")],
    )


if __name__ == "__main__":
    main()