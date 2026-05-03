"""
mintpy_processor.py  v1.4
==========================
Runs MintPy SBAS time-series processing on ISCE2 topsApp interferograms.

Key issues solved:
    - prep_isce.py cannot parse topsApp.xml -> .rsc files written manually
    - DATE12 not set automatically -> injected via .rsc files
    - Ascending/descending have different widths -> staged into separate dirs
    - Geometry XML must use imageFile format (not component-wrapped)
    - lat/lon XML must say FLOAT not DOUBLE (data is float32 after merge)
    - conncomp is uint8 (1 byte) not int16
    - Auto reference point lands outside conncomp -> maskFile = maskConnComp.h5
    - unwrapError bridging fragile over vegetated areas -> default is now "no"
    - smallbaselineApp.py named .EXE on Windows -> multi-variant discovery
    - Stale stage dirs / templates / cfg cause wrong paths -> purged on re-run
    - Path truncation (interferograms -> interferog) -> realpath() at entry
    - _detect_width/_detect_length use brute-force range, not hardcoded list

Changelog:
    v1.1 - write_rsc_files writes P_BASELINE_TOP/BOTTOM_HDR
    v1.2 - reference point uses maskConnComp.h5; unwrapError=no when pairs<3
    v1.3 - lat/lon XML dtype DOUBLE->FLOAT; _detect_length() added;
           fix_geometry_xmls runs AFTER stage_track_pairs
    v1.4 - unwrapError.method default "bridging"->"no" (fragile, see note)
           run_smallbaseline discovery rewritten; FileNotFoundError path removed
           _detect_width/_detect_length use range(1200,1600) not fixed list
"""

import glob
import json
import logging
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET

import pandas as pd

logger = logging.getLogger(__name__)

IFG_DIR    = os.path.join("outputs", "interferograms")
PAIRS_CSV  = os.path.join("outputs", "sbas_pairs.csv")
OUTPUT_DIR = os.path.join("outputs", "mintpy")
TRACK      = "all"

MIN_COHERENCE     = 0.4
MIN_TEMP_COH      = 0.7
REF_MIN_COHERENCE = 0.85
TROPO_METHOD      = "no"
GEOCODE_STEP      = "-0.0001, 0.0001"


# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

def _load_bperp_map(pairs_csv):
    try:
        df = pd.read_csv(pairs_csv)
        return {
            f"{str(r['reference_date']).replace('-','')}_{str(r['secondary_date']).replace('-','')}":
            float(r.get("spatial_baseline_m", 0.0))
            for _, r in df.iterrows()
        }
    except Exception:
        logger.warning("Cannot read pairs CSV for bperp — using 0.0")
        return {}


def _detect_width(pair_dir, pair_id):
    unw_file = os.path.join(pair_dir, f"{pair_id}_fine.unw")
    unw_xml  = unw_file + ".xml"
    if os.path.isfile(unw_xml):
        try:
            root = ET.parse(unw_xml).getroot()
            for prop in root.iter("property"):
                if prop.get("name") == "width":
                    val = prop.find("value")
                    if val is not None and val.text:
                        return int(val.text)
        except Exception:
            pass
    if os.path.isfile(unw_file):
        size = os.path.getsize(unw_file)
        for w in range(1200, 1600):
            if size % (w * 2 * 4) == 0:
                length = size // (w * 2 * 4)
                if 1000 < length < 5000:
                    return w
    return None


def _detect_length(pair_dir, pair_id):
    unw_file = os.path.join(pair_dir, f"{pair_id}_fine.unw")
    unw_xml  = unw_file + ".xml"
    if os.path.isfile(unw_xml):
        try:
            root = ET.parse(unw_xml).getroot()
            for prop in root.iter("property"):
                if prop.get("name") == "length":
                    val = prop.find("value")
                    if val is not None and val.text:
                        return int(val.text)
        except Exception:
            pass
    if os.path.isfile(unw_file):
        size = os.path.getsize(unw_file)
        for w in range(1200, 1600):
            if size % (w * 2 * 4) == 0:
                return size // (w * 2 * 4)
    return None


def _write_imageFile_xml(file_path, width, length, dtype, bands, scheme="BIL"):
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


# ─────────────────────────────────────────────
# 1. DETECT TRACKS
# ─────────────────────────────────────────────

def detect_tracks(ifg_dir, pairs_csv):
    tracks = {"ascending": [], "descending": []}
    try:
        df = pd.read_csv(pairs_csv)
    except Exception:
        logger.warning("Cannot read pairs CSV — treating all as ascending")
        pair_dirs = sorted([
            d for d in os.listdir(ifg_dir)
            if os.path.isdir(os.path.join(ifg_dir, d))
            and len(d) == 17 and d[8] == "_"
        ])
        tracks["ascending"] = pair_dirs
        return {k: v for k, v in tracks.items() if v}

    for _, row in df.iterrows():
        ref8     = str(row.get("reference_date", "")).replace("-", "")
        sec8     = str(row.get("secondary_date", "")).replace("-", "")
        pid      = f"{ref8}_{sec8}"
        pair_dir = os.path.join(ifg_dir, pid)
        if not os.path.isdir(pair_dir):
            continue
        flight = str(row.get("flight_dir", "")).lower()
        if not flight:
            try:
                path_num = int(row.get("path", 0))
                flight = "ascending" if path_num >= 100 else "descending"
            except (ValueError, TypeError):
                flight = "ascending"
        if "asc" in flight:
            tracks["ascending"].append(pid)
        else:
            tracks["descending"].append(pid)
    return {k: v for k, v in tracks.items() if v}


# ─────────────────────────────────────────────
# 2. FIX XML FILES
# ─────────────────────────────────────────────

def fix_geometry_xmls(pair_dirs):
    """
    Unconditionally rewrite geometry XMLs to minimal imageFile format.

    WHY UNCONDITIONAL: native topsApp XMLs contain doc-only <property>
    elements (no <value> child). MintPy's read_isce_xml() crashes with
    AttributeError: 'NoneType'.text on those properties.
    Always emitting our clean format eliminates this entirely.

    Called AFTER stage_track_pairs() so realpath() writes through symlinks
    to the actual files — readable via both original and stage paths.
    """
    specs = {
        "hgt.rdr":            ("FLOAT", 1),
        "los.rdr":            ("FLOAT", 2),
        "lat.rdr":            ("FLOAT", 1),
        "lon.rdr":            ("FLOAT", 1),
        "incidenceAngle.rdr": ("FLOAT", 1),
        "azimuthAngle.rdr":   ("FLOAT", 1),
    }
    rewritten = 0
    for pair_dir in pair_dirs:
        pair_id = os.path.basename(pair_dir)
        width   = _detect_width(pair_dir, pair_id)
        length  = _detect_length(pair_dir, pair_id)
        if not width:
            logger.warning("Cannot detect width for %s — using 1372", pair_id)
            width = 1372
        for fname, (dtype, bands) in specs.items():
            fpath = os.path.join(pair_dir, fname)
            if not os.path.isfile(fpath):
                continue
            real_path = os.path.realpath(fpath)
            _write_imageFile_xml(real_path, width, length, dtype, bands)
            rewritten += 1
    logger.info("Rewrote %d geometry XML files", rewritten)


def fix_conncomp_xmls(pair_dirs):
    for pair_dir in pair_dirs:
        pair_id = os.path.basename(pair_dir)
        cc_file = os.path.join(pair_dir, f"{pair_id}_fine.unw.conncomp")
        if not os.path.isfile(cc_file):
            continue
        width  = _detect_width(pair_dir, pair_id) or 1372
        length = _detect_length(pair_dir, pair_id)
        real_path = os.path.realpath(cc_file)
        _write_imageFile_xml(real_path, width, length, "BYTE", 1)


def fix_unw_xmls(pair_dirs):
    for pair_dir in pair_dirs:
        pair_id = os.path.basename(pair_dir)
        ref8    = pair_id[:8]
        sec8    = pair_id[9:]
        date12  = f"{ref8[2:]}-{sec8[2:]}"
        width   = _detect_width(pair_dir, pair_id) or 1372
        length  = _detect_length(pair_dir, pair_id)
        for suffix, bands in [("unw", 2), ("cor", 1)]:
            fpath = os.path.join(pair_dir, f"{pair_id}_fine.{suffix}")
            if not os.path.isfile(fpath):
                continue
            real_path = os.path.realpath(fpath)
            xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<imageFile>
    <property name="width"><value>{width}</value></property>
    <property name="length"><value>{length}</value></property>
    <property name="number_bands"><value>{bands}</value></property>
    <property name="data_type"><value>FLOAT</value></property>
    <property name="scheme"><value>BIL</value></property>
    <property name="byte_order"><value>l</value></property>
    <property name="filename"><value>{real_path}</value></property>
    <property name="date12"><value>{date12}</value></property>
    <property name="ref_date"><value>{ref8}</value></property>
    <property name="sec_date"><value>{sec8}</value></property>
</imageFile>
"""
            with open(real_path + ".xml", "w") as f:
                f.write(xml)


# ─────────────────────────────────────────────
# 3. WRITE .rsc FILES
# ─────────────────────────────────────────────

def write_rsc_files(pair_dirs, pairs_csv, orbit_direction="DESCENDING"):
    bperp_map = _load_bperp_map(pairs_csv)
    written   = 0
    for pair_dir in pair_dirs:
        pair_id = os.path.basename(pair_dir)
        ref8    = pair_id[:8]
        sec8    = pair_id[9:]
        date12  = f"{ref8}_{sec8}"
        bp      = bperp_map.get(pair_id, 0.0)
        width   = _detect_width(pair_dir, pair_id) or 1372
        length  = _detect_length(pair_dir, pair_id) or 0
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
            f"WAVELENGTH               0.05546576\n"
            f"PLATFORM                 SENTINEL1\n"
            f"SENSOR                   SENTINEL1\n"
            f"ORBIT_DIRECTION          {orbit_direction}\n"
            f"RANGE_PIXEL_SIZE         2.329562\n"
            f"AZIMUTH_PIXEL_SIZE       14.1\n"
            f"EARTH_RADIUS             6371000.0\n"
            f"HEIGHT                   693000.0\n"
            f"STARTING_RANGE           845960.0\n"
            f"RLOOKS                   19\n"
            f"ALOOKS                   7\n"
            f"P_BASELINE_TOP_HDR       {bp:.4f}\n"
            f"P_BASELINE_BOTTOM_HDR    {bp:.4f}\n"
        )
        for fname in [
            f"{pair_id}_fine.unw",
            f"{pair_id}_fine.cor",
            f"{pair_id}_fine.unw.conncomp",
        ]:
            fpath    = os.path.join(pair_dir, fname)
            rsc_path = (
                os.path.realpath(fpath) + ".rsc"
                if os.path.islink(fpath)
                else fpath + ".rsc"
            )
            with open(rsc_path, "w") as f:
                f.write(rsc)
            written += 1
    logger.info("Written %d .rsc files", written)


# ─────────────────────────────────────────────
# 4. CREATE BASELINES FOLDER
# ─────────────────────────────────────────────

def create_baselines_dir(pairs_csv, ifg_dir, pair_ids):
    baselines_dir = os.path.join(ifg_dir, "baselines")
    os.makedirs(baselines_dir, exist_ok=True)
    bperp_map = _load_bperp_map(pairs_csv)
    created   = 0
    for pair_id in pair_ids:
        bp     = bperp_map.get(pair_id, 0.0)
        sec8   = pair_id[9:]
        folder = os.path.join(baselines_dir, pair_id)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, f"{sec8}.txt"), "w") as f:
            f.write(f"P_BASELINE_TOP_HDR    {bp:.4f}\n")
            f.write(f"P_BASELINE_BOTTOM_HDR {bp:.4f}\n")
        created += 1
    logger.info("Created %d baseline entries in %s", created, baselines_dir)
    return baselines_dir


# ─────────────────────────────────────────────
# 5. STAGE TRACK PAIRS
# ─────────────────────────────────────────────

def stage_track_pairs(ifg_dir, pair_ids, track_name):
    stage_dir = os.path.join(ifg_dir, f"_stage_{track_name}")
    os.makedirs(stage_dir, exist_ok=True)
    for pid in pair_ids:
        src  = os.path.abspath(os.path.join(ifg_dir, pid))
        dest = os.path.join(stage_dir, pid)
        if os.path.islink(dest):
            os.unlink(dest)
        if not os.path.exists(dest):
            os.symlink(src, dest)
    logger.info("Staged %d pairs for track '%s' -> %s",
                len(pair_ids), track_name, stage_dir)
    return stage_dir


# ─────────────────────────────────────────────
# 6. WRITE MintPy TEMPLATE
# ─────────────────────────────────────────────

def write_mintpy_template(
    ifg_dir, baselines_dir, output_dir, pair_ids, track_label,
    orbit_direction="DESCENDING",
):
    os.makedirs(output_dir, exist_ok=True)
    first_pair    = os.path.join(ifg_dir, pair_ids[0])
    template_path = os.path.join(output_dir, f"mintpy_{track_label}.txt")

    # Default "no" — bridging requires the auto-selected reference point
    # to fall inside a connected component, which frequently fails over
    # vegetated areas or short stacks. Set manually in the template if needed:
    #   mintpy.unwrapError.method = bridging     (needs ref in conncomp)
    #   mintpy.unwrapError.method = phase_closure (needs closed triplets)
    unwrap_method = "no"
    logger.info("Track '%s': unwrapError.method = no "
                "(edit template manually to enable bridging/phase_closure)",
                track_label)

    with open(template_path, "w") as f:
        f.write(f"""##----------- MintPy Template — {track_label.upper()} track -----------##
## Generated by mintpy_processor.py v1.4
## Run with: smallbaselineApp.py {template_path}
##
## Unwrapping error correction is disabled by default ("no").
## Enable if your data supports it:
##   mintpy.unwrapError.method = bridging      (ref point must be in conncomp)
##   mintpy.unwrapError.method = phase_closure (needs closed triplets)
##----------------------------------------------------------------------##

########## 1. Load Data
mintpy.load.processor        = isce

mintpy.load.unwFile          = {ifg_dir}/*/*_fine.unw
mintpy.load.corFile          = {ifg_dir}/*/*_fine.cor
mintpy.load.connCompFile     = {ifg_dir}/*/*_fine.unw.conncomp

mintpy.load.demFile          = {first_pair}/hgt.rdr
mintpy.load.incAngleFile     = {first_pair}/incidenceAngle.rdr
mintpy.load.azAngleFile      = {first_pair}/azimuthAngle.rdr
mintpy.load.lookupYFile      = {first_pair}/lat.rdr
mintpy.load.lookupXFile      = {first_pair}/lon.rdr
mintpy.load.baselineDir      = {baselines_dir}
mintpy.load.waterMaskFile    = auto
mintpy.load.shadowMaskFile   = auto

########## 2. Modify Network
mintpy.network.coherenceBased    = yes
mintpy.network.minCoherence      = {MIN_COHERENCE}
mintpy.network.startDate         = auto
mintpy.network.endDate           = auto
mintpy.network.excludeDate       = no

########## 3. Reference Point
mintpy.reference.maskFile        = maskConnComp.h5
mintpy.reference.minCoherence    = {REF_MIN_COHERENCE}

########## 4. Time Series Inversion
mintpy.networkInversion.weightFunc          = no
mintpy.networkInversion.minNormVelocity     = yes
mintpy.networkInversion.minTempCoh          = {MIN_TEMP_COH}
mintpy.networkInversion.minNumIfgram        = 3
mintpy.networkInversion.maskDataset         = coherence
mintpy.networkInversion.maskThreshold       = {MIN_COHERENCE}

########## 5. Tropospheric Delay Correction
mintpy.troposphericDelay.method             = {TROPO_METHOD}

########## 6. Topographic (DEM) Error Correction
mintpy.topographicResidual                  = yes
mintpy.topographicResidual.pixelwiseGeometry = yes

########## 7. Unwrapping Error Correction
mintpy.unwrapError.method                   = {unwrap_method}
mintpy.unwrapError.connCompMinArea          = 2.5e3

########## 8. Phase Residual RMS
mintpy.residualPhase.maskFile               = maskTempCoh.h5
mintpy.residualPhase.ramp                   = quadratic

########## 9. Velocity Estimation
mintpy.velocity.startDate                   = auto
mintpy.velocity.endDate                     = auto
mintpy.velocity.excludeDate                 = no

########## 10. Geocoding
mintpy.geocode.laloStep                     = {GEOCODE_STEP}
mintpy.geocode.interpMethod                 = nearest
mintpy.geocode.fillValue                    = np.nan

########## 11. Google Earth
mintpy.save.kmz                             = yes

########## 12. HDF-EOS5
mintpy.save.hdfEos5                         = no
""")

    logger.info("Template written -> %s", template_path)
    return template_path


# ─────────────────────────────────────────────
# 7. RUN smallbaselineApp.py
# ─────────────────────────────────────────────

def run_smallbaseline(template_path, work_dir):
    """
    Run MintPy smallbaselineApp.py for one track.

    Discovery order (cross-platform):
      1. shutil.which() — smallbaselineApp.py / .exe / .EXE / bare name
      2. Scripts\\ folder — Windows Anaconda install (handles .EXE naming)
      3. python -m mintpy.cli.smallbaselineApp — always works if pip-installed

    The module fallback (#3) is unconditional — never raises FileNotFoundError
    if mintpy is installed, regardless of PATH configuration.
    """
    import shutil as _shutil
    import sys    as _sys

    os.makedirs(work_dir, exist_ok=True)

    # Purge stale cfg so MintPy always copies a fresh default from its package
    stale_cfg = os.path.join(work_dir, "smallbaselineApp.cfg")
    if os.path.isfile(stale_cfg):
        os.remove(stale_cfg)
        logger.info("Removed stale smallbaselineApp.cfg (will be regenerated)")

    log_file = os.path.join(work_dir, "mintpy.log")

    # 1. PATH
    _cmd = None
    for _cand in [
        "smallbaselineApp.py", "smallbaselineApp.py.exe",
        "smallbaselineApp.py.EXE", "smallbaselineApp",
    ]:
        _found = _shutil.which(_cand)
        if _found:
            _cmd = [_found]
            logger.info("Found smallbaselineApp: %s", _found)
            break

    # 2. Windows Scripts\ folder
    if _cmd is None:
        _scripts_dir = os.path.join(os.path.dirname(_sys.executable), "Scripts")
        for _name in [
            "smallbaselineApp.py", "smallbaselineApp.py.exe",
            "smallbaselineApp.py.EXE", "smallbaselineApp.exe",
        ]:
            _p = os.path.join(_scripts_dir, _name)
            if os.path.isfile(_p):
                _cmd = [_p]
                logger.info("Found smallbaselineApp in Scripts: %s", _p)
                break

    # 3. Module fallback — always works if mintpy is installed
    if _cmd is None:
        _cmd = [_sys.executable, "-m", "mintpy.cli.smallbaselineApp"]
        logger.info("Using module fallback: python -m mintpy.cli.smallbaselineApp")

    try:
        with open(log_file, "w") as log:
            result = subprocess.run(
                _cmd + [template_path],
                cwd=work_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=7200,
            )
        if result.returncode == 0:
            logger.info("MintPy completed successfully -> %s", work_dir)
            return {"status": "success", "work_dir": work_dir}
        else:
            logger.error("MintPy failed — check %s", log_file)
            return {"status": "failed", "log": log_file}
    except subprocess.TimeoutExpired:
        logger.error("MintPy timeout after 2 hours")
        return {"status": "failed", "reason": "timeout"}
    except Exception as e:
        logger.error("MintPy run error: %s", e)
        return {"status": "failed", "reason": str(e)}


# ─────────────────────────────────────────────
# 8. PRINT RESULTS
# ─────────────────────────────────────────────

def print_summary(track_results):
    print("\n" + "=" * 50)
    print("  MINTPY PROCESSING SUMMARY")
    print("=" * 50)
    for track, result in track_results.items():
        status = result.get("status", "unknown")
        wdir   = result.get("work_dir", "")
        print(f"  {track:15s} : {status}")
        if status == "success" and wdir:
            vel = os.path.join(wdir, "velocity.h5")
            kmz = os.path.join(wdir, "geo", "geo_velocity.kmz")
            if os.path.isfile(vel): print(f"    velocity   : {vel}")
            if os.path.isfile(kmz): print(f"    Google KMZ : {kmz}")
    print("=" * 50 + "\n")


# ─────────────────────────────────────────────
# 9. MAIN ENTRY POINT
# ─────────────────────────────────────────────

def run(
    ifg_dir    = IFG_DIR,
    pairs_csv  = PAIRS_CSV,
    output_dir = OUTPUT_DIR,
    track      = TRACK,
):
    """Full MintPy processing pipeline. Called by main_v3.py Step 4."""
    # Resolve all paths to absolute — prevents path truncation bugs
    ifg_dir    = os.path.realpath(ifg_dir)
    pairs_csv  = os.path.realpath(pairs_csv)
    output_dir = os.path.realpath(output_dir)

    logger.info("ifg_dir resolved    : %s", ifg_dir)
    logger.info("output_dir resolved : %s", output_dir)

    if not os.path.isdir(ifg_dir):
        logger.error("ifg_dir does not exist: %s", ifg_dir)
        return {}

    all_tracks = detect_tracks(ifg_dir, pairs_csv)
    if not all_tracks:
        logger.error("No valid pair directories found in %s", ifg_dir)
        return {}

    if track != "all":
        all_tracks = {k: v for k, v in all_tracks.items() if k == track}
        if not all_tracks:
            logger.error("No pairs found for track: %s", track)
            return {}

    logger.info("Tracks to process: %s", list(all_tracks.keys()))

    track_results = {}

    for track_name, pair_ids in all_tracks.items():
        logger.info("Processing track: %s  (%d pairs)", track_name, len(pair_ids))

        orbit_dir    = "ASCENDING" if track_name == "ascending" else "DESCENDING"
        pair_dirs    = [os.path.join(ifg_dir, p) for p in pair_ids]
        track_output = os.path.join(output_dir, f"mintpy_{track_name}")

        # Purge stale artifacts before writing anything new
        stale_stage = os.path.join(ifg_dir, f"_stage_{track_name}")
        if os.path.isdir(stale_stage):
            shutil.rmtree(stale_stage)
            logger.info("Removed stale stage dir: %s", stale_stage)

        stale_tpl = os.path.join(track_output, f"mintpy_{track_name}.txt")
        if os.path.isfile(stale_tpl):
            os.remove(stale_tpl)
            logger.info("Removed stale template: %s", stale_tpl)

        # Write .rsc and baselines BEFORE staging
        write_rsc_files(pair_dirs, pairs_csv=pairs_csv, orbit_direction=orbit_dir)
        baselines_dir = create_baselines_dir(pairs_csv, ifg_dir, pair_ids)

        # Stage symlinks
        stage_dir = stage_track_pairs(ifg_dir, pair_ids, track_name)

        # Fix XMLs AFTER staging (realpath writes through symlinks)
        stage_pair_dirs = [os.path.join(stage_dir, p) for p in pair_ids]
        fix_geometry_xmls(stage_pair_dirs)
        fix_conncomp_xmls(stage_pair_dirs)
        fix_unw_xmls(stage_pair_dirs)

        # Validate stage
        first_stage_pair = os.path.join(stage_dir, pair_ids[0])
        missing_geom = [
            g for g in ("hgt.rdr", "lat.rdr", "lon.rdr",
                         "incidenceAngle.rdr", "azimuthAngle.rdr")
            if not os.path.isfile(os.path.join(first_stage_pair, g))
        ]
        if missing_geom:
            logger.error(
                "Track '%s': geometry files missing — symlinks may be dangling.\n"
                "  Stage: %s\n  Missing: %s\n"
                "  Fix: delete _stage_%s/ and rerun.",
                track_name, first_stage_pair, missing_geom, track_name,
            )
            track_results[track_name] = {
                "status": "failed",
                "reason": f"geometry files missing: {missing_geom}",
            }
            continue

        logger.info("Stage dir validated — geometry files present in %s",
                    first_stage_pair)

        template_path = write_mintpy_template(
            ifg_dir=stage_dir, baselines_dir=baselines_dir,
            output_dir=track_output, pair_ids=pair_ids,
            track_label=track_name, orbit_direction=orbit_dir,
        )

        result = run_smallbaseline(template_path, track_output)
        track_results[track_name] = result

    print_summary(track_results)

    log_path = os.path.join(output_dir, "mintpy_results.json")
    os.makedirs(output_dir, exist_ok=True)
    with open(log_path, "w") as f:
        json.dump(track_results, f, indent=2)
    logger.info("Results log -> %s", log_path)

    return track_results


# ─────────────────────────────────────────────
# 10. CLI
# ─────────────────────────────────────────────

def main():
    import argparse
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%H:%M:%S",
    )
    p = argparse.ArgumentParser(
        description="MintPy SBAS Processor v1.4 — Sentinel-1 / ISCE2 topsApp",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--ifg-dir",    default=IFG_DIR)
    p.add_argument("--pairs-csv",  default=PAIRS_CSV)
    p.add_argument("--output-dir", default=OUTPUT_DIR)
    p.add_argument("--track", default=TRACK,
                   choices=["ascending", "descending", "all"])
    args = p.parse_args()
    run(ifg_dir=args.ifg_dir, pairs_csv=args.pairs_csv,
        output_dir=args.output_dir, track=args.track)


if __name__ == "__main__":
    main()
