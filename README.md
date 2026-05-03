# InSAR SBAS Pipeline

End-to-end Sentinel-1 InSAR surface deformation monitoring pipeline. Automates the full workflow from satellite data discovery through to MintPy time-series analysis and velocity maps.

```
ASF Search → Copernicus Download → ISCE2 Interferograms → MintPy Time-Series
```

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Credentials Setup](#credentials-setup)
5. [Quick Start](#quick-start)
6. [Step-by-Step Guide](#step-by-step-guide)
7. [Configuration](#configuration)
8. [CLI Reference](#cli-reference)
9. [Output Structure](#output-structure)
10. [Visualization](#visualization)
11. [Known Issues](#known-issues)
12. [Troubleshooting](#troubleshooting)

---

## Overview

The pipeline consists of five Python modules:

| Script | Step | Environment | Description |
|---|---|---|---|
| `main_v3.py` | — | Any | Orchestrator — runs all steps via CLI |
| `sbas_search.py` | 1 | Windows / Linux | ASF query + SBAS pair selection |
| `download_products.py` | 2 | Windows / Linux | Copernicus S3 + ASF download |
| `isce2_processor.py` | 3 | **WSL2 / Linux only** | topsApp + SNAPHU interferograms |
| `mintpy_processor.py` | 4 | Windows / Linux | MintPy SBAS time-series |

> **Note:** ISCE2 (Step 3) requires Linux. On Windows, use WSL2 for Step 3 and run all other steps in your normal Anaconda environment.

---

## Requirements

### Two Conda Environments

#### `mintpy` — Windows Anaconda (Steps 1, 2, 4)

```bash
conda create -n mintpy python=3.10 -y
conda activate mintpy
conda install -c conda-forge mintpy -y
pip install asf-search boto3 python-dotenv tqdm geopandas shapely contextily pandas numpy
```

#### `isce2` — WSL2 Ubuntu (Step 3)

```bash
# Inside WSL2
conda create -n isce2 python=3.9 -y
conda activate isce2
conda install -c conda-forge isce2 snaphu -y
pip install asf-search boto3 python-dotenv tqdm pandas numpy
```

> ISCE2 takes 20–40 minutes to install.

---

## Installation

```bash
git clone https://github.com/your-username/insar-sbas-pipeline.git
cd insar-sbas-pipeline
```

Place your `.env` credentials file (see below) in the project folder.

---

## Credentials Setup

Create a `.env` file in your project/work directory:

```ini
# Copernicus Data Space (for SLC download)
# Register at: https://dataspace.copernicus.eu
COPERNICUS_ACCESS_KEY=your_access_key
COPERNICUS_SECRET_KEY=your_secret_key

# NASA Earthdata (fallback download + DEM)
# Register at: https://urs.earthdata.nasa.gov
NASA_EARTHDATA_TOKEN=your_jwt_token
```

For WSL2 DEM download, also create `~/.netrc`:

```bash
cat > ~/.netrc << 'EOF'
machine urs.earthdata.nasa.gov login your_username password your_password
EOF
chmod 600 ~/.netrc
```

---

## Quick Start

### Galápagos / Sierra Negra Example (9 pairs, 6 acquisitions)

**Step 1 & 2 — Search and Download (Windows CMD, `mintpy` env):**

```cmd
conda activate mintpy
cd "F:\Projects\my_project"

python main_v3.py ^
  --aoi "{\"type\":\"Polygon\",\"coordinates\":[[[-91.61,-0.43],[-91.46,-0.43],[-91.46,-0.32],[-91.61,-0.32],[-91.61,-0.43]]]}" ^
  --start 2023-06-01 --end 2023-08-04 ^
  --skip-isce2 --skip-mintpy
```

**Step 3 — ISCE2 (WSL2, `isce2` env):**

```bash
conda activate isce2
cd "/mnt/f/Projects/my_project"

python main_v3.py \
  --work-dir "/mnt/f/Projects/my_project" \
  --skip-search --skip-download --skip-mintpy \
  --swaths 2 --isce2-workers 1 \
  --dem-path "/mnt/f/Projects/my_project/outputs/dem/dem_merged.hgt"
```

**Step 4 — MintPy (Windows CMD, `mintpy` env):**

```cmd
python main_v3.py ^
  --work-dir "F:\Projects\my_project" ^
  --skip-search --skip-download --skip-isce2 ^
  --mintpy-track ascending
```

---

## Step-by-Step Guide

### Step 1 — ASF Search & Pair Selection

Queries ASF for Sentinel-1 SLC products over your AOI and selects SBAS pairs within your temporal and spatial baseline constraints.

Configure `sbas_search.py`:

```python
AOI = {
    "min_lat": -0.433956,
    "max_lat": -0.318945,
    "min_lon": -91.613616,
    "max_lon": -91.462896,
}

SEARCH_PARAMS = {
    "start":           "2023-06-01T00:00:00Z",
    "end":             "2023-08-04T23:59:59Z",
    "platform":        asf.PLATFORM.SENTINEL1,
    "processingLevel": asf.PRODUCT_TYPE.SLC,
    "beamMode":        asf.BEAMMODE.IW,
    "maxResults":      200,
}

SBAS_CONSTRAINTS = {
    "max_temporal_baseline": 24,   # days
    "max_spatial_baseline":  300,  # metres
}

FLIGHT_DIR = "ASCENDING"   # "" = both directions
TRACK_PATH = None          # int to filter one path, None = all paths
```

Run:

```cmd
python main_v3.py --work-dir "F:\my_project" --skip-download --skip-isce2 --skip-mintpy
```

**Outputs:** `outputs/sbas_pairs.csv`, `outputs/sbas_network.png`, `outputs/aoi_map.png`

---

### Step 2 — Copernicus Download

Downloads the unique SLC scenes identified in `sbas_pairs.csv` from Copernicus Data Space S3. Falls back to ASF if a scene is not available on S3 (common for data older than ~12 months).

```cmd
python main_v3.py --work-dir "F:\my_project" --skip-search --skip-isce2 --skip-mintpy
```

**Outputs:** `outputs/slc_products/ascending/*.SAFE` (one per scene, ~4–8 GB each)

---

### Step 3 — ISCE2 Interferogram Generation

> **WSL2 required.** Run inside Ubuntu terminal with `isce2` conda env.

Processes each SBAS pair through:
1. `topsApp.py --end=filter` (coregistration, filtering)
2. `snaphu` (phase unwrapping)
3. Geometry merge (lat/lon/hgt/los)
4. XML finalization and cleanup

```bash
# Map Windows drive — replace f with your drive letter
ls /mnt/f/Projects/my_project

# Run ISCE2
conda activate isce2
python main_v3.py \
  --work-dir "/mnt/f/Projects/my_project" \
  --skip-search --skip-download --skip-mintpy \
  --swaths 2 --isce2-workers 1 \
  --dem-path "/path/to/dem_merged.hgt"
```

**DEM:** If you don't have a DEM, omit `--dem-path` and the pipeline will attempt to auto-download via ISCE2's `dem.py`. If that fails, download tiles manually:

```bash
# AWS (no auth required)
wget "https://s3.amazonaws.com/elevation-tiles-prod/skadi/N00/N00W092.hgt.gz"
gunzip N00W092.hgt.gz
gdal_merge.py -o dem_merged.hgt -of ENVI N00W092.hgt S01W092.hgt
```

**Outputs:** `outputs/interferograms/YYYYMMDD_YYYYMMDD/` per pair, containing:
- `*_fine.unw` — unwrapped phase
- `*_fine.cor` — coherence
- `*_fine.unw.conncomp` — connected components
- `hgt.rdr`, `lat.rdr`, `lon.rdr`, `incidenceAngle.rdr`, `azimuthAngle.rdr`

**Estimated runtime:** ~2–3 hours per pair, ~20 hours for 9 pairs (single worker, WSL2)

---

### Step 4 — MintPy SBAS Time-Series

> **Windows CMD, `mintpy` env.**

```cmd
python main_v3.py ^
  --work-dir "F:\my_project" ^
  --skip-search --skip-download --skip-isce2 ^
  --mintpy-track ascending
```

**Outputs:** `outputs/mintpy/mintpy_ascending/`

| File | Description |
|---|---|
| `inputs/ifgramStack.h5` | Full interferogram stack |
| `inputs/geometryRadar.h5` | Geometry (lat, lon, height, angles) |
| `timeseries.h5` | Displacement time series (radar coords) |
| `timeseries_demErr.h5` | Time series after DEM error correction |
| `velocity.h5` | Linear LOS velocity map |
| `temporalCoherence.h5` | Per-pixel temporal coherence |
| `geo/geo_velocity.h5` | Geocoded velocity |
| `geo/geo_velocity.kmz` | Google Earth KMZ |

---

## Configuration

### SBAS Network Planning

| Acquisitions | 12-day pairs | 24-day pairs | Total |
|:---:|:---:|:---:|:---:|
| 4 | 3 | 2 | 5 |
| 5 | 4 | 3 | 7 |
| **6** | **5** | **4** | **9** |
| 8 | 7 | 6 | 13 |

For a well-connected network: use `max_temporal_baseline = 24` days with 6+ acquisitions.

### Unwrapping Error Correction

The template defaults to `mintpy.unwrapError.method = no` because `bridging` requires the auto-selected reference point to fall inside a connected component — this fails frequently over vegetated areas. To enable manually:

```ini
# In mintpy_ascending.txt — edit after pipeline writes the template
mintpy.unwrapError.method = bridging      # ref point must be inside conncomp
mintpy.unwrapError.method = phase_closure # needs closed triplets (3+ pairs)
```

### Temporal Coherence Threshold

For areas with high vegetation (e.g. Galápagos), lower the threshold:

```ini
mintpy.networkInversion.minTempCoh = 0.3   # default 0.7
mintpy.networkInversion.maskDataset = no   # disable coherence masking
```

Edit `smallbaselineApp.cfg` inside the `mintpy_ascending/` work directory.

---

## CLI Reference

```
main_v3.py [options]

Pipeline control:
  --skip-search        Skip Step 1 (use existing sbas_pairs.csv)
  --skip-download      Skip Step 2 (use existing SLC files)
  --skip-isce2         Skip Step 3 (use existing interferograms)
  --skip-mintpy        Skip Step 4

AOI & dates:
  --aoi STR            GeoJSON Polygon string or .geojson file path
  --start YYYY-MM-DD   Search start date
  --end   YYYY-MM-DD   Search end date

Paths:
  --work-dir DIR       Root directory (all outputs under here)
  --dem-path FILE      DEM .hgt file (auto-downloaded if omitted)

ISCE2:
  --swaths STR         IW subswaths: "2" (fast) or "1,2,3" (full)
  --isce2-workers N    Parallel workers (keep at 1 for WSL2)
  --no-cleanup         Keep intermediate files (~60 GB/pair)

MintPy:
  --mintpy-track STR   ascending / descending / all
  --flight-dir STR     ascending / descending / "" (both)
```

---

## Output Structure

```
my_project/
├── .env                          ← credentials
├── main_v3.py
├── sbas_search.py
├── download_products.py
├── isce2_processor.py
├── mintpy_processor.py
└── outputs/
    ├── sbas_pairs.csv            ← Step 1 output
    ├── sbas_network.png
    ├── aoi_map.png
    ├── products_metadata.json
    ├── slc_products/
    │   └── ascending/
    │       └── S1A_IW_SLC__*.SAFE/   ← Step 2 output
    ├── interferograms/
    │   ├── YYYYMMDD_YYYYMMDD/    ← Step 3 output (one per pair)
    │   │   ├── *_fine.unw
    │   │   ├── *_fine.cor
    │   │   ├── *_fine.unw.conncomp
    │   │   ├── hgt.rdr
    │   │   ├── lat.rdr / lon.rdr
    │   │   ├── incidenceAngle.rdr
    │   │   └── azimuthAngle.rdr
    │   ├── baselines/
    │   └── _stage_ascending/     ← MintPy staging (auto-created)
    └── mintpy/
        └── mintpy_ascending/     ← Step 4 output
            ├── inputs/
            │   ├── ifgramStack.h5
            │   └── geometryRadar.h5
            ├── timeseries.h5
            ├── timeseries_demErr.h5
            ├── velocity.h5
            ├── temporalCoherence.h5
            ├── geo/
            │   ├── geo_velocity.h5
            │   └── geo_velocity.kmz
            └── pic/              ← auto-saved plots
```

---

## Visualization

Run from the `mintpy_ascending/` work directory:

```cmd
cd "F:\my_project\outputs\mintpy\mintpy_ascending"

REM Interactive time series — click any pixel to plot its displacement history
tsview.py timeseries_demErr.h5 --dem inputs\geometryRadar.h5 --mask maskTempCoh.h5 -u cm

REM Save velocity map
view.py geo\geo_velocity.h5 velocity -u cm/yr --colormap RdYlBu_r --figsize 12 8 --dpi 150 --nodisplay --save --outfile velocity_map.png

REM Save displacement time series panels (all epochs)
view.py timeseries_demErr.h5 -u cm --figsize 16 5 --dpi 150 --noaxis --nodisplay --save --outfile displacement_timeseries.png

REM Network connectivity plot
plot_network.py inputs\ifgramStack.h5 --nodisplay --save --outfile network.png

REM Google Earth time series KMZ
save_kmz_timeseries.py geo\geo_timeseries_demErr.h5 --mask geo\geo_maskTempCoh.h5 -u cm -o timeseries.kmz
```

---

## Known Issues

### WSL2 → NTFS Symlinks (Access Denied)

WSL2 symlinks into Windows NTFS paths may fail with `[Errno 1] Operation not permitted` when Windows Developer Mode is not enabled. The pipeline handles this gracefully — `ensure_safe_symlinks()` skips on failure, and `find_slc_path()` locates `.SAFE` directories directly.

### ISCE2 2.6.3 — `filename` Attribute Regression

ISCE2 2.6.3 removed the `filename` attribute from the Image class but still writes it to XML. `Image.load()` logs a non-fatal error that causes `check=True` subprocess calls to fail even on success. The pipeline avoids the ISCE2 Python API entirely — `_read_image_dims()` parses XML directly using stdlib `xml.etree`.

### MintPy — Doc-Only XML Properties

Native topsApp XMLs contain `<property>` elements with only a `<doc>` child and no `<value>` child. MintPy's `read_isce_xml()` calls `child.find('value').text` which crashes with `AttributeError: 'NoneType'.text`. `fix_geometry_xmls()` unconditionally rewrites all geometry XMLs to a clean minimal format before MintPy runs.

### Dimension Mismatch (1555 vs 1556 rows)

Sentinel-1 burst coverage can vary by 1 row across acquisitions with different precise orbit files. If MintPy raises `TypeError: Can't broadcast (1555, 1291) -> (1556, 1291)`, use `pad_mismatched_pairs.py` to pad the affected pairs.

### ASF Baseline API Mismatch

`stack()` returns baselines indexed by internal scene IDs that differ from scene names returned by `search()`. The pipeline uses a deterministic orbit-based estimate for all pairs instead. Real B⊥ values are extracted from `topsProc.xml` by ISCE2 and written to `baselines/` for MintPy.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `filt_topophase.flat not found` | topsApp.xml property names wrong | Ensure `<component name="reference">` nesting |
| `No annotation xml file found` | Path with spaces fed to ISCE2 | Use symlink with no spaces: `sudo ln -s "/path with spaces" /mnt/x/nospace` |
| `NoneType has no attribute text` | Doc-only XML properties | Upgrade to latest `mintpy_processor.py` — fix_geometry_xmls() handles this |
| `reference point NOT in connectComponent` | bridging enabled, ref outside conncomp | Set `mintpy.unwrapError.method = no` in template |
| `Not enough reliable pixels` | minTempCoh too high | Lower to 0.3 or 0.0; set maskDataset = no |
| `smallbaselineApp.py not found` | Windows PATH issue | Update to v1.4 — uses module fallback automatically |
| `Operation not permitted` (DEM copy) | WSL2 → NTFS metadata | Fixed in v1.4 — uses `copyfile` not `copy2` |
| Pairs `already_exists` skipped | Previous run left partial outputs | Delete the pair directory and rerun |

---

## License

MIT License — see `LICENSE` for details.
