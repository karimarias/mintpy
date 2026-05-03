"""
main_v3.py  —  InSAR SBAS Pipeline Orchestrator
================================================
Runs the full end-to-end Sentinel-1 InSAR SBAS workflow:

    STEP 1 — ASF Search & SBAS Pairing      (sbas_search.py)
    STEP 2 — Copernicus S3 Product Download  (download_products.py)
    STEP 3 — ISCE2 Interferogram Generation  (isce2_processor.py)  ← WSL2
    STEP 4 — MintPy SBAS Time-Series         (mintpy_processor.py)

Quick-start (Windows — Steps 1 & 2):
    python main_v3.py \
        --aoi '{"type":"Polygon","coordinates":[[[-91.61,-0.43],[-91.46,-0.43],[-91.46,-0.32],[-91.61,-0.32],[-91.61,-0.43]]]}' \
        --start 2023-06-01 --end 2023-08-04 \
        --skip-isce2 --skip-mintpy

ISCE2 only (WSL2 — Step 3):
    python main_v3.py \
        --work-dir "/mnt/f/Projects/my_project" \
        --skip-search --skip-download --skip-mintpy \
        --swaths 2 --isce2-workers 1 \
        --dem-path "/path/to/dem_merged.hgt"

MintPy only (Windows — Step 4):
    python main_v3.py \
        --work-dir "F:\\Projects\\my_project" \
        --skip-search --skip-download --skip-isce2 \
        --mintpy-track ascending

Credentials are loaded from a .env file in --work-dir:
    COPERNICUS_ACCESS_KEY=your_key
    COPERNICUS_SECRET_KEY=your_secret
    NASA_EARTHDATA_TOKEN=your_token
"""

import sys
import os
import logging
import time
import argparse
import json

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="InSAR SBAS Pipeline v3 — Sentinel-1 / ISCE2 / MintPy",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--aoi",   type=str, default=None,
                        help="GeoJSON Polygon string or path to .geojson file.")
    parser.add_argument("--start", type=str, default=None,
                        help="Start date (YYYY-MM-DD).")
    parser.add_argument("--end",   type=str, default=None,
                        help="End date (YYYY-MM-DD).")
    parser.add_argument("--skip-search",   action="store_true")
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-isce2",    action="store_true")
    parser.add_argument("--skip-mintpy",   action="store_true")
    parser.add_argument("--work-dir",      type=str, default=None)
    parser.add_argument("--pairs-csv",     type=str, default=None)
    parser.add_argument("--slc-dir",       type=str, default=None)
    parser.add_argument("--isce2-output-dir", type=str, default=None)
    parser.add_argument("--dem-path",      type=str, default="")
    parser.add_argument("--isce2-workers", type=int, default=1)
    parser.add_argument("--isce2-batch-size", type=int, default=1)
    parser.add_argument("--swaths",        type=str, default="2")
    parser.add_argument("--no-cleanup",    action="store_true")
    parser.add_argument("--mintpy-output-dir", type=str, default=None)
    parser.add_argument("--mintpy-track",  type=str, default="all",
                        choices=["ascending", "descending", "all"])
    parser.add_argument("--flight-dir",    type=str, default="",
                        choices=["ascending", "descending", ""])
    return parser.parse_args()


def parse_aoi(aoi_str: str) -> dict:
    clean = aoi_str.strip("'\"")
    if os.path.isfile(clean):
        with open(clean) as f:
            aoi_data = json.load(f)
    else:
        aoi_data = json.loads(clean)
    if aoi_data.get("type") == "Feature":
        aoi_data = aoi_data["geometry"]
    elif aoi_data.get("type") == "FeatureCollection":
        aoi_data = aoi_data["features"][0]["geometry"]
    coords = aoi_data["coordinates"][0]
    lons = [c[0] for c in coords]
    lats = [c[1] for c in coords]
    return {
        "min_lon": min(lons), "max_lon": max(lons),
        "min_lat": min(lats), "max_lat": max(lats),
    }


def update_search_config(args, sbas_search) -> dict:
    if args.start:
        if hasattr(sbas_search, "SEARCH_PARAMS"):
            sbas_search.SEARCH_PARAMS["start"] = f"{args.start}T00:00:00Z"
        for attr in ["start_date", "START", "start"]:
            if hasattr(sbas_search, attr):
                setattr(sbas_search, attr, args.start); break
    if args.end:
        if hasattr(sbas_search, "SEARCH_PARAMS"):
            sbas_search.SEARCH_PARAMS["end"] = f"{args.end}T23:59:59Z"
        for attr in ["end_date", "END", "end"]:
            if hasattr(sbas_search, attr):
                setattr(sbas_search, attr, args.end); break
    if hasattr(args, "flight_dir") and args.flight_dir:
        sbas_search.FLIGHT_DIR = args.flight_dir.upper()
        logger.info("Flight direction filter: %s", args.flight_dir.upper())
    aoi = getattr(sbas_search, "AOI", None)
    if args.aoi:
        try:
            aoi = parse_aoi(args.aoi)
            sbas_search.AOI = aoi
            if hasattr(sbas_search, "SITE_NAME"):
                sbas_search.SITE_NAME = "Custom CLI AOI"
            logger.info("AOI: lat [%.4f -> %.4f]  lon [%.4f -> %.4f]",
                        aoi["min_lat"], aoi["max_lat"],
                        aoi["min_lon"], aoi["max_lon"])
        except Exception as e:
            logger.error("Failed to parse --aoi: %s", e); sys.exit(1)
    return aoi


def run_step1(args):
    logger.info("STEP 1: ASF Search & SBAS Pairing ...")
    try:
        import sbas_search
        update_search_config(args, sbas_search)
        sbas_search.main()
        logger.info("STEP 1 completed successfully.")
    except ImportError:
        logger.error("sbas_search.py not found in project directory."); sys.exit(1)
    except Exception as e:
        logger.error("STEP 1 failed: %s", e); sys.exit(1)


def run_step2(args, aoi=None):
    logger.info("STEP 2: Copernicus S3 Product Download ...")
    logger.info("  Download dir: %s", args.slc_dir)
    try:
        import download_products
        download_products.PAIRS_CSV     = args.pairs_csv
        download_products.PRODUCTS_JSON = os.path.join(
            os.path.dirname(args.pairs_csv), "products_metadata.json")
        download_products.DOWNLOAD_LOG  = os.path.join(
            os.path.dirname(args.pairs_csv), "download_log.json")
        dem_output_dir = None
        if not args.dem_path and aoi:
            root = args.work_dir or "."
            dem_output_dir = os.path.join(root, "dem")
            logger.info("  DEM will be auto-downloaded to: %s", dem_output_dir)
        download_products.main(
            download_dir=args.slc_dir, aoi=aoi, dem_output_dir=dem_output_dir)
        if dem_output_dir and not args.dem_path:
            dem_candidate = os.path.join(dem_output_dir, "dem_merged.hgt")
            if os.path.isfile(dem_candidate):
                args.dem_path = dem_candidate
                logger.info("  DEM path set: %s", args.dem_path)
        logger.info("STEP 2 completed successfully.")
    except ImportError:
        logger.error("download_products.py not found."); sys.exit(1)
    except Exception as e:
        logger.error("STEP 2 failed: %s", e)
        logger.error("Check credentials in .env: COPERNICUS_ACCESS_KEY, COPERNICUS_SECRET_KEY")
        sys.exit(1)


def run_step3(args, aoi):
    logger.info("STEP 3: ISCE2 Interferogram Generation ...")
    logger.info("  SLC dir    : %s", args.slc_dir)
    logger.info("  Output dir : %s", args.isce2_output_dir)
    logger.info("  DEM        : %s", args.dem_path or "(auto-download)")
    logger.info("  Swaths     : %s", args.swaths)
    logger.info("  Workers    : %d (WSL2: always 1)", args.isce2_workers)
    try:
        import isce2_processor
        if args.no_cleanup:
            isce2_processor.CLEANUP_INTERMEDIATE = False
        swaths = [int(s.strip()) for s in args.swaths.split(",")]
        isce2_processor.run(
            pairs_csv=args.pairs_csv, slc_dir=args.slc_dir,
            output_dir=args.isce2_output_dir, dem_path=args.dem_path,
            num_workers=args.isce2_workers, batch_size=args.isce2_batch_size,
            aoi=aoi, swaths=swaths)
        logger.info("STEP 3 completed successfully.")
    except ImportError:
        logger.error("isce2_processor.py not found."); sys.exit(1)
    except Exception as e:
        logger.error("STEP 3 failed: %s", e); sys.exit(1)


def run_step4(args):
    logger.info("STEP 4: MintPy SBAS Time-Series ...")
    logger.info("  IFG dir    : %s", args.isce2_output_dir)
    logger.info("  Pairs CSV  : %s", args.pairs_csv)
    logger.info("  Output dir : %s", args.mintpy_output_dir)
    logger.info("  Track      : %s", args.mintpy_track)
    try:
        import mintpy_processor
        results = mintpy_processor.run(
            ifg_dir=args.isce2_output_dir, pairs_csv=args.pairs_csv,
            output_dir=args.mintpy_output_dir, track=args.mintpy_track)
        if not results:
            logger.error("STEP 4 failed — mintpy_processor returned no results.")
            sys.exit(1)
        failed  = [t for t, r in results.items() if r.get("status") != "success"]
        success = [t for t, r in results.items() if r.get("status") == "success"]
        if success:
            logger.info("STEP 4 completed for track(s): %s", success)
        if failed:
            logger.warning("STEP 4 failed for track(s): %s — check mintpy.log", failed)
    except ImportError:
        logger.error("mintpy_processor.py not found."); sys.exit(1)
    except Exception as e:
        logger.error("STEP 4 failed: %s", e); sys.exit(1)


def resolve_paths(args):
    root = args.work_dir or "."
    if args.pairs_csv is None:
        args.pairs_csv = os.path.join(root, "outputs", "sbas_pairs.csv")
    if args.slc_dir is None:
        args.slc_dir = os.path.join(root, "outputs", "slc_products")
    if args.isce2_output_dir is None:
        args.isce2_output_dir = os.path.join(root, "outputs", "interferograms")
    if args.mintpy_output_dir is None:
        args.mintpy_output_dir = os.path.join(root, "outputs", "mintpy")
    os.makedirs(os.path.join(root, "outputs"), exist_ok=True)
    if args.work_dir:
        os.chdir(args.work_dir)
    env_file = os.path.join(root, ".env")
    if os.path.isfile(env_file):
        try:
            from dotenv import load_dotenv
            load_dotenv(env_file, override=True)
            logger.info("Loaded credentials from %s", env_file)
        except ImportError:
            with open(env_file) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#") and "=" in line:
                        k, v = line.split("=", 1)
                        os.environ[k.strip()] = v.strip()
            logger.info("Loaded credentials from %s (manual)", env_file)
    return args


def main():
    args = parse_args()
    args = resolve_paths(args)

    logger.info("=" * 50)
    logger.info("  InSAR SBAS Pipeline  —  v3")
    logger.info("=" * 50)
    steps = []
    steps.append("Search"   if not args.skip_search   else "[SKIP] Search")
    steps.append("Download" if not args.skip_download else "[SKIP] Download")
    steps.append("ISCE2"    if not args.skip_isce2    else "[SKIP] ISCE2")
    steps.append("MintPy"   if not args.skip_mintpy   else "[SKIP] MintPy")
    logger.info("  Steps: %s", " -> ".join(steps))

    # AOI resolved from sbas_search.AOI or --aoi flag — no hardcoded default
    aoi = None
    start_time = time.time()

    if not args.skip_search:
        run_step1(args)
        try:
            import sbas_search
            aoi = sbas_search.AOI
        except Exception:
            pass
        print("-" * 50)
    else:
        logger.info("STEP 1 skipped — using existing %s", args.pairs_csv)
        if args.aoi:
            try: aoi = parse_aoi(args.aoi)
            except Exception: pass

    if not args.skip_download:
        run_step2(args, aoi=aoi); print("-" * 50)
    else:
        logger.info("STEP 2 skipped — using existing SLC files in %s", args.slc_dir)

    if not args.skip_isce2:
        run_step3(args, aoi or {}); print("-" * 50)
    else:
        logger.info("STEP 3 skipped")

    if not args.skip_mintpy:
        run_step4(args); print("-" * 50)
    else:
        logger.info("STEP 4 skipped")

    elapsed = time.time() - start_time
    logger.info("=" * 50)
    logger.info("  PIPELINE FINISHED in %.1f minutes.", elapsed / 60)
    logger.info("=" * 50)


if __name__ == "__main__":
    main()
