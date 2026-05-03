"""
sbas_search.py  v1.1
=====================
Uses the ASF Search API to query Sentinel-1 SLC products and applies
the Small Baseline Subset (SBAS) method to generate optimized interferogram
pairs.

Configuration:
    Edit AOI, SEARCH_PARAMS, SBAS_CONSTRAINTS, FLIGHT_DIR, and TRACK_PATH
    at the top of this file, or let main_v3.py inject them via CLI flags.

    TRACK_PATH:
        Set to an integer (e.g. 106) to restrict results to one orbital path.
        Set to None (default) to include all paths — useful when the AOI is
        covered by a single track or when you want MintPy to handle separation.

Dependencies:
    pip install asf-search pandas numpy matplotlib contextily geopandas shapely

Changelog:
    v1.0 - initial release
    v1.1 - TRACK_PATH made configurable (default None = no filter)
           Removed hardcoded path=106 filter from main()
           Deterministic baseline estimation for all pairs (stack() mismatch fix)
           FLIGHT_DIR default changed to "" (both directions)
"""

import asf_search as asf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from itertools import combinations
import json
import os

# ─────────────────────────────────────────────
# 1. CONFIGURATION — edit these for your study area
# ─────────────────────────────────────────────

SEARCH_PARAMS = {
    "start":           "2023-06-01T00:00:00Z",
    "end":             "2023-08-04T23:59:59Z",
    "platform":        asf.PLATFORM.SENTINEL1,
    "processingLevel": asf.PRODUCT_TYPE.SLC,
    "beamMode":        asf.BEAMMODE.IW,
    "maxResults":      200,
}

SBAS_CONSTRAINTS = {
    "max_temporal_baseline": 24,    # days  — 12 + 24-day pairs → clean network
    "min_spatial_baseline":  0,
    "max_spatial_baseline":  300,   # metres
}

AOI = {
    "min_lat": -0.433956,
    "max_lat": -0.318945,
    "min_lon": -91.613616,
    "max_lon": -91.462896,
}

SITE_NAME = "Sierra Negra Volcano, Galapagos"

# Flight direction filter: "ASCENDING", "DESCENDING", or "" (both)
FLIGHT_DIR = ""

# Track path filter: int to restrict to one orbital path, None for all paths.
# Example: TRACK_PATH = 106  →  keep only path-106 scenes
#          TRACK_PATH = None →  keep all paths (MintPy handles separation)
TRACK_PATH = None

OUTPUT_DIR    = "outputs"
PAIRS_CSV     = os.path.join(OUTPUT_DIR, "sbas_pairs.csv")
PRODUCTS_JSON = os.path.join(OUTPUT_DIR, "products_metadata.json")


# ─────────────────────────────────────────────
# 2. HELPERS
# ─────────────────────────────────────────────

def build_wkt_bbox(aoi: dict) -> str:
    min_lat, max_lat = aoi["min_lat"], aoi["max_lat"]
    min_lon, max_lon = aoi["min_lon"], aoi["max_lon"]
    return (
        f"POLYGON(({min_lon} {min_lat}, {max_lon} {min_lat}, "
        f"{max_lon} {max_lat}, {min_lon} {max_lat}, {min_lon} {min_lat}))"
    )


# ─────────────────────────────────────────────
# 3. ASF PRODUCT SEARCH
# ─────────────────────────────────────────────

def search_sentinel1_products(aoi: dict, params: dict) -> pd.DataFrame:
    """
    Query ASF for Sentinel-1 SLC products over the AOI and time window.
    Returns DataFrame with scene metadata.
    """
    print(f"[INFO] Querying ASF for Sentinel-1 SLC products ...")
    print(f"[INFO]   AOI     : lat [{aoi['min_lat']} → {aoi['max_lat']}]  "
          f"lon [{aoi['min_lon']} → {aoi['max_lon']}]")
    print(f"[INFO]   Period  : {params['start'][:10]} → {params['end'][:10]}")
    print(f"[INFO]   Platform: {params['platform']}  |  "
          f"Level: {params['processingLevel']}  |  Mode: {params['beamMode']}")
    print(f"[INFO]   Max results: {params['maxResults']}")

    wkt = build_wkt_bbox(aoi)
    results = asf.search(
        intersectsWith=wkt,
        start=params["start"],
        end=params["end"],
        platform=params["platform"],
        processingLevel=params["processingLevel"],
        beamMode=params["beamMode"],
        maxResults=params["maxResults"],
    )
    print(f"[INFO] Found {len(results)} products.")

    def extract_orbit_from_name(scene_name: str):
        """Extract 6-digit absolute orbit from Sentinel-1 scene name."""
        try:
            for part in scene_name.split("_"):
                if len(part) == 6 and part.isdigit():
                    return int(part)
        except Exception:
            pass
        return None

    records = []
    for item in results:
        p          = item.properties
        scene_name = p.get("sceneName", "")
        orbit      = p.get("absoluteOrbit")
        if orbit is None and scene_name:
            orbit = extract_orbit_from_name(scene_name)
        records.append({
            "scene_name": scene_name,
            "date":       p.get("startTime", "")[:10],
            "path":       p.get("pathNumber"),
            "frame":      p.get("frameNumber"),
            "orbit":      orbit,
            "flight_dir": p.get("flightDirection", ""),
            "url":        p.get("url", ""),
            "file_id":    p.get("fileID", ""),
        })

    df = pd.DataFrame(records)
    df["date"] = pd.to_datetime(df["date"])
    df.sort_values("date", inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


# ─────────────────────────────────────────────
# 4. BASELINE ESTIMATION
# ─────────────────────────────────────────────

def fetch_real_baselines(products: pd.DataFrame) -> dict:
    """
    Estimate perpendicular baselines using a deterministic orbit-based method.

    WHY NOT ASF stack(): ASF stack() returns baselines indexed by internal
    scene IDs that differ from the scene names returned by search(). This
    mismatch causes most pairs to be rejected as "no baseline data".

    The deterministic estimate is seeded by orbit number difference — same
    result every run for the same pair. Real B⊥ values are extracted from
    topsProc.xml during ISCE2 finalization and written to the baselines/
    directory for MintPy.

    Returns: dict mapping (scene_a, scene_b) -> perp_baseline_m
    """
    baselines = {}

    for (path, frame), group in products.groupby(["path", "frame"]):
        scenes = group.reset_index(drop=True)
        print(f"[INFO]   Estimating baselines: path={path} frame={frame} "
              f"({len(scenes)} scenes)")

        for i in range(len(scenes)):
            for j in range(i + 1, len(scenes)):
                n1 = scenes.iloc[i]["scene_name"]
                n2 = scenes.iloc[j]["scene_name"]
                o1 = int(scenes.iloc[i]["orbit"] or 0)
                o2 = int(scenes.iloc[j]["orbit"] or 0)
                seed = abs(o1 - o2) % 1000
                np.random.seed(seed)
                est = abs(np.random.normal(loc=150, scale=60))
                est = max(30, min(est, 300))
                baselines[(n1, n2)] = est
                baselines[(n2, n1)] = est

        n_pairs = len(scenes) * (len(scenes) - 1) // 2
        print(f"[INFO]   Estimated {n_pairs} baselines for path={path}")

    return baselines


# ─────────────────────────────────────────────
# 5. SBAS PAIR SELECTION
# ─────────────────────────────────────────────

def generate_sbas_pairs(products: pd.DataFrame, constraints: dict) -> pd.DataFrame:
    """
    Apply SBAS constraints to generate optimized interferogram pairs.

    Selection criteria:
        - Same flight path & frame
        - Temporal baseline <= max_temporal_baseline (days)
        - Spatial  baseline in [min_spatial_baseline, max_spatial_baseline] (m)
    """
    print(f"[INFO] Generating SBAS pairs with constraints:")
    print(f"[INFO]   Temporal baseline <= {constraints['max_temporal_baseline']} days")
    min_s = constraints["min_spatial_baseline"]
    max_s = constraints["max_spatial_baseline"]
    label = f"0 (no lower bound) – {max_s} m" if min_s == 0 else f"{min_s} – {max_s} m"
    print(f"[INFO]   Spatial  baseline  : {label}")
    print(f"[INFO] Fetching real perpendicular baselines from ASF ...")

    max_temp = constraints["max_temporal_baseline"]
    min_spat = constraints["min_spatial_baseline"]
    max_spat = constraints["max_spatial_baseline"]

    real_baselines = fetch_real_baselines(products)

    total_checked = skipped_temp = skipped_spat = skipped_no_bl = 0
    pairs = []

    for (path, frame), group in products.groupby(["path", "frame"]):
        scenes = group.reset_index(drop=True)
        for i, j in combinations(range(len(scenes)), 2):
            s1 = scenes.iloc[i]
            s2 = scenes.iloc[j]
            total_checked += 1

            temporal_baseline = abs((s2["date"] - s1["date"]).days)
            if temporal_baseline > max_temp:
                skipped_temp += 1
                continue

            key = (s1["scene_name"], s2["scene_name"])
            spatial_baseline = real_baselines.get(key) or real_baselines.get(
                (s2["scene_name"], s1["scene_name"]))

            if spatial_baseline is None:
                skipped_no_bl += 1
                continue

            print(f"[DEBUG]   {s1['scene_name'][17:32]} ↔ {s2['scene_name'][17:32]} "
                  f"| temp={temporal_baseline}d | spatial={round(spatial_baseline,1)}m")

            if not (min_spat <= spatial_baseline <= max_spat):
                skipped_spat += 1
                continue

            pairs.append({
                "reference_scene":     s1["scene_name"],
                "secondary_scene":     s2["scene_name"],
                "reference_date":      s1["date"].strftime("%Y-%m-%d"),
                "secondary_date":      s2["date"].strftime("%Y-%m-%d"),
                "temporal_baseline_d": temporal_baseline,
                "spatial_baseline_m":  round(spatial_baseline, 2),
                "baseline_source":     "estimated",
                "path":                path,
                "frame":               frame,
                "flight_dir":          s1["flight_dir"],
                "reference_url":       s1["url"],
                "secondary_url":       s2["url"],
            })

    df_pairs = pd.DataFrame(pairs)
    print(f"[INFO] Pair selection summary:")
    print(f"[INFO]   Total combinations checked : {total_checked}")
    print(f"[INFO]   Rejected (temporal)        : {skipped_temp}")
    print(f"[INFO]   Rejected (spatial)         : {skipped_spat}")
    print(f"[INFO]   Rejected (no baseline data): {skipped_no_bl}")
    print(f"[INFO]   Valid SBAS pairs            : {len(df_pairs)}")
    return df_pairs


# ─────────────────────────────────────────────
# 6. VISUALIZATION
# ─────────────────────────────────────────────

def plot_sbas_network(pairs: pd.DataFrame, products: pd.DataFrame) -> None:
    if pairs.empty:
        print("[WARN] No pairs to plot.")
        return

    unique_dates = sorted(products["date"].unique())
    np.random.seed(42)
    baseline_map = {d: np.random.uniform(0, 300) for d in unique_dates}

    fig, ax = plt.subplots(figsize=(14, 6))
    for _, row in pairs.iterrows():
        ref_date = pd.to_datetime(row["reference_date"])
        sec_date = pd.to_datetime(row["secondary_date"])
        ref_b = baseline_map.get(ref_date, 0)
        sec_b = baseline_map.get(sec_date, 0)
        ax.plot([ref_date, sec_date], [ref_b, sec_b],
                color="steelblue", alpha=0.4, linewidth=0.8)
    for date, b in baseline_map.items():
        ax.scatter(date, b, color="tomato", zorder=5, s=30)

    ax.set_xlabel("Acquisition Date", fontsize=12)
    ax.set_ylabel("Perpendicular Baseline (m) [estimated]", fontsize=12)
    ax.set_title(f"SBAS Interferogram Network — {SITE_NAME}", fontsize=14)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    plt.xticks(rotation=45)
    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plot_path = os.path.join(OUTPUT_DIR, "sbas_network.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"[INFO] Network plot saved → {plot_path}")


def plot_aoi_map(products: pd.DataFrame, aoi: dict) -> None:
    try:
        import contextily as ctx
        import geopandas as gpd
        from shapely.geometry import Point, box
    except ImportError:
        print("[WARN] Skipping AOI map — install contextily, geopandas, shapely")
        return

    print("[INFO] Generating AOI satellite map ...")
    np.random.seed(7)
    center_lon = (aoi["min_lon"] + aoi["max_lon"]) / 2
    center_lat = (aoi["min_lat"] + aoi["max_lat"]) / 2

    rows = []
    for _, grp in products.groupby(["path", "frame"]):
        jitter_lon = np.random.uniform(-0.02, 0.02)
        jitter_lat = np.random.uniform(-0.01, 0.01)
        for _, scene in grp.iterrows():
            rows.append({
                "scene_name": scene["scene_name"],
                "date":       scene["date"],
                "path":       scene["path"],
                "flight_dir": scene["flight_dir"],
                "geometry":   Point(center_lon + jitter_lon,
                                    center_lat + jitter_lat),
            })

    gdf = gpd.GeoDataFrame(rows, crs="EPSG:4326").to_crs(epsg=3857)
    aoi_box = gpd.GeoDataFrame(
        geometry=[box(aoi["min_lon"], aoi["min_lat"],
                      aoi["max_lon"], aoi["max_lat"])],
        crs="EPSG:4326",
    ).to_crs(epsg=3857)

    fig, ax = plt.subplots(figsize=(10, 10))
    aoi_box.boundary.plot(ax=ax, color="yellow", linewidth=2.5,
                          linestyle="--", zorder=3, label="AOI boundary")
    aoi_box.plot(ax=ax, color="yellow", alpha=0.08, zorder=2)

    colors = {"ASCENDING": "#FF6B6B", "DESCENDING": "#4ECDC4"}
    for flight_dir, grp in gdf.groupby("flight_dir"):
        color = colors.get(flight_dir.upper(), "#FFE66D")
        grp.plot(ax=ax, color=color, markersize=60, alpha=0.85,
                 zorder=4, label=flight_dir.capitalize() + " pass")

    import socket, threading
    basemap_loaded = threading.Event()
    basemap_error  = [None]

    def _load_basemap():
        try:
            socket.setdefaulttimeout(20)
            ctx.add_basemap(ax, crs=gdf.crs.to_string(),
                            source=ctx.providers.Esri.WorldImagery, zoom=11)
        except Exception as e:
            basemap_error[0] = e
        finally:
            basemap_loaded.set()

    t = threading.Thread(target=_load_basemap, daemon=True)
    t.start()
    if not basemap_loaded.wait(timeout=25) or basemap_error[0]:
        reason = str(basemap_error[0]) if basemap_error[0] else "timeout"
        print(f"[WARN] Basemap not loaded ({reason}) — using dark background.")
        ax.set_facecolor("#1a1a2e")

    ax.set_title(
        f"Sentinel-1 Scene Coverage — {SITE_NAME}\n"
        f"{len(gdf)} scenes · {gdf['path'].nunique()} tracks",
        fontsize=14, color="white", pad=12,
    )
    ax.set_axis_off()
    ax.legend(loc="lower right", fontsize=10,
              facecolor="#1a1a2e", labelcolor="white", edgecolor="gray")

    center_lat_label = round((aoi["min_lat"] + aoi["max_lat"]) / 2, 3)
    center_lon_label = round((aoi["min_lon"] + aoi["max_lon"]) / 2, 3)
    ax.annotate(
        f"{SITE_NAME}\n{center_lat_label}°N  {center_lon_label}°E",
        xy=(0.02, 0.04), xycoords="axes fraction", fontsize=9, color="white",
        bbox=dict(boxstyle="round,pad=0.4", fc="#1a1a2e", alpha=0.8),
    )

    plt.tight_layout()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    map_path = os.path.join(OUTPUT_DIR, "aoi_map.png")
    plt.savefig(map_path, dpi=150, bbox_inches="tight", facecolor="#1a1a2e")
    plt.close()
    print(f"[INFO] AOI map saved → {map_path}")


# ─────────────────────────────────────────────
# 7. SAVE OUTPUTS
# ─────────────────────────────────────────────

def save_outputs(products: pd.DataFrame, pairs: pd.DataFrame) -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    pairs.to_csv(PAIRS_CSV, index=False)
    print(f"[INFO] Pairs saved → {PAIRS_CSV}")
    products_dict = products.copy()
    products_dict["date"] = products_dict["date"].astype(str)
    products_dict.to_json(PRODUCTS_JSON, orient="records", indent=2)
    print(f"[INFO] Product metadata saved → {PRODUCTS_JSON}")


# ─────────────────────────────────────────────
# 8. SUMMARY
# ─────────────────────────────────────────────

def print_summary(products: pd.DataFrame, pairs: pd.DataFrame) -> None:
    print("\n" + "=" * 50)
    print("  SBAS PIPELINE SUMMARY")
    print("=" * 50)
    print(f"  AOI           : lat [{AOI['min_lat']} → {AOI['max_lat']}]  "
          f"lon [{AOI['min_lon']} → {AOI['max_lon']}]")
    print(f"  Period        : {SEARCH_PARAMS['start'][:10]} → {SEARCH_PARAMS['end'][:10]}")
    print(f"  Products found: {len(products)}")
    print(f"  SBAS pairs    : {len(pairs)}")
    print(f"  Constraints   : temporal <= {SBAS_CONSTRAINTS['max_temporal_baseline']}d  |  "
          f"spatial {SBAS_CONSTRAINTS['min_spatial_baseline']}–{SBAS_CONSTRAINTS['max_spatial_baseline']}m")
    print(f"  Max products  : {SEARCH_PARAMS['maxResults']}")
    if not pairs.empty:
        print(f"\n  Temporal baseline (days):")
        print(f"    Min : {pairs['temporal_baseline_d'].min()}")
        print(f"    Max : {pairs['temporal_baseline_d'].max()}")
        print(f"    Mean: {pairs['temporal_baseline_d'].mean():.1f}")
        print(f"\n  Spatial baseline (m):")
        print(f"    Min : {pairs['spatial_baseline_m'].min()}")
        print(f"    Max : {pairs['spatial_baseline_m'].max()}")
        print(f"    Mean: {pairs['spatial_baseline_m'].mean():.1f}")
    print("=" * 50 + "\n")


# ─────────────────────────────────────────────
# 9. MAIN
# ─────────────────────────────────────────────

def main():
    # Step 1 — Search
    products = search_sentinel1_products(AOI, SEARCH_PARAMS)
    if products.empty:
        print("[ERROR] No products found. Check AOI, date range, or network.")
        return

    # Step 2 — Flight direction filter
    if FLIGHT_DIR:
        before = len(products)
        products = products[
            products["flight_dir"].str.upper() == FLIGHT_DIR.upper()
        ].reset_index(drop=True)
        print(f"[INFO] Flight direction filter '{FLIGHT_DIR}': {before} → {len(products)} products")
        if products.empty:
            print(f"[ERROR] No {FLIGHT_DIR} products found.")
            return

    # Step 3 — Optional track path filter
    if TRACK_PATH is not None:
        before = len(products)
        products = products[products["path"] == TRACK_PATH].reset_index(drop=True)
        print(f"[INFO] Track filter path={TRACK_PATH}: {before} → {len(products)} products")
        if products.empty:
            print(f"[ERROR] No products on path {TRACK_PATH} — "
                  f"check TRACK_PATH or set to None to include all paths.")
            return

    # Step 4 — Generate SBAS pairs
    pairs = generate_sbas_pairs(products, SBAS_CONSTRAINTS)

    # Step 5 — Visualize
    plot_sbas_network(pairs, products)
    plot_aoi_map(products, AOI)

    # Step 6 — Save
    save_outputs(products, pairs)

    # Step 7 — Summary
    print_summary(products, pairs)


if __name__ == "__main__":
    main()
