"""14_dt_export.py <root> — inputs for the Doñana Digital Twin (Evers et al. 2026, doi:10.32942/x27q4j;
code CC0 at gitlab.lifewatch.dev/lfw024-donana-digital-twin), produced from the HEALPix products.

What the DT reads (read-only inspection of donanadt-shinyapp, 22 Oct 2026):
  (a) pond module vegetation_donana/R/shrubs/functions/shrub_invasion_risk.R: one yearly hydroperiod raster
      cortehydroperiod_nor_<y-1>_<y>.tif, EPSG:32629, projected bilinearly to the shrub grid (EPSG:25829),
      min-max normalised; "future risk" where hydroperiod < 30 days; scenarios = +/-20 % of hydroperiod.
  (b) rabbit / dung-beetle models: monthly flood maps LAST_monthly_floodmaps/<y>_<m>.asc = per-month MAX of the
      LAST masks, resampled by max to a 500 m template (DT border + 1 km).

Products written here:
  1. data/products/monthly_flood.zarr (GRID4EARTH, depth 19): per month Sep 2017 – Aug 2026
     flooded_any (0/1, 255 no obs; the DT's monthly-max rule) and flooded_pct (share of valid S2 observations
     flooded, 0–100, 255 no obs), n_obs.
  2. data/products/dt_export/hydroperiod_<y>_<y+1>.tif — days, Float32, EPSG:32629, 20 m (Sentinel-2 lattice),
     single band like the DT input (nodata -9999); hydroperiod_unc_<y>_<y+1>.tif — half-sample mean absolute
     deviation (days); hydroperiod_nvalid_<y>_<y+1>.tif. 2016/17 is NOT from S2 (no S2 coverage before
     Apr 2017 in this run): use LAST H1617 (Zenodo 3516929) for that year.
  3. data/products/dt_export/monthly_flood_<y>_<y+1>.tif — 12 bands (Sep..Aug), 500 m, EPSG:32629, value =
     flooded area share of the 500 m pixel under the monthly-max rule (0–1), -9999 where < 50 % observed.
  4. data/products/ponds_hydroperiod_series.csv — per OSM pond/lagoon and hydrological year 2016/17–2025/26:
     LAST (2016/17, 2017/18) and S2 (2017/18–2025/26) hydroperiod, n_valid, max gap, half-sample uncertainty.
Raster export = exact transpose of the cell/pixel overlap operator (pixel value = overlap-area-weighted mean
of the cells it intersects); no reprojection kernel is involved.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import zarr
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).parent))
from common import dggs_attrs, load_grid

ROOT = Path(sys.argv[1])
OUT = ROOT / "data/products/dt_export"; OUT.mkdir(parents=True, exist_ok=True)
W, win, cells = load_grid(ROOT)
n = len(cells["cell_ids"]); A = cells["area_utm_m2"].astype(np.float64)
WT = W[20].T.tocsr()
pix_tot = np.asarray(W[20].sum(0)).ravel()
nx, ny = win["nx20"], win["ny20"]
T20 = from_origin(win["x0"], win["y1"], 20, 20)


def to_pixels(v, valid):
    m = valid.astype(np.float32)
    num = WT @ (np.where(valid, v, 0).astype(np.float32) * 1.0)
    den = WT @ m
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where((den > 0) & (den >= 0.5 * pix_tot), num / den, -9999).astype(np.float32)
    return out.reshape(ny, nx)


def write(path, arrs, transform, desc):
    prof = {"driver": "GTiff", "dtype": "float32", "count": len(arrs), "height": arrs[0].shape[0], "width": arrs[0].shape[1],
            "crs": "EPSG:32629", "transform": transform, "nodata": -9999, "compress": "deflate", "predictor": 3, "tiled": True}
    with rasterio.open(path, "w", **prof) as d:
        for i, a in enumerate(arrs, 1):
            d.write(a, i)
        d.update_tags(**desc)


# ---- 2. yearly hydroperiod exports
hz = zarr.open_group(ROOT / "data/products/hydroperiod.zarr", mode="r")
years = hz.attrs["hydrological_years"]
for k, yy in enumerate(years):
    y = int(yy[:4])
    hp = hz["hydroperiod/19/hydroperiod_days"][k].astype(np.float32)
    ok = hp >= 0
    tag = {"source": "Sentinel-2 L2A, Kordelas et al. 2018 open-water rule on HEALPix depth-19 cells (WGS84, NESTED)",
           "hydrological_year": yy, "units": "days", "rule": "nearest valid observation in time",
           "provenance": "LifeWatch ERIC prototype; Copernicus Sentinel data (modified)"}
    write(OUT / f"hydroperiod_{y}_{y + 1}.tif", [to_pixels(hp, ok)], T20, tag)
    unc = hz["hydroperiod/19/halfsample_mad_days"][k].astype(np.float32)
    write(OUT / f"hydroperiod_unc_{y}_{y + 1}.tif", [to_pixels(unc, np.isfinite(unc))], T20, tag | {"units": "days (half-sample mean abs. deviation)"})
    nv = hz["hydroperiod/19/n_valid"][k].astype(np.float32)
    write(OUT / f"hydroperiod_nvalid_{y}_{y + 1}.tif", [to_pixels(nv, nv >= 0)], T20, tag | {"units": "valid S2 observations"})

# ---- 1. monthly flood on HEALPix
s2 = zarr.open_group(ROOT / "data/s2/codes.zarr", mode="r")
sc = pd.read_csv(ROOT / "data/s2/scenes.csv").set_index("date")
dates = pd.to_datetime(pd.Series(s2.attrs["dates"]))
months = pd.period_range("2017-09", "2026-08", freq="M")
g = zarr.open_group(ROOT / "data/products/monthly_flood.zarr", mode="w")
G = g.create_group("19", attributes=dggs_attrs(19) | {"months": [str(m) for m in months],
                   "flooded_any": "1 if any valid S2 observation in the month is open water (DT monthly-max rule); 0 none; 255 no obs",
                   "flooded_pct": "percent of valid S2 observations in the month that are open water; 255 no obs"})
G.create_array("cell_ids", data=cells["cell_ids"], chunks=(1 << 20,), dimension_names=["cells"])
kw = dict(chunks=(1, 1 << 20), dtype="uint8", fill_value=255, dimension_names=["month", "cells"],
          compressors=zarr.codecs.BloscCodec(cname="zstd", clevel=7, shuffle="bitshuffle"))
fa = G.create_array("flooded_any", shape=(len(months), n), **kw)
fp = G.create_array("flooded_pct", shape=(len(months), n), **kw)
no = G.create_array("n_obs", shape=(len(months), n), **kw)
st = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")["stratum"][:]
# 500 m blocks of the 20 m lattice (25 x 25 pixels, origin = window origin)
B = 25; nbx, nby = int(np.ceil(nx / B)), int(np.ceil(ny / B))
T500 = from_origin(win["x0"], win["y1"], 500, 500)


def block(a20, valid20):
    a = np.where(valid20, a20, 0.0); v = valid20.astype(float)
    pa = np.zeros((nby * B, nbx * B)); pv = np.zeros_like(pa)
    pa[:ny, :nx] = a; pv[:ny, :nx] = v
    s = pa.reshape(nby, B, nbx, B).sum((1, 3)); c = pv.reshape(nby, B, nbx, B).sum((1, 3))
    return np.where(c >= 0.5 * B * B, s / np.maximum(c, 1), -9999).astype(np.float32)


yearly = {}
for k, mth in enumerate(months):
    idx = [i for i, d in enumerate(dates) if d.to_period("M") == mth and s2.attrs["dates"][i] in sc.index
           and sc.loc[s2.attrs["dates"][i], "status"] == "ok"]
    nobs = np.zeros(n, np.int16); nw = np.zeros(n, np.int16)
    for i in idx:
        c = s2["codes"][i]; v = c != 255
        nobs += v; nw += v & ((c & 3) == 1)
    has = nobs > 0
    fa[k] = np.where(has, (nw > 0).astype(np.uint8), 255)
    fp[k] = np.where(has, np.round(100.0 * nw / np.maximum(nobs, 1)), 255).astype(np.uint8)
    no[k] = np.minimum(nobs, 254).astype(np.uint8)
    pix = to_pixels((nw > 0).astype(np.float32), has & (st != 0))
    pv = pix >= 0
    yearly.setdefault(mth.year - (mth.month < 9), []).append(block(np.where(pv, pix, 0), pv))
for y, bands in yearly.items():
    write(OUT / f"monthly_flood_{y}_{y + 1}.tif", bands, T500,
          {"bands": "Sep..Aug", "units": "flooded area share of 500 m pixel (monthly-max rule)", "hydrological_year": f"{y}/{str(y + 1)[2:]}",
           "grid": "500 m blocks of the S2 20 m lattice, origin x0/y1 of the processing window, EPSG:32629"})

# ---- 4. per-pond series (LAST 2016/17 & 2017/18 + S2)
stg = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")
lag_id = stg["lagoon_id"][:]
lh = zarr.open_group(ROOT / "data/ref/last_hydroperiod.zarr", mode="r")
cyc = list(lh.attrs["cycles"])
rows = []
for lg in stg.attrs["lagoons"]:
    sel = lag_id == lg["id"]
    if not sel.any():
        continue
    for tag in ("2016/17", "2017/18"):
        if tag in cyc:
            H = lh["hydroperiod_days"][cyc.index(tag)]; ok = sel & (H >= 0)
            if ok.any():
                rows.append({"pond_id": lg["id"], "name": lg["name"], "polygon_area_ha": round(lg["polygon_area_ha"], 3),
                             "n_cells": int(sel.sum()), "hydro_year": tag, "source": "LAST-EBD Landsat (Zenodo 3516929)",
                             "hydroperiod_days": float(np.average(H[ok], weights=A[ok]))})
lp = pd.read_csv(ROOT / "data/products/lagoons_hydroperiod.csv")
for _, r in lp.iterrows():
    rows.append({"pond_id": r.lagoon_id, "name": r["name"], "polygon_area_ha": r.polygon_area_ha, "n_cells": r.n_cells,
                 "hydro_year": r.hydro_year, "source": "Sentinel-2 (this pipeline)", "hydroperiod_days": r.mean_hp_days,
                 "max_cell_hydroperiod_days": r.max_hp_days, "flooded_share_ge1d": r.flooded_share_ge1d,
                 "median_n_valid": r.median_n_valid, "median_max_gap_days": r.get("median_max_gap_days"),
                 "halfsample_mad_days": r.get("halfsample_mad_days")})
pd.DataFrame(rows).sort_values(["name", "hydro_year", "source"]).to_csv(ROOT / "data/products/ponds_hydroperiod_series.csv", index=False)
print("exported", len(list(OUT.glob("*.tif"))), "GeoTIFFs;", len(rows), "pond-year rows")
