"""13_wcd_validate.py <root> — Copernicus Land "Water Cover Duration" (WCD, HR-WSI, 10 m, days of water
per hydrological year Sep–Aug, Sentinel-1 + Sentinel-2) on the same HEALPix d19 cells, compared with
(a) LAST-EBD Landsat masks turned into hydroperiod with our rule (2017/18, 2018/19) and LAST annual map H1718,
(b) the Sentinel-2 Kordelas-2018 replication (this pipeline, 2017/18 – 2025/26), by stratum
(seasonal marsh = turbid marsh; lagoons; other).

Access (checked 22 Oct 2026): product metadata doi:10.2909/39e1ff33-01d9-45c2-9edc-00924383f34f; files on the
Copernicus Data Space Ecosystem (CDSE) only, e.g. /eodata/CLMS/bio-geophysical/water_extent/
clms_wsi_water-cover-duration_europe_utm_10m_yearly_v1/<yyyy>/09/01/CLMS_WSI_WCD_010m_T29SQB_<yyyy>0901P1Y_COMB_V100.
OData catalogue search is public; the download (https://download.dataspace.copernicus.eu/odata/v1/Products(<Id>)/$value)
needs a CDSE access token (free registration) -> HTTP 401 without it. CDSE STAC lists the collection but no items yet.
Licence: Copernicus full, free and open data policy (Reg. (EU) 1159/2013): credit the source, state modifications.
Set CDSE_TOKEN (a valid access token) in the environment; nothing is stored in the repository.

Resampling: WCD days are an intensive per-area quantity -> area-weighted (exact overlap, same operator as S2,
script 01) onto cells; values > 366 treated as no data (flag values per product user manual).
Outputs: data/raw/clms_wcd/*, data/products/wcd.zarr, results_wcd.json
"""
import json, os, sys, zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
import zarr
from rasterio.windows import from_bounds

sys.path.insert(0, str(Path(__file__).parent))
from common import conservative, dggs_attrs, jdump, load_grid

ROOT = Path(sys.argv[1])
RAW = ROOT / "data/raw/clms_wcd"; RAW.mkdir(parents=True, exist_ok=True)
CAT = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DL = "https://download.dataspace.copernicus.eu/odata/v1/Products({})/$value"
YEARS = list(range(2016, 2026))


def fetch():
    tok = os.environ.get("CDSE_TOKEN")
    found = []
    for y in YEARS:
        for t in ("T29SQA", "T29SQB"):
            name = f"CLMS_WSI_WCD_010m_{t}_{y}0901P1Y_COMB_V100"
            dest = RAW / f"{name}.zip"
            r = requests.get(CAT, params={"$filter": f"Name eq '{name}'"}, timeout=60); r.raise_for_status()
            v = r.json()["value"]
            found.append({"name": name, "in_catalogue": bool(v), "id": v[0]["Id"] if v else None})
            if not v or dest.exists():
                continue
            if not tok:
                continue
            with requests.get(DL.format(v[0]["Id"]), headers={"Authorization": f"Bearer {tok}"}, stream=True, timeout=600) as d:
                d.raise_for_status()
                with open(dest, "wb") as f:
                    for b in d.iter_content(1 << 20):
                        f.write(b)
    return found


found = fetch()
zips = sorted(RAW.glob("*.zip"))
res = {"catalogue": found, "downloaded": [z.name for z in zips]}
if not zips:
    res["status"] = "not downloaded: CDSE_TOKEN not set (download returns HTTP 401 without a CDSE account token)"
    jdump(res, ROOT / "results_wcd.json"); print(res["status"]); sys.exit(0)

W, win, cells = load_grid(ROOT)
n = len(cells["cell_ids"]); A = cells["area_utm_m2"].astype(float)
st = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")["stratum"][:]
STR = {"marsh": st == 1, "lagoon": st == 2, "other": np.isin(st, (3, 4))}
g = zarr.open_group(ROOT / "data/products/wcd.zarr", mode="w")
G = g.create_group("19", attributes=dggs_attrs(19) | {"source": "CLMS HR-WSI Water Cover Duration 10 m (Copernicus), resampled area-weighted"})
G.create_array("cell_ids", data=cells["cell_ids"], chunks=(1 << 20,), dimension_names=["cells"])
wcd = G.create_array("water_cover_duration_days", shape=(len(YEARS), n), chunks=(1, 1 << 20), dtype="float32",
                     fill_value=np.nan, dimension_names=["year", "cells"])
for k, y in enumerate(YEARS):
    val = np.full(n, np.nan, np.float32); cov = np.zeros(n, np.float32)
    for t in ("T29SQB", "T29SQA"):
        z = RAW / f"CLMS_WSI_WCD_010m_{t}_{y}0901P1Y_COMB_V100.zip"
        if not z.exists():
            continue
        with zipfile.ZipFile(z) as zz:
            tif = [m for m in zz.namelist() if m.endswith(".tif") and "WCD" in m.upper() and "QA" not in m.upper().split("_")[-1]][0]
            zz.extract(tif, RAW / "tmp")
        with rasterio.open(RAW / "tmp" / tif) as r:
            assert r.crs.to_epsg() == 32629 and abs(r.res[0] - 10) < 1e-6
            img = np.zeros((win["ny10"], win["nx10"]), np.float32)
            w = from_bounds(win["x0"], win["y0"], win["x1"], win["y1"], r.transform).round_offsets().round_lengths()
            a = r.read(1, window=w, boundless=True, fill_value=65535).astype(np.float32)
            img[:a.shape[0], :a.shape[1]] = a
        (RAW / "tmp" / tif).unlink()
        ok = img <= 366
        v, f = conservative(W[10], img, ok)
        take = (f > cov) & (f >= 0.999)
        val[take] = v[take]; cov[take] = f[take]
    wcd[k] = val

hz = zarr.open_group(ROOT / "data/products/hydroperiod.zarr", mode="r")
hyears = [int(t[:4]) for t in hz.attrs["hydrological_years"]]


def compare(a, b, sel):
    ok = sel & np.isfinite(a) & (b >= 0)
    if ok.sum() == 0:
        return {}
    d = a[ok] - b[ok]; w = A[ok]
    fa, fb = a[ok] >= 15, b[ok] >= 15   # 'flooded at least ~2 weeks' class for omission/commission
    tp = w[fa & fb].sum(); fp = w[fa & ~fb].sum(); fn = w[~fa & fb].sum(); tn = w[~fa & ~fb].sum()
    N = tp + fp + fn + tn; po = (tp + tn) / N; pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / N ** 2
    return {"n_cells": int(ok.sum()), "mean_wcd": float(np.average(a[ok], weights=w)), "mean_ref": float(np.average(b[ok], weights=w)),
            "bias_days": float(np.average(d, weights=w)), "mae_days": float(np.average(np.abs(d), weights=w)),
            "r": float(np.corrcoef(a[ok], b[ok])[0, 1]),
            "ge15d_kappa": float((po - pe) / (1 - pe)) if pe < 1 else None,
            "ge15d_omission": float(fn / (tp + fn)) if tp + fn else None, "ge15d_commission": float(fp / (tp + fp)) if tp + fp else None}


res["vs_s2_kordelas"] = {}
for k, y in enumerate(YEARS):
    if y in hyears:
        a = wcd[k]; b = hz["hydroperiod/19/hydroperiod_days"][hyears.index(y)]
        res["vs_s2_kordelas"][f"{y}/{str(y + 1)[2:]}"] = {s: compare(a, b, m) for s, m in STR.items()}
lh = zarr.open_group(ROOT / "data/ref/last_hydroperiod.zarr", mode="r")
if "2017/18" in lh.attrs["cycles"] and 2017 in YEARS:
    H = lh["hydroperiod_days"][list(lh.attrs["cycles"]).index("2017/18")]
    res["vs_LAST_H1718"] = {s: compare(wcd[YEARS.index(2017)], H, m) for s, m in STR.items()}
res["status"] = "ok"
jdump(res, ROOT / "results_wcd.json")
print(json.dumps({k: v for k, v in res.items() if k != "catalogue"}, indent=1)[:3000])
