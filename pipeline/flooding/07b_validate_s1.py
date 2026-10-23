"""07b_validate_s1.py <root> — Sentinel-1 open water (and flooded-vegetation flag) vs LAST masks (<= 1 day,
Sep 2017 – Aug 2019) and vs S2 open water (<= 1 day, all years), by stratum, and S1 omission vs ERA5 10 m wind
speed at the overpass hour (wind roughening). Same metrics/block bootstrap as 07_validate.py.
Output: results_validation_s1.json, data/validation_pairs_s1.csv
"""
import glob, importlib.util, sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import zarr

sys.path.insert(0, str(Path(__file__).parent))
from common import jdump

ROOT = Path(sys.argv[1])
sys.argv = [sys.argv[0], str(ROOT)]
src = open(Path(__file__).parent / "07_validate.py").read().split("# ---------- R1")[0]
ns = {"__file__": str(Path(__file__).parent / "07_validate.py")}
exec(compile(src, "07_validate_head", "exec"), ns)
A, st, NB, STRATA = ns["A"], ns["st"], ns["NB"], ns["STRATA"]
conf_blocks, summarise, metrics, boundary = ns["conf_blocks"], ns["summarise"], ns["metrics"], ns["boundary"]
s2, scenes, lm, ldates = ns["s2"], ns["scenes"], ns["lm"], ns["ldates"]
s2dates = pd.to_datetime(pd.Series(s2.attrs["dates"]))
s1 = zarr.open_group(ROOT / "data/s1/codes.zarr", mode="r")
s1d = pd.to_datetime(pd.Series(s1.attrs["dates"]))
s1sc = pd.read_csv(ROOT / "data/s1/scenes.csv")
okk = set(zip(s1sc.date[s1sc.status == "ok"], s1sc.rel_orbit[s1sc.status == "ok"]))
orb = s1.attrs["rel_orbit"]

# wind at overpass: ascending orbits (74, 147) ~18 UTC; descending (154) ~06 UTC
wind = None
fs = sorted(glob.glob(str(ROOT / "data/raw/era5/wind_*.nc")))
if fs:
    parts = []
    for f in fs:
        w = xr.open_dataset(f)
        tv = "valid_time" if "valid_time" in w.dims else "time"
        parts.append(np.sqrt(w.u10 ** 2 + w.v10 ** 2).mean([d for d in w.u10.dims if d != tv]).to_series())
    wind = pd.concat(parts).sort_index()


def wind_at(d, o):
    if wind is None:
        return np.nan
    h = 18 if o in (74, 147) else 6
    t = pd.Timestamp(d) + pd.Timedelta(hours=h)
    return float(wind.get(t, np.nan))


rows = []
accL = {s: np.zeros((4, NB)) for s in STRATA}; accS = {s: np.zeros((4, NB)) for s in STRATA}
fvL = []
for i, d in enumerate(s1d):
    if (s1.attrs["dates"][i], orb[i]) not in okk:
        continue
    code = s1["codes"][i]
    ok1 = code != 255
    p = code == 1
    row = {"s1_date": s1.attrs["dates"][i], "rel_orbit": orb[i], "wind_ms": wind_at(d, orb[i])}
    # vs LAST
    dl = (ldates - d).abs(); j = int(dl.idxmin())
    if dl[j] <= pd.Timedelta(days=1) and pd.Timestamp("2017-09-01") <= d <= pd.Timestamp("2019-08-31"):
        ref = lm["masks"][j]; ok = ok1 & (ref != 255); rb = ref == 1; bnd = boundary(ref)
        for s, sm in STRATA.items():
            cb = conf_blocks(p, rb, ok & sm & ~bnd); accL[s] += cb
            mt = metrics(cb.sum(1))
            row[f"last_{s}_kappa"] = mt.get("kappa"); row[f"last_{s}_PA"] = mt.get("PA_water"); row[f"last_{s}_UA"] = mt.get("UA_water")
        fv = ok & (code == 2)
        row["fv_flag_ha"] = float(A[fv].sum()); row["fv_flag_last_water_share"] = float(A[fv & rb].sum() / max(A[fv].sum(), 1e-9))
        fvL.append((float(A[fv].sum()), float(A[fv & rb].sum())))
    # vs S2 open water
    ds = (s2dates - d).abs(); k = int(ds.idxmin())
    dk = s2.attrs["dates"][k]
    if ds[k] <= pd.Timedelta(days=1) and dk in scenes.index and scenes.loc[dk, "status"] == "ok":
        c2 = s2["codes"][k]; ok = ok1 & (c2 != 255); r2 = (c2 & 3) == 1
        for s, sm in STRATA.items():
            cb = conf_blocks(p, r2, ok & sm); accS[s] += cb
            mt = metrics(cb.sum(1))
            row[f"s2_{s}_kappa"] = mt.get("kappa"); row[f"s2_{s}_PA"] = mt.get("PA_water")
    if len(row) > 3:
        rows.append(row)
df = pd.DataFrame(rows); df.to_csv(ROOT / "data/validation_pairs_s1.csv", index=False)
res = {"n_pairs_last": int(df.filter(like="last_complete_kappa").notna().sum().sum()) if len(df) else 0,
       "n_pairs_s2": int(df.filter(like="s2_complete_kappa").notna().sum().sum()) if len(df) else 0,
       "vs_LAST_EAP": {s: summarise(accL[s]) for s in STRATA},
       "vs_S2_open_water": {s: summarise(accS[s]) for s in STRATA},
       "flooded_veg_flag_vs_LAST": {"flag_ha_total": sum(a for a, b in fvL), "share_LAST_water": (sum(b for a, b in fvL) / max(sum(a for a, b in fvL), 1e-9))}}
if len(df) and df.wind_ms.notna().any():
    for ref in ("last", "s2"):
        col = f"{ref}_marsh_PA"
        x = df[[col, "wind_ms"]].dropna()
        if len(x) > 4:
            res[f"marsh_PA_by_wind_{ref}"] = {f"<{t}": float(x[x.wind_ms < t][col].mean()) for t in (4, 6)} | \
                                             {f">={t}": float(x[x.wind_ms >= t][col].mean()) for t in (4, 6)} | \
                                             {"spearman": float(x[col].corr(x.wind_ms, method="spearman")), "n": int(len(x))}
jdump(res, ROOT / "results_validation_s1.json")
print(res)
