"""15_checks.py <root> — robustness checks.

C1 Deviation test: segmentation RGB from 2x COG overviews (20 m) vs native 10 m (120 dates run both ways):
   difference in Tfinal and in mapped open water.
C2 Year-to-year stability on cells whose true state should not change: (a) land never flooded in any of the 492
   LAST masks 1984–2019 and outside marsh/lagoons ("stable dry") -> share of area given >= 30 flooded days per
   year (false flooding); (b) cells flooded in >= 98 % of LAST masks outside the sea/estuary ("stable wet") ->
   share given < 300 days (false drying).
C3 Minimum mappable lagoon size: for each OSM lagoon polygon and each S2–LAST date pair (<= 1 day) on which LAST
   maps >= 50 % of the lagoon cells as water, is the lagoon detected by S2 (>= 1 open-water cell; >= 50 % of
   cells)? Detection rate by polygon-area class.
C4 External check: EBD-CSIC press release (2 Jul 2025): marsh 7.9 % flooded on 29 Dec 2024, 97.4 % on
   19 Mar 2025 -> our seasonal-marsh flooded share on the nearest usable S2 dates (different marsh outline).
Output: results_checks.json
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

sys.path.insert(0, str(Path(__file__).parent))
from common import jdump

ROOT = Path(sys.argv[1])
res = {}
A = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")["area_utm_m2"][:].astype(float)
stg = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")
st = stg["stratum"][:]; fr = stg["last_freq_all"][:]; lag = stg["lagoon_id"][:]
# C1
new = pd.read_csv(ROOT / "data/s2/scenes.csv").set_index("date")
oldp = ROOT / "data/s2_scenes_rgb10m_sensitivity.csv"
if oldp.exists():
    old = pd.read_csv(oldp).set_index("date")
    j = old[old.status == "ok"].join(new[new.status == "ok"], lsuffix="_10m", rsuffix="_ov", how="inner")
    d = j.Tfinal_ov - j.Tfinal_10m
    dw = (j.open_water_ha_ov - j.open_water_ha_10m) / j.open_water_ha_10m.clip(lower=1)
    res["C1_rgb_overview_vs_10m"] = {"n_dates": int(len(j)), "Tfinal_diff_median": float(d.median()),
                                     "Tfinal_absdiff_p90": float(d.abs().quantile(0.9)),
                                     "open_water_rel_diff_median_pct": float(100 * dw.median()),
                                     "open_water_rel_absdiff_p90_pct": float(100 * dw.abs().quantile(0.9)),
                                     "note": "Tfinal in stretched units (0-255)"}
# C2
hz = zarr.open_group(ROOT / "data/products/hydroperiod.zarr", mode="r")
years = hz.attrs["hydrological_years"]
dry = (fr == 0) & np.isin(st, (4,))
wet = (fr >= 0.98) & (st != 0)
c2 = {}
for k, y in enumerate(years):
    hp = hz["hydroperiod/19/hydroperiod_days"][k]
    okd = dry & (hp >= 0); okw = wet & (hp >= 0)
    c2[y] = {"stable_dry_false_flood_ge30d_pct": float(100 * A[okd & (hp >= 30)].sum() / A[okd].sum()),
             "stable_dry_any_flood_pct": float(100 * A[okd & (hp >= 1)].sum() / A[okd].sum()),
             "stable_wet_false_dry_lt300d_pct": float(100 * A[okw & (hp < 300)].sum() / max(A[okw].sum(), 1))}
res["C2_stability"] = c2
res["C2_areas_ha"] = {"stable_dry": float(A[dry].sum() / 1e4), "stable_wet": float(A[wet].sum() / 1e4)}
# C3
s2 = zarr.open_group(ROOT / "data/s2/codes.zarr", mode="r")
lm = zarr.open_group(ROOT / "data/ref/last_masks.zarr", mode="r")
pairs = pd.read_csv(ROOT / "data/validation_pairs.csv")
s2d = list(s2.attrs["dates"]); ld = list(lm.attrs["dates"])
info = {l["id"]: l for l in stg.attrs["lagoons"]}
ids = [i for i in np.unique(lag) if i > 0]
order = np.argsort(lag); sl = lag[order]
cellsets = {i: order[np.searchsorted(sl, i):np.searchsorted(sl, i, side="right")] for i in ids}
rows = []
for _, p in pairs.iterrows():
    c = s2["codes"][s2d.index(p.s2_date)]; r = lm["masks"][ld.index(p.last_date.replace("-", ""))]
    for i, cs in cellsets.items():
        v = (c[cs] != 255) & (r[cs] != 255)
        if v.mean() < 0.8:
            continue
        lw = (r[cs][v] == 1).mean()
        if lw < 0.5:
            continue
        sw = ((c[cs][v] & 3) == 1)
        rows.append({"id": i, "area_ha": info[i]["polygon_area_ha"], "date": p.s2_date, "any": bool(sw.any()), "half": bool(sw.mean() >= 0.5)})
dd = pd.DataFrame(rows)
if len(dd):
    dd["size"] = pd.cut(dd.area_ha, [0, 0.5, 1, 2, 5, 10, 1000], labels=["<0.5", "0.5-1", "1-2", "2-5", "5-10", ">10"])
    g = dd.groupby("size", observed=True).agg(n=("any", "size"), lagoons=("id", "nunique"), detected_any=("any", "mean"), detected_half=("half", "mean"))
    res["C3_detection_by_size"] = g.round(3).reset_index().to_dict(orient="records")
# C4
m = pd.read_csv(ROOT / "data/products/marsh_area_by_date.csv", parse_dates=["date"])
m = m[m.source == "S2"]
c4 = {}
for target, ebd in (("2024-12-29", 7.9), ("2025-03-19", 97.4)):
    t = pd.Timestamp(target); k = (m.date - t).abs().idxmin()
    c4[target] = {"ebd_pct": ebd, "s2_date": str(m.date[k].date()), "s2_open_water_pct": round(100 * m.flooded_share[k], 1),
                  "s2_valid_marsh_pct": round(100 * m.valid_marsh_share[k], 1)}
res["C4_ebd_press_release"] = c4
jdump(res, ROOT / "results_checks.json")
import json; print(json.dumps(res, indent=1)[:4000])
