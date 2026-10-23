"""10b_anomaly.py <root> — hydroperiod anomaly against the LAST-EBD Landsat baseline.

Baseline: LAST annual hydroperiod maps 1984/85 – 2016/17 (33 cycles; the Landsat-TM era, same years as the
LAST flood masks), per HEALPix d19 cell: mean and standard deviation; marsh-mean series.
S2 years 2017/18 – 2025/26 from 10_hydroperiod.py. Offset between sensors measured on the one overlap year
(2017/18: S2 marsh mean minus LAST H1718 marsh mean) and reported; anomalies are given both raw and
offset-corrected. Marsh-mean trend: Theil–Sen slope with 95 % CI (scipy), LAST 1984/85–2016/17 and
S2-era separately (9 points: indicative only).
Outputs: results_anomaly.json, data/products/anomaly.zarr (GRID4EARTH, level 19: baseline mean/sd, anomaly per S2 year),
data/products/marsh_series.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from common import dggs_attrs, jdump

ROOT = Path(sys.argv[1])
c = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")
ids = c["cell_ids"][:]; A = c["area_utm_m2"][:].astype(float); n = len(ids)
st = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")["stratum"][:]
marsh = st == 1
lh = zarr.open_group(ROOT / "data/ref/last_hydroperiod.zarr", mode="r")
cyc = list(lh.attrs["cycles"])
order = sorted(range(len(cyc)), key=lambda k: cyc[k])
series = []
s = np.zeros(n); s2 = np.zeros(n); cnt = np.zeros(n)
for k in order:
    y = int(cyc[k][:4]); H = lh["hydroperiod_days"][k].astype(float)
    ok = marsh & (H >= 0)
    series.append({"hy": y, "source": "LAST", "marsh_mean_hp": float(np.average(H[ok], weights=A[ok])),
                   "marsh_ge90_ha": float(A[ok & (H >= 90)].sum() / 1e4)})
    if 1984 <= y <= 2016:
        v = H >= 0
        s[v] += H[v]; s2[v] += H[v] ** 2; cnt[v] += 1
base_mean = np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)
base_sd = np.sqrt(np.maximum(s2 / np.maximum(cnt, 1) - base_mean ** 2, 0))
hz = zarr.open_group(ROOT / "data/products/hydroperiod.zarr", mode="r")
yrs = [int(t[:4]) for t in hz.attrs["hydrological_years"]]
HP = hz["hydroperiod/19/hydroperiod_days"]
for k, y in enumerate(yrs):
    H = HP[k].astype(float); ok = marsh & (H >= 0)
    if ok.sum():
        series.append({"hy": y, "source": "S2", "marsh_mean_hp": float(np.average(H[ok], weights=A[ok])),
                       "marsh_ge90_ha": float(A[ok & (H >= 90)].sum() / 1e4)})
ser = pd.DataFrame(series); ser.to_csv(ROOT / "data/products/marsh_series.csv", index=False)
L = ser[(ser.source == "LAST") & ser.hy.between(1984, 2016)]
bm, bsd = L.marsh_mean_hp.mean(), L.marsh_mean_hp.std()
S = ser[ser.source == "S2"].set_index("hy")
lastov = ser[(ser.source == "LAST") & (ser.hy == 2017)].marsh_mean_hp
offset = float(S.marsh_mean_hp.get(2017) - lastov.iloc[0]) if len(lastov) and 2017 in S.index else np.nan
res = {"baseline_years": "1984/85-2016/17", "baseline_marsh_mean_hp": bm, "baseline_marsh_sd": bsd,
       "baseline_min": float(L.marsh_mean_hp.min()), "baseline_max": float(L.marsh_mean_hp.max()),
       "offset_S2_minus_LAST_2017_18_days": offset, "years": {}}
for y, r in S.iterrows():
    a = r.marsh_mean_hp - bm
    res["years"][f"{y}/{str(y + 1)[2:]}"] = {"marsh_mean_hp": r.marsh_mean_hp, "anomaly_days": a,
                                              "anomaly_days_offset_corrected": a - offset if np.isfinite(offset) else None,
                                              "z": a / bsd, "z_offset_corrected": (a - offset) / bsd if np.isfinite(offset) else None,
                                              "rank_in_baseline_pct": float((L.marsh_mean_hp < r.marsh_mean_hp).mean() * 100)}
ts = stats.theilslopes(L.marsh_mean_hp, L.hy)
res["theil_sen_LAST_1984_2016_days_per_year"] = [ts.slope, ts.low_slope, ts.high_slope]
if len(S) >= 5:
    t2 = stats.theilslopes(S.marsh_mean_hp, S.index)
    res["theil_sen_S2_era_days_per_year"] = [t2.slope, t2.low_slope, t2.high_slope]
g = zarr.open_group(ROOT / "data/products/anomaly.zarr", mode="w")
G = g.create_group("19", attributes=dggs_attrs(19) | {"baseline": "LAST-EBD annual hydroperiod 1984/85-2016/17 (CC BY 4.0)"})
G.create_array("cell_ids", data=ids, chunks=(1 << 20,), dimension_names=["cells"])
G.create_array("baseline_mean_days", data=base_mean.astype("float32"), chunks=(1 << 20,), dimension_names=["cells"])
G.create_array("baseline_sd_days", data=base_sd.astype("float32"), chunks=(1 << 20,), dimension_names=["cells"])
an = G.create_array("anomaly_days", shape=(len(yrs), n), chunks=(1, 1 << 20), dtype="float32", fill_value=np.nan,
                    dimension_names=["year", "cells"])
for k in range(len(yrs)):
    H = HP[k].astype(float)
    an[k] = np.where(H >= 0, H - base_mean, np.nan).astype("float32")
G.attrs["hydrological_years"] = hz.attrs["hydrological_years"]
jdump(res, ROOT / "results_anomaly.json")
print(res)
