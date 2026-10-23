"""10_hydroperiod.py <root> — hydroperiod per HEALPix d19 cell per hydrological year (1 Sep – 31 Aug).

Interpolation rule (stated, simple, reproducible): every day of the year takes the state (flooded / not)
of the valid observation nearest in time for that cell (Voronoi in time; ties go to the later observation;
days before the first / after the last valid observation take that observation's state).
Hydroperiod = number of flooded days. Reported with it, per cell: n_valid (valid S2 observations) and
max_gap_days (longest interval between consecutive valid observations, including year edges).
Flooded = Kordelas open water (class 1). Variants: open + water-vegetation (class 1|2); S2 + S1 open water
(S1 dates fill in where the cell has no valid S2 observation that day).
Uncertainty: half-sample test — recompute with a random 50 % of dates (5 seeded draws); per-cell mean
absolute difference and the change in the marsh-mean hydroperiod.
Same rule applied to the LAST Landsat masks for 2017/18 and 2018/19, and checked against the LAST annual
hydroperiod map H1718 (Zenodo 3516929) to verify that our rule reproduces theirs.
Outputs: data/products/hydroperiod.zarr (GRID4EARTH DGGS-Zarr, levels 19/17/15), results_hydroperiod.json,
data/products/hydroperiod_summary.csv, data/products/lagoons_hydroperiod.csv
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

sys.path.insert(0, str(Path(__file__).parent))
from common import dggs_attrs, jdump

ROOT = Path(sys.argv[1])
c = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")
cells = {k: c[k][:] for k in c.array_keys()}
n = len(cells["cell_ids"]); A = cells["area_utm_m2"].astype(np.float64) / 1e4
stg = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")
st = stg["stratum"][:]; lag_id = stg["lagoon_id"][:]; lagoons = stg.attrs["lagoons"]
marsh = st == 1
s2 = zarr.open_group(ROOT / "data/s2/codes.zarr", mode="r")
dates = pd.to_datetime(pd.Series(s2.attrs["dates"]))
scenes = pd.read_csv(ROOT / "data/s2/scenes.csv").set_index("date")
okdates = set(scenes.index[scenes.status == "ok"])
YEARS = list(range(2017, 2026))
s1p = ROOT / "data/s1/codes.zarr"
s1 = zarr.open_group(s1p, mode="r") if s1p.exists() else None


def hydroperiod(obs, L):
    """obs: list of (day_index, valid bool[n], flooded bool[n]) sorted by day. Returns hp, n_valid, max_gap."""
    acc = np.zeros(n, np.float32); start = np.zeros(n, np.float32); prev_t = np.full(n, -1.0, np.float32)
    prev_w = np.zeros(n, bool); nv = np.zeros(n, np.int16); gap = np.zeros(n, np.float32)
    for t, v, w in obs:
        first = v & (prev_t < 0)
        later = v & (prev_t >= 0)
        m = (prev_t + t) / 2.0
        acc[later] += np.where(prev_w[later], m[later] - start[later], 0)
        gap[later] = np.maximum(gap[later], t - prev_t[later])
        start[later] = m[later]
        start[first] = 0.0; gap[first] = t
        prev_t[v] = t; prev_w[v] = w[v]; nv[v] += 1
    seen = prev_t >= 0
    acc[seen] += np.where(prev_w[seen], L - start[seen], 0)
    gap[seen] = np.maximum(gap[seen], L - prev_t[seen])
    hp = np.where(seen, np.round(acc), -1).astype(np.int16)
    return hp, nv, np.where(seen, gap, -1).astype(np.int16)


def year_obs(y, variant, keep=None):
    y0 = pd.Timestamp(f"{y}-09-01"); L = (pd.Timestamp(f"{y + 1}-09-01") - y0).days
    obs = []
    for i, d in enumerate(dates):
        if not (y0 <= d < y0 + pd.Timedelta(days=L)) or s2.attrs["dates"][i] not in okdates:
            continue
        if keep is not None and not keep(i):
            continue
        code = s2["codes"][i]
        v = code != 255
        w = ((code & 3) == 1) if variant == "open" else ((code & 3) >= 1)
        obs.append(((d - y0).days, v, w))
    if variant == "s2s1" and s1 is not None:
        s1d = pd.to_datetime(pd.Series(s1.attrs["dates"]))
        have = {o[0]: o for o in obs}
        for i, d in enumerate(s1d):
            if not (y0 <= d < y0 + pd.Timedelta(days=L)):
                continue
            code = s1["codes"][i]
            v = code != 255
            t = (d - y0).days
            if t in have:
                v = v & ~have[t][1]
            obs.append((t, v, code == 1))
        obs.sort(key=lambda o: o[0])
    return obs, L


OUT = ROOT / "data/products"; OUT.mkdir(parents=True, exist_ok=True)
g = zarr.open_group(OUT / "hydroperiod.zarr", mode="w")
g.attrs.update({"title": "Doñana hydroperiod on HEALPix (prototype, LifeWatch ERIC Spain brief)",
                "hydrological_years": [f"{y}/{str(y + 1)[2:]}" for y in YEARS],
                "rule": "nearest valid observation in time; flooded = Kordelas 2018 open water (S2)",
                "source": "Sentinel-2 L2A (Copernicus), Earth Search + Planetary Computer; pipeline scripts 01-10"})
lev = g.create_group("hydroperiod")
L19 = lev.create_group("19", attributes=dggs_attrs(19))
L19.create_array("cell_ids", data=cells["cell_ids"], chunks=(1 << 20,), dimension_names=["cells"])
specs = {"hydroperiod_days": "int16", "n_valid": "int16", "max_gap_days": "int16", "hydroperiod_days_open_plus_veg": "int16",
         "halfsample_mad_days": "float32"}
arrs = {k: L19.create_array(k, shape=(len(YEARS), n), chunks=(1, 1 << 20), dtype=t, fill_value=-1,
                            dimension_names=["year", "cells"]) for k, t in specs.items()}
if s1 is not None:
    arrs["hydroperiod_days_s2s1"] = L19.create_array("hydroperiod_days_s2s1", shape=(len(YEARS), n), chunks=(1, 1 << 20),
                                                     dtype="int16", fill_value=-1, dimension_names=["year", "cells"])
rows, lrows = [], []
rng = np.random.default_rng(7)
for k, y in enumerate(YEARS):
    obs, L = year_obs(y, "open")
    if not obs:
        continue
    hp, nv, gp = hydroperiod(obs, L)
    arrs["hydroperiod_days"][k] = hp; arrs["n_valid"][k] = nv; arrs["max_gap_days"][k] = gp
    hpv, _, _ = hydroperiod(*year_obs(y, "veg")[:1], L)
    arrs["hydroperiod_days_open_plus_veg"][k] = hpv
    diffs, mean_shift = [], []
    for r in range(5):
        keep_idx = set(rng.choice(len(dates), len(dates) // 2, replace=False).tolist())
        o2, _ = year_obs(y, "open", keep=lambda i: i in keep_idx)
        if len(o2) < 3:
            continue
        h2, _, _ = hydroperiod(o2, L)
        okc = (hp >= 0) & (h2 >= 0)
        diffs.append(np.where(okc, np.abs(h2.astype(np.float32) - hp), np.nan))
        mean_shift.append(float(np.average(h2[okc & marsh], weights=A[okc & marsh]) - np.average(hp[okc & marsh], weights=A[okc & marsh])))
    mad = np.nanmean(np.stack(diffs), 0) if diffs else np.full(n, np.nan)
    arrs["halfsample_mad_days"][k] = mad.astype(np.float32)
    row = {"hydro_year": f"{y}/{str(y + 1)[2:]}", "n_s2_dates_used": len(obs)}
    if s1 is not None:
        o3, _ = year_obs(y, "s2s1"); h3, nv3, _ = hydroperiod(o3, L)
        arrs["hydroperiod_days_s2s1"][k] = h3
        okm3 = marsh & (h3 >= 0)
        row["marsh_mean_hp_s2s1"] = float(np.average(h3[okm3], weights=A[okm3]))
    okm = marsh & (hp >= 0)
    row.update({"marsh_mean_hp_days": float(np.average(hp[okm], weights=A[okm])),
                "marsh_mean_hp_open_plus_veg": float(np.average(hpv[okm], weights=A[okm])),
                "marsh_flooded_ge1d_ha": float(A[okm & (hp >= 1)].sum()),
                "marsh_flooded_ge90d_ha": float(A[okm & (hp >= 90)].sum()),
                "marsh_median_n_valid": float(np.median(nv[okm])), "marsh_median_max_gap": float(np.median(gp[okm])),
                "marsh_p90_max_gap": float(np.quantile(gp[okm], 0.9)),
                "marsh_halfsample_mad_days_median": float(np.nanmedian(mad[okm])),
                "marsh_halfsample_shift_of_mean_days": [round(x, 2) for x in mean_shift]})
    rows.append(row); print(row, flush=True)
    for lg in lagoons:
        sel = lag_id == lg["id"]
        if sel.sum() == 0:
            continue
        okl = sel & (hp >= 0)
        lrows.append({"hydro_year": row["hydro_year"], "lagoon_id": lg["id"], "name": lg["name"],
                      "polygon_area_ha": round(lg["polygon_area_ha"], 2), "n_cells": int(sel.sum()),
                      "mean_hp_days": float(np.average(hp[okl], weights=A[okl])) if okl.any() else np.nan,
                      "max_hp_days": int(hp[okl].max()) if okl.any() else -1,
                      "flooded_share_ge1d": float(A[okl & (hp >= 1)].sum() / A[okl].sum()) if okl.any() else np.nan,
                      "median_n_valid": float(np.median(nv[okl])) if okl.any() else 0,
                      "median_max_gap_days": float(np.median(gp[okl])) if okl.any() else np.nan,
                      "halfsample_mad_days": float(np.nanmean(mad[okl])) if okl.any() else np.nan})
pd.DataFrame(rows).to_csv(OUT / "hydroperiod_summary.csv", index=False)
pd.DataFrame(lrows).to_csv(OUT / "lagoons_hydroperiod.csv", index=False)

# multiscale coarsening (mean of available children; NESTED parent = id >> 2*(19-l))
layout = [{"asset": "19", "dggs": dggs_attrs(19)["dggs"]}]
for lvl in (17, 15):
    par = (cells["cell_ids"] >> np.uint64(2 * (19 - lvl)))
    up, inv = np.unique(par, return_inverse=True)
    G = lev.create_group(str(lvl), attributes=dggs_attrs(lvl))
    G.create_array("cell_ids", data=up, dimension_names=["cells"])
    out = np.full((len(YEARS), len(up)), np.nan, np.float32)
    for k in range(len(YEARS)):
        h = arrs["hydroperiod_days"][k].astype(np.float32); ok = h >= 0
        s = np.bincount(inv[ok], weights=h[ok] * A[ok], minlength=len(up)); w = np.bincount(inv[ok], weights=A[ok], minlength=len(up))
        out[k] = np.where(w > 0, s / np.maximum(w, 1e-12), np.nan)
    G.create_array("hydroperiod_days", data=out, dimension_names=["year", "cells"])
    layout.append({"asset": str(lvl), "dggs": dggs_attrs(lvl)["dggs"]})
lev.attrs.update({"zarr_conventions": [{"uuid": "d35379db-88df-4056-af3a-620245f8e347", "name": "multiscales",
                                        "schema_url": "https://raw.githubusercontent.com/zarr-conventions/multiscales/refs/tags/v1/schema.json",
                                        "spec_url": "https://github.com/zarr-conventions/multiscales/blob/v1/README.md",
                                        "description": "Multiscale layout of zarr datasets"}] + dggs_attrs(19)["zarr_conventions"],
                  "multiscales": {"layout": layout}})

# ---- LAST masks with the same rule; check against LAST annual map H1718 ----
lm = zarr.open_group(ROOT / "data/ref/last_masks.zarr", mode="r")
ld = pd.to_datetime(pd.Series(lm.attrs["dates"]), format="%Y%m%d")
lh = zarr.open_group(ROOT / "data/ref/last_hydroperiod.zarr", mode="r")
cyc = list(lh.attrs["cycles"])
chk = {}
for y in (2017, 2018):
    y0 = pd.Timestamp(f"{y}-09-01"); L = (pd.Timestamp(f"{y + 1}-09-01") - y0).days
    obs = []
    for j, d in enumerate(ld):
        if y0 <= d < y0 + pd.Timedelta(days=L):
            m = lm["masks"][j]; obs.append(((d - y0).days, m != 255, m == 1))
    hpl, nvl, _ = hydroperiod(obs, L)
    k = YEARS.index(y); hps = arrs["hydroperiod_days"][k]
    okm = marsh & (hpl >= 0) & (hps >= 0)
    e = {"n_last_masks": len(obs),
         "marsh_mean_hp_last_rule": float(np.average(hpl[okm], weights=A[okm])),
         "marsh_mean_hp_s2": float(np.average(hps[okm], weights=A[okm])),
         "s2_minus_last_mean_days": float(np.average(hps[okm].astype(float) - hpl[okm], weights=A[okm])),
         "mae_days": float(np.average(np.abs(hps[okm].astype(float) - hpl[okm]), weights=A[okm])),
         "pearson_r": float(np.corrcoef(hps[okm], hpl[okm])[0, 1])}
    tag = f"{y}/{str(y + 1)[2:]}"
    if tag in cyc:
        H = lh["hydroperiod_days"][cyc.index(tag)]
        ok2 = marsh & (H >= 0) & (hpl >= 0)
        e["our_rule_vs_LAST_H1718_mae_days"] = float(np.average(np.abs(hpl[ok2].astype(float) - H[ok2]), weights=A[ok2]))
        e["our_rule_vs_LAST_H1718_bias_days"] = float(np.average(hpl[ok2].astype(float) - H[ok2], weights=A[ok2]))
        e["LAST_H1718_marsh_mean"] = float(np.average(H[ok2], weights=A[ok2]))
        ok3 = marsh & (H >= 0) & (hps >= 0)
        e["s2_vs_LAST_H1718_mae_days"] = float(np.average(np.abs(hps[ok3].astype(float) - H[ok3]), weights=A[ok3]))
        e["s2_vs_LAST_H1718_bias_days"] = float(np.average(hps[ok3].astype(float) - H[ok3], weights=A[ok3]))
        e["s2_vs_LAST_H1718_r"] = float(np.corrcoef(hps[ok3], H[ok3])[0, 1])
    chk[tag] = e
jdump({"years": rows, "last_overlap": chk}, ROOT / "results_hydroperiod.json")
print(json.dumps(chk, indent=1))

# ---- flooded share of the seasonal marsh per usable date (S2; S1 fill dates) ----
rows = []
for i, d in enumerate(s2.attrs["dates"]):
    if d not in okdates:
        continue
    code = s2["codes"][i]; v = marsh & (code != 255)
    if A[v].sum() >= 0.5 * A[marsh].sum():
        rows.append({"date": d, "source": "S2", "valid_marsh_share": float(A[v].sum() / A[marsh].sum()),
                     "flooded_share": float(A[v & ((code & 3) == 1)].sum() / A[v].sum()),
                     "flooded_share_incl_waterveg": float(A[v & ((code & 3) >= 1)].sum() / A[v].sum())})
if s1 is not None:
    for i, d in enumerate(s1.attrs["dates"]):
        if "F" not in s1.attrs["why"][i]:
            continue
        code = s1["codes"][i]; v = marsh & (code != 255)
        if A[v].sum() >= 0.5 * A[marsh].sum():
            rows.append({"date": d, "source": "S1", "valid_marsh_share": float(A[v].sum() / A[marsh].sum()),
                         "flooded_share": float(A[v & (code == 1)].sum() / A[v].sum()),
                         "flooded_veg_flag_share": float(A[v & (code == 2)].sum() / A[v].sum())})
pd.DataFrame(rows).sort_values("date").to_csv(OUT / "marsh_area_by_date.csv", index=False)
