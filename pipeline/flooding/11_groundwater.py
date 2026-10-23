"""11_groundwater.py <root> — dune-lagoon hydroperiod vs groundwater level (CHG piezometers) and rainfall.

Lagoons (OSM polygons -> HEALPix cells, 06_strata.py) paired with the CHG piezometer beside them:
  Laguna de Santa Olalla <-> 05.51.108 "Santa Olalla - S1"; Laguna Dulce <-> 05.51.109 "Laguna Dulce - S1";
  Laguna del Zahillo <-> 05.51.110 "Zahillo".
Annual lagoon hydroperiod (mean over lagoon cells): LAST-EBD annual maps (Zenodo 3516929) to 2016/17,
S2 (this pipeline) from 2017/18; 2017/18 is available from both and is used to report the sensor offset.
Groundwater: hydrological-year (Sep–Aug) mean of "Cota piezométrica" (m a.s.l.).
Rain: ERA5 monthly total precipitation, hydrological-year sum (mm), mean of the 0.25° cells over the area.
Model: OLS  HP ~ GW + P + S2_dummy  (standardised GW, P), years with all three; 95 % CI by bootstrap over
years (5000 draws). Also Spearman rho of HP vs GW, and of first differences (removes shared trends).
Per S2 date (2017–2026): lagoon flooded share vs piezometric level (nearest reading within 20 days).
Outputs: results_groundwater.json, data/products/lagoon_gw_annual.csv, data/products/lagoon_gw_dates.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import zarr
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from common import jdump

ROOT = Path(sys.argv[1])
PAIRS = {"Laguna de Santa Olalla": "05.51.108", "Laguna Dulce": "05.51.109", "Laguna del Zahillo": "05.51.110"}
c = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")
A = c["area_utm_m2"][:].astype(float)
stg = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")
lag_id = stg["lagoon_id"][:]
ids = {}
for lg in stg.attrs["lagoons"]:
    if lg["name"] in PAIRS and lg["n_cells"] > 0:
        ids.setdefault(lg["name"], []).append(lg["id"])
sel = {k: np.isin(lag_id, v) for k, v in ids.items()}

# groundwater
pz = pd.read_csv(ROOT / "data/raw/chg/piezometria_within10km_santa_olalla.csv", dtype={"cod_estsub": str})
pz = pz[pz.parameter.str.startswith("Cota")].copy()
pz["value"] = pd.to_numeric(pz.value.astype(str).str.replace(",", "."), errors="coerce")
pz["date"] = pd.to_datetime(pz.date)
pz["hy"] = pz.date.dt.year - (pz.date.dt.month < 9)
# rain
ds = xr.open_dataset(ROOT / "data/raw/era5/tp_monthly.nc")
tv = "valid_time" if "valid_time" in ds.dims else "time"
tp = ds["tp"].mean([d for d in ds["tp"].dims if d != tv]).to_series()
tp.index = pd.to_datetime(tp.index)
mm = tp * 1000 * tp.index.days_in_month
rain = mm.groupby(tp.index.year - (tp.index.month < 9)).sum()
rain_n = mm.groupby(tp.index.year - (tp.index.month < 9)).size()
rain = rain[rain_n == 12]

# lagoon hydroperiod: LAST annual maps + S2
lh = zarr.open_group(ROOT / "data/ref/last_hydroperiod.zarr", mode="r")
hp_rows = []
for k, cyc in enumerate(lh.attrs["cycles"]):
    y = int(cyc[:4]); H = lh["hydroperiod_days"][k]
    for name, s in sel.items():
        ok = s & (H >= 0)
        if ok.sum():
            hp_rows.append({"lagoon": name, "hy": y, "source": "LAST", "hp": float(np.average(H[ok], weights=A[ok]))})
hz = zarr.open_group(ROOT / "data/products/hydroperiod.zarr", mode="r")
years = [int(t[:4]) for t in hz.attrs["hydrological_years"]]
for k, y in enumerate(years):
    H = hz["hydroperiod/19/hydroperiod_days"][k]
    for name, s in sel.items():
        ok = s & (H >= 0)
        if ok.sum():
            hp_rows.append({"lagoon": name, "hy": y, "source": "S2", "hp": float(np.average(H[ok], weights=A[ok]))})
hp = pd.DataFrame(hp_rows)
res = {"pairs": PAIRS, "lagoon_cells": {k: int(v.sum()) for k, v in sel.items()}, "lagoons": {}}
annual = []
rng = np.random.default_rng(3)
for name, code in PAIRS.items():
    if name not in sel:
        continue
    g = pz[pz.cod_estsub == code]
    gw = g.groupby("hy").value.mean(); gwn = g.groupby("hy").value.size()
    gw = gw[gwn >= 6]
    h = hp[hp.lagoon == name]
    last = h[h.source == "LAST"].set_index("hy").hp; s2 = h[h.source == "S2"].set_index("hy").hp
    comb = pd.concat([last[last.index <= 2016], s2]).rename("hp")
    df = pd.DataFrame({"hp": comb, "gw": gw, "rain": rain})
    df["s2"] = (df.index >= 2017).astype(int)
    df = df.dropna()
    df.insert(0, "lagoon", name)
    annual.append(df.reset_index(names="hy"))
    r = {"n_years": int(len(df)), "years": [int(df.index.min()), int(df.index.max())],
         "overlap_2017_18": {"LAST": float(last.get(2017, np.nan)), "S2": float(s2.get(2017, np.nan))},
         "gw_station": code, "gw_mean_first5y": float(df.gw.iloc[:5].mean()), "gw_mean_last5y": float(df.gw.iloc[-5:].mean())}
    if len(df) >= 8:
        rho, p = stats.spearmanr(df.hp, df.gw); r["spearman_hp_gw"] = [float(rho), float(p)]
        rho2, p2 = stats.spearmanr(df.hp.diff().dropna(), df.gw.diff().dropna()); r["spearman_diff_hp_gw"] = [float(rho2), float(p2)]
        rho3, p3 = stats.spearmanr(df.hp, df.rain); r["spearman_hp_rain"] = [float(rho3), float(p3)]
        X = np.column_stack([np.ones(len(df)), (df.gw - df.gw.mean()) / df.gw.std(), (df.rain - df.rain.mean()) / df.rain.std(), df.s2])
        beta = np.linalg.lstsq(X, df.hp.values, rcond=None)[0]
        bs = []
        for _ in range(5000):
            ix = rng.integers(0, len(df), len(df))
            bs.append(np.linalg.lstsq(X[ix], df.hp.values[ix], rcond=None)[0])
        bs = np.array(bs)
        r["ols_days_per_sd"] = {nm: [round(float(beta[i]), 1), round(float(np.quantile(bs[:, i], 0.025)), 1),
                                     round(float(np.quantile(bs[:, i], 0.975)), 1)]
                                for i, nm in enumerate(["intercept", "groundwater", "rain", "s2_dummy"])}
        r["gw_sd_m"] = float(df.gw.std()); r["rain_sd_mm"] = float(df.rain.std())
        pred = X @ beta
        r["ols_r2"] = float(1 - ((df.hp - pred) ** 2).sum() / ((df.hp - df.hp.mean()) ** 2).sum())
    res["lagoons"][name] = r
pd.concat(annual).to_csv(ROOT / "data/products/lagoon_gw_annual.csv", index=False)

# per S2 date
s2c = zarr.open_group(ROOT / "data/s2/codes.zarr", mode="r")
sc = pd.read_csv(ROOT / "data/s2/scenes.csv").set_index("date")
rows = []
for i, d in enumerate(s2c.attrs["dates"]):
    if d not in sc.index or sc.loc[d, "status"] != "ok":
        continue
    code = s2c["codes"][i]
    for name, s in sel.items():
        v = s & (code != 255)
        if A[v].sum() < 0.8 * A[s].sum():
            continue
        g = pz[pz.cod_estsub == PAIRS[name]]
        dt = (g.date - pd.Timestamp(d)).abs()
        if len(g) == 0 or dt.min() > pd.Timedelta(days=20):
            gwv = np.nan
        else:
            gwv = float(g.value.values[int(np.argmin(dt.values))])
        rows.append({"date": d, "lagoon": name, "flooded_share": float(A[v & ((code & 3) == 1)].sum() / A[v].sum()), "gw": gwv})
dd = pd.DataFrame(rows); dd.to_csv(ROOT / "data/products/lagoon_gw_dates.csv", index=False)
for name in sel:
    x = dd[(dd.lagoon == name)].dropna()
    if len(x) > 10:
        rho, p = stats.spearmanr(x.flooded_share, x.gw)
        res["lagoons"][name]["per_date_spearman_flooded_vs_gw"] = [float(rho), float(p), int(len(x))]
        wet = x[x.flooded_share >= 0.1].gw; dry = x[x.flooded_share < 0.1].gw
        res["lagoons"][name]["gw_when_wet_vs_dry_median"] = [float(wet.median()) if len(wet) else None,
                                                             float(dry.median()) if len(dry) else None]
res["licence_note"] = ("CHG piezometer data: no explicit licence; CHG aviso legal authorises reproduction citing the source, "
                       "data 'no tienen carácter oficial'. Results kept local until CHG confirms terms.")
jdump(res, ROOT / "results_groundwater.json")
print(res)
