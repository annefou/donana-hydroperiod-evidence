"""Is the Climate DT (gen 2, IFS-NEMO and IFS-FESOM) good enough to drive Donana pond-hydroperiod scenarios?

Inputs (from scripts 01-03): ROOT/data/raw/dt/<param>/*.covjson, ROOT/data/raw/era5/*.zip, ROOT/data/raw/eobs/*.nc
Outputs: ROOT/results.json, ROOT/data/derived/site_daily.nc (site-mean daily series: obs + DT statistics inputs;
DT-derived, not raw DT), ROOT/data/derived/hy_table.csv (per-hydrological-year metrics).

Site series: precipitation = mean of the 5 DT cells (marsh centre/north/east, aquifer at El Rocio, dunes of the
Biological Reserve); temperature = mean of 3 cells. Observations: E-OBS v31.0e 0.1 deg (nearest land cell to each
DT cell centre, then averaged; primary reference for precipitation) and ERA5 0.25 deg (nearest grid point; primary
for temperature, Tmax and Hargreaves PET because E-OBS tg is daily mean only here).
Hydrological year (HY) y = 1 Sep y .. 31 Aug y+1; HY 1990..2023 (34 years) used, i.e. Sep 1990 - Aug 2024.
The DT is a free-running coupled simulation: its years are not in phase with real years. Only statistics are
compared (distributions, quantiles, CV, autocorrelation, seasonal cycle, trends).

Bias correction (documented; fitted on HY 1990..2013 = Sep 1990-Aug 2014, all DT "hist"; evaluated out of sample on
HY 2014..2023 = Sep 2014-Aug 2024, i.e. 4 months of "hist" then SSP3-7.0):
  BC1 monthly scaling: P_bc = P_dt * mean_obs(month)/mean_dt(month); T_bc = T_dt + [mean_obs(month)-mean_dt(month)],
      applied separately to daily mean, max and min temperature (each against ERA5).
  BC2 empirical quantile mapping of daily precipitation per calendar month (100 quantiles, wet-day frequency
      adjusted implicitly; values beyond the fitted range scaled by the top-quantile ratio); temperature as BC1.
Reference observation for fitting: E-OBS (precipitation), ERA5 (temperature).
Usage: python -I scripts/04_analyse.py <ROOT>
"""
import sys, os, json, glob, zipfile, io
import numpy as np, pandas as pd, xarray as xr
from scipy import stats

ROOT = sys.argv[1]
RAW = os.path.join(ROOT, "data", "raw")
DER = os.path.join(ROOT, "data", "derived")
os.makedirs(DER, exist_ok=True)
CELLS = {3490355: (36.9352, -6.3281), 3490408: (37.1221, -6.2402), 3490451: (37.1221, -6.5039),
         3490362: (36.9819, -6.4600), 3490357: (36.9819, -6.1963)}
MODELS = ("IFS-NEMO", "IFS-FESOM")
HY0, HY1 = 1990, 2023
CAL, VAL = (1990, 2013), (2014, 2023)       # HY ranges: calibration (in hist) / validation (mostly SSP3-7.0)


# ---------------- loading ----------------
def load_dt(model, param):
    """Hourly site-mean series over all cells available for this param (returns series, cell list)."""
    per = {}
    for f in sorted(glob.glob(os.path.join(RAW, "dt", param, f"{model}_*.covjson"))):
        cell = int(os.path.basename(f).split("_")[1])
        j = json.load(open(f))
        for cov in j["coverages"]:
            t = pd.to_datetime([s.replace("Z", "") for s in cov["domain"]["axes"]["t"]["values"]])
            key = list(cov["ranges"])[0]
            v = np.array([np.nan if x is None else x for x in cov["ranges"][key]["values"]], float)
            per.setdefault(cell, []).append(pd.Series(v, index=t))
    ser = {}
    for c, parts in per.items():
        s = pd.concat(parts).sort_index()
        ser[c] = s[~s.index.duplicated()]
    if not ser:
        return None, []
    df = pd.DataFrame(ser)
    return df.mean(axis=1, skipna=False), sorted(ser)


def dt_daily(model):
    out = {}
    p, pc = load_dt(model, "avg_tprate")
    t, tc = load_dt(model, "2t")
    e, ec = load_dt(model, "avg_ie")
    if p is not None:
        out["pr"] = p.resample("D").mean() * 86400.0                       # kg m-2 s-1 -> mm/day
        out["pr_n"] = p.resample("D").count()
    if t is not None:
        t = t - 273.15
        out["tas"] = t.resample("D").mean()
        out["tmax"] = t.resample("D").max()
        out["tmin"] = t.resample("D").min()
        out["t_n"] = t.resample("D").count()
    if e is not None:
        out["evap"] = -e.resample("D").mean() * 86400.0                     # mm/day, positive = evaporation
    df = pd.DataFrame(out)
    if "pr_n" in df:
        df.loc[df.pr_n < 24, "pr"] = np.nan
    if "t_n" in df:
        df.loc[df.t_n < 24, ["tas", "tmax", "tmin"]] = np.nan
    return df.drop(columns=[c for c in ("pr_n", "t_n") if c in df]), {"pr": pc, "tas": tc, "evap": ec}


def era5_daily():
    series, pts = {}, {}
    for f in sorted(glob.glob(os.path.join(RAW, "era5", "era5_*.zip"))):
        z = zipfile.ZipFile(f)
        name = [n for n in z.namelist() if n.endswith(".nc")][0]
        ds = xr.open_dataset(io.BytesIO(z.read(name)), engine="h5netcdf")
        key = (round(float(ds.latitude), 3), round(float(ds.longitude), 3))
        if key in pts:
            continue
        pts[key] = os.path.basename(f)
        df = pd.DataFrame({"tp": ds.tp.values * 1000.0, "t2m": ds.t2m.values - 273.15},
                          index=pd.to_datetime(ds.valid_time.values))
        # tp at valid_time = accumulation over the preceding hour -> shift by 1 h before daily sums
        tp = df.tp.copy(); tp.index = tp.index - pd.Timedelta("1h")
        series[key] = pd.DataFrame({"pr": tp.resample("D").sum(min_count=24), "tas": df.t2m.resample("D").mean(),
                                    "tmax": df.t2m.resample("D").max(), "tmin": df.t2m.resample("D").min()})
    allp = pd.concat(series.values(), keys=range(len(series)))
    return allp.groupby(level=1).mean(), [list(k) for k in pts]


def eobs_daily():
    out, used = {}, {}
    for var in ("rr", "tg"):
        ds = xr.open_dataset(os.path.join(RAW, "eobs", f"eobs_{var}_donana.nc"))
        land = ds[var].notnull().mean("time") > 0.9
        lat2, lon2 = xr.broadcast(ds.latitude, ds.longitude)
        cols = []
        for c, (la, lo) in CELLS.items():
            d = ((lat2 - la) ** 2 + ((lon2 - lo) * np.cos(np.deg2rad(la))) ** 2).where(land)
            iy, ix = np.unravel_index(int(np.nanargmin(d.values)), d.shape)
            used.setdefault(var, []).append([round(float(ds.latitude[iy]), 3), round(float(ds.longitude[ix]), 3)])
            cols.append(ds[var][:, iy, ix].to_series())
        out[var] = pd.concat(cols, axis=1).mean(axis=1, skipna=False)
    return pd.DataFrame({"pr": out["rr"], "tas": out["tg"]}), used


# ---------------- metrics ----------------
def hy_index(idx):
    return np.where(idx.month >= 9, idx.year, idx.year - 1)


def ra_mj(doy, lat_deg):
    """Extraterrestrial radiation (FAO-56 eq. 21), MJ m-2 day-1."""
    phi = np.deg2rad(lat_deg)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)
    dl = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)
    ws = np.arccos(-np.tan(phi) * np.tan(dl))
    return 24 * 60 / np.pi * 0.0820 * dr * (ws * np.sin(phi) * np.sin(dl) + np.cos(phi) * np.cos(dl) * np.sin(ws))


def hargreaves(d):
    """Hargreaves-Samani reference ET (mm/day) from daily Tmean/Tmax/Tmin of hourly data (FAO-56 eq. 52)."""
    ra = ra_mj(d.index.day_of_year.values, 37.03)
    return pd.Series(0.0023 * 0.408 * ra * (d.tas + 17.8) * np.sqrt((d.tmax - d.tmin).clip(lower=0)), index=d.index)


def hy_table(d):
    d = d[(d.index >= f"{HY0}-09-01") & (d.index <= f"{HY1 + 1}-08-31")].copy()
    d["hy"] = hy_index(d.index)
    g = d.groupby("hy")
    out = pd.DataFrame(index=sorted(d.hy.unique()))
    if "pr" in d:
        n = g.pr.count()
        out["P"] = g.pr.sum(min_count=1).where(n >= 360)
        om = d.index.month.isin([10, 11, 12, 1, 2, 3])
        out["P_OctMar_share"] = d[om].groupby("hy").pr.sum() / out.P
        out["wetdays_ge1mm"] = g.pr.apply(lambda x: int((x >= 1).sum()))
        out["days_ge20mm"] = g.pr.apply(lambda x: int((x >= 20).sum()))
        out["P_Sep_Nov"] = d[d.index.month.isin([9, 10, 11])].groupby("hy").pr.sum()
    if "tas" in d:
        out["T_mean"] = g.tas.mean()
        jja = d[d.index.month.isin([6, 7, 8])]
        if "tmax" in d and d.tmax.notna().any():
            out["Tmax_JJA_mean"] = jja.groupby("hy").tmax.mean()
            out["Tmax_annual_max"] = g.tmax.max()
            d["pet"] = hargreaves(d)
            out["PET_hargreaves"] = d.groupby("hy").pet.sum(min_count=360)
        else:
            out["T_JJA_mean"] = jja.groupby("hy").tas.mean()
    if "evap" in d and d.evap.notna().any():
        out["E_model"] = g.evap.sum(min_count=360)
    return out


def trend(y):
    y = y.dropna()
    x = y.index.values.astype(float)
    r = stats.linregress(x, y.values)
    tcrit = stats.t.ppf(0.975, len(y) - 2)
    # lag-1 autocorrelation-adjusted effective n for the CI
    res = y.values - (r.intercept + r.slope * x)
    r1 = np.corrcoef(res[:-1], res[1:])[0, 1]
    neff = len(y) * (1 - r1) / (1 + r1) if r1 > 0 else len(y)
    se_adj = r.stderr * np.sqrt((len(y) - 2) / max(neff - 2, 1))
    ts = stats.theilslopes(y.values, x)
    return {"ols_per_decade": 10 * r.slope, "ci95_per_decade": [10 * (r.slope - tcrit * r.stderr), 10 * (r.slope + tcrit * r.stderr)],
            "ci95_ar1_adjusted_per_decade": [10 * (r.slope - tcrit * se_adj), 10 * (r.slope + tcrit * se_adj)],
            "p_value": r.pvalue, "theil_sen_per_decade": 10 * ts.slope, "n": int(len(y))}


def lag1(y):
    y = y.dropna().values
    return float(np.corrcoef(y[:-1], y[1:])[0, 1])


def summary(tab, ref_P=None):
    s = {}
    for col in tab.columns:
        y = tab[col].dropna()
        if len(y) < 5:
            continue
        s[col] = {"mean": y.mean(), "sd": y.std(ddof=1), "cv": y.std(ddof=1) / y.mean() if col != "T_mean" else None,
                  "q10": y.quantile(.1), "q20": y.quantile(.2), "q50": y.quantile(.5), "q80": y.quantile(.8),
                  "q90": y.quantile(.9), "min": y.min(), "max": y.max(), "lag1_autocorr": lag1(y), "n": int(len(y))}
    if ref_P is not None and "P" in tab:
        p20, p80 = ref_P.quantile(.2), ref_P.quantile(.8)
        P = tab.P.dropna()
        dry = (P < p20).values
        s["dry_wet_years"] = {"obs_p20_mm": p20, "obs_p80_mm": p80, "n_years": int(len(P)),
                              "n_dry_below_obs_p20": int(dry.sum()), "n_wet_above_obs_p80": int((P > p80).sum()),
                              "longest_dry_run": int(max((len(r) for r in "".join("1" if v else "0" for v in dry).split("0")), default=0)),
                              "n_consecutive_dry_pairs": int(np.sum(dry[1:] & dry[:-1]))}
    return s


def compare(a, b):
    a, b = a.dropna(), b.dropna()
    ks = stats.ks_2samp(a, b)
    return {"ks_stat": ks.statistic, "ks_p": ks.pvalue, "mean_diff": a.mean() - b.mean(),
            "mean_ratio": a.mean() / b.mean() if b.mean() != 0 else None,
            "quantile_diff": {q: a.quantile(q) - b.quantile(q) for q in (0.1, 0.2, 0.5, 0.8, 0.9)},
            "sd_ratio": a.std(ddof=1) / b.std(ddof=1)}


def seasonal(d):
    d = d[(d.index >= f"{HY0}-09-01") & (d.index <= f"{HY1 + 1}-08-31")]
    m = {}
    if "pr" in d:
        mon = d.pr.resample("MS").sum(min_count=25)
        m["P_mm_month"] = mon.groupby(mon.index.month).mean().to_dict()
    if "tas" in d:
        m["T_C"] = d.tas.groupby(d.index.month).mean().to_dict()
    return m


# ---------------- bias correction ----------------
def cal_mask(idx, rng):
    hy = hy_index(idx)
    return (hy >= rng[0]) & (hy <= rng[1])


def bc_scaling(dt, obs_p, obs_t):
    out = dt.copy()
    mdt, mo = cal_mask(dt.index, CAL), cal_mask(obs_p.index, CAL)
    fp = {}
    for m in range(1, 13):
        a = obs_p[mo & (obs_p.index.month == m)].mean()
        b = dt.pr[mdt & (dt.index.month == m)].mean()
        fp[m] = a / b if b > 0 else 1.0
    out["pr"] = dt.pr * dt.index.month.map(fp).values
    ft = {}
    if "tas" in dt and obs_t is not None:
        mt = cal_mask(obs_t.index, CAL)
        for c in ("tas", "tmax", "tmin"):          # separate monthly offsets: keeps the observed diurnal range
            ft[c] = {m: obs_t[c][mt & (obs_t.index.month == m)].mean() - dt[c][mdt & (dt.index.month == m)].mean()
                     for m in range(1, 13)}
            out[c] = dt[c] + dt.index.month.map(ft[c]).values
    return out, {"pr_factor": fp, "t_offset_K": ft}


def bc_qm(dt, obs_p):
    out = dt.copy()
    mdt, mo = cal_mask(dt.index, CAL), cal_mask(obs_p.index, CAL)
    qs = np.linspace(0.005, 0.995, 100)
    pr = dt.pr.copy()
    for m in range(1, 13):
        xd = dt.pr[mdt & (dt.index.month == m)].dropna().values
        xo = obs_p[mo & (obs_p.index.month == m)].dropna().values
        qd, qo = np.quantile(xd, qs), np.quantile(xo, qs)
        sel = dt.index.month == m
        v = dt.pr[sel].values
        mapped = np.interp(v, qd, qo)
        top = v > qd[-1]
        mapped[top] = v[top] * (qo[-1] / qd[-1] if qd[-1] > 0 else 1.0)
        mapped[v <= 0] = 0.0
        pr[sel] = mapped
    out["pr"] = pr
    return out


# ---------------- main ----------------
obs_e, eobs_cells = eobs_daily()
obs_r, era5_pts = era5_daily()
DT, DTcells = {}, {}
for m in MODELS:
    d, c = dt_daily(m)
    if len(d):
        DT[m], DTcells[m] = d, c

tabs = {"E-OBS": hy_table(obs_e), "ERA5": hy_table(obs_r)}
for m, d in list(DT.items()):
    tabs[m] = hy_table(d)
    sc, factors = bc_scaling(d, obs_e.pr, obs_r)
    tabs[m + " BC-scaling"] = hy_table(sc)
    DT[m + " BC-scaling"] = sc
    DT[m + " BC-scaling"].attrs["factors"] = factors
    qm = bc_qm(d, obs_e.pr)
    for c in ("tas", "tmax", "tmin"):
        if c in sc:
            qm[c] = sc[c]                       # temperature: same monthly offset as BC1
    tabs[m + " BC-QM"] = hy_table(qm)
    DT[m + " BC-QM"] = qm

refP = tabs["E-OBS"].P
res = {"question": "Is DestinE Climate DT gen-2 usable to drive Donana pond-hydroperiod scenarios?",
       "site": "Donana, Natura 2000 ES0000024", "period_hydrological_years": [HY0, HY1],
       "hy_definition": "HY y = 1 Sep y - 31 Aug y+1",
       "dt": {"generation": "2", "models": list(DT_k for DT_k in DT if " BC" not in DT_k),
              "experiments": "hist (activity baseline) 1990-2014 + SSP3-7.0 (activity projections) 2015-2024",
              "cells_spherical_healpix_l10_nested": {m: DTcells[m] for m in DTcells},
              "cell_centres_latlon": {str(k): v for k, v in CELLS.items()}},
       "obs": {"E-OBS_v31.0e_0.1deg_cells_used": eobs_cells, "ERA5_gridpoints_used": era5_pts},
       "calibration_hy": list(CAL), "validation_hy": list(VAL),
       "summary": {}, "compare_vs_EOBS_P": {}, "compare_vs_ERA5_T": {}, "trend": {}, "seasonal": {},
       "validation_2014_2023": {}, "bc_factors": {}}
for k, t in tabs.items():
    res["summary"][k] = summary(t, refP)
    res["seasonal"][k] = seasonal(DT[k] if k in DT else (obs_e if k == "E-OBS" else obs_r))
    res["trend"][k] = {c: trend(t[c]) for c in ("P", "T_mean", "Tmax_JJA_mean", "PET_hargreaves") if c in t and t[c].notna().sum() > 10}
    if k not in ("E-OBS",):
        res["compare_vs_EOBS_P"][k] = {c: compare(t[c], tabs["E-OBS"][c]) for c in ("P", "P_OctMar_share", "wetdays_ge1mm", "days_ge20mm")}
    if k not in ("E-OBS", "ERA5") and "T_mean" in t:
        res["compare_vs_ERA5_T"][k] = {c: compare(t[c], tabs["ERA5"][c]) for c in ("T_mean", "Tmax_JJA_mean", "Tmax_annual_max", "PET_hargreaves") if c in t}
    v = t.loc[VAL[0]:VAL[1]]
    o = tabs["E-OBS"].loc[VAL[0]:VAL[1]]
    orr = tabs["ERA5"].loc[VAL[0]:VAL[1]]
    res["validation_2014_2023"][k] = {
        "P_mean": v.P.mean(), "P_mean_bias_pct_vs_EOBS": 100 * (v.P.mean() / o.P.mean() - 1),
        "P_cv": v.P.std() / v.P.mean(), "OctMar_share_mean": v.P_OctMar_share.mean(),
        "wetdays_mean": v.wetdays_ge1mm.mean(),
        "T_mean": v.T_mean.mean() if "T_mean" in v else None,
        "T_bias_K_vs_ERA5": (v.T_mean.mean() - orr.T_mean.mean()) if "T_mean" in v else None,
        "Tmax_JJA_bias_K_vs_ERA5": (v.Tmax_JJA_mean.mean() - orr.Tmax_JJA_mean.mean()) if "Tmax_JJA_mean" in v else None,
        "PET_bias_pct_vs_ERA5": 100 * (v.PET_hargreaves.mean() / orr.PET_hargreaves.mean() - 1) if "PET_hargreaves" in v else None,
        "ks_P_vs_EOBS_p": stats.ks_2samp(v.P.dropna(), o.P.dropna()).pvalue, "n_years": int(v.P.notna().sum())}
    if "factors" in getattr(DT.get(k), "attrs", {}):
        res["bc_factors"][k] = DT[k].attrs["factors"]

for k, t in tabs.items():
    t.add_prefix(k + "|").to_csv(os.path.join(DER, f"hy_{k.replace(' ', '_')}.csv"), index_label="hy")


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 4)
    if isinstance(o, np.integer):
        return int(o)
    return o


res["attribution"] = ("DT-derived statistics: this data is created based on data of the European Union, using the "
                      "Destination Earth Platform, but has been modified (point extraction, daily/annual aggregation, "
                      "bias correction). ERA5: generated using Copernicus Climate Change Service information. "
                      "E-OBS: EU-FP6 UERRA, C3S and the ECA&D data providers.")
json.dump(clean(res), open(os.path.join(ROOT, "results.json"), "w"), indent=1)
print("written results.json; tables:", ", ".join(tabs))
