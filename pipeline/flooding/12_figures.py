"""12_figures.py <root> — figures (PNG, 200 dpi) drawn directly from HEALPix cells (cell polygons from
healpix_geo.nested.vertices, WGS84), never from a UTM raster.

fig1_hydroperiod_years.png   hydroperiod per hydrological year, depth-15 cells (~200 m), 2017/18–2025/26
fig2_hydroperiod_latest.png  latest year at depth 17 (~50 m) + valid observations per cell
fig3_marsh_series.png        marsh-mean hydroperiod 1974/75–2025/26 (LAST vs S2) and anomaly map vs 1984–2017 baseline
fig4_validation.png          validation table (this study vs published)
fig5_lagoon_groundwater.png  Santa Olalla / Dulce hydroperiod vs groundwater level and rain
fig6_flooded_area.png        flooded marsh area per S2 date (and S1 fill dates)
"""
import json, sys
from pathlib import Path

import healpix_geo.nested as hn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm, LinearSegmentedColormap

ROOT = Path(sys.argv[1]); FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
S1C, S2C, S3C = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlecolor": INK, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.facecolor": SURF})
BLUES = ["#f0efec", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
HP_BOUNDS = [0, 1, 30, 60, 90, 120, 180, 270, 367]
HP_CMAP = ListedColormap(BLUES); HP_NORM = BoundaryNorm(HP_BOUNDS, HP_CMAP.N)
DIV = LinearSegmentedColormap.from_list("div", ["#a32b2b", "#e34948", "#f6c9c8", "#f0efec", "#b7d3f6", "#3987e5", "#104281"])


def polys(cids, level):
    lon, lat = hn.vertices(np.asarray(cids, dtype=np.uint64), level, ellipsoid="WGS84")
    lon = (lon + 180) % 360 - 180
    return np.stack([lon, lat], -1)


def draw(ax, cids, level, vals, cmap, norm):
    ok = np.isfinite(vals)
    pc = PolyCollection(polys(cids[ok], level), array=vals[ok], cmap=cmap, norm=norm, edgecolors="face",
                        linewidths=0.05, rasterized=True)
    ax.add_collection(pc); ax.set_xlim(-6.72, -6.12); ax.set_ylim(36.78, 37.36)
    ax.set_aspect(1 / np.cos(np.radians(37))); ax.grid(False)
    ax.tick_params(labelsize=7)
    return pc


hz = zarr.open_group(ROOT / "data/products/hydroperiod.zarr", mode="r")
years = hz.attrs["hydrological_years"]
c15 = hz["hydroperiod/15/cell_ids"][:]; h15 = hz["hydroperiod/15/hydroperiod_days"][:]
_st = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")["stratum"][:]
_c19 = hz["hydroperiod/19/cell_ids"][:]


def sea_share(level):
    par = _c19 >> np.uint64(2 * (19 - level)); u, iv = np.unique(par, return_inverse=True)
    return u, np.bincount(iv, weights=(_st == 0).astype(float), minlength=len(u)) / np.bincount(iv, minlength=len(u))


_u15, _s15 = sea_share(15)
h15 = np.where(_s15[None, :] > 0.5, np.nan, h15)
# ---- fig1
fig, axs = plt.subplots(3, 3, figsize=(11, 11.5), constrained_layout=True)
for k, ax in enumerate(axs.ravel()):
    if k >= len(years):
        ax.set_visible(False); continue
    pc = draw(ax, c15, 15, h15[k], HP_CMAP, HP_NORM)
    ax.set_title(years[k], fontsize=11)
cb = fig.colorbar(pc, ax=axs, shrink=0.6, ticks=HP_BOUNDS[:-1] + [365], label="days flooded (open water, Sentinel-2)")
fig.suptitle("Doñana hydroperiod per hydrological year (Sep–Aug), Sentinel-2 open water; HEALPix depth 15 (~200 m mean of depth-19 cells); sea and estuary masked",
             color=INK, fontsize=12)
fig.savefig(FIG / "fig1_hydroperiod_years.png", dpi=200); plt.close(fig)

# ---- fig2 latest year at depth 17 + n_valid
c19 = hz["hydroperiod/19/cell_ids"][:]
k = len(years) - 1
c17 = hz["hydroperiod/17/cell_ids"][:]; h17 = hz["hydroperiod/17/hydroperiod_days"][k]
_u17, _s17 = sea_share(17); h17 = np.where(_s17 > 0.5, np.nan, h17)
nv = hz["hydroperiod/19/n_valid"][k].astype(float); nv[nv < 0] = np.nan
par = c19 >> np.uint64(4); up, inv = np.unique(par, return_inverse=True)
nv17 = np.bincount(inv, weights=np.nan_to_num(nv), minlength=len(up)) / np.bincount(inv, minlength=len(up))
fig, axs = plt.subplots(1, 2, figsize=(12, 6.2), constrained_layout=True)
pc = draw(axs[0], c17, 17, h17, HP_CMAP, HP_NORM)
fig.colorbar(pc, ax=axs[0], shrink=0.75, ticks=HP_BOUNDS[:-1] + [365], label="days flooded")
axs[0].set_title(f"Hydroperiod {years[k]} (depth 17, ~50 m)")
vmin, vmax = np.nanpercentile(nv17, [1, 99])
pc2 = draw(axs[1], up, 17, nv17, LinearSegmentedColormap.from_list("b", BLUES[1:]), plt.Normalize(vmin, vmax))
fig.colorbar(pc2, ax=axs[1], shrink=0.75, label="valid Sentinel-2 observations per cell")
axs[1].set_title(f"Valid observations {years[k]}")
fig.savefig(FIG / "fig2_hydroperiod_latest.png", dpi=200); plt.close(fig)

# ---- fig3 marsh series + anomaly map
ser = pd.read_csv(ROOT / "data/products/marsh_series.csv")
an = json.load(open(ROOT / "results_anomaly.json"))
fig = plt.figure(figsize=(13, 5.6), constrained_layout=True)
gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1])
ax = fig.add_subplot(gs[0])
L = ser[ser.source == "LAST"]; S = ser[ser.source == "S2"]
ax.axhspan(an["baseline_marsh_mean_hp"] - an["baseline_marsh_sd"], an["baseline_marsh_mean_hp"] + an["baseline_marsh_sd"],
           color="#e1e0d9", alpha=0.6, lw=0)
ax.axhline(an["baseline_marsh_mean_hp"], color=MUTED, lw=1, ls="--")
ax.plot(L.hy, L.marsh_mean_hp, "-o", color=S1C, lw=2, ms=4, label="LAST-EBD Landsat maps (Zenodo 3516929)")
ax.plot(S.hy, S.marsh_mean_hp, "-s", color=S2C, lw=2, ms=5, label="Sentinel-2, this pipeline (HEALPix d19)")
ax.text(1975, an["baseline_marsh_mean_hp"] + an["baseline_marsh_sd"] - 9, "1984/85–2016/17 mean ± 1 sd", color=INK2, fontsize=8)
off = an.get("offset_S2_minus_LAST_2017_18_days")
ax.text(2001, 8, f"S2 reads {abs(off):.0f} days below Landsat in the 2017/18 overlap\n(23–44 d across checks): compare S2 years with each other", color=INK2, fontsize=8)
ax.set_xlabel("hydrological year (start)"); ax.set_ylabel("marsh-mean hydroperiod (days)")
ax.set_title("Seasonal marsh (30,900 ha): mean days flooded per year"); ax.legend(frameon=False, fontsize=8, loc="upper right")
ax2 = fig.add_subplot(gs[1])
ag = zarr.open_group(ROOT / "data/products/anomaly.zarr", mode="r")
a19 = ag["19/anomaly_days"][len(years) - 1]
ok = np.isfinite(a19)
s_ = np.bincount(inv[ok], weights=a19[ok], minlength=len(up)); w_ = np.bincount(inv[ok], minlength=len(up))
a17 = np.where(w_ > 0, s_ / np.maximum(w_, 1), np.nan)
st = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")["stratum"][:]
m17 = np.bincount(inv, weights=np.isin(st, (1, 2)).astype(float), minlength=len(up)) / np.bincount(inv, minlength=len(up))
a17[m17 < 0.5] = np.nan
pc = draw(ax2, up, 17, a17, DIV, TwoSlopeNorm(0, -150, 150))
fig.colorbar(pc, ax=ax2, shrink=0.8, label="days vs 1984/85–2016/17 mean (blue wetter, red drier)")
ax2.set_title(f"Raw anomaly {years[-1]} (S2 minus Landsat baseline), marsh + lagoons", fontsize=10); ax2.set_xlim(-6.55, -6.15); ax2.set_ylim(36.85, 37.20)
fig.savefig(FIG / "fig3_marsh_series_anomaly.png", dpi=200); plt.close(fig)

# ---- fig4 validation table
v = json.load(open(ROOT / "results_validation.json"))
rows = []
for sname, lab in (("complete", "Complete (no sea)"), ("marsh", "Seasonal marsh"), ("lagoon", "Lagoons (OSM)"), ("other", "Other land")):
    for var, vl in (("open_water", "open water (Kordelas)"), ("fixed_B11_0.186", "B11<0.186 (Landsat rule)")):
        e = v["combined"][var][sname]["EAP"]
        if not e:
            continue
        rows.append([lab, vl, f"{e['kappa']:.2f} [{e['kappa_ci95'][0]:.2f}–{e['kappa_ci95'][1]:.2f}]",
                     f"{100 * e['PA_water']:.0f}", f"{100 * e['UA_water']:.0f}", f"{100 * e['omission_water']:.0f}",
                     f"{100 * e['commission_water']:.0f}", f"{v['combined'][var][sname]['IAP'].get('kappa', float('nan')):.2f}"])
fp = v["field_points_20170424"]; fk = list(fp)[0] if fp else None
if fk:
    t_ = fp[fk]["Tfinal"]
    rows.append(["Field points 24 Apr 2017 (936)", f"open water, scene {fk}", "—", f"{100 * t_['In']['share_mapped_water']:.0f} (In)",
                 "—", "—", f"{100 * t_['Ss']['share_mapped_water']:.1f} (Ss)", "—"])
pub = [["Kordelas 2018, marsh (Table 8)", "unsupervised", "0.94", "91", "100", "9", "0", "—"],
       ["Kordelas 2018, complete (Table 8)", "unsupervised", "0.88", "90", "89", "10", "11", "—"],
       ["Kordelas 2018, temp. ponds (Table 8)", "unsupervised", "0.66", "90", "53", "10", "47", "—"],
       ["Kordelas 2019, marsh / Reserve (Fig.10)", "best variant", "0.94 / 0.85", "—", "—", "—", "—", "—"],
       ["Díaz-Delgado 2016 Landsat vs field (§3.1)", "band-5 threshold", "0.65 (0.61–0.66)", "—", "—", "—", "—", "—"]]
fig, ax = plt.subplots(figsize=(15, 0.33 * (len(rows) + len(pub)) + 0.8))
ax.axis("off")
cols = ["Stratum", "S2 water class", "Kappa [95% CI]", "PA %", "UA %", "Omission %", "Commission %", "Kappa incl. edges"]
t = ax.table(cellText=rows + pub, colLabels=cols, loc="upper left", cellLoc="left", colLoc="left",
             colWidths=[0.24, 0.17, 0.15, 0.07, 0.07, 0.08, 0.09, 0.09])
t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.45)
for (r, cidx), cell in t.get_celld().items():
    cell.set_edgecolor(GRID)
    if r == 0:
        cell.set_text_props(weight="bold", color=INK)
    if r > len(rows):
        cell.set_text_props(color=INK2); cell.set_facecolor("#f6f5f2")
ax.set_title(f"Sentinel-2 (HEALPix d19) vs LAST-EBD Landsat masks, {v['n_pairs']} date pairs (≤1 day), boundary cells excluded "
             "(Kordelas 2018 protocol). Grey rows: published values. The water-vegetation subclass never fired, so open+veg = open.",
             fontsize=9.5, color=INK, loc="left")
fig.savefig(FIG / "fig4_validation.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---- fig5 lagoon vs groundwater
p = ROOT / "data/products/lagoon_gw_annual.csv"
if p.exists():
    d = pd.read_csv(p)
    lags = [l for l in ("Laguna de Santa Olalla", "Laguna Dulce") if l in set(d.lagoon)]
    fig, axs = plt.subplots(3, len(lags), figsize=(6 * len(lags), 8.5), sharex=True, constrained_layout=True, squeeze=False)
    for j, lg in enumerate(lags):
        x = d[d.lagoon == lg]
        for src, col, mk in ((0, S1C, "o"), (1, S2C, "s")):
            y = x[x.s2 == src]
            axs[0, j].plot(y.hy, y.hp, mk + "-", color=col, ms=4, lw=1.6, label="LAST Landsat" if src == 0 else "Sentinel-2")
        axs[0, j].set_ylabel("lagoon hydroperiod (days)"); axs[0, j].set_title(lg); axs[0, j].legend(frameon=False, fontsize=8)
        axs[1, j].plot(x.hy, x.gw, "o-", color=S3C, ms=4, lw=1.6)
        axs[1, j].set_ylabel("piezometric level (m a.s.l.)\nCHG annual mean")
        axs[2, j].bar(x.hy, x.rain, color="#9ec5f4", width=0.8)
        axs[2, j].set_ylabel("rain Sep–Aug (mm, ERA5)"); axs[2, j].set_xlabel("hydrological year (start)")
    fig.savefig(FIG / "fig5_lagoon_groundwater.png", dpi=200); plt.close(fig)

# ---- fig6 flooded marsh area per date
sc = pd.read_csv(ROOT / "data/s2/scenes.csv"); sc = sc[sc.status == "ok"].copy(); sc["date"] = pd.to_datetime(sc.date)
f6 = ROOT / "data/products/marsh_area_by_date.csv"
if f6.exists():
    m = pd.read_csv(f6, parse_dates=["date"])
    fig, ax = plt.subplots(figsize=(13, 4.2), constrained_layout=True)
    for src, col, mk in (("S2", S1C, "o"), ("S1", S2C, "^")):
        y = m[m.source == src]
        if len(y):
            ax.plot(y.date, y.flooded_share * 100, mk, color=col, ms=3.5, label=f"{src} open water" + (" (cloudy-season fill)" if src == "S1" else ""))
    ax.set_ylabel("seasonal marsh flooded (% of valid cells)"); ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title("Flooded share of the seasonal marsh on each usable date, Sep 2017 – Aug 2026")
    fig.savefig(FIG / "fig6_flooded_area.png", dpi=200); plt.close(fig)
print("figures written")
