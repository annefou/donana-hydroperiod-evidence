"""Two figures from scripts/04_analyse.py outputs (ROOT/data/derived/hy_*.csv, ROOT/results.json).

fig1_annual_rainfall.png/.pdf: hydrological-year (Sep-Aug) rainfall and Oct-Mar share, HY 1990-2023, per source.
fig2_seasonal_cycle.png/.pdf: mean monthly precipitation and 2 m temperature (Sep..Aug), HY 1990-2023.
Categorical colours: validated default palette slots 1-4 (E-OBS, ERA5, IFS-NEMO, IFS-FESOM).
Usage: python -I scripts/05_figures.py <ROOT>
"""
import sys, os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = sys.argv[1]
DER = os.path.join(ROOT, "data", "derived")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
r = json.load(open(os.path.join(ROOT, "results.json")))
COL = {"E-OBS": "#2a78d6", "ERA5": "#eb6834", "IFS-NEMO": "#1baf7a", "IFS-FESOM": "#eda100"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})


def tab(k):
    f = os.path.join(DER, f"hy_{k.replace(' ', '_')}.csv")
    if not os.path.exists(f):
        return None
    t = pd.read_csv(f, index_col="hy")
    t.columns = [c.split("|", 1)[1] for c in t.columns]
    return t


order = [("E-OBS", "E-OBS", "obs"), ("ERA5", "ERA5", "obs"), ("IFS-NEMO", "IFS-NEMO", "raw"),
         ("IFS-NEMO BC-scaling", "IFS-NEMO", "bias-corr."), ("IFS-FESOM", "IFS-FESOM", "raw"),
         ("IFS-FESOM BC-scaling", "IFS-FESOM", "bias-corr.")]
order = [o for o in order if tab(o[0]) is not None]
p20 = r["summary"]["E-OBS"]["dry_wet_years"]["obs_p20_mm"]
p80 = r["summary"]["E-OBS"]["dry_wet_years"]["obs_p80_mm"]

fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.35, 1]})
rng = np.random.default_rng(0)
for ax, col, lab in ((axs[0], "P", "Hydrological-year rainfall (mm, Sep-Aug)"),
                     (axs[1], "P_OctMar_share", "Share of rainfall falling Oct-Mar")):
    for i, (k, fam, kind) in enumerate(order):
        y = tab(k)[col].dropna().values
        c = COL[fam]
        bp = ax.boxplot(y, positions=[i], widths=0.5, showfliers=False, patch_artist=True,
                        medianprops={"color": INK, "lw": 2}, whiskerprops={"color": INK2},
                        capprops={"color": INK2}, boxprops={"facecolor": "none", "edgecolor": INK2})
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(y)), y, s=22, color=c, alpha=0.85,
                   edgecolor="#fcfcfb", lw=0.8, zorder=3, marker="o" if kind != "bias-corr." else "D")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{fam.replace('IFS-', '')}\n{kind}" for _, fam, kind in order], fontsize=8.5)
    ax.set_ylabel(lab)
    ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for v, t in ((p20, f"obs P20 = {p20:.0f} mm (dry-year threshold)"), (p80, f"obs P80 = {p80:.0f} mm (wet-year threshold)")):
    axs[0].axhline(v, color=INK2, lw=1, ls="--", zorder=1)
    axs[0].annotate(t, (-0.45, v), xytext=(0, 3), textcoords="offset points", ha="left", va="bottom",
                    fontsize=7.5, color=INK2, bbox={"fc": "#fcfcfb", "ec": "none", "pad": 0.5})
axs[0].set_title("a  Annual rainfall, HY 1990-2023 (n = %d years per source)" % len(tab("E-OBS").P.dropna()), loc="left", fontsize=10.5, color=INK)
axs[1].set_title("b  Rainfall timing", loc="left", fontsize=10.5, color=INK)
fig.text(0.01, 0.005, "Donana ES0000024, site mean of DT HEALPix L10 cells (5 for rainfall). DT years are not in phase with real years: "
         "compare spreads, not points.\nDiamonds: bias-corrected (monthly scaling fitted on HY 1990-2013 vs E-OBS; HY 2014-2023 out of sample).",
         fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.06, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig1_annual_rainfall.{ext}"), dpi=200)

# ---- fig 2 ----
months = [9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8]
mlab = ["S", "O", "N", "D", "J", "F", "M", "A", "M", "J", "J", "A"]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for k in ("E-OBS", "ERA5", "IFS-NEMO", "IFS-FESOM"):
    s = r["seasonal"].get(k)
    if not s:
        continue
    for ax, key in ((axs[0], "P_mm_month"), (axs[1], "T_C")):
        if key not in s:
            continue
        y = [s[key][str(m)] for m in months]
        ax.plot(range(12), y, color=COL[k], lw=2, marker="o", ms=4,
                ls="-" if k.startswith("IFS") else (0, (4, 2)), label=k)
for ax, lab, t in ((axs[0], "Mean monthly precipitation (mm)", "a  Precipitation"),
                   (axs[1], "Mean 2 m temperature (deg C)", "b  Temperature")):
    ax.set_xticks(range(12)); ax.set_xticklabels(mlab); ax.set_xlim(-0.3, 11.3)
    ax.set_ylabel(lab); ax.set_title(t, loc="left", fontsize=10.5, color=INK)
    ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
axs[0].legend(frameon=False, fontsize=8.5, loc="upper right")
fig.suptitle("Mean seasonal cycle, HY 1990-2023 (Sep-Aug), raw Climate DT without bias correction", x=0.01, ha="left",
             fontsize=11, color=INK)
fig.text(0.01, 0.005, "Observations dashed (E-OBS 0.1 deg; ERA5 0.25 deg), Climate DT gen-2 solid (hist 1990-2014 + SSP3-7.0 2015-2024).",
         fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig2_seasonal_cycle.{ext}"), dpi=200)
print("figures written")
