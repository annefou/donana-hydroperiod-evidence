# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 03 — Hydroperiod, anomaly and groundwater
#
# Hydroperiod = days flooded per HEALPix cell and hydrological year (1 September to 31 August). The Sentinel-2 era is
# compared with the LAST-EBD Landsat record 1984/85–2016/17; the two overlap in 2017/18, where the offset between
# them is measured rather than assumed.

# %%
import json
from pathlib import Path

import pandas as pd

R = Path("../results/flooding")
h = json.loads((R / "results_hydroperiod.json").read_text())
cols = ["hydro_year", "n_s2_dates_used", "marsh_mean_hp_days", "marsh_flooded_ge1d_ha", "marsh_flooded_ge90d_ha",
        "marsh_halfsample_mad_days_median"]
pd.DataFrame(h["years"])[cols].round(1)

# %%
print("Overlap with LAST-EBD (2017/18):", json.dumps(h["last_overlap"], indent=1)[:1500])

# %%
a = json.loads((R / "results_anomaly.json").read_text())
pd.Series({k: v for k, v in a.items() if k != "years"}).to_frame("value")

# %% [markdown]
# Lagoon hydroperiod against groundwater (CHG piezometers; data reproduced citing the source, not official data):

# %%
g = json.loads((R / "results_groundwater.json").read_text())
print("lagoon -> piezometer:", g["pairs"])
print(json.dumps(g["lagoons"], indent=1)[:2000])

# %% [markdown]
# Internal checks (`pipeline/flooding/15_checks.py`): 10 m against the RGB overview, stability, detection by size,
# and the EBD press release.

# %%
print(json.dumps(json.loads((R / "results_checks.json").read_text()), indent=1)[:2500])
