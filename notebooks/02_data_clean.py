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
# # 02 — Flood mapping and its validation
#
# Flooding is mapped with a re-implementation of Kordelas et al. 2018 on every Sentinel-2 date, and checked against
# the LAST-EBD Landsat flood masks (same-day pairs) and the 936 field points of 24 April 2017. Kappa and overall
# accuracy (OA) by stratum; three variants: open water, open water plus flooded vegetation, a fixed B11 threshold.

# %%
import json
from pathlib import Path

import pandas as pd

R = Path("../results/flooding")
v = json.loads((R / "results_validation.json").read_text())
print("same-day pairs:", v["n_pairs"], v["pairs_by_sensor"])
pd.DataFrame({(var, s): {k: x["EAP"].get(k) for k in ("kappa", "OA")} for var, d in v["combined"].items()
              for s, x in d.items()}).T.round(3)

# %% [markdown]
# Sentinel-1 for cloudy dates, against LAST-EBD and against Sentinel-2 open water:

# %%
print(json.dumps(json.loads((R / "results_validation_s1.json").read_text()), indent=1)[:2500])
