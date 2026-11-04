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
# # 01 — Data
#
# **DRAFT, not reviewed; prototype (TRL 3: the scripts ran end to end on 22 October 2026, but not yet from a clean
# clone).** LifeWatch ERIC evidence for Doñana (Natura 2000 ES0000024) and the technical side of the FAIR2Adapt
# Doñana transfer case. Reference method: Kordelas et al. 2018, *Remote Sensing* 10(6):910,
# [doi:10.3390/rs10060910](https://doi.org/10.3390/rs10060910).
#
# The pipeline (`pipeline/`, `Snakefile`) is **not** rerun when this book is built: it reads every Sentinel-2 date
# 2017–2026 and needs a CDS account (ERA5) and DestinE upgraded access (Climate DT). These pages read the summary
# results committed in `results/` and the figures in `figures/`.
#
# | Input | Source | Licence | Account |
# |---|---|---|---|
# | Sentinel-2 L2A, Sentinel-1 | Earth Search, Microsoft Planetary Computer | Copernicus terms | none |
# | LAST-EBD Landsat flood masks and hydroperiod maps; field points 24 Apr 2017 | Díaz-Delgado et al., Zenodo | CC BY 4.0 | none |
# | ICTS-Doñana marsh water levels | ICTS-Doñana | CC BY 4.0 | none |
# | CHG piezometers | Confederación Hidrográfica del Guadalquivir | reproduction citing the source; "not official" data | none |
# | Natura 2000 boundary; OpenStreetMap water | EEA; OSM | EEA terms; ODbL | none |
# | ERA5 | Copernicus Climate Data Store | Copernicus licence | **CDS** |
# | E-OBS v31.0e | ECA&D | non-commercial research only | none |
# | Climate DT generation 2 (IFS-NEMO, IFS-FESOM) | DestinE Platform | DestinE terms: derived statistics only | **DestinE upgraded access** |
# | Copernicus Land Water Cover Duration (check not yet run) | Copernicus Data Space | Copernicus terms | **CDSE token** |
#
# Grid: HEALPix depth 19 (~12 m), healpix-geo, WGS84 ellipsoid, NESTED; reflectance resampled conservatively and the
# scene classification by majority, before classification.

# %%
import json
from pathlib import Path

import pandas as pd

R = Path("../results")
pd.DataFrame([{"file": str(f.relative_to(R)), "bytes": f.stat().st_size} for f in sorted(R.rglob("*.json"))])

# %% [markdown]
# Resampler checks on the HEALPix grid (`pipeline/flooding/01b_test_resampler.py`): a constant field stays constant,
# the integral is preserved, and the overlap weights match a Monte-Carlo estimate.

# %%
pd.Series(json.loads((R / "flooding/results_resampler_tests.json").read_text())).to_frame("value")
