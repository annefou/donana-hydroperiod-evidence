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
# # 04 — Figures and the Climate DT fitness test
#
# Written by `pipeline/flooding/12_figures.py` and `pipeline/climate_dt/05_figures.py`, committed in `figures/`.
#
# ![Hydroperiod by year](../figures/flooding/fig1_hydroperiod_years.png)
#
# ![Latest hydroperiod](../figures/flooding/fig2_hydroperiod_latest.png)
#
# ![Marsh series and anomaly](../figures/flooding/fig3_marsh_series_anomaly.png)
#
# ![Validation](../figures/flooding/fig4_validation.png)
#
# ![Lagoons and groundwater](../figures/flooding/fig5_lagoon_groundwater.png)
#
# ![Flooded area](../figures/flooding/fig6_flooded_area.png)
#
# ## Climate DT fitness test
# Does the DestinE Climate Change Adaptation Digital Twin (generation 2, IFS-NEMO and IFS-FESOM) reproduce Doñana
# rainfall and temperature statistics? Prototype test: one site, one realisation per model, 34 hydrological years.
# Verdict (`pipeline/climate_dt/REPORT.txt`): usable after bias correction, not as-is; both models as a two-member
# envelope, with E-OBS and ERA5 as reference.
#
# ![Annual rainfall](../figures/climate_dt/fig1_annual_rainfall.png)
#
# ![Seasonal cycle](../figures/climate_dt/fig2_seasonal_cycle.png)

# %%
print(open("../pipeline/climate_dt/REPORT.txt").read()[:3000])
