# Doñana (Natura 2000 ES0000024): flooding regime from Sentinel-2/-1 on HEALPix, and a Climate DT fitness test.
#
#   pixi run snakemake --cores 1              flooding regime (open data; ERA5 needs a CDS account). Long: every
#                                             Sentinel-2 date 2017-2026
#   pixi run snakemake --cores 1 climate_dt   Climate DT fitness test (DestinE upgraded access and a CDS account)
#   pixi run snakemake --cores 1 flood_wcd    Copernicus Land Water Cover Duration check (CDSE_TOKEN, free
#                                             Copernicus Data Space account); not run yet
#   pixi run snakemake --cores 1 --dry-run    what would run, without running it
#
# The scripts write into work/ (gitignored); the last step of each rule copies the summary results and figures
# into results/ and figures/, which are committed. The same steps are the pixi tasks in pixi.toml.

W = "work"
FLOOD_RESULTS = ["anomaly", "checks", "groundwater", "hydroperiod", "resampler_tests", "strata", "validation", "validation_s1"]
FLOOD_FIGS = ["fig1_hydroperiod_years", "fig2_hydroperiod_latest", "fig3_marsh_series_anomaly", "fig4_validation",
              "fig5_lagoon_groundwater", "fig6_flooded_area"]
DT_FIGS = ["fig1_annual_rainfall", "fig2_seasonal_cycle"]


rule all:
    input:
        expand("results/flooding/results_{r}.json", r=FLOOD_RESULTS),
        expand("figures/flooding/{f}.png", f=FLOOD_FIGS),


rule climate_dt:
    input:
        "results/climate_dt/results.json",
        expand("figures/climate_dt/{f}.png", f=DT_FIGS),


rule flooding:
    output:
        expand("results/flooding/results_{r}.json", r=FLOOD_RESULTS),
        expand("figures/flooding/{f}.png", f=FLOOD_FIGS),
    params: R=f"{W}/flooding"
    shell:
        """
        python -I pipeline/flooding/00_fetch_inputs.py {params.R} && python -I pipeline/flooding/00b_fetch_era5.py {params.R}
        python -I pipeline/flooding/01_grid.py {params.R} && python -I pipeline/flooding/01b_test_resampler.py {params.R}
        python -I pipeline/flooding/04_last_to_healpix.py {params.R} && python -I pipeline/flooding/04b_last_hydroperiod.py {params.R}
        python -I pipeline/flooding/06_strata.py {params.R}
        python -I pipeline/flooding/02_catalog.py {params.R} && python -I pipeline/flooding/02b_catalog_pc.py {params.R}
        python -I pipeline/flooding/05_run_s2.py {params.R} 2017-09-01 2026-08-31
        python -I pipeline/flooding/08_s1_catalog.py {params.R} && python -I pipeline/flooding/09_s1_process.py {params.R}
        python -I pipeline/flooding/07_validate.py {params.R} && python -I pipeline/flooding/07b_validate_s1.py {params.R}
        python -I pipeline/flooding/10_hydroperiod.py {params.R} && python -I pipeline/flooding/10b_anomaly.py {params.R}
        python -I pipeline/flooding/11_groundwater.py {params.R} && python -I pipeline/flooding/12_figures.py {params.R}
        python -I pipeline/flooding/14_dt_export.py {params.R} && python -I pipeline/flooding/15_checks.py {params.R}
        mkdir -p results/flooding figures/flooding
        cp {params.R}/results_*.json results/flooding/ && cp {params.R}/figures/*.png figures/flooding/
        """


rule climate_dt_run:
    output:
        "results/climate_dt/results.json",
        expand("figures/climate_dt/{f}.png", f=DT_FIGS),
    params: R=f"{W}/climate_dt"
    shell:
        """
        sh pipeline/climate_dt/01_fetch_dt_all.sh {params.R} && python -I pipeline/climate_dt/02_fetch_era5.py {params.R}
        python -I pipeline/climate_dt/03_fetch_eobs.py {params.R} && python -I pipeline/climate_dt/04_analyse.py {params.R}
        python -I pipeline/climate_dt/05_figures.py {params.R}
        mkdir -p results/climate_dt figures/climate_dt
        cp {params.R}/results.json results/climate_dt/ && cp {params.R}/figures/*.png figures/climate_dt/
        """


rule flood_wcd:
    input: "results/flooding/results_validation.json"
    output: "results/flooding/results_wcd.json"
    params: R=f"{W}/flooding"
    shell: "python -I pipeline/flooding/13_wcd_validate.py {params.R} && cp {params.R}/results_wcd.json results/flooding/"
