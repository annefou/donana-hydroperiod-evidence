"""00b_fetch_era5.py <root> — ERA5 monthly-mean total precipitation over Doñana, 1983–2026.
Source: Copernicus Climate Change Service (C3S) Climate Data Store, dataset
'reanalysis-era5-single-levels-monthly-means' (doi:10.24381/cds.f17050d7), licence CC BY 4.0.
Needs a CDS account: the cdsapi client reads ~/.cdsapirc (never stored in this repository).
Area 37.10N -6.60E 36.85N -6.25E (2x2 0.25-degree cells).                -> data/raw/era5/tp_monthly.nc
Plus hourly 10 m wind at the S1 overpass hours (see below).
"""
import sys
from pathlib import Path

import cdsapi

ROOT = Path(sys.argv[1]); out = ROOT / "data/raw/era5/tp_monthly.nc"; out.parent.mkdir(parents=True, exist_ok=True)
if not out.exists():
    cdsapi.Client(quiet=True).retrieve("reanalysis-era5-single-levels-monthly-means", {
        "product_type": ["monthly_averaged_reanalysis"], "variable": ["total_precipitation"],
        "year": [str(y) for y in range(1983, 2027)], "month": [f"{m:02d}" for m in range(1, 13)], "time": ["00:00"],
        "area": [37.10, -6.60, 36.85, -6.25], "data_format": "netcdf", "download_format": "unarchived"}, str(out))
print(out, out.stat().st_size)
# 10 m wind at the Sentinel-1 overpass hours (descending ~06:30 UTC, ascending ~18:20 UTC), 2017-09 .. 2026-08,
# dataset 'reanalysis-era5-single-levels' (doi:10.24381/cds.adbb2d47), CC BY 4.0   -> data/raw/era5/wind_<year>.nc
for y in range(2017, 2027):
    o = ROOT / f"data/raw/era5/wind_{y}.nc"
    if o.exists():
        continue
    cdsapi.Client(quiet=True).retrieve("reanalysis-era5-single-levels", {
        "product_type": ["reanalysis"], "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
        "year": [str(y)], "month": [f"{m:02d}" for m in range(1, 13)], "day": [f"{d:02d}" for d in range(1, 32)],
        "time": ["06:00", "07:00", "18:00", "19:00"], "area": [37.05, -6.50, 36.95, -6.30],
        "data_format": "netcdf", "download_format": "unarchived"}, str(o))
    print(o, o.stat().st_size)
