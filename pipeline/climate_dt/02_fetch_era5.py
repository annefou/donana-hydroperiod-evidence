"""Fetch ERA5 hourly point time series (total precipitation, 2 m temperature, evaporation) at the Donana points.

Source: Copernicus Climate Change Service (C3S) Climate Data Store, dataset
"reanalysis-era5-single-levels-timeseries" (ERA5 hourly time-series data on single levels from 1940 to present),
https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-timeseries ,
DOI 10.24381/cds.adbb2d47 (parent ERA5 single levels). ERA5 is on a 0.25 deg grid; the time-series service
returns the nearest grid point (no interpolation).
Licence: CC-BY 4.0 (Copernicus licence, https://cds.climate.copernicus.eu/licences/cc-by). Cite Hersbach et al.
(2020) QJRMS 146:1999-2049 and "Generated using Copernicus Climate Change Service information".
Credential: ~/.cdsapirc (read by cdsapi; never copied or printed).
Points: the same DT cell centres as scripts/01_fetch_dt.py.

Usage: python -I scripts/02_fetch_era5.py <ROOT>
"""
import sys, os
import cdsapi

ROOT = sys.argv[1]
OUT = os.path.join(ROOT, "data", "raw", "era5")
os.makedirs(OUT, exist_ok=True)
CELLS = {3490355: (36.9352, -6.3281), 3490408: (37.1221, -6.2402), 3490451: (37.1221, -6.5039),
         3490362: (36.9819, -6.4600), 3490357: (36.9819, -6.1963)}
c = cdsapi.Client(quiet=True)
for cell, (la, lo) in CELLS.items():
    fn = os.path.join(OUT, f"era5_{cell}.zip")
    if os.path.exists(fn) and os.path.getsize(fn) > 10000:
        continue
    c.retrieve("reanalysis-era5-single-levels-timeseries",
               {"variable": ["2m_temperature", "total_precipitation", "evaporation"],
                "location": {"latitude": la, "longitude": lo},
                "date": ["1989-09-01/2024-12-31"], "data_format": "netcdf"}, fn)
    print("OK", cell, flush=True)
