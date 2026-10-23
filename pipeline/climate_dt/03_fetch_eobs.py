"""Subset E-OBS v31.0e (0.1 deg regular grid, ensemble mean) daily precipitation (rr) and mean temperature (tg)
over Donana, 1989-09-01 to 2024-12-31, by HTTP range reads (the full European files are 1.6 / 4.6 GB; only the
needed days are streamed, nothing large is written to disk).

Source: E-OBS gridded dataset, Copernicus Climate Change Service / ECA&D (KNMI), version 31.0e,
https://knmi-ecad-assets-prd.s3.amazonaws.com/ensembles/data/Grid_0.1deg_reg_ensemble/rr_ens_mean_0.1deg_reg_v31.0e.nc
https://knmi-ecad-assets-prd.s3.amazonaws.com/ensembles/data/Grid_0.1deg_reg_ensemble/tg_ens_mean_0.1deg_reg_v31.0e.nc
Landing page: https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php (no registration needed).
Licence / terms: E-OBS data are free for non-commercial research and education; acknowledge "the E-OBS
dataset from the EU-FP6 project UERRA (https://www.uerra.eu) and the Copernicus Climate Change Service, and the
data providers in the ECA&D project (https://www.ecad.eu)" and cite Cornes et al. (2018) JGR-Atmos 123:9391-9409.
Output: ROOT/data/raw/eobs/eobs_{rr,tg}_donana.nc (box 36.7-37.3 N, 6.7-6.0 W).

Usage: python -I scripts/03_fetch_eobs.py <ROOT>
"""
import sys, os
import numpy as np, pandas as pd, xarray as xr, fsspec, h5py

ROOT = sys.argv[1]
OUT = os.path.join(ROOT, "data", "raw", "eobs")
os.makedirs(OUT, exist_ok=True)
BASE = "https://knmi-ecad-assets-prd.s3.amazonaws.com/ensembles/data/Grid_0.1deg_reg_ensemble/"
for var in ("rr", "tg"):
    fn = os.path.join(OUT, f"eobs_{var}_donana.nc")
    if os.path.exists(fn):
        continue
    f = fsspec.open(BASE + f"{var}_ens_mean_0.1deg_reg_v31.0e.nc", block_size=64 * 2**20, cache_type="readahead").open()
    h = h5py.File(f, "r")
    lat, lon = h["latitude"][:], h["longitude"][:]
    t = pd.Timestamp("1950-01-01") + pd.to_timedelta(h["time"][:], "D")
    iy = np.where((lat > 36.7) & (lat < 37.3))[0]
    ix = np.where((lon > -6.7) & (lon < -6.0))[0]
    it = np.where((t >= "1989-09-01") & (t <= "2024-12-31"))[0]
    d = h[var]
    raw = d[it[0]:it[-1] + 1, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1].astype("float32")
    a1 = lambda k, dflt: float(np.ravel(d.attrs.get(k, dflt))[0])
    fill = a1("_FillValue", -9999)
    raw[raw == fill] = np.nan
    val = raw * a1("scale_factor", 1.0) + a1("add_offset", 0.0)
    xr.Dataset({var: (("time", "latitude", "longitude"), val)},
               coords={"time": t[it[0]:it[-1] + 1], "latitude": lat[iy[0]:iy[-1] + 1],
                       "longitude": lon[ix[0]:ix[-1] + 1]},
               attrs={"source": BASE + f"{var}_ens_mean_0.1deg_reg_v31.0e.nc", "units": np.ravel(d.attrs["units"])[0].decode() if isinstance(np.ravel(d.attrs["units"])[0], bytes) else str(np.ravel(d.attrs["units"])[0])}).to_netcdf(fn)
    print("OK", var, val.shape, flush=True)
