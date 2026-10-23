"""Fetch hourly point time series from the DestinE Climate Change Adaptation Digital Twin (Climate DT),
generation 2, at one HEALPix level-10 cell over Donana (Natura 2000 ES0000024).

Source: Destination Earth Climate DT, generation 2 (1990-2049, ~5 km, hourly, HEALPix nside 1024),
DOI 10.21957/79c6af3105, served by Polytope:
  IFS-NEMO  -> https://polytope.mn5.apps.dte.destination-earth.eu
  IFS-FESOM -> https://polytope.lumi.apps.dte.destination-earth.eu
Access: DestinE Platform account with upgraded access; credential read by polytope-client from
~/.polytopeapirc (never copied or printed here).
Licence: DestinE Terms and Conditions v2.0 (https://platform.destine.eu/terms-and-conditions/).
Raw extracts are NOT redistributable: they stay in ROOT/data/raw/dt/ (git-ignored). Only derived statistics
are published, with the attribution "This data is created based on data of the European Union, using the
Destination Earth Platform, but has been modified."

Request keys follow destination-earth-digital-twins/polytope-examples climate-dt/ (explorer/polytope_zarr.py,
explorer/destine_portfolio.py, commit 3a13f05):
  hist 1990-2014: activity "baseline", experiment "hist"; 2015+: activity "projections", experiment "SSP3-7.0";
  stream "clte" (hourly), levtype "sfc", resolution "high"; "time" lists all 24 hours (a /to/by range returns
  only the first hour); server-side "timeseries" feature, ONE point per request, <=10-year chunks.
Params: avg_tprate (hourly-mean precipitation rate, kg m-2 s-1), 2t (K), avg_ie (mean moisture-flux /
evaporation rate, kg m-2 s-1, negative = upward). Checked live 2026-10-22 (avg_pevr, mx2t, mn2t refused:
"The datacube does not contain a date axis").
Points are spherical HEALPix (nest) level-10 cell centres (the DT grid as delivered), NOT WGS84 centres;
see CELLS below (computed with healpix_geo.nested.healpix_to_lonlat(..., ellipsoid="sphere")). No regridding.

Usage: python -I scripts/01_fetch_dt.py <ROOT> <MODEL> <PARAM> <CELL>
"""
import sys, os, json, time
from polytope.api import Client

ROOT, MODEL, PARAM, CELL = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
# cell id -> (name, lat, lon) ; spherical HEALPix nested level 10 centres
CELLS = {3490355: ("marsh_centre", 36.9352, -6.3281), 3490408: ("marsh_north", 37.1221, -6.2402),
         3490451: ("aquifer_rocio", 37.1221, -6.5039), 3490362: ("reserve_dunes", 36.9819, -6.4600),
         3490357: ("marsh_east", 36.9819, -6.1963)}
CHUNKS = [("hist", 1990, 1999), ("hist", 2000, 2009), ("hist", 2010, 2014), ("SSP3-7.0", 2015, 2024)]
address = ("polytope.mn5.apps.dte.destination-earth.eu" if MODEL == "IFS-NEMO"
           else "polytope.lumi.apps.dte.destination-earth.eu")
OUT = os.path.join(ROOT, "data", "raw", "dt", PARAM)
os.makedirs(OUT, exist_ok=True)
_, la, lo = CELLS[CELL]
client = Client(address=address, quiet=True)
for exp, y0, y1 in CHUNKS:
    fn = os.path.join(OUT, f"{MODEL}_{CELL}_{exp}_{y0}_{y1}.covjson")
    if os.path.exists(fn) and os.path.getsize(fn) > 10000:
        continue
    req = {"class": "d1", "dataset": "climate-dt", "type": "fc", "expver": "0001", "generation": "2",
           "realization": "1", "activity": "baseline" if exp == "hist" else "projections",
           "experiment": exp, "levtype": "sfc", "resolution": "high", "stream": "clte", "model": MODEL,
           "param": PARAM, "date": f"{y0}0101/to/{y1}1231", "time": "/".join(f"{h:02d}00" for h in range(24)),
           "feature": {"type": "timeseries", "points": [[la, lo]], "axes": ["latitude", "longitude"],
                       "time_axis": "date"}}
    json.dump({"address": address, "request": req}, open(fn.replace(".covjson", ".request.json"), "w"))
    for attempt in range(3):
        t0 = time.time()
        try:
            client.retrieve("destination-earth", req, fn)
            print("OK", MODEL, PARAM, CELL, exp, y0, y1, f"{time.time() - t0:.0f}s", flush=True)
            break
        except Exception as e:
            msg = str(e).split("Details:")[-1][:300].replace("\n", " ")
            print("FAIL", MODEL, PARAM, CELL, exp, y0, y1, attempt, type(e).__name__, msg, flush=True)
            time.sleep(30)
