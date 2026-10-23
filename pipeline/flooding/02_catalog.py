"""02_catalog.py — Sentinel-2 L2A items for tiles 29SQA/29SQB, 1 Sep 2017 – 31 Aug 2026 (Earth Search v1).

Prefers sentinel-2-c1-l2a (Collection-1 reprocessing, harmonised offsets, has CLD_20m cloud probability);
falls back to sentinel-2-l2a where c1 has no item for that date/tile (mostly 2017–2018 and 2022).
Within sentinel-2-l2a, duplicate items for the same date/tile keep the highest processing version.
Output: data/catalog_s2.csv
"""
import sys
import warnings
from pathlib import Path

import pandas as pd
from pystac_client import Client

warnings.filterwarnings("ignore")
ROOT = Path(sys.argv[1])
YEARS = [int(a) for a in sys.argv[2].split(",")] if len(sys.argv) > 2 else list(range(2017, 2027))
BANDS = {"scl": "scl", "cloud": "cloud", "B02": "blue", "B03": "green", "B04": "red",
         "B05": "rededge1", "B07": "rededge3", "B11": "swir16", "B12": "swir22"}
c = Client.open("https://earth-search.aws.element84.com/v1")
rows = []
for col in ("sentinel-2-c1-l2a", "sentinel-2-l2a"):
    for y in YEARS:
        s = c.search(collections=[col], bbox=[-6.55, 36.85, -6.30, 37.10],
                     datetime=(f"{y}-04-01/{y}-08-31" if y == 2017 else f"{y}-01-01/{y}-12-31"), max_items=5000)
        for it in s.items():
            p = it.properties
            tile = (p.get("grid:code") or p.get("s2:mgrs_tile") or "").replace("MGRS-", "")
            if tile not in ("29SQA", "29SQB"):
                continue
            r = {"collection": col, "id": it.id, "tile": tile, "datetime": p["datetime"],
                 "date": p["datetime"][:10], "cloud_cover": p.get("eo:cloud_cover"),
                 "baseline": p.get("s2:processing_baseline"), "offset_applied": p.get("earthsearch:boa_offset_applied"),
                 "platform": p.get("platform"), "rel_orbit": p.get("sat:relative_orbit")}
            for k, a in BANDS.items():
                r[k] = it.assets[a].href if a in it.assets else ""
            if col == "sentinel-2-c1-l2a":
                rb = it.assets["swir16"].extra_fields.get("raster:bands", [{}])[0]
                r["scale"], r["offset"] = rb.get("scale", 1e-4), rb.get("offset", 0.0)
            else:
                needs = (str(p.get("s2:processing_baseline", "0")) >= "04.00") and not p.get("earthsearch:boa_offset_applied", False)
                r["scale"], r["offset"] = 1e-4, (-0.1 if needs else 0.0)
            rows.append(r)
        print(col, y, len(rows), flush=True)
df = pd.DataFrame(rows)
df["version"] = df["id"].str.extract(r"_(\d)_L2A$").fillna("9").astype(int)
df["pref"] = (df["collection"] == "sentinel-2-c1-l2a").astype(int)
df = df.sort_values(["date", "tile", "pref", "version"]).drop_duplicates(["date", "tile"], keep="last")
df = df[(df.date >= "2017-04-01") & (df.date <= "2026-08-31")]
df["hydro_year"] = df["date"].apply(lambda d: f"{int(d[:4]) - (int(d[5:7]) < 9)}/{str(int(d[:4]) - (int(d[5:7]) < 9) + 1)[2:]}")
df = df.drop(columns=["pref"])
out = ROOT / "data/catalog_s2.csv"
if len(sys.argv) > 2 and out.exists():
    old = pd.read_csv(out); df = pd.concat([old, df[~df.date.isin(old.date)]]).sort_values(["date", "tile"])
df.to_csv(out, index=False)
print(df.groupby(["hydro_year", "collection"]).size().unstack(fill_value=0))
print("dates:", df.date.nunique())
