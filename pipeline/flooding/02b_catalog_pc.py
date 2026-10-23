"""02b_catalog_pc.py — fill date/tile gaps of the Earth Search catalogue with Microsoft Planetary Computer
'sentinel-2-l2a' items (same ESA L2A product family; no cloud-probability asset, SCL only).
Also replaces Earth Search rows whose assets are JP2 in the requester-pays bucket.
Earth Search lacks e.g. 15 Dec 2017 – 27 Mar 2018 and every other 2017 date for tiles 29SQA/29SQB.
"""
import sys, warnings
from pathlib import Path

import pandas as pd
from pystac_client import Client

warnings.filterwarnings("ignore")
ROOT = Path(sys.argv[1])
out = ROOT / "data/catalog_s2.csv"
old = pd.read_csv(out, dtype=str, keep_default_na=False)
BANDS = ["scl", "cloud", "B02", "B03", "B04", "B05", "B07", "B11", "B12"]
# rows pointing at JP2 files in the requester-pays bucket (s3://sentinel-s2-l2a) are treated as missing
bad = old[BANDS].apply(lambda r: any(x.startswith("s3://") or x.endswith(".jp2") for x in r), axis=1)
print("replacing requester-pays rows:", list(zip(old.date[bad], old.tile[bad])))
old = old[~bad]
have = set(zip(old.date, old.tile))
c = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
rows = []
for y in range(2017, 2027):
    rng = f"{y}-04-01/{y}-12-31" if y == 2017 else (f"{y}-01-01/{y}-08-31" if y == 2026 else f"{y}-01-01/{y}-12-31")
    for it in c.search(collections=["sentinel-2-l2a"], bbox=[-6.55, 36.85, -6.30, 37.10], datetime=rng).items():
        p = it.properties
        tile = p.get("s2:mgrs_tile", "")
        if tile not in ("29SQA", "29SQB"):
            continue
        d = p["datetime"][:10]
        if (d, tile) in have:
            continue
        base = str(p.get("s2:processing_baseline", "0"))
        r = {"collection": "pc:sentinel-2-l2a", "id": it.id, "tile": tile, "datetime": p["datetime"], "date": d,
             "cloud_cover": p.get("eo:cloud_cover"), "baseline": base, "offset_applied": False,
             "platform": p.get("platform"), "rel_orbit": p.get("sat:relative_orbit"), "scale": 1e-4,
             "offset": -0.1 if base >= "04.00" else 0.0, "cloud": ""}
        for k in ("B02", "B03", "B04", "B05", "B07", "B11", "B12"):
            r[k] = it.assets[k].href
        r["scl"] = it.assets["SCL"].href
        rows.append(r); have.add((d, tile))
    print(y, len(rows), flush=True)
new = pd.DataFrame(rows)
new["version"] = 9
new["hydro_year"] = new["date"].apply(lambda d: f"{int(d[:4]) - (int(d[5:7]) < 9)}/{str(int(d[:4]) - (int(d[5:7]) < 9) + 1)[2:]}")
df = pd.concat([old, new]).sort_values(["date", "tile"])
df.to_csv(out, index=False)
print(df.groupby(["hydro_year", "collection"]).size().unstack(fill_value=0)); print("dates", df.date.nunique())
