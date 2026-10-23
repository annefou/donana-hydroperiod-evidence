"""08_s1_catalog.py — Sentinel-1 RTC items (Microsoft Planetary Computer 'sentinel-1-rtc', gamma0 terrain-
corrected from ESA S1 IW GRD; Copernicus data, free, full and open; PC hosting terms apply) covering the
Doñana marsh, 1 Sep 2017 – 31 Aug 2026. Earth Search 'sentinel-1-grd' is requester-pays and was not used.
Output: data/catalog_s1.csv
"""
import sys, warnings
from pathlib import Path

import pandas as pd
from pystac_client import Client

warnings.filterwarnings("ignore")
ROOT = Path(sys.argv[1])
c = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
rows = []
for y in range(2017, 2027):
    for it in c.search(collections=["sentinel-1-rtc"], bbox=[-6.50, 36.90, -6.25, 37.10],
                       datetime=f"{y}-01-01/{y}-12-31").items():
        p = it.properties
        rows.append({"id": it.id, "datetime": p["datetime"], "date": p["datetime"][:10],
                     "orbit_state": p.get("sat:orbit_state"), "rel_orbit": p.get("sat:relative_orbit"),
                     "platform": p.get("platform"), "vv": it.assets["vv"].href, "vh": it.assets["vh"].href,
                     "bbox": ",".join(f"{v:.4f}" for v in it.bbox)})
    print(y, len(rows), flush=True)
df = pd.DataFrame(rows).drop_duplicates("id")
df = df[(df.date >= "2017-09-01") & (df.date <= "2026-08-31")].sort_values("datetime")
df["hydro_year"] = df["date"].apply(lambda d: f"{int(d[:4]) - (int(d[5:7]) < 9)}/{str(int(d[:4]) - (int(d[5:7]) < 9) + 1)[2:]}")
df.to_csv(ROOT / "data/catalog_s1.csv", index=False)
print(df.groupby(["hydro_year", "rel_orbit"]).size().unstack(fill_value=0))
