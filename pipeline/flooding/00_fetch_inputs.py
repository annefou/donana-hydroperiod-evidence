"""00_fetch_inputs.py <root> [--only name,...] — download every non-satellite input (satellite scenes are
read on the fly from cloud-optimised GeoTIFFs; see 02_catalog.py / 02b_catalog_pc.py / 08_s1_catalog.py).

Sources, licences (checked 22 Oct 2026) and destinations:
 n2k      EEA Natura 2000 sites ES0000024, ES6150009 (ArcGIS REST, Natura2000Sites MapServer layer 2)
          https://bio.discomap.eea.europa.eu/arcgis/rest/services/ProtectedSites/Natura2000Sites/MapServer/2/query
          EEA standard re-use policy (CC BY 4.0)                                 -> data/raw/n2k/<code>.geojson
 masks    LAST-EBD "Flood Masks Doñana 1984/2019", Díaz-Delgado et al., doi:10.5281/zenodo.3518182, CC BY 4.0
          https://zenodo.org/api/records/3518182/files/inunda.rar/content (122.7 MB)  -> data/raw/zenodo/inunda.rar
 hydro    LAST-EBD "Hydroperiod Doñana 1974/2019", doi:10.5281/zenodo.3516929, CC BY 4.0
          https://zenodo.org/api/records/3516929/files/hidroperiodos.zip/content (40.6 MB) -> data/raw/zenodo/
 field    Díaz-Delgado (2025) ground-truth points 24 Apr 2017, doi:10.5281/zenodo.17158473, CC BY 4.0 -> data/raw/zenodo/
 osm      OpenStreetMap natural=water in the site bbox via Overpass API, © OSM contributors, ODbL 1.0
          https://overpass-api.de/api/interpreter                                -> data/raw/osm/water.json
 icts     ICTS-Doñana Hidromet API, daily marsh water level, 5 stations, CC BY 4.0 (API metadata and site footer)
          https://datos-automaticos.icts-donana.es/api/v1/                         -> data/raw/icts/*.csv
 chg      CHG (Confederación Hidrográfica del Guadalquivir) piezometers: WFS station list + IDE API series
          https://idechg.chguadalquivir.es/geoserver/ggiscloud_root/ows (WFS piezometros)
          https://idechg.chguadalquivir.es/nodo-areas/api/redes-control/getDownload/PZ (POST)
          No explicit licence; CHG aviso legal: reproduction authorised citing the source; data "not official".
          Stations within 10 km of 36.98 N, 6.47 W (Santa Olalla / Dulce)         -> data/raw/chg/*.csv
All requests retry 6 times with exponential back-off (EEA/EMODnet endpoints returned 502s on 22 Oct 2026).
"""
import csv, json, math, sys, time
from pathlib import Path

import requests
from pyproj import Transformer
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(sys.argv[1])
ONLY = sys.argv[sys.argv.index("--only") + 1].split(",") if "--only" in sys.argv else None
RAW = ROOT / "data/raw"
S = requests.Session()
S.mount("https://", HTTPAdapter(max_retries=Retry(total=6, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504),
                                                  allowed_methods=None)))
S.headers["User-Agent"] = "lifewatch-donana-flooding/0.1 (research; reproducible pipeline)"


def dl(url, dest, **kw):
    dest = Path(dest); dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    with S.get(url, stream=True, timeout=600, **kw) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for b in r.iter_content(1 << 20):
                f.write(b)
        tmp.rename(dest)
    print("ok", dest, dest.stat().st_size)
    return dest


def n2k():
    q = "https://bio.discomap.eea.europa.eu/arcgis/rest/services/ProtectedSites/Natura2000Sites/MapServer/2/query"
    for c in ("ES0000024", "ES6150009"):
        dl(q, RAW / f"n2k/{c}.geojson", params={"where": f"SITECODE='{c}'", "returnGeometry": "true", "outSR": "4326",
                                               "outFields": "SITECODE,SITENAME,RELEASE_DATE", "f": "geojson"})


def zenodo():
    dl("https://zenodo.org/api/records/3518182/files/inunda.rar/content", RAW / "zenodo/inunda.rar")
    dl("https://zenodo.org/api/records/3516929/files/hidroperiodos.zip/content", RAW / "zenodo/hidroperiodos.zip")
    rec = S.get("https://zenodo.org/api/records/17158473", timeout=60).json()
    dl(rec["files"][0]["links"]["self"], RAW / "zenodo/groundtruth_20170424.xlsx")


def osm():
    q = ('[out:json][timeout:120];(way["natural"="water"](36.78,-6.80,37.38,-6.09);'
         'relation["natural"="water"](36.78,-6.80,37.38,-6.09););out geom;')
    dest = RAW / "osm/water.json"
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        r = S.post("https://overpass-api.de/api/interpreter", data={"data": q}, timeout=300); r.raise_for_status()
        dest.write_bytes(r.content); print("ok", dest)


def icts():
    H = "https://datos-automaticos.icts-donana.es/api/v1"
    out = RAW / "icts"; out.mkdir(parents=True, exist_ok=True)
    st = {s["station_id"]: s for s in S.get(H + "/catalog/stations/", timeout=60).json()["stations"]}
    json.dump(st, open(out / "stations.json", "w"), indent=1)
    for sid in (8437, 8466, 1072, 8433, 8436):
        for ins in S.get(f"{H}/catalog/instruments/?station={sid}&variable=waterLevel", timeout=60).json()["instruments"]:
            q = S.get(f"{H}/catalog/query-params/?station={sid}&variable=waterLevel&instrument={ins['instrument']}",
                      timeout=60).json()
            for comb in q["combinations"]:
                url, page, rows = comb["example_urls"]["daily"], 1, []
                while True:
                    d = S.get(url + f"&start_date=2015-01-01&end_date=2026-10-22&page={page}&page_size=5000", timeout=120).json()
                    rows += d.get("results") or d.get("data") or []
                    if not d.get("next"):
                        break
                    page += 1
                fn = out / f"waterLevel_daily_{sid}_{st[sid]['name'].replace(' ', '_')}_{ins['instrument']}.csv"
                if rows:
                    with open(fn, "w", newline="") as f:
                        w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)
                print("ok", fn, len(rows))


def chg():
    out = RAW / "chg"; out.mkdir(parents=True, exist_ok=True)
    wfs = dl("https://idechg.chguadalquivir.es/geoserver/ggiscloud_root/ows", out / "piezometros_wfs.geojson",
             params={"service": "WFS", "version": "1.0.0", "request": "GetFeature",
                     "typeName": "ggiscloud_root:piezometros", "outputFormat": "application/json"})
    A = "https://idechg.chguadalquivir.es/nodo-areas/api/redes-control"
    to_ll = Transformer.from_crs(25830, 4326, always_xy=True)
    X, Y = Transformer.from_crs(4326, 25830, always_xy=True).transform(-6.47, 36.98)
    rows = []
    for f in json.load(open(wfs))["features"]:
        x, y = f["geometry"]["coordinates"][:2]
        d = math.hypot(x - X, y - Y) / 1000
        if d > 10:
            continue
        p = f["properties"]; lon, lat = to_ll.transform(x, y)
        r = S.post(A + "/getDownload/PZ", json={"estacion": p["cod_estsub"], "fechaIni": "01/01/1970",
                                                "fechaFin": "22/10/2026"}, timeout=120)
        r.raise_for_status()
        for v in r.json():
            dd, mm, yy = v["fecha"].split("/")
            rows.append({"cod_estsub": p["cod_estsub"], "nom_estsub": p["nom_estsub"], "masa_subterranea": p.get("massub"),
                         "lat": round(lat, 6), "lon": round(lon, 6), "cota_boca_msnm": p.get("cota_estsub"),
                         "dist_km": round(d, 2), "date": f"{yy}-{mm}-{dd}", "parameter": v["nomParam"],
                         "unit": v["unidadParam"], "value": v["valor"]})
        time.sleep(0.5)
    with open(out / "piezometria_within10km_santa_olalla.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0])); w.writeheader(); w.writerows(rows)
    print("ok chg", len(rows))


for name, fn in (("n2k", n2k), ("zenodo", zenodo), ("osm", osm), ("icts", icts), ("chg", chg)):
    if ONLY is None or name in ONLY:
        fn()
