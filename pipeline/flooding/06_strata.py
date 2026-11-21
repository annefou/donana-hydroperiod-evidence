"""06_strata.py — evaluation strata on HEALPix d19 cells.

sea_estuary : cells flooded in >= 95 % of LAST masks 1984–2019 forming the largest connected component (base-face lattice,
              4-connectivity) — the Atlantic strip and the Guadalquivir channel.
              Excluded from accuracy and hydroperiod statistics (Kordelas 2018 §3.1 excludes coastal sea).
lagoon      : OpenStreetMap natural=water polygons named "Laguna…" or "Charco…" (dune ponds and lagoons;
              © OpenStreetMap contributors, ODbL 1.0), polygon -> cells with healpix-geo polygon_coverage
              (WGS84) + cell-centre-inside test.
artificial  : OSM natural=water with water in {basin, reservoir, wastewater, canal} or names
              "Piscifactoría/Balsa/Embalse/Colector" (fish farms, rafts, reservoirs).
marsh       : remaining cells flooded in >= 25 % of the LAST Dec–Mar masks 1984–2019 (seasonal marsh envelope).
other       : everything else in the domain (dunes, pine/shrub, farmland, marsh edges).
Also stores LAST flood frequencies (all masks; Dec–Mar masks) per cell.
"""
import json, sys
from pathlib import Path

import healpix_geo.nested as hn
import numpy as np
import shapely
import zarr
from scipy import ndimage
from shapely.geometry import Polygon

sys.path.insert(0, str(Path(__file__).parent))
from common import DEPTH, face_index, jdump

ROOT = Path(sys.argv[1])
c = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")
cells = {k: c[k][:] for k in c.array_keys()}
ids = cells["cell_ids"]; n = len(ids)
m = zarr.open_group(ROOT / "data/ref/last_masks.zarr", mode="r")
dates = list(m.attrs["dates"])
months = np.array([int(d[4:6]) for d in dates])
cnt = np.zeros(n, np.int32); cntw = np.zeros(n, np.int32); nv = np.zeros(n, np.int32); nvw = np.zeros(n, np.int32)
for i in range(len(dates)):
    a = m["masks"][i]
    ok = a != 255
    cnt += (a == 1); nv += ok
    if months[i] in (12, 1, 2, 3):
        cntw += (a == 1); nvw += ok
freq = cnt / np.maximum(nv, 1); freqw = cntw / np.maximum(nvw, 1)

ri, ci, shape = face_index(cells)
perm = np.zeros(shape, bool); perm[ri, ci] = freq >= 0.95
lab, _ = ndimage.label(perm)
sizes = np.bincount(lab[ri, ci])
sizes[0] = 0
sea_label = int(np.argmax(sizes))  # largest permanent-water component = Atlantic strip + Guadalquivir channel
sea = lab[ri, ci] == sea_label

osm = json.load(open(ROOT / "data/raw/osm/water.json"))["elements"]


def rings(el):
    if el["type"] == "way":
        g = el.get("geometry", [])
        return [Polygon([(p["lon"], p["lat"]) for p in g])] if len(g) >= 4 else []
    out = []
    for mem in el.get("members", []):
        if mem.get("role") == "outer" and len(mem.get("geometry", [])) >= 4:
            out.append(Polygon([(p["lon"], p["lat"]) for p in mem["geometry"]]))
    return out


def to_cells(polys, ids_out=None):
    sel = np.zeros(n, bool)
    lon = cells["lon"]; lat = cells["lat"]
    order = np.argsort(ids)
    for p in polys:
        if ids_out is not None:  # one id per input polygon, even if it yields no cell
            ids_out[1] += 1
            ids_out[2].append(p)
        if not p.is_valid:
            p = p.buffer(0)
        if p.geom_type == "MultiPolygon":
            p = max(p.geoms, key=lambda q: q.area)
        if p.is_empty or p.geom_type != "Polygon":
            continue
        cid, _, _ = hn.polygon_coverage(np.asarray(p.exterior.coords)[:-1], DEPTH, ellipsoid="WGS84", flat=True)
        pos = np.searchsorted(ids[order], cid)
        pos = pos[(pos < n)]
        hit = order[pos][ids[order][pos] == cid[:len(pos)]] if len(pos) == len(cid) else order[pos][np.isin(ids[order][pos], cid)]
        if len(hit) == 0:
            continue
        inside = shapely.contains_xy(p, lon[hit], lat[hit])
        sel[hit[inside]] = True
        if ids_out is not None:
            ids_out[0][hit[inside]] = ids_out[1] - 1
    return sel


lag_p, art_p, lag_names = [], [], []
for el in osm:
    t = el.get("tags", {}); name = t.get("name", "")
    if name.startswith(("Laguna", "Charco")):
        lag_p += rings(el); lag_names.append(name)
    elif t.get("water") in ("basin", "reservoir", "wastewater", "canal") or name.startswith(("Piscifactor", "Balsa", "Embalse", "Colector")):
        art_p += rings(el)
lag_id = np.zeros(n, np.int32)
lag_info = [lag_id, 1, []]
lag_named = []
for el in osm:
    t = el.get("tags", {}); name = t.get("name", "")
    if name.startswith(("Laguna", "Charco")):
        for p in rings(el):
            lag_named.append(name)
lagoon = to_cells(lag_p, lag_info) & ~sea
lag_id[~lagoon] = 0
artificial = to_cells(art_p) & ~sea & ~lagoon
marsh = (freqw >= 0.25) & ~sea & ~lagoon & ~artificial
stratum = np.full(n, 4, np.uint8)  # 4 other
stratum[marsh] = 1; stratum[lagoon] = 2; stratum[artificial] = 3; stratum[sea] = 0
g = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="w")
g.create_array("stratum", data=stratum, chunks=(1 << 20,))
g.create_array("lagoon_id", data=lag_id, chunks=(1 << 20,))
from pyproj import Geod
geod = Geod(ellps="WGS84")
g.attrs["lagoons"] = [{"id": k + 1, "name": nm, "polygon_area_ha": abs(geod.geometry_area_perimeter(pp)[0]) / 1e4,
                       "n_cells": int((lag_id == k + 1).sum())} for k, (nm, pp) in enumerate(zip(lag_named, lag_info[2]))]
g.create_array("last_freq_all", data=freq.astype("float32"), chunks=(1 << 20,))
g.create_array("last_freq_decmar", data=freqw.astype("float32"), chunks=(1 << 20,))
g.attrs.update({"stratum": {"0": "sea_estuary", "1": "marsh", "2": "lagoon", "3": "artificial", "4": "other"}})
A = cells["area_utm_m2"] / 1e4
res = {nm: round(float(A[stratum == v].sum()), 1) for v, nm in enumerate(["sea_estuary", "marsh", "lagoon", "artificial", "other"])}
res["n_lagoon_polygons"] = len(lag_names); res["lagoon_names"] = sorted(set(lag_names))
res["n_last_masks"] = len(dates)
jdump(res, ROOT / "data/grid/strata_info.json")
jdump(res, ROOT / "results_strata.json")  # stratum areas (ha), published to results/flooding/
print(res)
