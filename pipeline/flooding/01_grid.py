"""01_grid.py — HEALPix depth-19 domain for Doñana and exact area-weighted UTM->HEALPix operators.

Domain: cells (healpix-geo, WGS84 ellipsoid, NESTED, depth 19) whose centre lies inside
ES0000024 ∪ ES6150009 (EEA Natura 2000 release, data/raw/n2k), restricted to the part covered
by Sentinel-2 tiles 29SQA/29SQB (UTM 29N, x >= 699960 m; the 2.2 % west strip in 29SPB is excluded).

Weights: for every domain cell, the cell boundary (4 vertices from healpix_geo.nested.vertices,
WGS84) is projected to UTM 29N and intersected exactly (shapely) with each candidate pixel square
of the S2 10 m and 20 m lattices (both tiles share one lattice). W[cell, pixel] = overlap area (m²).
A reflectance cell value is then sum_p W v_p / sum_p W over valid pixels (first-order conservative,
'covered' normalisation; no kernel, no PSF). Categorical layers use the class with the largest
summed overlap (majority), never an average.

Usage: .venv/bin/python -I scripts/01_grid.py <flooding_dir>
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path

import geopandas as gpd
import healpix_geo.nested as hn
import numpy as np
import pandas as pd
import shapely
import zarr
from pyproj import Transformer
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parent))
from common import weights  # exact overlap operator

ROOT = Path(sys.argv[1])
OUT = ROOT / "data" / "grid"
OUT.mkdir(parents=True, exist_ok=True)
DEPTH = 19
X0_TILE = 699960.0  # west edge of 29SQA/29SQB
TO_UTM = Transformer.from_crs(4326, 32629, always_xy=True)


def domain_polygon():
    g = pd.concat([gpd.read_file(ROOT / f"data/raw/n2k/{c}.geojson") for c in ("ES0000024", "ES6150009")])
    u = shapely.union_all(shapely.make_valid(g.geometry.values))  # lon/lat
    gu = gpd.GeoSeries([u], crs=4326).to_crs(32629).iloc[0]
    clip = shapely.box(X0_TILE + 40, 0, 1e7, 1e8)
    du = gu.intersection(clip)
    return gpd.GeoSeries([du], crs=32629).to_crs(4326).iloc[0], gu.area, du.area


def cover(poly):
    parts = getattr(poly, "geoms", [poly])
    ids = []
    for p in parts:
        ext = np.asarray(p.exterior.coords)[:-1]
        c, _, _ = hn.polygon_coverage(ext, DEPTH, ellipsoid="WGS84", flat=True)
        ids.append(c)
    ids = np.unique(np.concatenate(ids))
    lon, lat = hn.healpix_to_lonlat(ids, DEPTH, ellipsoid="WGS84")
    lon = (lon + 180.0) % 360.0 - 180.0
    inside = shapely.contains_xy(poly, lon, lat)
    return ids[inside], lon[inside], lat[inside]


def main():
    t0 = time.time()
    poly, site_area, dom_area = domain_polygon()
    ids, lon, lat = cover(poly)
    base, fx, fy = hn.healpix_to_base_cell_coordinates(ids, DEPTH)
    assert np.all(base == base[0]), "domain spans more than one HEALPix base face"
    print(f"cells {len(ids):,} base face {base[0]}  ({time.time()-t0:.0f}s)")

    # exact cell boundaries (WGS84) in UTM 29N
    vlon, vlat = hn.vertices(ids, DEPTH, ellipsoid="WGS84")
    vlon = (vlon + 180.0) % 360.0 - 180.0
    qx, qy = TO_UTM.transform(vlon.ravel(), vlat.ravel())
    qx = qx.reshape(-1, 4); qy = qy.reshape(-1, 4)
    qarea = shapely.area(shapely.polygons(np.stack([qx, qy], -1)))

    # common UTM window (both tiles share the lattice; origin multiple of 20 m)
    x0 = X0_TILE + 20 * np.floor((qx.min() - X0_TILE) / 20)
    x1 = X0_TILE + 20 * np.ceil((qx.max() - X0_TILE) / 20)
    y1 = 4100040.0 + 20 * np.ceil((qy.max() - 4100040.0) / 20)
    y0 = 4100040.0 + 20 * np.floor((qy.min() - 4100040.0) / 20)
    win = {"crs": "EPSG:32629", "x0": x0, "x1": x1, "y0": y0, "y1": y1}
    store = zarr.open_group(OUT / "weights.zarr", mode="w")
    for res in (10, 20):
        nx = int(round((x1 - x0) / res)); ny = int(round((y1 - y0) / res))
        win[f"nx{res}"] = nx; win[f"ny{res}"] = ny
        t = time.time()
        W = weights(qx, qy, res, x0, y1, nx, ny)
        rs = np.asarray(W.sum(1)).ravel()
        rel = np.abs(rs / qarea - 1)
        print(f"W{res}: nnz {W.nnz:,}  max |rowsum/cellarea-1| {rel.max():.2e}  ({time.time()-t:.0f}s)")
        g = store.create_group(f"W{res}")
        for name, arr in (("data", W.data.astype("float32")), ("indices", W.indices.astype("int32")),
                          ("indptr", W.indptr.astype("int64"))):
            g.create_array(name, data=arr, chunks=(min(len(arr), 4_000_000),))
        g.attrs.update({"shape": list(W.shape), "units": "m2 overlap (UTM 29N plane)",
                        "row": "domain cell index (cells.zarr order)", "col": "row-major pixel index in window"})
    store.attrs.update({"window": win, "depth": DEPTH})

    cells = pd.DataFrame({"cell_ids": ids.astype("uint64"), "fx": fx.astype("int32"), "fy": fy.astype("int32"),
                          "lon": lon, "lat": lat, "area_utm_m2": qarea.astype("float32")})
    cg = zarr.open_group(OUT / "cells.zarr", mode="w")
    for k in cells.columns:
        cg.create_array(k, data=cells[k].values, chunks=(1_000_000,))
    gpd.GeoSeries([poly], crs=4326).to_file(OUT / "domain.geojson")
    info = {"depth": DEPTH, "n_cells": int(len(ids)), "base_face": int(base[0]),
            "cell_area_m2_mean_utm_plane": float(qarea.mean()),
            "site_union_area_ha": site_area / 1e4, "domain_area_ha": dom_area / 1e4,
            "excluded_west_strip_ha": (site_area - dom_area) / 1e4, "window": win}
    (OUT / "grid_info.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1), f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
