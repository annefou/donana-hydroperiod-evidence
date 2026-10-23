"""Shared helpers: grid loading, conservative / majority resampling onto HEALPix, I/O."""
import json
from pathlib import Path

import numpy as np
import zarr
from scipy import sparse
from multiprocessing import Pool
import shapely

DEPTH = 19
CLEAR_SCL = (2, 4, 5, 6, 7)        # dark area, vegetation, bare, water, unclassified
CLOUD_SCL = (3, 8, 9, 10)          # cloud shadow, cloud medium, cloud high, cirrus


def load_grid(root: Path):
    g = zarr.open_group(Path(root) / "data/grid/weights.zarr", mode="r")
    win = g.attrs["window"]
    W = {}
    for res in (10, 20):
        sub = g[f"W{res}"]
        W[res] = sparse.csr_matrix((sub["data"][:], sub["indices"][:], sub["indptr"][:]),
                                   shape=tuple(sub.attrs["shape"]))
    c = zarr.open_group(Path(root) / "data/grid/cells.zarr", mode="r")
    cells = {k: c[k][:] for k in c.array_keys()}
    return W, win, cells


def conservative(W, img, valid):
    """First-order conservative (area-weighted) value per cell from a pixel image on the window lattice.

    Returns (value, covered_fraction): value = sum W*v*m / sum W*m; fraction = sum W*m / sum W.
    """
    m = valid.ravel().astype(np.float32)
    v = np.where(valid, img, 0).ravel().astype(np.float32)
    num = W @ (v * m)
    den = W @ m
    tot = np.asarray(W.sum(1)).ravel() if not hasattr(W, "_rowsum") else W._rowsum
    with np.errstate(invalid="ignore", divide="ignore"):
        val = np.where(den > 0, num / den, np.nan).astype(np.float32)
    return val, (den / tot).astype(np.float32)


def class_fractions(W, cls_img, classes, valid):
    """Overlap-area fraction of each class per cell (used for majority / purity; never averaged)."""
    tot = W @ valid.ravel().astype(np.float32)
    out = {}
    for c in classes:
        m = ((cls_img == c) & valid).ravel().astype(np.float32)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[c] = np.where(tot > 0, (W @ m) / tot, 0).astype(np.float32)
    return out, tot


def majority(W, cls_img, classes, valid):
    fr, tot = class_fractions(W, cls_img, classes, valid)
    stack = np.stack([fr[c] for c in classes])
    maj = np.asarray(classes, dtype=np.uint8)[np.argmax(stack, 0)]
    maj[tot <= 0] = 0
    return maj, fr


def face_raster(cells, values, fill=np.nan, dtype=np.float32):
    """Place per-cell values on the HEALPix base-face (fx, fy) lattice (native HEALPix topology)."""
    fx, fy = cells["fx"], cells["fy"]
    x0, y0 = fx.min(), fy.min()
    a = np.full((fy.max() - y0 + 1, fx.max() - x0 + 1), fill, dtype=dtype)
    a[fy - y0, fx - x0] = values
    return a


def face_index(cells):
    fx, fy = cells["fx"], cells["fy"]
    return fy - fy.min(), fx - fx.min(), (fy.max() - fy.min() + 1, fx.max() - fx.min() + 1)


def dggs_attrs(level=DEPTH):
    return {
        "zarr_conventions": [{
            "uuid": "7b255807-140c-42ca-97f6-7a1cfecdbc38", "name": "dggs",
            "schema_url": "https://raw.githubusercontent.com/zarr-conventions/dggs/refs/tags/v1/schema.json",
            "spec_url": "https://github.com/zarr-conventions/dggs/blob/v1/README.md",
            "description": "Discrete Global Grid Systems convention for zarr"}],
        "dggs": {"name": "healpix", "refinement_level": level, "indexing_scheme": "nested",
                 "ellipsoid": {"name": "wgs84", "semi_major_axis": 6378137.0,
                               "inverse_flattening": 298.257223563},
                 "spatial_dimension": "cells", "coordinate": "cell_ids", "compression": "none"}}


def jdump(obj, path):
    Path(path).write_text(json.dumps(obj, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


# ---- exact cell/pixel overlap operator (used by 01_grid.py and s1scene.py) ----
def _overlap_chunk(args):
    qx, qy, res, x0, y1, nx, ny = args
    n = qx.shape[0]
    ix0 = np.floor((qx.min(1) - x0) / res).astype(np.int64)
    ix1 = np.floor((qx.max(1) - x0) / res).astype(np.int64)
    iy0 = np.floor((y1 - qy.max(1)) / res).astype(np.int64)  # row index from top
    iy1 = np.floor((y1 - qy.min(1)) / res).astype(np.int64)
    span = int(max((ix1 - ix0).max(), (iy1 - iy0).max())) + 1
    quads = shapely.polygons(np.stack([qx, qy], -1))
    rows, cols, vals = [], [], []
    for dx in range(span):
        for dy in range(span):
            ix = ix0 + dx
            iy = iy0 + dy
            ok = (ix <= ix1) & (iy <= iy1) & (ix >= 0) & (iy >= 0) & (ix < nx) & (iy < ny)
            if not ok.any():
                continue
            k = np.nonzero(ok)[0]
            bx = x0 + ix[k] * res
            by = y1 - (iy[k] + 1) * res
            a = shapely.area(shapely.intersection(quads[k], shapely.box(bx, by, bx + res, by + res)))
            m = a > 1e-9
            rows.append(k[m]); cols.append(iy[k][m] * nx + ix[k][m]); vals.append(a[m])
    return np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)


def weights(qx, qy, res, x0, y1, nx, ny, nproc=16, chunk=200_000):
    jobs = [(qx[i:i + chunk], qy[i:i + chunk], res, x0, y1, nx, ny) for i in range(0, len(qx), chunk)]
    with Pool(nproc) as pool:
        out = pool.map(_overlap_chunk, jobs)
    R, C, V = [], [], []
    for j, (r, c, v) in enumerate(out):
        R.append(r + j * chunk); C.append(c); V.append(v)
    W = sparse.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(len(qx), nx * ny))
    W.sum_duplicates()
    return W


