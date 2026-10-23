"""Sentinel-1 RTC (gamma0, linear power, 10 m; Microsoft Planetary Computer 'sentinel-1-rtc', derived from
ESA S1 GRD) onto HEALPix d19 cells. Linear gamma0 is resampled area-weighted (exact cell/pixel overlap,
same operator as for S2) and converted to dB afterwards. Scenes come in EPSG:32629 or EPSG:32630; the
exact overlap operator is built once per CRS lattice (10 m pixel edges at multiples of 10 m).
"""
import importlib.util
import os
from pathlib import Path

import healpix_geo.nested as hn
import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import Window

os.environ.update(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="3")
import sys
sys.path.insert(0, str(Path(__file__).parent))
from common import weights
_CACHE = {}


def operator(epsg, cells, W10=None, win=None):
    if epsg == 32629 and W10 is not None:
        return W10, {"x0": win["x0"], "y1": win["y1"], "nx": win["nx10"], "ny": win["ny10"]}
    if epsg in _CACHE:
        return _CACHE[epsg]
    vlon, vlat = hn.vertices(cells["cell_ids"], 19, ellipsoid="WGS84")
    vlon = (vlon + 180.0) % 360.0 - 180.0
    qx, qy = Transformer.from_crs(4326, epsg, always_xy=True).transform(vlon.ravel(), vlat.ravel())
    qx = qx.reshape(-1, 4); qy = qy.reshape(-1, 4)
    x0 = 10 * np.floor(qx.min() / 10); y1 = 10 * np.ceil(qy.max() / 10)
    nx = int(np.ceil((qx.max() - x0) / 10)) + 1; ny = int(np.ceil((y1 - qy.min()) / 10)) + 1
    W = weights(qx, qy, 10, x0, y1, nx, ny)
    _CACHE[epsg] = (W, {"x0": x0, "y1": y1, "nx": nx, "ny": ny})
    return _CACHE[epsg]


def read(href, lat):
    with rasterio.open(href) as r:
        epsg = r.crs.to_epsg()
        return epsg, r


def resample(href_vv, href_vh, cells, W10, win, decim=2):
    """VV (and VH if href_vh) -> cells. decim=2 reads the 2x COG overview (20 m) and replicates it 2x2 onto
    the 10 m lattice (bandwidth; S1 IW GRD resolution is ~20 x 22 m, so little information is lost)."""
    out = {}
    with rasterio.open(href_vv) as r:
        epsg = r.crs.to_epsg()
    W, lw = operator(epsg, cells, W10, win)
    A = np.asarray(W.sum(1)).ravel()
    for pol, href in (("vv", href_vv), ("vh", href_vh)):
        if not href:
            continue
        with rasterio.open(href) as r:
            c0 = int(round((lw["x0"] - r.transform.c) / 10)); r0 = int(round((r.transform.f - lw["y1"]) / 10))
            win_ = Window(c0, r0, lw["nx"], lw["ny"])
            if decim > 1:
                a = r.read(1, window=win_, boundless=True, fill_value=r.nodata,
                           out_shape=(int(np.ceil(lw["ny"] / decim)), int(np.ceil(lw["nx"] / decim))))
                a = np.repeat(np.repeat(a, decim, 0), decim, 1)[:lw["ny"], :lw["nx"]]
            else:
                a = r.read(1, window=win_, boundless=True, fill_value=r.nodata)
        valid = np.isfinite(a) & (a > 0) & (a != r.nodata)
        m = valid.ravel().astype(np.float32)
        num = W @ (np.where(valid, a, 0).ravel().astype(np.float32))
        den = W @ m
        with np.errstate(invalid="ignore", divide="ignore"):
            lin = np.where(den / A >= 0.999, num / den, np.nan)
            out[pol] = (10 * np.log10(lin)).astype(np.float32)
    return out, epsg
