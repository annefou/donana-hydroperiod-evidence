"""04_last_to_healpix.py — LAST-EBD Landsat flood masks (Zenodo 10.5281/zenodo.3518182, CC BY 4.0) onto the HEALPix d19 cells.

Categorical -> nearest: each cell takes the 30 m pixel containing its centre (healpix-geo WGS84 centre,
transformed to ED50 / UTM 30N, EPSG:23030, the grid of the LAST products; the .tif files carry no CRS,
the .TFW world files give the transform). Masks: 1 = flooded, 0 = not flooded (no nodata value).
Each mask is extracted from the RAR one at a time and deleted after use (disk).
Output: data/ref/last_masks.zarr (masks[n_masks, n_cells] uint8; attrs dates, sensors)
"""
import sys, zipfile
from pathlib import Path

import libarchive
import numpy as np
import rasterio
import zarr
from affine import Affine
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).parent))
from common import load_grid

ROOT = Path(sys.argv[1])
RAR = ROOT / "data/raw/zenodo/inunda.rar"
TMP = ROOT / "data/raw/last_tmp"; TMP.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "data/ref"; OUT.mkdir(parents=True, exist_ok=True)

c = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")
lon, lat = c["lon"][:], c["lat"][:]
X, Y = Transformer.from_crs(4326, 23030, always_xy=True).transform(lon, lat)


def tfw_affine(txt):
    a, d, b, e, c0, f = [float(v) for v in txt.split()]
    return Affine(a, b, c0 - a / 2 - b / 2, d, e, f - d / 2 - e / 2)  # TFW gives pixel-centre of UL


def sample(arr, aff):
    col = np.floor((X - aff.c) / aff.a).astype(np.int64)
    row = np.floor((Y - aff.f) / aff.e).astype(np.int64)
    ok = (row >= 0) & (row < arr.shape[0]) & (col >= 0) & (col < arr.shape[1])
    out = np.full(X.shape, 255, np.uint8)
    out[ok] = arr[row[ok], col[ok]]
    return out


n = len(lon)
tfw = {}
with libarchive.file_reader(str(RAR)) as a:
    for e in a:
        if e.pathname.upper().endswith(".TFW"):
            tfw[e.pathname[:-4]] = b"".join(e.get_blocks()).decode()
names = sorted(tfw)
store = zarr.open_group(OUT / "last_masks.zarr", mode="w")
arr = store.create_array("masks", shape=(len(names), n), chunks=(1, 1 << 20), dtype="uint8",
                         compressors=zarr.codecs.BloscCodec(cname="zstd", clevel=7, shuffle="bitshuffle"))
done = 0
with libarchive.file_reader(str(RAR)) as a:
    for e in a:
        stem = e.pathname[:-4]
        if not e.pathname.lower().endswith(".tif") or stem not in tfw:
            continue
        p = TMP / "m.tif"
        with open(p, "wb") as f:
            for b in e.get_blocks():
                f.write(b)
        with rasterio.open(p) as r:
            m = r.read(1)
        aff = tfw_affine(tfw[stem])
        arr[names.index(stem)] = sample(m, aff)
        p.unlink(); done += 1
        if done % 50 == 0:
            print(done, stem, flush=True)
store.attrs.update({"names": names, "dates": [s[:8] for s in names], "sensors": [s[8:12] for s in names],
                    "source": "Díaz-Delgado et al., Flood Masks Doñana 1984/2019, Zenodo 10.5281/zenodo.3518182, CC BY 4.0",
                    "values": "1 flooded, 0 not flooded, 255 outside mask extent", "resampling": "nearest (cell centre)"})
print("masks", done)
