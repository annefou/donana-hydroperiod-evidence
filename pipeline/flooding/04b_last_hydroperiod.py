"""04b_last_hydroperiod.py — LAST-EBD annual hydroperiod maps (Hydroperiod Doñana 1974/2019,
doi:10.5281/zenodo.3516929, CC BY 4.0; 43 cycles 1974-75 .. 2017-18, 1979-80 missing; 30 m, int16 days,
ED50 / UTM 30N EPSG:23030) onto HEALPix d19 cells by nearest (cell centre). One file extracted at a time.
Output: data/ref/last_hydroperiod.zarr (hydroperiod_days[n_years, n_cells] int16; attrs cycles)
"""
import sys, zipfile
from pathlib import Path

import numpy as np
import rasterio
import zarr
from pyproj import Transformer

ROOT = Path(sys.argv[1])
TMP = ROOT / "data/raw/last_tmp"; TMP.mkdir(parents=True, exist_ok=True)
c = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")
X, Y = Transformer.from_crs(4326, 23030, always_xy=True).transform(c["lon"][:], c["lat"][:])
with zipfile.ZipFile(ROOT / "data/raw/zenodo/hidroperiodos.zip") as z:
    names = sorted(x for x in z.namelist() if x.endswith("_v01.tif"))
    g = zarr.open_group(ROOT / "data/ref/last_hydroperiod.zarr", mode="w")
    arr = g.create_array("hydroperiod_days", shape=(len(names), len(X)), chunks=(1, 1 << 20), dtype="int16",
                         fill_value=-1, compressors=zarr.codecs.BloscCodec(cname="zstd", clevel=7))
    cycles = []
    for k, nm in enumerate(names):
        z.extract(nm, TMP)
        with rasterio.open(TMP / nm) as r:
            h = r.read(1); a = r.transform; nod = r.nodata
        col = np.floor((X - a.c) / a.a).astype(np.int64); row = np.floor((Y - a.f) / a.e).astype(np.int64)
        ok = (row >= 0) & (row < h.shape[0]) & (col >= 0) & (col < h.shape[1])
        v = np.full(len(X), -1, np.int16); v[ok] = h[row[ok], col[ok]]
        if nod is not None:
            v[v == nod] = -1
        arr[k] = v
        (TMP / nm).unlink()
        yy = Path(nm).stem.split("_")[0][1:]  # e.g. H1718 -> 1718
        cycles.append(f"{1900 + int(yy[:2]) if int(yy[:2]) > 50 else 2000 + int(yy[:2])}/{yy[2:]}")
        print(nm, cycles[-1], int(v.max()), flush=True)
g.attrs.update({"cycles": cycles, "source": "Hydroperiod Doñana 1974/2019, Zenodo 10.5281/zenodo.3516929, CC BY 4.0",
                "resampling": "nearest (cell centre)", "units": "days flooded per hydrological cycle (Sep-Aug)"})
