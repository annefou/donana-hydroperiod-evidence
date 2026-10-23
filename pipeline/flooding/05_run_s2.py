"""05_run_s2.py <root> <start> <end> [nproc] — classify every S2 date in [start, end] (resumable).

Writes data/s2/codes.zarr (codes[date_index, cell] uint8, one row per catalogue date; 255 = no obs)
and appends per-date metadata/thresholds to data/s2/scenes.csv.
"""
import csv, json, sys, time, traceback
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

sys.path.insert(0, str(Path(__file__).parent))
from common import dggs_attrs, face_index, load_grid
import s2scene

ROOT = Path(sys.argv[1]); START, END = sys.argv[2], sys.argv[3]
NPROC = int(sys.argv[4]) if len(sys.argv) > 4 else 10
cat = pd.read_csv(ROOT / "data/catalog_s2.csv", dtype=str, keep_default_na=False)
dates = sorted(cat.date.unique())
W, win, cells = load_grid(ROOT)
FI = face_index(cells)
OUT = ROOT / "data/s2"; OUT.mkdir(parents=True, exist_ok=True)
zpath = OUT / "codes.zarr"
if not zpath.exists():
    g = zarr.open_group(zpath, mode="w", attributes=dggs_attrs())
    g.create_array("codes", shape=(len(dates), len(cells["cell_ids"])), chunks=(1, 1 << 20), dtype="uint8",
                   fill_value=255, compressors=zarr.codecs.BloscCodec(cname="zstd", clevel=7, shuffle="bitshuffle"),
                   dimension_names=["time", "cells"])
    g.create_array("cell_ids", data=cells["cell_ids"], chunks=(1 << 20,), dimension_names=["cells"])
    g.create_array("time", data=np.array(dates, dtype="datetime64[D]").astype("int64"), dimension_names=["time"],
                   attributes={"units": "days since 1970-01-01", "calendar": "proleptic_gregorian"})
    g.attrs.update({"codes": "bits0-1: Kordelas 2018 class 0 dry/1 open water/2 water-vegetation; bit2: B11<0.186; 255 no valid obs",
                    "dates": dates})
g = zarr.open_group(zpath, mode="r+")
assert list(g.attrs["dates"]) == dates
meta = OUT / "scenes.csv"
done = set(pd.read_csv(meta).date) if meta.exists() else set()
todo = [d for d in dates if START <= d <= END and d not in done]
print(f"{len(todo)} dates to do", flush=True)


def work(d):
    rows = cat[cat.date == d].to_dict("records")
    t = time.time()
    try:
        code, info = s2scene.process(rows, W, win, cells, FI)
    except Exception as e:
        return d, None, {"date": d, "status": f"error: {type(e).__name__}: {str(e)[:150]}"}
    info["seconds"] = round(time.time() - t, 1)
    return d, code, info


FIELDS = ["date", "status", "items", "valid_share", "Tinit_fallback", "cloud_share", "maj_snow_share", "Tinit", "Mopt", "Tfinal",
          "Tupper", "TMNDVI", "n_segments", "n_inundated_segments", "n_patch_segments", "swir_p1", "swir_p99",
          "open_water_ha", "water_veg_ha", "fixed_water_ha", "valid_ha", "seconds"]
new = not meta.exists()
with open(meta, "a", newline="") as fh, Pool(NPROC) as pool:
    w = csv.DictWriter(fh, FIELDS, extrasaction="ignore")
    if new:
        w.writeheader()
    for k, (d, code, info) in enumerate(pool.imap_unordered(work, todo)):
        if code is not None:
            g["codes"][dates.index(d)] = code
        if not str(info.get("status", "")).startswith("error"):
            w.writerow(info); fh.flush()
        print(k, d, info.get("status"), round(info.get("valid_share", -1), 2), info.get("Tfinal"),
              round(info.get("open_water_ha", 0)), round(info.get("water_veg_ha", 0)), info.get("seconds"), flush=True)
