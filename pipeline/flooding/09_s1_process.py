"""09_s1_process.py <root> [nproc] — Sentinel-1 open water and flagged flooded vegetation on HEALPix d19.

Scene selection (per date and relative orbit; all slices of the datatake are combined per cell):
  V  validation: within 1 day of a LAST mask (Sep 2017 – Aug 2019) or of an S2 scene with >= 80 % valid cells
  F  flood-season fill: Oct–May dates with no S2 scene having >= 50 % valid cells within 2 days
  R  dry reference: Jul–Aug, up to 2 dates per relative orbit per hydrological year
Open water: VV(dB) < T, T = Otsu threshold of the VV histogram over marsh + lagoon + artificial-water cells
  (sea/estuary excluded), accepted if -24 <= T <= -12 dB, else the scene is not classified.
Flooded vegetation (FLAG, kept separate, not merged): marsh cells not open water with
  VV - VV_ref >= 3 dB, VV_ref = per-orbit median VV over all R scenes (double-bounce increase).
Codes: 0 dry, 1 open water, 2 flooded-vegetation flag, 255 no data.
Known error: wind roughening raises VV over open water -> omission (see validation by wind speed).
"""
import csv, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

sys.path.insert(0, str(Path(__file__).parent))
from common import dggs_attrs, load_grid
import s1scene

ROOT = Path(sys.argv[1]); NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 6
W, win, cells = load_grid(ROOT)
s1scene.operator(32630, cells)  # exact overlap operator for UTM 30N scenes, built once before forking workers
st = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")["stratum"][:]
wet = np.isin(st, (1, 2, 3)); marsh = st == 1
s1 = pd.read_csv(ROOT / "data/catalog_s1.csv")
s2 = pd.read_csv(ROOT / "data/s2/scenes.csv")
s2d = pd.to_datetime(s2.date)
last = pd.to_datetime(pd.Series(zarr.open_group(ROOT / "data/ref/last_masks.zarr", mode="r").attrs["dates"]), format="%Y%m%d")
s1["d"] = pd.to_datetime(s1.date)


def near(d, series, days):
    return (abs(series - d) <= pd.Timedelta(days=days)).any()


keys = s1.groupby(["date", "rel_orbit"]).size().reset_index()[["date", "rel_orbit"]]
sel = []
for _, k in keys.iterrows():
    d = pd.Timestamp(k.date); why = []
    if (pd.Timestamp("2017-09-01") <= d <= pd.Timestamp("2019-08-31") and near(d, last, 1)) or \
            near(d, s2d[s2.valid_share >= 0.8], 1):
        why.append("V")
    if d.month in (10, 11, 12, 1, 2, 3, 4, 5) and not near(d, s2d[s2.valid_share >= 0.5], 2):
        why.append("F")
    if d.month in (7, 8):
        why.append("R")
    if why:
        sel.append((k.date, int(k.rel_orbit), "".join(why)))
sel = pd.DataFrame(sel, columns=["date", "rel_orbit", "why"])
sel["hy"] = sel.date.apply(lambda d: int(d[:4]) - (int(d[5:7]) < 9))
r = sel[sel.why.str.contains("R")]
r = r.groupby(["hy", "rel_orbit"]).head(2)
sel = sel[sel.why.str.contains("V|F") | sel.index.isin(r.index)].reset_index(drop=True)
print(sel.why.value_counts().to_dict(), len(sel), flush=True)


def work(i):
    d, o = sel.date[i], sel.rel_orbit[i]
    its = s1[(s1.date == d) & (s1.rel_orbit == o)]
    vv = np.full(len(st), np.nan, np.float32); vh = vv.copy()
    t = time.time()
    try:
        import planetary_computer as pc
        for _, it in its.iterrows():
            out, epsg = s1scene.resample(pc.sign(it.vv), None, cells, W[10], win)
            f = np.isnan(vv) & np.isfinite(out["vv"])
            vv[f] = out["vv"][f]
    except Exception as e:
        return i, None, None, f"error {type(e).__name__}: {str(e)[:120]}", 0
    return i, vv.astype(np.float16), vh.astype(np.float16), "ok", time.time() - t


def otsu(x, lo=-35, hi=5, nb=200):
    h, e = np.histogram(x, bins=nb, range=(lo, hi))
    c = 0.5 * (e[1:] + e[:-1]); w0 = np.cumsum(h); w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1); m1 = (np.sum(h * c) - np.cumsum(h * c)) / np.maximum(w1, 1)
    return float(c[np.argmax(w0 * w1 * (m0 - m1) ** 2)])


OUT = ROOT / "data/s1"; OUT.mkdir(parents=True, exist_ok=True)
A = cells["area_utm_m2"] / 1e4
isref = sel.why.str.contains("R").values
# pass 1: dry references (per relative orbit median VV over the R scenes)
refstack = {}
with Pool(NPROC) as pool:
    for i, vv, vh, status, sec in pool.imap_unordered(work, list(np.nonzero(isref)[0])):
        print("ref", sel.date[i], sel.rel_orbit[i], status, round(sec), flush=True)
        if vv is not None:
            refstack.setdefault(sel.rel_orbit[i], []).append(vv)
ref = {o: np.nanmedian(np.stack(v).astype(np.float32), 0) for o, v in refstack.items()}
del refstack
g = zarr.open_group(OUT / "codes.zarr", mode="w", attributes=dggs_attrs())
codes = g.create_array("codes", shape=(len(sel), len(st)), chunks=(1, 1 << 20), dtype="uint8", fill_value=255,
                       compressors=zarr.codecs.BloscCodec(cname="zstd", clevel=7, shuffle="bitshuffle"),
                       dimension_names=["time", "cells"])
g.create_array("cell_ids", data=cells["cell_ids"], chunks=(1 << 20,), dimension_names=["cells"])
g.attrs.update({"dates": list(sel.date), "rel_orbit": [int(x) for x in sel.rel_orbit], "why": list(sel.why),
                "codes": "0 dry, 1 open water (VV < Otsu T), 2 flooded-vegetation FLAG (VV - dry ref >= 3 dB, marsh only), 255 no data"})
rows = []
todo = [i for i in range(len(sel)) if "V" in sel.why[i] or "F" in sel.why[i]]
with Pool(NPROC) as pool:
    for k, (i, vv, vh, status, sec) in enumerate(pool.imap_unordered(work, todo)):
        rec = {"date": sel.date[i], "rel_orbit": sel.rel_orbit[i], "why": sel.why[i], "status": status}
        if vv is None:
            rows.append(rec); print(k, rec, flush=True); continue
        vv = vv.astype(np.float32)
        ok = np.isfinite(vv) & (st != 0)
        rec["valid_share"] = float(ok.mean())
        x = vv[ok & wet]
        if x.size < 10000:
            rec["status"] = "too few cells"; rows.append(rec); continue
        T = otsu(x); rec["T_vv"] = T
        if not -24 <= T <= -12:
            rec["status"] = "threshold out of range"; rows.append(rec); continue
        c = np.where(ok, 0, 255).astype(np.uint8)
        c[ok & (vv < T)] = 1
        o = sel.rel_orbit[i]
        if o in ref:
            c[ok & marsh & (c == 0) & (vv - ref[o] >= 3.0)] = 2
        codes[i] = c
        rec.update({"open_water_ha": float(A[c == 1].sum()), "flooded_veg_flag_ha": float(A[c == 2].sum()),
                    "seconds": round(sec)})
        rows.append(rec)
        print(k, sel.date[i], o, sel.why[i], round(T, 1), round(rec["open_water_ha"]), round(rec["flooded_veg_flag_ha"]), round(sec), flush=True)
        if k % 20 == 0:
            pd.DataFrame(rows).to_csv(OUT / "scenes.csv", index=False)
pd.DataFrame(rows).to_csv(OUT / "scenes.csv", index=False)
print(pd.DataFrame(rows).status.value_counts())
