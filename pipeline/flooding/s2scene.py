"""Read one Sentinel-2 date (tiles 29SQA + 29SQB), resample EACH tile onto HEALPix d19 cells
(reflectance: conservative; SCL: majority; CLD: conservative mean of a probability), then classify
with kordelas.classify on the HEALPix base-face lattice.

Per-cell code (uint8): bits 0-1 Kordelas class (0 dry, 1 open water, 2 water-vegetation);
bit 2 = fixed-threshold water (B11 reflectance < FIXED_B11, Díaz-Delgado et al. 2016 rule applied to S2);
255 = no valid observation (cloud, shadow, no data, or scene rejected).
"""
import os
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window

sys.path.insert(0, str(Path(__file__).parent))
from common import CLEAR_SCL, CLOUD_SCL, conservative, majority, class_fractions
import kordelas

os.environ.update(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", AWS_NO_SIGN_REQUEST="YES",
                  GDAL_HTTP_MULTIPLEX="YES", GDAL_HTTP_VERSION="2", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                  GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="3", VSI_CACHE="FALSE")
FIXED_B11 = 0.186
TINIT_FALLBACK_B11 = 0.088   # IMPL: median Tfinal (B11 reflectance) of the first 32 scenes with a valid valley (Apr–Dec 2017)
MIN_VALID_SHARE = 0.20       # scenes with < 20 % valid domain cells are not classified
CLD_MAX = 20.0               # mean cloud probability (%) allowed in a cell (c1 items only)
# rows of each tile needed (QB covers the domain north of its south edge; QA only the rest)
TILE_YRANGE = {"29SQB": (4090200.0, 1e9), "29SQA": (0.0, 4090200.0 + 40.0)}


RGB_OVERVIEW = True  # DEVIATION (bandwidth): B02-B04 (segmentation input only, §2.4.1) read from the 2x COG
                     # overview (20 m) and replicated 2x2 onto the 10 m lattice; SWIR/red-edge bands untouched.


def read_band(href, res, win, tile, decim=1):
    nx, ny = win[f"nx{res}"], win[f"ny{res}"]
    out = np.zeros((ny, nx), np.uint16)
    if "blob.core.windows.net" in href:
        import planetary_computer
        href = planetary_computer.sign(href)
    with rasterio.open(href) as r:
        tb = r.bounds
        ylo, yhi = TILE_YRANGE.get(tile, (0, 1e9))
        bx0, bx1 = max(win["x0"], tb.left), min(win["x1"], tb.right)
        by0, by1 = max(win["y0"], tb.bottom, ylo), min(win["y1"], tb.top, yhi)
        if bx1 <= bx0 or by1 <= by0:
            return out
        c0 = int(round((bx0 - tb.left) / res)); r0 = int(round((tb.top - by1) / res))
        w = int(round((bx1 - bx0) / res)); h = int(round((by1 - by0) / res))
        if decim > 1:
            a = r.read(1, window=Window(c0, r0, w, h), out_shape=(int(np.ceil(h / decim)), int(np.ceil(w / decim))))
            a = np.repeat(np.repeat(a, decim, 0), decim, 1)[:h, :w]
        else:
            a = r.read(1, window=Window(c0, r0, w, h))
        oc = int(round((bx0 - win["x0"]) / res)); orow = int(round((win["y1"] - by1) / res))
        out[orow:orow + a.shape[0], oc:oc + a.shape[1]] = a
    return out


def _fetch(jobs, win, nthreads=8):
    """read several (href, res, tile) windows concurrently (latency-bound HTTP range requests)."""
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(nthreads) as ex:
        fut = {k: ex.submit(read_band, h, r, win, t, 2 if (RGB_OVERVIEW and r == 10) else 1) for k, (h, r, t) in jobs.items()}
        return {k: f.result() for k, f in fut.items()}


def process(rows, W, win, cells, face_idx, return_arrays=False):
    """rows: catalogue rows (dicts) of one date. Returns (code uint8[n_cells], info dict[, arrays])."""
    n = len(cells["cell_ids"])
    A20 = np.asarray(W[20].sum(1)).ravel()
    tiles = sorted(rows, key=lambda r: r["tile"] != "29SQB")  # QB first: wins ties
    jobs = {(r["tile"], "scl"): (r["scl"], 20, r["tile"]) for r in tiles}
    jobs.update({(r["tile"], "cloud"): (r["cloud"], 20, r["tile"]) for r in tiles if r.get("cloud")})
    pre = _fetch(jobs, win)
    per = []
    for r in tiles:
        scl = pre[(r["tile"], "scl")]
        cov = (W[20] @ (scl > 0).ravel().astype(np.float32)) / A20
        per.append({"row": r, "scl": scl, "cov": cov})
    covs = np.stack([p["cov"] for p in per])
    best = np.argmax(covs, 0)
    cov = covs[best, np.arange(n)]
    maj = np.zeros(n, np.uint8); cloudf = np.zeros(n, np.float32); cld = np.full(n, np.nan, np.float32)
    for t, p in enumerate(per):
        sel = best == t
        if not sel.any():
            continue
        m, fr = majority(W[20], p["scl"], list(range(1, 12)), p["scl"] > 0)
        maj[sel] = m[sel]
        cloudf[sel] = sum(fr[c] for c in CLOUD_SCL)[sel]
        if p["row"].get("cloud"):
            c_img = pre[(p["row"]["tile"], "cloud")].astype(np.float32)
            v, _ = conservative(W[20], c_img, p["scl"] > 0)
            cld[sel] = v[sel]
    valid = (cov >= 0.999) & np.isin(maj, CLEAR_SCL) & (cloudf <= 0.0) & ~(cld > CLD_MAX)
    info = {"date": rows[0]["date"], "items": ";".join(r["id"] for r in tiles), "valid_share": float(valid.mean()),
            "cloud_share": float((cloudf > 0).mean()), "maj_snow_share": float((maj == 11).mean())}
    code = np.full(n, 255, np.uint8)
    if valid.mean() < MIN_VALID_SHARE:
        info["status"] = "rejected: valid share"
        return (code, info, None) if return_arrays else (code, info)
    need = [(b, res) for b, res in (("B11", 20), ("B05", 20), ("B07", 20), ("B02", 10), ("B03", 10), ("B04", 10))]
    used = [t for t, p in enumerate(per) if ((best == t) & valid).any()]
    pre = _fetch({(per[t]["row"]["tile"], b): (per[t]["row"][b], res, per[t]["row"]["tile"]) for t in used for b, res in need}, win)
    bands = {}
    for b, res in (("B11", 20), ("B05", 20), ("B07", 20), ("B02", 10), ("B03", 10), ("B04", 10)):
        val = np.full(n, np.nan, np.float32)
        for t, p in enumerate(per):
            sel = (best == t) & valid
            if not sel.any():
                continue
            dn = pre[(p["row"]["tile"], b)]
            v, _ = conservative(W[res], dn.astype(np.float32), dn > 0)
            val[sel] = v[sel] * float(p["row"]["scale"]) + float(p["row"]["offset"])
        bands[b] = val
    valid &= np.all([np.isfinite(bands[b]) for b in bands], axis=0)
    ri, ci, shape = face_idx

    def fr(v, fill=np.nan):
        a = np.full(shape, fill, np.float32); a[ri, ci] = v; return a

    vface = np.zeros(shape, bool); vface[ri, ci] = valid
    swir_n, (p1, p99) = kordelas.stretch(fr(bands["B11"]), vface)
    rgb = [kordelas.stretch(fr(bands[b]), vface)[0] for b in ("B02", "B03", "B04")]
    with np.errstate(invalid="ignore", divide="ignore"):
        mndvi = (bands["B07"] - bands["B05"]) / (bands["B07"] + bands["B05"])
    fb = (TINIT_FALLBACK_B11 - p1) / max(p99 - p1, 1e-6) * 255.0
    cls, kinfo = kordelas.classify(None, vface, swir_n, rgb, fr(mndvi), tinit_fallback=fb)
    info.update(kinfo); info.update({"swir_p1": p1, "swir_p99": p99})
    if cls is None:
        info["status"] = "rejected: " + kinfo.get("status", "")
        return (code, info, None) if return_arrays else (code, info)
    c = cls[ri, ci]
    fixed = (bands["B11"] < FIXED_B11).astype(np.uint8)
    code = np.where(valid, c | (fixed << 2), 255).astype(np.uint8)
    A = A20 / 1e4
    ok = code != 255
    info.update({"open_water_ha": float(A[ok & ((code & 3) == 1)].sum()),
                 "water_veg_ha": float(A[ok & ((code & 3) == 2)].sum()),
                 "fixed_water_ha": float(A[ok & (((code >> 2) & 1) == 1)].sum()),
                 "valid_ha": float(A[ok].sum())})
    if return_arrays:
        return code, info, {"bands": bands, "valid": valid, "maj": maj, "cld": cld, "cloudf": cloudf}
    return code, info
