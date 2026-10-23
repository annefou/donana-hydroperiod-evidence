"""Kordelas et al. 2018 (Remote Sens. 10:910, doi:10.3390/rs10060910) unsupervised inundation mapping,
re-implemented on HEALPix depth-19 cells (base-face (fx, fy) lattice). Section numbers refer to the paper.

§2.4.1 mean-shift segmentation of B2/B3/B4, each stretched 1–99 % -> [0,255], hs = 3, hr = 3
§2.4.2 SWIR (B11) stretched to [0,255]; Tinit = first deep valley of the SWIR histogram; segments with
       > 70 % of cells below Tinit -> centroids; square patches of 20k x 20k 10-m pixels (k = 1..20;
       here 16k x 16k 12.4-m cells, same ground size), kept if bimodal, MCET threshold per patch;
       median over k, then median over segments = Mopt; Tfinal = max(Mopt, Tinit); open water: SWIR < Tfinal
§2.4.5 water-vegetation: Tfinal <= SWIR < Tupper (next valley after the first) and MNDVI > TMNDVI,
       MNDVI = (B7 - B5)/(B7 + B5), TMNDVI = first valley of the MNDVI histogram above 0.4.
Choices not specified in the paper (marked IMPL): histogram smoothing, the 'deep valley' rule, the
bimodality test, the minimum segment size and the cap on the number of segments.
"""
import cv2
import numpy as np
from scipy.ndimage import gaussian_filter1d
from skimage.measure import label

SMOOTH_SIGMA = 1.5        # IMPL: Gaussian smoothing of the 256-bin histogram (bins)
VALLEY_DEPTH = 0.25       # IMPL: valley is 'deep' if h(valley) <= (1 - 0.25) * min(neighbouring peaks)
MIN_PEAK_SHARE = 0.002    # IMPL: a peak must hold >= 0.2 % of the cells within +-5 bins
MIN_SEG_CELLS = 20        # IMPL: segments smaller than 20 cells (~0.3 ha) are not used as patch centres
MAX_SEGMENTS = 300        # IMPL: random (seeded) subset of qualifying segments, for run time
PATCH_UNIT = 16           # 20 S2 10-m pixels ~= 16 HEALPix d19 cells (12.4 m)
FIRST_PEAK_MAX = 60      # IMPL: the first histogram peak must be a low-SWIR (water) mode
TUPPER_MAX = 100         # IMPL: Tupper searched only up to 100 (paper: observed Tupper 29-67, Section 3.2)


def stretch(v, valid):
    p1, p99 = np.nanpercentile(v[valid], [1, 99])
    return np.clip((v - p1) / max(p99 - p1, 1e-6) * 255.0, 0, 255), (float(p1), float(p99))


def _peaks_valleys(h):
    # IMPL: reflect-free zero padding so that a peak in the first/last bin (e.g. open water clipped
    # to 0 by the 1-99 % stretch) counts as a peak
    hs = gaussian_filter1d(np.concatenate([[0.0], h.astype(float), [0.0]]), SMOOTH_SIGMA, mode="constant")
    tot = hs.sum()
    if tot <= 0:
        return hs, [], []
    d = np.diff(hs)
    peaks = [i for i in range(1, len(hs) - 1) if d[i - 1] > 0 and d[i] <= 0]
    peaks = [p for p in peaks if hs[max(0, p - 5):p + 6].sum() / tot >= MIN_PEAK_SHARE]
    valleys = []
    for a, b in zip(peaks[:-1], peaks[1:]):
        v = a + int(np.argmin(hs[a:b + 1]))
        if hs[v] <= (1 - VALLEY_DEPTH) * min(hs[a], hs[b]):
            valleys.append(v)
    return hs[1:-1], [p - 1 for p in peaks], [v - 1 for v in valleys]


def valleys(values, nbins=256, rng=(0, 256)):
    h, edges = np.histogram(values, bins=nbins, range=rng)
    _, _, v = _peaks_valleys(h)
    return [edges[i] + 0.5 * (edges[1] - edges[0]) for i in v]


def bimodal(h):
    _, peaks, vals = _peaks_valleys(h)
    return len(vals) >= 1


def mcet(h):
    """Minimum cross-entropy threshold (Li & Lee 1993; iterative Li & Tam 1998) on a 256-bin histogram."""
    x = np.arange(len(h), dtype=float) + 1.0  # +1 avoids log(0)
    t = (h * x).sum() / max(h.sum(), 1)
    for _ in range(100):
        lo = x <= t
        hl, hh = h[lo].sum(), h[~lo].sum()
        if hl == 0 or hh == 0:
            return None
        ml = (h[lo] * x[lo]).sum() / hl
        mh = (h[~lo] * x[~lo]).sum() / hh
        tn = (mh - ml) / (np.log(mh) - np.log(ml))
        if abs(tn - t) < 0.5:
            t = tn
            break
        t = tn
    return float(t - 1.0)


def first_valley(values):
    """Tinit: first deep valley after a low (water) peak; IMPL guards: first peak <= 60, valley <= TUPPER_MAX."""
    h, _ = np.histogram(values, bins=256, range=(0, 256))
    _, peaks, vals = _peaks_valleys(h)
    if not peaks or peaks[0] > FIRST_PEAK_MAX:
        return None, []
    vals = [v + 0.5 for v in vals if v > peaks[0]]
    if not vals or vals[0] > TUPPER_MAX:
        return None, vals
    return vals[0], vals


def classify(face, valid, swir_n, rgb_n, mndvi, seed=0, tinit_fallback=None):
    """face-lattice arrays (2-D). Returns class array (0 dry, 1 open water, 2 water-vegetation), info dict."""
    info = {}
    sv = swir_n[valid]
    tinit, v = first_valley(sv)
    info["Tinit_fallback"] = tinit is None
    if tinit is None:
        if tinit_fallback is None:
            info["status"] = "no SWIR valley"
            return None, info
        tinit = float(tinit_fallback)
    info["Tinit"] = tinit
    # §2.4.1 segmentation
    img = np.zeros(swir_n.shape + (3,), np.uint8)
    for k in range(3):
        img[..., k] = np.where(valid, rgb_n[k], 0).astype(np.uint8)
    ms = cv2.pyrMeanShiftFiltering(img, 3, 3)
    key = ms[..., 0].astype(np.int32) * 65536 + ms[..., 1].astype(np.int32) * 256 + ms[..., 2]
    key[~valid] = -1
    seg = label(key, background=-1, connectivity=1)
    lab = seg[valid]
    inund = swir_n[valid] < tinit
    n = np.bincount(lab, minlength=seg.max() + 1)
    ni = np.bincount(lab, weights=inund, minlength=seg.max() + 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = ni / n
    cand = np.nonzero((frac > 0.7) & (n >= MIN_SEG_CELLS))[0]
    cand = cand[cand > 0]
    info["n_segments"] = int(seg.max()); info["n_inundated_segments"] = int(len(cand))
    rng = np.random.default_rng(seed)
    if len(cand) > MAX_SEGMENTS:
        cand = rng.choice(cand, MAX_SEGMENTS, replace=False)
    rr, cc = np.nonzero(valid)
    sr = np.bincount(lab, weights=rr, minlength=seg.max() + 1)
    sc = np.bincount(lab, weights=cc, minlength=seg.max() + 1)
    fm = []
    H, Wd = swir_n.shape
    for s in cand:
        r0, c0 = int(sr[s] / n[s]), int(sc[s] / n[s])
        ts = []
        for k in range(1, 21):
            hk = PATCH_UNIT * k // 2
            sl = (slice(max(0, r0 - hk), min(H, r0 + hk)), slice(max(0, c0 - hk), min(Wd, c0 + hk)))
            pv = swir_n[sl][valid[sl]]
            if pv.size < 50:
                continue
            h = np.bincount(np.clip(pv.astype(int), 0, 255), minlength=256)
            if not bimodal(h):
                continue
            t = mcet(h)
            if t is not None:
                ts.append(t)
        if ts:
            fm.append(np.median(ts))
    mopt = float(np.median(fm)) if fm else float("nan")
    tfinal = max(mopt, tinit) if fm else tinit
    info.update({"Mopt": mopt, "n_patch_segments": len(fm), "Tfinal": tfinal})
    cls = np.zeros(swir_n.shape, np.uint8)
    cls[valid & (swir_n < tfinal)] = 1
    # §2.4.5 water-vegetation subclass
    tup = [x for x in v if tinit < x <= TUPPER_MAX]
    tupper = tup[0] if tup else None
    tm = None
    mv = mndvi[valid]
    mv = mv[np.isfinite(mv) & (mv > 0.4)]
    if mv.size > 1000:
        vm = valleys(mv, nbins=120, rng=(0.4, 1.0))
        tm = vm[0] if vm else None
    info["Tupper"] = tupper; info["TMNDVI"] = tm
    if tupper is not None and tm is not None and tupper > tfinal:
        wv = valid & (swir_n >= tfinal) & (swir_n < tupper) & (mndvi > tm)
        cls[wv] = 2
    info["status"] = "ok"
    return cls, info
