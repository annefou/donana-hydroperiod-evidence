"""01b_test_resampler.py — tests of the area-weighted UTM->HEALPix operator.

T1 constant field reproduced exactly; T2 integral conserved; T3 overlap weights agree with an
independent Monte-Carlo estimate that uses healpix-geo's own point-in-cell assignment
(lonlat_to_healpix, WGS84) instead of the projected cell polygon; T4 majority of a categorical
field equals the class with the largest Monte-Carlo share.
"""
import sys
from pathlib import Path

import healpix_geo.nested as hn
import numpy as np
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).parent))
from common import DEPTH, conservative, jdump, load_grid, majority

ROOT = Path(sys.argv[1])
W, win, cells = load_grid(ROOT)
rng = np.random.default_rng(0)
res = {}
for r in (10, 20):
    nx, ny = win[f"nx{r}"], win[f"ny{r}"]
    one = np.ones((ny, nx), np.float32)
    v, f = conservative(W[r], one, np.ones_like(one, bool))
    res[f"T1_const_maxdev_{r}m"] = float(np.nanmax(np.abs(v - 1)))
    img = rng.gamma(2.0, 0.05, (ny, nx)).astype(np.float32)
    v, f = conservative(W[r], img, np.ones_like(img, bool))
    A = np.asarray(W[r].sum(1)).ravel()
    colw = np.asarray(W[r].sum(0)).ravel()
    res[f"T2_integral_rel_err_{r}m"] = float(abs((v * A).sum() / (img.ravel() * colw).sum() - 1))

# T3: Monte-Carlo with healpix-geo point membership, 300 random cells, 20 m and 10 m lattice
to_utm = Transformer.from_crs(4326, 32629, always_xy=True)
to_ll = Transformer.from_crs(32629, 4326, always_xy=True)
idx = rng.choice(len(cells["cell_ids"]), 300, replace=False)
errs = {10: [], 20: []}
maj_ok = 0
for i in idx:
    cx, cy = to_utm.transform(cells["lon"][i], cells["lat"][i])
    px = cx + rng.uniform(-12, 12, 40000); py = cy + rng.uniform(-12, 12, 40000)
    lon, lat = to_ll.transform(px, py)
    inside = hn.lonlat_to_healpix(lon, lat, np.uint8(DEPTH), ellipsoid="WGS84") == cells["cell_ids"][i]
    px, py = px[inside], py[inside]
    for r in (10, 20):
        nx = win[f"nx{r}"]
        col = (np.floor((win["y1"] - py) / r).astype(int) * nx + np.floor((px - win["x0"]) / r).astype(int))
        u, cnt = np.unique(col, return_counts=True)
        mc = dict(zip(u, cnt / cnt.sum()))
        row = W[r][i]
        wv = dict(zip(row.indices, row.data / row.data.sum()))
        keys = set(mc) | set(wv)
        errs[r].append(max(abs(mc.get(k, 0) - wv.get(k, 0)) for k in keys))
        if r == 20:
            top_mc = max(mc, key=mc.get); top_w = max(wv, key=wv.get)
            maj_ok += int(top_mc == top_w or abs(mc[top_mc] - mc.get(top_w, 0)) < 0.02)
for r in (10, 20):
    e = np.array(errs[r])
    res[f"T3_weight_maxabs_err_{r}m_median"] = float(np.median(e))
    res[f"T3_weight_maxabs_err_{r}m_p99"] = float(np.quantile(e, 0.99))
res["T3_note"] = "Monte-Carlo with ~25,800 points per cell: sampling s.e. of a share ~0.003"
res["T4_majority_agree_20m"] = maj_ok / len(idx)
jdump(res, ROOT / "results_resampler_tests.json")
print(res)
