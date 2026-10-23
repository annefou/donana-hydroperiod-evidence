"""03_test_scene.py <root> <date> — run one date and save a diagnostic figure (face lattice)."""
import sys, time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import face_index, load_grid
import s2scene

ROOT = Path(sys.argv[1]); date = sys.argv[2]
cat = pd.read_csv(ROOT / "data/catalog_s2.csv", dtype=str, keep_default_na=False)
rows = cat[cat.date == date].to_dict("records")
W, win, cells = load_grid(ROOT)
fi = face_index(cells)
t = time.time()
code, info, arr = s2scene.process(rows, W, win, cells, fi, return_arrays=True)
print(f"{time.time()-t:.0f}s", info)
ri, ci, shape = fi
def fr(v, fill=np.nan):
    a = np.full(shape, fill, np.float32); a[ri, ci] = v; return a
fig, ax = plt.subplots(1, 3, figsize=(18, 6.5))
if arr is not None:
    ax[0].imshow(np.clip(fr(arr["bands"]["B11"]), 0, 0.5), cmap="gray"); ax[0].set_title(f"B11 {date}")
c = fr(np.where(code == 255, np.nan, code & 3))
ax[1].imshow(c, cmap="viridis", vmin=0, vmax=2, interpolation="nearest"); ax[1].set_title("Kordelas: 0 dry 1 open 2 water-veg")
ax[2].imshow(fr(np.where(code == 255, np.nan, (code >> 2) & 1)), interpolation="nearest"); ax[2].set_title("B11<0.186")
for a in ax: a.set_axis_off()
(ROOT / "figures/diag").mkdir(parents=True, exist_ok=True)
fig.savefig(ROOT / f"figures/diag/scene_{date}.png", dpi=70, bbox_inches="tight")
