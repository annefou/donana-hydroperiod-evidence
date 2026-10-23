"""07_validate.py <root> — accuracy of the HEALPix S2 (and S1) flood maps against independent references.

R1 LAST-EBD Landsat flood masks (Zenodo 3518182; themselves a band-5 threshold product, Kappa ~0.65 vs field
   points per Díaz-Delgado et al. 2016) on S2 dates within <= 1 day (Sep 2017 – Aug 2019; also Apr–Aug 2017).
R2 Field points 24 Apr 2017 (Zenodo 17158473; class + S2 reflectance per point, no coordinates):
   point-level test of the scene thresholds.
R3 ICTS-Doñana marsh water-level gauges (2020–2026): gauge flooded (level > 0.05 m) vs S2 water at the gauge.

Kordelas 2018 §3.1 protocol: sea/coastal water excluded; results with boundary pixels excluded (EAP: cells within
~37 m (±3 cells on the face lattice) of a change of class in the LAST mask) and included (IAP); confusion matrices
summed over dates ("combined", Table 8). Strata: marsh, lagoon, other (+artificial), complete (all but sea).
Uncertainty: 95 % CI by bootstrap over spatial blocks (HEALPix depth-10 parent cells, ~6.4 km), 2000 draws.
Spatial block hold-out: 5-fold GroupKFold over the same blocks (no parameter of the S2 method is fitted to the
reference; the folds measure how much accuracy varies between independent areas).
Outputs: results_validation.json, data/validation_pairs.csv
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).parent))
from common import face_index, jdump

ROOT = Path(sys.argv[1])
c = zarr.open_group(ROOT / "data/grid/cells.zarr", mode="r")
cells = {k: c[k][:] for k in c.array_keys()}
n = len(cells["cell_ids"])
A = cells["area_utm_m2"].astype(np.float64) / 1e4
stg = zarr.open_group(ROOT / "data/grid/strata.zarr", mode="r")
st = stg["stratum"][:]
block = (cells["cell_ids"] >> np.uint64(2 * (19 - 10))).astype(np.int64)
ub, bidx = np.unique(block, return_inverse=True)
NB = len(ub)
ri, ci, shape = face_index(cells)
s2 = zarr.open_group(ROOT / "data/s2/codes.zarr", mode="r")
s2dates = pd.to_datetime(pd.Series(s2.attrs["dates"]))
scenes = pd.read_csv(ROOT / "data/s2/scenes.csv").set_index("date")
lm = zarr.open_group(ROOT / "data/ref/last_masks.zarr", mode="r")
ldates = pd.to_datetime(pd.Series(lm.attrs["dates"]), format="%Y%m%d")
lsens = list(lm.attrs["sensors"])
STRATA = {"complete": st != 0, "marsh": st == 1, "lagoon": st == 2, "other": np.isin(st, (3, 4))}
VARIANTS = {"open_water": lambda code: (code & 3) == 1,
            "open_plus_waterveg": lambda code: (code & 3) >= 1,
            "fixed_B11_0.186": lambda code: ((code >> 2) & 1) == 1}


def boundary(ref):
    a = np.zeros(shape, np.int8) - 1
    r16 = ref.astype(np.int16)
    a[ri, ci] = np.where(r16 == 255, -1, r16).astype(np.int8)
    mx = ndimage.maximum_filter(a, size=7); mn = ndimage.minimum_filter(np.where(a < 0, 2, a), size=7)
    return (mx[ri, ci] != mn[ri, ci])


def conf_blocks(pred, ref, m):
    """per-block confusion counts in ha: [TP, FP, FN, TN] x NB"""
    out = np.zeros((4, NB))
    for k, sel in enumerate((pred & ref, pred & ~ref, ~pred & ref, ~pred & ~ref)):
        out[k] = np.bincount(bidx[m & sel], weights=A[m & sel], minlength=NB)
    return out


def metrics(cm):
    tp, fp, fn, tn = cm
    N = tp + fp + fn + tn
    if N <= 0:
        return {}
    po = (tp + tn) / N
    pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / N ** 2
    return {"kappa": (po - pe) / (1 - pe) if pe < 1 else np.nan, "OA": po,
            "PA_water": tp / (tp + fn) if tp + fn else np.nan, "UA_water": tp / (tp + fp) if tp + fp else np.nan,
            "omission_water": fn / (tp + fn) if tp + fn else np.nan, "commission_water": fp / (tp + fp) if tp + fp else np.nan,
            "ref_water_ha": tp + fn, "map_water_ha": tp + fp, "area_ha": N}


rng = np.random.default_rng(42)
BOOT = rng.integers(0, NB, size=(2000, NB))
folds = rng.permutation(NB) % 5


def summarise(cmb):
    """cmb: [4, NB] summed over dates -> metrics + bootstrap CI + block-fold spread"""
    m = metrics(cmb.sum(1))
    if not m:
        return m
    ks = []
    for b in BOOT:
        k = metrics(cmb[:, b].sum(1)).get("kappa")
        ks.append(k)
    ks = np.array(ks, float)
    m["kappa_ci95"] = [float(np.nanquantile(ks, 0.025)), float(np.nanquantile(ks, 0.975))]
    fk = [metrics(cmb[:, folds == f].sum(1)).get("kappa", np.nan) for f in range(5)]
    m["kappa_block_folds"] = [round(float(x), 3) for x in fk]
    return {k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in m.items()}


# ---------- R1: S2 vs LAST ----------
pairs = []
for j, ld in enumerate(ldates):
    if ld < pd.Timestamp("2017-04-01"):
        continue
    dd = (s2dates - ld).abs()
    i = int(dd.idxmin())
    if dd[i] > pd.Timedelta(days=1):
        continue
    d = s2.attrs["dates"][i]
    if d not in scenes.index or str(scenes.loc[d, "status"]) != "ok":
        continue
    pairs.append((i, j, d, str(ld.date()), lsens[j], int(dd[i].days)))
print(len(pairs), "S2-LAST pairs", flush=True)
acc = {v: {s: {"EAP": np.zeros((4, NB)), "IAP": np.zeros((4, NB))} for s in STRATA} for v in VARIANTS}
per_date = []
for i, j, d, ldd, sens, dd in pairs:
    code = s2["codes"][i]; ref = lm["masks"][j]
    ok = (code != 255) & (ref != 255)
    rb = ref == 1
    bnd = boundary(ref)
    row = {"s2_date": d, "last_date": ldd, "sensor": sens, "dt_days": dd, "valid_share": float(ok.mean())}
    for v, f in VARIANTS.items():
        p = f(code)
        for s, sm in STRATA.items():
            for mode, mm in (("EAP", ok & sm & ~bnd), ("IAP", ok & sm)):
                cb = conf_blocks(p, rb, mm)
                acc[v][s][mode] += cb
                if v == "open_water" and mode == "EAP":
                    mt = metrics(cb.sum(1))
                    row[f"{s}_kappa"] = round(mt.get("kappa", np.nan), 3)
                    row[f"{s}_PA"] = round(mt.get("PA_water", np.nan), 3); row[f"{s}_UA"] = round(mt.get("UA_water", np.nan), 3)
                    row[f"{s}_ref_ha"] = round(mt.get("ref_water_ha", 0)); row[f"{s}_map_ha"] = round(mt.get("map_water_ha", 0))
    # water-vegetation subclass: how much of it is LAST water (commission check)
    wv = ok & ((code & 3) == 2) & (st != 0)
    row["waterveg_ha"] = float(A[wv].sum()); row["waterveg_last_water_share"] = float(A[wv & rb].sum() / max(A[wv].sum(), 1e-9))
    per_date.append(row)
pd.DataFrame(per_date).to_csv(ROOT / "data/validation_pairs.csv", index=False)
res = {"n_pairs": len(pairs), "pairs_by_sensor": pd.Series([p[4] for p in pairs]).value_counts().to_dict(),
       "combined": {v: {s: {m: summarise(acc[v][s][m]) for m in ("EAP", "IAP")} for s in STRATA} for v in VARIANTS}}
# L8-only (Landsat 7 SLC-off masks may carry gap artefacts)
l8 = [k for k, p in enumerate(pairs) if p[4].startswith("l8")]
res["n_pairs_l8"] = len(l8)

# ---------- R2: field points ----------
gt = pd.read_excel(ROOT / "data/raw/zenodo/groundtruth_20170424.xlsx", header=1)
fld = {}
near = scenes.loc[[d for d in scenes.index if "2017-04-17" <= d <= "2017-04-27"]]
for d, r in near.iterrows():
    if r.status != "ok":
        continue
    T = r.swir_p1 + r.Tfinal / 255 * (r.swir_p99 - r.swir_p1)
    Ti = r.swir_p1 + r.Tinit / 255 * (r.swir_p99 - r.swir_p1)
    out = {"scene": d, "T_final_B11": round(T, 4), "T_init_B11": round(Ti, 4)}
    for name, thr in (("Tfinal", T), ("Tinit", Ti), ("fixed_0.186", 0.186)):
        w = gt.B11 < thr
        out[name] = {cls: {"n": int((gt["Inundation class"] == cls).sum()),
                           "share_mapped_water": round(float(w[gt["Inundation class"] == cls].mean()), 3)}
                     for cls in gt["Inundation class"].unique()}
    fld[d] = out
res["field_points_20170424"] = fld
res["field_points_note"] = ("936 GNSS points, 24 Apr 2017 (Díaz-Delgado 2025, Zenodo 17158473). Classes as coded in the file: "
                            "In (inundated), Ss, Sh, Emp, En (codes not expanded in the file). Reflectance is the "
                            "file's own S2 pixel value; threshold from our scene of the nearest date.")
jdump(res, ROOT / "results_validation.json")
print(json.dumps({v: {s: res["combined"][v][s]["EAP"].get("kappa") for s in STRATA} for v in VARIANTS}, indent=1))
print(json.dumps(res["combined"]["open_water"]["complete"], indent=1))
