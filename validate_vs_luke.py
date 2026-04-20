"""Validate herbs_map outputs against Luke's originals for every Sess with outputs.

For each Sess folder under MICROSCOPE/Jul30 herbs to ccf:
    1. Parse session label + recording side from the folder name.
    2. Find Luke's matching xpos/ypos .mat files (by unit count).
    3. Run herbs_map.run_session in-memory (no saving) and compare to
       Luke's unit_sites/unit_3d_coords/unit_voxels .npy.

Prints a one-line PASS/FAIL per session and a summary.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, "/home/thw051/CODE/PYTHON/HERBS")
import herbs_map as hm

JUL30 = Path("/n/files/Neurobio/MICROSCOPE/Tom/CODE/LC code/processed python files/Jul30 herbs to ccf")
MAT_ROOT = Path("/n/files/Neurobio/MICROSCOPE/Tom/CODE/LC code/processed mat + csv files")

sess_num_re = re.compile(r"Sess\s+(\d+)\s")
label_re = re.compile(r"(\d{8}_M?\d+_g\d+)")
record_side_re = re.compile(r"hand\s+([LR])\s+record", re.IGNORECASE)


def parse_sess(sess_dir: Path):
    name = sess_dir.name
    sess_m = sess_num_re.search(name)
    label_m = label_re.search(name)
    side_m = record_side_re.search(name)
    if not (sess_m and label_m and side_m):
        return None
    return dict(
        sess_num=int(sess_m.group(1)),
        label=label_m.group(1),
        recording_side=side_m.group(1).upper(),
    )


def find_xpos_ypos_for_sess(sess_num: int, n_units: int):
    """Find the best-matching xpos/ypos .mat files for a Sess by unit count.

    Returns (xpos_path, ypos_path) or (None, None) if no good match.
    Prefers newest mtime among candidates that match n_units.
    """
    pat_x = f"xpos*sess{sess_num}.mat"
    pat_y = f"ypos*sess{sess_num}.mat"
    cands_x = sorted(MAT_ROOT.rglob(pat_x), key=lambda p: p.stat().st_mtime, reverse=True)
    cands_y = sorted(MAT_ROOT.rglob(pat_y), key=lambda p: p.stat().st_mtime, reverse=True)

    for xp in cands_x:
        try:
            xv = loadmat(str(xp))["xpos"].flatten()
        except Exception:
            continue
        if len(xv) != n_units:
            continue
        # pick a ypos in the same date folder if possible, else matching length
        same_dir = [y for y in cands_y if y.parent == xp.parent]
        for yp in same_dir + cands_y:
            try:
                yv = loadmat(str(yp))["ypos"].flatten()
            except Exception:
                continue
            if len(yv) == n_units:
                return xp, yp
    return None, None


def load_luke_outputs(sess_dir: Path):
    out = {}
    for name in ("unit_sites", "unit_3d_coords", "unit_voxels"):
        p = sess_dir / f"{name}.npy"
        if not p.exists():
            return None
        out[name] = np.load(p)
    return out


def run_and_compare(sess_dir: Path, info: dict):
    pkls = hm.find_probe_pkls(sess_dir)
    if len(pkls) == 0:
        return dict(status="SKIP_NO_PKLS", **info)

    luke = load_luke_outputs(sess_dir)
    if luke is None:
        return dict(status="SKIP_NO_LUKE_OUTPUTS", **info)
    n_units = len(luke["unit_sites"])

    xp, yp = find_xpos_ypos_for_sess(info["sess_num"], n_units)
    if xp is None:
        return dict(status="SKIP_NO_XPOS_MATCH", n_units=n_units, **info)

    xpos = loadmat(str(xp))["xpos"].flatten()
    ypos = loadmat(str(yp))["ypos"].flatten()

    sl, rsl, vx, cs = hm.load_pkls(pkls)
    try:
        hb = hm.infer_hemisphere_backwards(pkls, info["recording_side"])
    except ValueError as e:
        return dict(status="SKIP_NO_PKL_SIDE", err=str(e), n_units=n_units, **info)
    if hb:
        hm.flip_hemisphere(sl, vx)
    hm.normalize_rel_site_locs(rsl, n_shanks=len(rsl))

    u3d, uvx_raw, ust = hm.map_units(xpos, ypos, sl, rsl, vx, cs)
    uvx = hm.voxels_to_ccf_axes(uvx_raw)

    results = {}
    for name, mine in [("unit_sites", ust), ("unit_3d_coords", u3d), ("unit_voxels", uvx)]:
        l = luke[name]
        results[name] = dict(equal=bool(np.array_equal(mine, l)),
                             max_diff=float(np.abs(mine.astype(float) - l.astype(float)).max()))

    all_equal = all(r["equal"] for r in results.values())
    return dict(
        status="PASS" if all_equal else "FAIL",
        n_units=n_units,
        hemisphere_backwards=hb,
        xpos=str(xp),
        results=results,
        **info,
    )


def main():
    sess_dirs = sorted([p for p in JUL30.iterdir() if p.is_dir() and p.name.startswith("Sess")])
    results = []
    for sd in sess_dirs:
        info = parse_sess(sd)
        if info is None:
            print(f"SKIP_PARSE  {sd.name}")
            continue
        r = run_and_compare(sd, info)
        results.append((sd, r))
        if r["status"] == "PASS":
            print(f"PASS  Sess {info['sess_num']:>2}  {info['label']:<22} hb={r['hemisphere_backwards']}  n={r['n_units']}")
        elif r["status"] == "FAIL":
            diffs = {k: v['max_diff'] for k,v in r['results'].items()}
            print(f"FAIL  Sess {info['sess_num']:>2}  {info['label']:<22} hb={r['hemisphere_backwards']}  n={r['n_units']}  {diffs}")
        else:
            print(f"{r['status']:<22} Sess {info['sess_num']:>2}  {info['label']}")

    n_pass = sum(1 for _, r in results if r["status"] == "PASS")
    n_fail = sum(1 for _, r in results if r["status"] == "FAIL")
    n_skip = len(results) - n_pass - n_fail
    print(f"\nSummary: {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP  (of {len(results)})")


if __name__ == "__main__":
    main()
