"""Map Neuropixels-sorted units to CCF coordinates from HERBS probe pkl exports.

The notebook mapToUnits.ipynb is the canonical reference; this module factors
its pipeline into reusable functions so a CLI driver and the notebook can share
the same code path.

Inputs per session:
    * 4 HERBS `probe_*.pkl` files (one per NP2 shank, ordered 1..4 A..P)
      Each pkl contains (under `data`):
        - sites_loc_b:        list[2] of ndarray[n_sites, 3]   CCF coords per site
        - sites_loc_relative: list[2] of ndarray[n_sites, 3]   shank-local (y, x-col)
        - sites_vox:          list[2] of ndarray[n_sites, 3]   Allen 10 um voxel idx
        - sites_label:        list[2] of ndarray[n_sites]      CCF region IDs
    * xpos.mat, ypos.mat — 1D arrays of per-unit probe coords from spike sorting

Outputs:
    * unit_3d_coords  (n_units, 3)  HERBS xyz (+x right, +y anterior, +z dorsal)
    * unit_voxels     (n_units, 3)  (AP, DV, LR) in 10 um voxels
    * unit_sites      (n_units,)    Allen CCF region ID per unit
"""
from __future__ import annotations

import pickle
import re
import sys
from math import floor
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.io import loadmat, savemat


# CCF constants (Allen 10 um)
CCF_WIDTH_X = 1140
CCF_WIDTH_Y = 1320
CCF_WIDTH_Z = 800
CCF_RESOLUTION = 0.1  # mm per voxel (10 um)

# NP2 probe geometry (as used in the notebook)
SHANK_PITCH_UM = 250
COL_PITCH_UM = 32

_REL_KEY_CANDIDATES = (
    "sites_loc_relative",
    "site_loc_relative",
    "sites_loc_rel",
    "site_loc_rel",
    "rel_sites_loc",
)


# ---------- NumPy 1 <-> 2 pickle compatibility ----------

def _ensure_numpy_core_alias() -> None:
    """Alias numpy._core -> numpy.core so NumPy-2 pickles load on NumPy-1."""
    try:
        import numpy._core  # noqa: F401
    except ModuleNotFoundError:
        import numpy.core as _np_core
        sys.modules["numpy._core"] = _np_core
        sys.modules.setdefault("numpy._core.multiarray", _np_core.multiarray)
        sys.modules.setdefault("numpy._core._multiarray_umath", _np_core.multiarray)
        sys.modules.setdefault("numpy._core.umath", _np_core.umath)


class _NPCoreFixUnpickler(pickle.Unpickler):
    """Fallback unpickler: rewrites numpy._core.* -> numpy.core.* on the fly."""

    def find_class(self, module, name):
        if module.startswith("numpy._core"):
            module = module.replace("numpy._core", "numpy.core", 1)
        return super().find_class(module, name)


# ---------- Core pipeline ----------

def _get_rel(d: dict):
    for k in _REL_KEY_CANDIDATES:
        if k in d:
            return d[k], k
    return None, None


def _nskey(p: Path):
    return [int(s) if s.isdigit() else s.lower() for s in re.split(r"(\d+)", p.name)]


def find_probe_pkls(session_dir: Path) -> list[Path]:
    """Return sorted probe_*.pkl files in a session directory (natural order)."""
    return sorted(Path(session_dir).glob("probe_*.pkl"), key=_nskey)


_PKL_SIDE_RE = re.compile(r"probe_\d+_([LR])[_.]", re.IGNORECASE)


def detect_pkl_side(pkl_paths: Sequence[Path]) -> str | None:
    """Return 'L' or 'R' if all pkls agree on side in their filenames, else None.

    Matches names like probe_1_L_tritc.pkl, probe_2_R_cy5.pkl.
    """
    sides = set()
    for p in pkl_paths:
        m = _PKL_SIDE_RE.search(Path(p).name)
        if m:
            sides.add(m.group(1).upper())
    if len(sides) == 1:
        return next(iter(sides))
    return None


def infer_hemisphere_backwards(pkl_paths: Sequence[Path], recording_side: str) -> bool:
    """Return True if the pkl filename side disagrees with the recording side.

    The HERBS pkls are labelled by which side they were *drawn on in HERBS*.
    If that disagrees with the physical recording side, the shanks were drawn in
    the wrong hemisphere and CCF x must be flipped.
    """
    side = recording_side.strip().upper()
    if side not in ("L", "R"):
        raise ValueError(f"recording_side must be 'L' or 'R', got {recording_side!r}")
    pkl_side = detect_pkl_side(pkl_paths)
    if pkl_side is None:
        raise ValueError(
            "Could not determine a single pkl hemisphere from filenames; "
            "pass hemisphere_backwards explicitly."
        )
    return pkl_side != side


def find_xy_mats(session_dir: Path) -> tuple[Path, Path]:
    """Return (xpos_mat, ypos_mat) for a session directory; raise if ambiguous."""
    session_dir = Path(session_dir)
    xpos = sorted(session_dir.glob("*xpos*.mat"))
    ypos = sorted(session_dir.glob("*ypos*.mat"))
    if len(xpos) != 1 or len(ypos) != 1:
        raise FileNotFoundError(
            f"Expected exactly one *xpos*.mat and one *ypos*.mat in {session_dir}, "
            f"got {len(xpos)} and {len(ypos)}"
        )
    return xpos[0], ypos[0]


def load_pkls(pkl_paths: Sequence[Path]) -> tuple[list, list, list, list]:
    """Load the 4 probe pkls and return (site_locs, rel_site_locs, voxels, ccf_sites).

    Each element of each list is itself a list[2] (columns) of ndarrays.
    Raises KeyError if a pkl is missing sites_loc_relative (or a known alias).
    """
    _ensure_numpy_core_alias()
    site_locs, rel_site_locs, voxels, ccf_sites = [], [], [], []
    for p in pkl_paths:
        with open(p, "rb") as f:
            try:
                loaded = pickle.load(f)
            except ModuleNotFoundError:
                f.seek(0)
                loaded = _NPCoreFixUnpickler(f).load()
        d = loaded["data"]
        rel, used_key = _get_rel(d)
        if rel is None:
            raise KeyError(
                f"{p.name}: no sites_loc_relative-like key; "
                f"available: {sorted(d.keys())[:10]}..."
            )
        if used_key != "sites_loc_relative":
            print(f"[{p.name}] using '{used_key}' for sites_loc_relative")
        site_locs.append(d["sites_loc_b"])
        rel_site_locs.append(rel)
        voxels.append(d["sites_vox"])
        ccf_sites.append(d["sites_label"])
    return site_locs, rel_site_locs, voxels, ccf_sites


def flip_hemisphere(site_locs, voxels) -> None:
    """In-place: negate CCF x in site_locs, mirror voxel x around CCF_WIDTH_X."""
    for shank in site_locs:
        for col in shank:
            col[:, 0] = -col[:, 0]
    for shank in voxels:
        for col in shank:
            col[:, 0] = CCF_WIDTH_X - col[:, 0]


def normalize_rel_site_locs(rel_site_locs, n_shanks: int = 4) -> None:
    """In-place: flip y so lowest = 0, rewrite col-x to absolute shank positions.

    After this:
      [:, 0] = y with 0 at probe tip (original was with 0 at top)
      [:, 1] = x in um, with shank `s` columns at {0+250*s, 32+250*s}
    """
    for shank in range(n_shanks):
        for col in range(2):
            arr = rel_site_locs[shank][col]
            arr[:, 0] -= np.min(arr[:, 0])
            arr[:, 0] = (arr[:, 0] - np.max(arr[:, 0])) * -1
            x_vals = arr[:, 1]
            x_vals[x_vals == 16] = COL_PITCH_UM + SHANK_PITCH_UM * shank
            x_vals[x_vals == -8] = 0 + SHANK_PITCH_UM * shank
            arr[:, 1] = x_vals


def map_units(
    xpos: np.ndarray,
    ypos: np.ndarray,
    site_locs,
    rel_site_locs,
    voxels,
    ccf_sites,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assign each (xpos, ypos) unit to its nearest electrode site.

    Returns (unit_3d_coords, unit_voxels, unit_sites) shaped (n_units, 3/3/1).
    Units with invalid/empty-shank mappings get zeros.
    """
    n_units = len(xpos)
    n_shanks = len(rel_site_locs)
    unit_3d = np.zeros((n_units, 3))
    unit_vx = np.zeros((n_units, 3))
    unit_st = np.zeros(n_units)

    for i in range(n_units):
        cur_xpos = xpos[i]
        cur_shank_idx = int(np.clip(round(cur_xpos / SHANK_PITCH_UM), 0, n_shanks - 1))
        cur_col_idx = int(np.clip(
            floor((cur_xpos - SHANK_PITCH_UM * cur_shank_idx) / COL_PITCH_UM), 0, 1
        ))
        # flip columns if more anterior (notebook heuristic)
        try:
            v = voxels[cur_shank_idx]
            if (
                len(v) >= 2
                and len(v[0]) > 0
                and len(v[1]) > 0
                and v[0][0][1] < v[1][0][1]
            ):
                cur_col_idx = 1 - cur_col_idx
        except Exception:
            pass

        cur_shank = rel_site_locs[cur_shank_idx]
        if cur_col_idx >= len(cur_shank) or cur_shank[cur_col_idx].size == 0:
            continue
        cur_col = cur_shank[cur_col_idx]

        # nearest by y
        y_col = cur_col[:, 0]
        closest = int(np.argmin(np.abs(y_col - ypos[i])))

        unit_3d[i] = site_locs[cur_shank_idx][cur_col_idx][closest]
        unit_vx[i] = voxels[cur_shank_idx][cur_col_idx][closest]
        unit_st[i] = int(ccf_sites[cur_shank_idx][cur_col_idx][closest])

    return unit_3d, unit_vx, unit_st


def voxels_to_ccf_axes(unit_voxels: np.ndarray) -> np.ndarray:
    """Transform (ML L-R, PA, VD) voxel coords -> (AP, DV, LR) in 10 um voxel units.

    Returns a NEW array (does not modify input).
    """
    out = unit_voxels.copy()
    out[:, 2] = CCF_WIDTH_Z - out[:, 2]
    out[:, 1] = CCF_WIDTH_Y - out[:, 1]
    out = out[:, [1, 2, 0]]
    out = out / CCF_RESOLUTION
    return out


def save_outputs(
    out_dir: Path,
    unit_3d: np.ndarray,
    unit_vx: np.ndarray,
    unit_st: np.ndarray,
    prefix: str = "",
) -> None:
    """Save unit_{sites,3d_coords,voxels} as both .npy and .mat."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stub = out_dir / (prefix + "unit_sites")
    np.save(str(stub) + ".npy", unit_st)
    savemat(str(stub) + ".mat", {"unit_sites": unit_st})
    stub = out_dir / (prefix + "unit_3d_coords")
    np.save(str(stub) + ".npy", unit_3d)
    savemat(str(stub) + ".mat", {"unit_3d_coords": unit_3d})
    stub = out_dir / (prefix + "unit_voxels")
    np.save(str(stub) + ".npy", unit_vx)
    savemat(str(stub) + ".mat", {"unit_voxels": unit_vx})


def run_session(
    session_dir: Path,
    out_dir: Path | None = None,
    hemisphere_backwards: bool | None = None,
    recording_side: str | None = None,
    reverse_xpos: bool = False,
    pkl_files: Sequence[Path] | None = None,
    xpos_mat: Path | None = None,
    ypos_mat: Path | None = None,
    output_prefix: str = "",
) -> dict:
    """End-to-end mapping for one session. Returns a summary dict.

    Parameters
    ----------
    hemisphere_backwards : bool or None
        If None, inferred from `recording_side` vs pkl filename sides.
        If True/False, used verbatim.
    recording_side : 'L' or 'R' or None
        Physical recording hemisphere (for inferring `hemisphere_backwards`).
    reverse_xpos : bool
        If True, replace xpos with max(xpos) - xpos before mapping. Needed
        when the xpos source (e.g. MATLAB's apFromFrontProbePerChannel) uses
        the opposite AP convention from the HERBS pkl shank numbering (1..4 A..P).
    """
    session_dir = Path(session_dir)
    if pkl_files is None:
        pkl_files = find_probe_pkls(session_dir)
    if len(pkl_files) == 0:
        raise FileNotFoundError(f"No probe_*.pkl in {session_dir}")
    if xpos_mat is None or ypos_mat is None:
        xpos_mat, ypos_mat = find_xy_mats(session_dir)
    if out_dir is None:
        out_dir = session_dir

    # Resolve hemisphere_backwards
    pkl_side = detect_pkl_side(pkl_files)
    if hemisphere_backwards is None:
        if recording_side is None:
            raise ValueError(
                "Must pass either hemisphere_backwards or recording_side "
                "(so it can be inferred from pkl filename sides)."
            )
        hemisphere_backwards = infer_hemisphere_backwards(pkl_files, recording_side)

    xpos = np.asarray(loadmat(str(xpos_mat))["xpos"]).flatten()
    ypos = np.asarray(loadmat(str(ypos_mat))["ypos"]).flatten()
    if reverse_xpos:
        xpos = xpos.max() - xpos

    site_locs, rel_site_locs, voxels, ccf_sites = load_pkls(pkl_files)
    if hemisphere_backwards:
        flip_hemisphere(site_locs, voxels)
    normalize_rel_site_locs(rel_site_locs, n_shanks=len(rel_site_locs))

    unit_3d, unit_vx_raw, unit_st = map_units(
        xpos, ypos, site_locs, rel_site_locs, voxels, ccf_sites
    )
    unit_vx = voxels_to_ccf_axes(unit_vx_raw)

    save_outputs(out_dir, unit_3d, unit_vx, unit_st, prefix=output_prefix)

    return dict(
        session_dir=str(session_dir),
        out_dir=str(out_dir),
        n_pkls=len(pkl_files),
        n_units=int(len(xpos)),
        n_mapped=int((unit_st != 0).sum()),
        xpos_mat=str(xpos_mat),
        ypos_mat=str(ypos_mat),
        pkl_side=pkl_side,
        recording_side=recording_side,
        hemisphere_backwards=hemisphere_backwards,
        reverse_xpos=reverse_xpos,
    )
