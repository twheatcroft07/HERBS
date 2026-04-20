"""Map every data1 session that has pkls + xpos + ypos, using dec24data.csv for recording_side.

Looks up recording_side from /home/thw051/CODE/MATLAB/BEHAVIOUR/dec24data.csv (Sabatini lab),
falling back to alreadyCuratedData.csv for ServerInterface/Neuropixels20240723 sessions.
Always runs with reverse_xpos=True because our MATLAB-derived xpos uses the opposite AP
convention from the HERBS pkl shank numbering.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, "/home/thw051/CODE/PYTHON/HERBS")
import herbs_map as hm

DATA1 = Path("/n/data1/hms/neurobio/sabatini/Tom")
CSVS = [
    Path("/home/thw051/CODE/MATLAB/BEHAVIOUR/dec24data.csv"),
    Path("/home/thw051/CODE/MATLAB/BEHAVIOUR/alreadyCuratedData.csv"),
]

label_re = re.compile(r"(\d{8})_M?(\d+)_g(\d+)")
path_label_re = re.compile(r"(\d{8}_M?\d+_g\d+)")


def load_recording_sides() -> dict:
    """Build {session_label: 'L' or 'R'} from the Sabatini lab CSVs."""
    sides = {}
    for csv_path in CSVS:
        if not csv_path.exists():
            continue
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sp = row.get("ServerPath") or row.get("serverPath") or ""
                m = path_label_re.search(sp)
                side = (row.get("RecordingSide") or row.get("recordingSide") or "").strip().upper()
                if m and side in ("L", "R"):
                    sides[m.group(1)] = side
    return sides


def find_ready_sessions() -> list[Path]:
    """Session dirs under data1 that have: 4 pkls AND xpos.mat AND ypos.mat."""
    candidates = []
    for root in (DATA1 / "ANIMALS", DATA1 / "ServerInterface" / "Neuropixels20240723"):
        for p in root.rglob("*_g*/"):
            if not p.is_dir():
                continue
            if not path_label_re.fullmatch(p.name):
                continue
            pkls = list(p.glob("probe_*.pkl"))
            xpos = list(p.glob("*_xpos.mat"))
            ypos = list(p.glob("*_ypos.mat"))
            if len(pkls) >= 4 and xpos and ypos:
                candidates.append(p)
    return sorted(set(candidates))


def main():
    sides = load_recording_sides()
    sessions = find_ready_sessions()
    print(f"Found {len(sessions)} ready sessions. Mapping with reverse_xpos=True.")
    n_ok, n_skip, n_fail = 0, 0, 0
    for s in sessions:
        side = sides.get(s.name)
        if side is None:
            print(f"[SKIP] {s.name}: no recording_side in CSVs")
            n_skip += 1
            continue
        try:
            summary = hm.run_session(
                session_dir=s,
                recording_side=side,
                reverse_xpos=True,
            )
            print(
                f"[OK]   {s.name}: {summary['n_units']} units, "
                f"hb={summary['hemisphere_backwards']}, pkl_side={summary['pkl_side']}, "
                f"side={side}"
            )
            n_ok += 1
        except Exception as e:
            print(f"[FAIL] {s.name}: {type(e).__name__}: {e}")
            n_fail += 1
    print(f"\nSummary: {n_ok} ok, {n_skip} skip, {n_fail} fail (of {len(sessions)})")


if __name__ == "__main__":
    main()
