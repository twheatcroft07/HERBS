"""CLI driver for herbs_map — map NP units to CCF coords for one or more sessions.

Usage:
    python run_map.py --session /path/to/session [--hemisphere-backwards] [--out /path]
    python run_map.py --sessions-glob '/n/data1/.../ANIMALS/M335/*/2025*_M335_g0'
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import herbs_map


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Map NP units to CCF coords (HERBS pkl inputs).")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--session", type=Path, help="One session directory.")
    src.add_argument(
        "--sessions-glob",
        type=str,
        help="Shell glob for multiple session directories (will be globbed here).",
    )
    p.add_argument(
        "--recording-side",
        choices=("L", "R"),
        default=None,
        help="Physical recording hemisphere. Used to auto-infer --hemisphere-backwards "
             "by comparing against the pkl filename side (probe_N_L/R_*.pkl).",
    )
    p.add_argument(
        "--hemisphere-backwards",
        choices=("auto", "true", "false"),
        default="auto",
        help="Override the CCF-X flip. 'auto' (default) derives it from --recording-side.",
    )
    p.add_argument(
        "--reverse-xpos",
        action="store_true",
        help="Mirror xpos before mapping (max-x - x). Set when xpos comes from MATLAB's "
             "apFromFrontProbePerChannel — it uses the opposite AP convention to HERBS.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Override output directory (default: save next to inputs).",
    )
    p.add_argument(
        "--prefix",
        type=str,
        default="",
        help="Prefix for output filenames (default: '' -> unit_sites.npy etc).",
    )
    p.add_argument(
        "--continue-on-error",
        action="store_true",
        help="In multi-session mode, log errors and keep going.",
    )
    return p.parse_args(argv)


def _resolve_hemisphere(arg: str) -> bool | None:
    return {"auto": None, "true": True, "false": False}[arg]


def run_one(session: Path, args) -> dict:
    summary = herbs_map.run_session(
        session_dir=session,
        out_dir=args.out,
        hemisphere_backwards=_resolve_hemisphere(args.hemisphere_backwards),
        recording_side=args.recording_side,
        reverse_xpos=args.reverse_xpos,
        output_prefix=args.prefix,
    )
    print(
        f"[OK] {session.name}: "
        f"{summary['n_units']} units, {summary['n_mapped']} mapped  "
        f"(hemisphere_backwards={summary['hemisphere_backwards']}, "
        f"pkl_side={summary['pkl_side']}, reverse_xpos={summary['reverse_xpos']}) "
        f"-> {summary['out_dir']}"
    )
    return summary


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.session is not None:
        sessions = [args.session]
    else:
        sessions = sorted(Path(p) for p in glob.glob(args.sessions_glob))
        if not sessions:
            print(f"No sessions matched: {args.sessions_glob}", file=sys.stderr)
            return 2
        print(f"Matched {len(sessions)} sessions")

    n_ok, n_fail = 0, 0
    for s in sessions:
        try:
            run_one(s, args)
            n_ok += 1
        except Exception as e:
            print(f"[FAIL] {s}: {type(e).__name__}: {e}", file=sys.stderr)
            n_fail += 1
            if not args.continue_on_error:
                return 1

    if n_fail:
        print(f"\nDone with errors: {n_ok} ok, {n_fail} failed", file=sys.stderr)
        return 1
    print(f"\nDone: {n_ok} session(s) mapped.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
