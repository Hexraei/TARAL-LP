"""Run the experimental CPU LP solver on one user-provided MPS file.

This is a smoke-test runner, not a reproduction of the checked-in result ledgers.
"""
import argparse
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "engine"))
sys.path.insert(0, str(ROOT / "parsers"))
from revised_simplex import solve_lp
from mps_free import read_mps_extended


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mps", type=Path, help="path to an LP MPS model")
    args = parser.parse_args()
    A, b, c, kinds, names, offset, lower = read_mps_extended(args.mps)
    start = time.process_time()
    result = solve_lp(A, b, c, kinds)
    elapsed = time.process_time() - start
    print({
        "model": args.mps.name,
        "rows": len(b),
        "columns": len(names),
        "objective": result["objective"] + offset,
        "max_residual": result["max_residual"],
        "phase1_iterations": result["phase1_iterations"],
        "phase2_iterations": result["phase2_iterations"],
        "cpu_seconds": elapsed,
    })


if __name__ == "__main__":
    main()
