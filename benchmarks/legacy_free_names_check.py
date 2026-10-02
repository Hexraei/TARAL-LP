"""Regression check for the historical Python parser's free-variable split.

Two free columns sharing a 7-character prefix (ABCDEFG1, ABCDEFG2) plus a real column named like a
split clone (ABCDEFGP). Original optimum: x1 = 1, x2 = 2, ABCDEFGP = 0, objective 5. The parser must
keep four distinct split columns and map the solution back exactly; the independent original-model
check must then confirm feasibility and the objective.

Usage (repository root, NumPy and SciPy installed): python benchmarks/legacy_free_names_check.py
"""
import sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "engine"), str(ROOT / "parsers"), str(ROOT / "benchmarks")]
from mps_free import read_mps_extended
from revised_simplex import solve_lp
import orig_check as oc


def fixed(kind="", f1="", f2="", f3="", f4="", f5=""):
    line = " " + kind.ljust(2) + " " + f1.ljust(8) + "  " + f2.ljust(8) + "  " + f3.rjust(12)
    if f4:
        line += "   " + f4.ljust(8) + "  " + f5.rjust(12)
    return line.rstrip()


MODEL = "\n".join([
    "NAME          FRNAMES", "ROWS", fixed("N", "OBJ"), fixed("E", "R1"), fixed("E", "R2"), fixed("L", "R3"),
    "COLUMNS",
    fixed("", "ABCDEFG1", "OBJ", "1", "R1", "1"),
    fixed("", "ABCDEFG2", "OBJ", "2", "R2", "1"),
    fixed("", "ABCDEFGP", "OBJ", "1", "R3", "1"),
    "RHS", fixed("", "RHS", "R1", "1", "R2", "2"), fixed("", "RHS", "R3", "3"),
    "BOUNDS", fixed("FR", "BND", "ABCDEFG1"), fixed("FR", "BND", "ABCDEFG2"),
    "ENDATA", ""])


def main():
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "frnames.mps"
        path.write_text(MODEL)
        A, b, c, kinds, names, offset, lower = read_mps_extended(path)
        assert len(names) == 5, f"expected 2+2 split columns and 1 plain column, got {names}"
        r = solve_lp(A, b, c, kinds)
        x = {}
        for n, v in zip(names, list(r["x"] + lower)):
            if n.endswith("|FR+"): n = n[:-4]
            elif n.endswith("|FR-"): n, v = n[:-4], -v
            x[n] = x.get(n, 0.0) + v
        assert set(x) == {"ABCDEFG1", "ABCDEFG2", "ABCDEFGP"}, x
        chk = oc.check(oc.parse_mps(path), x)
        obj = r["objective"] + offset
        print({"objective": obj, "x": x, "original_check": chk})
        assert abs(obj - 5) < 1e-9 and abs(chk["objective"] - 5) < 1e-9
        assert abs(x["ABCDEFG1"] - 1) < 1e-9 and abs(x["ABCDEFG2"] - 2) < 1e-9 and abs(x["ABCDEFGP"]) < 1e-9
        assert chk["row_violation_rel"] <= 1e-9 and chk["bound_violation"] <= 1e-9
    print("PASS legacy free-name split")


if __name__ == "__main__":
    main()
