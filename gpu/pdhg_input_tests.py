#!/usr/bin/env python3
"""CPU-test the PDHG input gate using the real MPS reader, without CUDA.

Run: python3 gpu/pdhg_input_tests.py [--pdhg /path/to/cuda-built/pdhg]
The optional binary runs rejection checks on both backends and precisions.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PROBE = r'''
#include "gpu/pdhg_input.hpp"
#include <cassert>
int main(int argc, char** argv) {
    Model m;
    assert(pdhg_input_error(m) == nullptr);
    m.is_int = {0, 0};
    assert(pdhg_input_error(m) == nullptr);
    m.is_int[1] = 1;
    assert(pdhg_input_error(m) != nullptr);
    m.qobj.push_back({0, 0, 1.0});
    assert(pdhg_input_error(m) != nullptr);
    m.is_int.clear();
    assert(pdhg_input_error(m) != nullptr);
    assert(argc == 2);
    return pdhg_reject_unsupported_input(read_mps(argv[1])) ? 2 : 0;
}
'''


def mps(name, bounds=" UP B X 2\n", quadratic="", marker=False, maximize=False):
    sense = "OBJSENSE\n MAX\n" if maximize else ""
    start = " M0 'MARKER' 'INTORG'\n" if marker else ""
    end = " M1 'MARKER' 'INTEND'\n" if marker else ""
    return (f"NAME {name}\n{ sense }ROWS\n N OBJ\n G R\nCOLUMNS\n"
            f"{start} X OBJ 1 R 1\n{end}RHS\n RHS R 0.5\n"
            f"BOUNDS\n{bounds}{quadratic}ENDATA\n")


def check_rejection(result, label, expected):
    if result.returncode != 2:
        raise AssertionError(f"{label}: expected exit 2, got {result.returncode}: {result.stderr}")
    if result.stdout.strip() != "status unsupported_model":
        raise AssertionError(f"{label}: unexpected stdout: {result.stdout}")
    if "continuous LPs only" not in result.stderr or any(s not in result.stderr for s in expected):
        raise AssertionError(f"{label}: missing diagnostic: {result.stderr}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdhg", type=Path)
    args = parser.parse_args()
    binary = args.pdhg.resolve() if args.pdhg else None
    cases = [
        ("continuous_lp", mps("LP"), ()),
        ("maximize_lp", mps("MAXLP", maximize=True), ()),
        ("zero_quadratic", mps("ZERO", quadratic="QUADOBJ\n X X 0\n"), ()),
        ("quadobj", mps("QP", quadratic="QUADOBJ\n X X 2\n"), ("quadratic",)),
        ("qmatrix", mps("QM", quadratic="QMATRIX\n X X 2\n"), ("quadratic",)),
        ("qsection", mps("QS", quadratic="QSECTION OBJ\n X X 2\n"), ("quadratic",)),
        ("maximize_qp", mps("MAXQP", quadratic="QUADOBJ\n X X -2\n", maximize=True), ("quadratic",)),
        ("intorg", mps("MILP", marker=True), ("integer",)),
        ("binary_bound", mps("BV", bounds=" BV B X\n"), ("integer",)),
        ("integer_lower", mps("LI", bounds=" LI B X 0\n UP B X 2\n"), ("integer",)),
        ("integer_upper", mps("UI", bounds=" UI B X 2\n"), ("integer",)),
        ("fixed_integer", mps("FIX", bounds=" FX B X 1\n", marker=True), ("integer",)),
        ("miqp", mps("MIQP", marker=True, quadratic="QUADOBJ\n X X 2\n"), ("quadratic", "integer")),
    ]
    with tempfile.TemporaryDirectory(prefix="pdhg-input-") as folder:
        folder = Path(folder)
        source, probe = folder / "probe.cpp", folder / "probe"
        source.write_text(PROBE)
        subprocess.run([os.environ.get("CXX", "g++"), "-O2", "-std=c++17", "-Wall", "-Wextra", "-pedantic",
                        "-I", str(ROOT), str(source), str(ROOT / "src/mps.cpp"), "-o", str(probe)], check=True)
        checked = 0
        for name, text, expected in cases:
            path = folder / (name + ".mps")
            path.write_text(text)
            result = subprocess.run([str(probe), str(path)], capture_output=True, text=True)
            if expected:
                check_rejection(result, name, expected)
            elif result.returncode != 0 or result.stdout or result.stderr:
                raise AssertionError(f"{name}: continuous LP rejected: {result}")
            checked += 1
            print(f"PASS host gate: {name}")
            if binary and expected:
                for device in ("cpu", "gpu"):
                    for fp32 in (False, True):
                        outputs = [folder / "answer.json", folder / "answer.sol", folder / "answer.dual"]
                        flags = ["--json", str(outputs[0]), "--sol", str(outputs[1]), "--dual", str(outputs[2])]
                        for existing in (False, True):
                            for output in outputs:
                                if existing:
                                    output.write_text("sentinel\n")
                                elif output.exists():
                                    output.unlink()
                            cmd = [str(binary), str(path), "--device", device, *flags]
                            if fp32:
                                cmd.append("--fp32")
                            result = subprocess.run(cmd, capture_output=True, text=True)
                            label = f"{name}/{device}/fp{32 if fp32 else 64}/existing={existing}"
                            check_rejection(result, label, expected)
                            if any((o.read_text() != "sentinel\n" if existing else o.exists()) for o in outputs):
                                raise AssertionError(f"{label}: output files modified")
                            checked += 1
                            print(f"PASS CLI: {label}")
        print(f"PASS: {checked} checks (plus direct Model assertions)")
        if not binary:
            print("CUDA-built CLI not supplied; CUDA compilation/runtime and CLI integration not tested.")


if __name__ == "__main__":
    main()
