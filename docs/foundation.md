# Inspectable solver foundation and independent validation

Measured scope: the evaluated CPU solver builds warning-free using C++17 and
its standard library. Source/header inspection and the build definition do not
link HiGHS, SciPy, NumPy, CUDA or an external numerical solver into that engine.
This is a dependency observation, not proof of originality, ownership or solver
correctness. [Inspection evidence](../results/core_inspection_20261003/).
Compiler/runtime libraries appear in the executable's dynamic dependencies;
"standard-library-only" does not mean a binary with no runtime dependencies.

## Source layout and extension boundaries

| Module | Responsibility |
| --- | --- |
| `src/taral.hpp`, `src/ipm.hpp` | Model/result types and solver entry points |
| `src/mps.cpp` | Model input |
| `src/simplex.cpp`, `src/dual.cpp` | Continuous LP methods |
| `src/ipm.cpp` | Continuous convex quadratic and interior-point calculations |
| `src/milp.cpp` | Mixed-integer linear search |
| `src/lu.cpp`, `src/crossover.cpp` | Basis factorization and previous-point initialization |
| `src/kkt_gate.cpp`, `src/nonoptimal_certificate.cpp` | Original-model checks and LP proof checks |
| `src/main.cpp` | Command-line dispatch and result export |

The [CMake build](../CMakeLists.txt) enumerates the CPU translation units.
[The C++ example](../examples/api_lp.cpp) demonstrates using the LP interface
without the command-line entry point. This is inspectable modular structure,
not a measured guarantee that a new problem class can be added without work.
Mixed-integer quadratic/nonlinear extensions and stable API versioning remain
unfinished. [Interface scope and tests](interface.md).

## Independent checks, separate from engine success flags

- [`benchmarks/orig_check.py`](../benchmarks/orig_check.py) parses the original
  MPS independently of the C++ reader, maps returned variable values, and
  recomputes original constraints, bounds and objective.
- [`benchmarks/netlib_gate.py`](../benchmarks/netlib_gate.py) applies the recorded
  acceptance rules using HiGHS reference answers. The evaluated CPU
  [Netlib record](../results/netlib_postmerge_58f77af/) contains 93/93 accepted
  cases in each method, 92/93 at the tighter check. The
  [post-quadratic-change Kaggle confirmation](../results/qp_postmerge_gate_20261003/kaggle_second_host/)
  has unchanged per-case outcomes on that host.
- [`tools/advqp/qpcheck.py`](../tools/advqp/qpcheck.py) checks quadratic results
  using original-data feasibility, recomputed objective and independently
  computed optimality conditions. [Raw stress results and costs](qp_ipm_diagnosis.md)
  retain failures and unresolved outcomes.
- [`benchmarks/simplex_certificate_check.py`](../benchmarks/simplex_certificate_check.py)
  and [`benchmarks/nonoptimal_certificate_check.py`](../benchmarks/nonoptimal_certificate_check.py)
  verify LP evidence separately. These are floating-point checks, not exact
  rational proof, and do not cover general MILP/QP proof export.
- [Published application measurements](../examples/literature_lp/) retain full
  returned points, model hashes and independent reference checks.

The offset and negative-curvature diagnostic audits are not part of this
published measurement claim yet; their commit-pinned records and final
writeups remain pending. Reference software is used by test tools, not by the
CPU optimization core. Neither passing this selection nor having checkers
establishes full industrial reliability. [Limitations](limitations.md).

## Access is not an open-source license

The source is inspectable only within the [proprietary evaluation terms](../LICENSE).
The license permits official competition judges/organizers to view the
submission for evaluation and restricts copying, modification, redistribution,
commercial use and AI/ML use. This document changes no permission, repository
visibility or website state. "Extensible" describes module structure, not a
license grant to extend or reuse it.
