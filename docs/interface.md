# Command-line and C++ interfaces

Measured scope: the [31-case command-line suite](../results/interface_20261003/cli_test.log)
passes on the evaluated solver, and [one C++ LP example](../examples/api_lp.cpp)
builds warning-free and returns the published steel-production optimum 192000.
[Hashes and test conditions](../results/interface_20261003/provenance.json).
This demonstrates a working basic interface, not a stable versioned API,
complete integration coverage or a graphical interface.

## Build and solve

```sh
g++ -O3 -std=c++17 -o taral src/*.cpp
./taral examples/literature_lp/steel_production.mps --time-limit 30 --sol /tmp/steel.sol --json /tmp/steel.json
./taral --help
```

Input is an MPS model file. JSON is the machine-readable result format; solution
files contain variable names and values. [Captured help](../results/interface_20261003/help.txt)
is the exact reference for this evaluated source snapshot.

| Option | Meaning |
| --- | --- |
| `--time-limit S` | Wall-time budget, positive finite seconds; default 60. Parsing/setup can overrun it; not a hard operating-system deadline. |
| `--node-limit N` | Nonnegative integer branch-and-bound node limit; default unlimited. |
| `--method simplex\|dual\|ipm` | Default primal simplex, dual simplex or interior-point method for continuous models. Integer models use mixed-integer search; convex QP uses its quadratic method. This flag is not a pure-method MILP comparison. |
| `--sol FILE`, `--json FILE` | Solution and result-summary destinations. Do not interpret a file's existence as success; inspect status and verify the point. |
| `--warm-sol FILE`, `--warm-dual FILE`, `--cross-tol T` | Experimental reuse of an approximate point and optional row multipliers to construct a simplex starting basis. |
| `--no-fallback` | Disable switching a stalled default LP solve from primal to dual simplex. |
| `--audit-prop`, `--no-prop-prune` | Mixed-integer search diagnostic controls. |

Exit codes: 0 definitive status (optimal, infeasible or unbounded), 2 usage
error, 3 parse error, 4 time/iteration/node limit, 5 other failure. A definitive
status is still subject to the stated floating-point validation rules. A
mixed-integer time limit can carry a feasible incumbent without proving it
optimal; inspect `has_solution`, `best_bound` and `gap`. LP certificate quality
is separately reported and can fail the tighter check while meeting ordinary
acceptance. See [quality definitions](simplex_certificate_quality.md).

The 31-case suite covers help, valid options/methods, missing model, invalid
limits/methods, missing output arguments, unknown/extra arguments and `--`.
It is not comprehensive parser fuzzing, output-I/O failure coverage or every
combination of interface options.

## C++ interface

[`src/taral.hpp`](../src/taral.hpp) defines the model, sparse coefficient entries,
LP results and mixed-integer options/results. `read_mps(path)` throws
`ParseError`; `solve(model, seconds)` solves an LP. `solve_lp` and
`solve_lp_dual` accept effective bounds and an optional previous basis.
`solve_milp` is the mixed-integer entry point. The continuous QP interface is
in [`src/ipm.hpp`](../src/ipm.hpp). Select the entry point for the model class;
the minimal example below deliberately rejects integer/quadratic models.

```sh
g++ -O3 -std=c++17 -Isrc examples/api_lp.cpp $(find src -name '*.cpp' ! -name main.cpp | sort) -o /tmp/api-lp
/tmp/api-lp examples/literature_lp/steel_production.mps
```

The measured output is `optimal 192000`. Model/result ownership is by C++ value
and standard-library containers; no binary-compatibility guarantee is made.
