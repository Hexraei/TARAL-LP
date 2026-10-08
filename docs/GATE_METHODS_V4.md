# Gate methods: near-integrality repair candidate, pinned paired gate v4 (EXECUTED)

Refresh of the stale v3-era GATE_METHODS text ("NOT EXECUTED/cleared", candidate hash
9043600e). This text describes the v4 gate as actually executed; it supersedes the v3
provenance text for this gate. The v3 candidate and its partial gate are separate records.

## Builds and provenance
- Baseline source: c4800e9944587f221cfe25b9aebfa2a31750ff8b.
- Candidate: trial-milp.cpp sha256 532f722dec4a456bfe57ac4bc7c98791823b6d8ec2a0a3623724846a7da5023c,
  replacing src/milp.cpp ONLY atop that baseline.
- Binaries (auditor-verified rebuild targets):
  build-base/taral sha256 8ce643231c7683a31d0b4cc8c6375e55ffbb1969cb2a36ab07ff5c5b2ed1d451;
  build-grid/taral sha256 982470cebbf24564b340253f67cdfaa8dc9a8a333dee9d620b434c60c3f9a460.
  "build-grid" is a compatibility directory name; the actual arm is the near-integrality
  candidate with the objective lattice disabled (TARAL_GENERIC_GRID unset).
- Final ledger: fullpair/ledger_300s.json sha256 ca84574ee8c2cdece57228ad50aba6a77367225992fa87ebefdce9fb3ef53221.
- Producer audit: nearint-v4-original-audit.json sha256 331d51bfe8f65d3882275b287f03b7c2194d9960fed923e84a34179f697bc57d
  (separate from the raw archive). Status: auditor-inspected; the reviewer confirmed its
  hash and reconciled all 130 rows and all 106 recomputed point objectives against its
  own independent audit - summaries match.
- Raw archive root: nearint-v4-gate/ sha256 4a7dac8752075f1b1dc2f7f62d1e9222d09bc2d6a8e30926ccc22e4813247c70
  (contains compiler/flags.make/host/environment/run.sh records).
- Completion markers: COMPLETE=PAIRED300_COMPLETE, RUN_EXIT=0.

## Host and toolchain
Azure VM (SKU producer-reported as Standard_D4s_v4; not established by the archive host
record): 4 vCPU Intel Xeon Platinum 8272CL 2.60GHz, Ubuntu 22.04,
kernel 6.8.0-1064-azure, G++ 11.4.0, CMake 3.22. Both arms built with flags.make
CXX_FLAGS = -O3 -DNDEBUG and NO -march=native (portable build). FLAGS.txt declares C++17
but flags.make carries no explicit -std=c++17; the compiler default applied - this
declaration is not promoted to a command flag here. environment.txt: LANG=C.UTF-8,
LC_ALL empty.

## Protocol
Paired arms on identical models: TARAL_NO_FJ=1 in every subprocess (FJ off both arms),
TARAL_GENERIC_GRID unset, --persistent-nodes, --node-limit 1000000, --time-limit 300
serial per case, per-case --json/--sol retained. All 65 MIPLIB 3 cases, both arms,
130/130 rows. Resume validates binary/model SHA before skipping. Feasibility 1e-6
numerical gate, not exact-arithmetic proof. Polishing LP work is counted in
iterations/wall, not the B&B node limit.

## References and audit
Cached reference provenance retained: 43 Optimal / 22 timecap references; unknown
references are never called established; misc06 uses the tight mip_rel_gap=0 reference
(reference highspy 1.15.1). Independent original-model point/bound checks per row
(runner audit highspy 1.12.0, numpy 2.2.6; producer independent reader 1.15.1).
Objective/bound missing or nonfinite is classified; absent solution, unknown reference
or missing bound produces a null check, never a certified pass. Summary statuses report
valid/invalid/absent points, valid/invalid/unknown bounds, and optimum agreement only
on Optimal rows with known references.

## Verified headline facts (from the final ledger, sha ca84574e)
130/130 rows, 65 cases x 2 arms; 0 runner errors; 0 numerical_failure rows on either
arm; status tallies identical per arm (34 optimal / 28 time_limit / 3 node_limit);
audit totals across both arms: 106 valid points / 24 absent; 86 reference-checked
bounds / 44 unestablished; 56 optimum-agreement rows (the remaining 74 rows have no
established optimum-agreement comparison); binary shas and source pin c4800e9944587f221cfe25b9aebfa2a31750ff8b consistent across
all rows. Independent review: accepted for the pinned portable-Release, FJ-off gate only, with
the caveats below. Acceptance is not a universal optimum proof: optimum agreement is
established on 56 rows only; 74 rows carry no established optimum-agreement comparison.

## Caveats (preserved verbatim in scope)
- stein45 agrees with its cached incumbent only (reference time-limited, bound 27,
  gap 0.10): no independent optimum agreement; stein45 stays out of optimum-agreement
  claims.
- Capped outputs are not byte-comparable.
- Native-vs-portable floating-point sensitivity is explicit: these figures come from a
  portable (-O3 -DNDEBUG, no -march=native) build; native-flag builds may differ.
- This gate ran FJ off; it says nothing about FJ-on behavior.
- The pinned binary does not test current main (main has since diverged in
  src/presolve.cpp and src/taral.hpp); the integrated current-main candidate requires
  its own post-integration verification before any merge of the fix into main.
- No speed claims. No current-main claims. Known-defect caveat on mixed-integer
  correctness wording stays until the integrated candidate lands and re-passes.
