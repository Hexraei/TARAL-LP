# Unit commitment

Deterministic textbook-scale MILP illustration, not an industrial deployment.

The supplied measurement was produced at solver commit
`c4800e9944587f221cfe25b9aebfa2a31750ff8b`. The stored result records TARAL
optimal status and exact objective agreement with HiGHS. These files do not
include a solution vector or independent row-check ledger. Runtime is a
single small-case measurement, not a speed comparison.

To reproduce, build TARAL, install `highspy`, then run from this directory:

```sh
TARAL_ENGINE=/absolute/path/to/taral python3 unit_commitment.py
```

The script writes the MPS, TARAL JSON and solution, and a HiGHS comparison
result. It regenerates the checked-in result JSON; keep separate copies when
comparing revisions.
