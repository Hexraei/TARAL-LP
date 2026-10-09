# Exact VALUES nonconvexity witness

VALUES is refused by design; an exact feasible negative-curvature witness verifies nonconvexity of the encoded benchmark as published.

This proof concerns the pinned decimal model, not an inferred intended model. It does not establish an optimum, a reference objective, or a benchmark pass.

## Pinned model

Source: https://raw.githubusercontent.com/YimingYAN/QP-Test-Problems/871cc300154cfb438acc20f52b2c585db3f11e31/QPS_Files/VALUES.QPS

Model SHA256: `3589091f3d73ed1680d2c06ed79b3c8f6d5959667bf07bb182c9246f90c80398`.

`VALUES.mps` contains exactly those source bytes. The model has 202 columns, one equality, and bounds [0,10]. QUADOBJ stores one symmetric Hessian triangle, with objective `c'x + 0.5 x'Qx`. The equality is `sum(as1..as101) - sum(a1..a101) = 0`.

## Exact proof

`exact_certificate.json` gives the column-ordered integer direction d. Exact rational arithmetic verifies:

- `Ad = 0`.
- `d'Qd = -795901283767/62500 < 0`.
- The midpoint x has every coordinate 5. For `t = 1/205859`, both `x+t*d` and `x-t*d` satisfy the equality and have every coordinate in [4,6], hence satisfy the original bounds.
- `f(x+t*d) + f(x-t*d) - 2*f(x) = -795901283767/2648620492562500 < 0`.

The negative second difference between feasible endpoints contradicts convexity on the original feasible set. The result uses the literal decimal coefficients. A possible coefficient-rounding or truncation explanation for the source data remains unverified.

## Replay

From this directory, using Python 3.10 or newer:

```sh
sha256sum -c SHA256SUMS
python3 replay_exact.py
```

The replay uses only the Python standard library, including `fractions.Fraction`. It checks the model hash, column order, equality direction, negative curvature and feasible endpoints. No eigensolver, HiGHS, NumPy, SciPy or solver run is needed.

TARAL's recorded behavior on main `4cd9f93844e3530c641b2cdeb1d9b032323a8626` was a nonconvex refusal with zero iterations. That behavior is separate from the exact proof and does not supply an optimality claim.
