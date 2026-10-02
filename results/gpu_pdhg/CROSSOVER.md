# PDHG-to-simplex crossover prototype (local, measured Oct 2; prototype, not a headline claim)

Handoff: `gpu/pdhg.cu` writes an approximate primal point and row multipliers; `taral --warm-sol X.sol --warm-dual X.dual
--cross-tol T` (src/crossover.cpp) turns them into a simplex basis and the existing warm-start path finishes exactly.
Basis rule: the m variables (structurals and row logicals) farthest from their bounds are basic; slots left over when fewer
than m are farther than T from a bound go to the smallest reduced costs |c - A'y|; everything else sits at its nearest bound.
Hardware: RTX 5050 Laptop (PDHG), Intel Core 7 240H (simplex). Default engine path is unchanged without the flags.

## dfl001 (6,071 rows; the exact primal simplex does not finish in 60 s)
| Route | PDHG | Simplex after handoff | Total | Simplex iterations |
|---|---|---|---|---|
| Cold dual simplex (no handoff) | - | 11.8 s | **11.8 s** | 21,256 |
| PDHG tol 1e-4 + dual simplex (cross-tol 1e-5) | 0.47 s (1.9 s on a first run) | 12.7 s | 13.2 s | 23,653 |
| PDHG tol 1e-6 + dual simplex | 0.68 s | 17.8 s | 18.5 s | 29,869 |
| PDHG tol 1e-8 + dual simplex | 3.7 s | 18.7 s | 22.4 s | 30,532 |
| PDHG tol 1e-6 + primal simplex | 0.68 s | not done in 60 s | - | - |
| PDHG tol 1e-6, primal point only (no multipliers), dual or primal | 0.68 s | not done in 60 s | - | - |
All handoff routes that finished reached the same objective as the cold dual, 11266396.0467 (HiGHS agrees).

Finding: the handoff does **not** pay on dfl001. The PDHG point is not near a vertex: only about 3,500 of the 6,071 basis slots
are determined by variables clearly away from a bound, and the rest is guesswork from reduced costs, so the simplex still
needs about as many iterations as from the slack basis (23.7k against 21.3k) and the PDHG time is added on top. A tighter
PDHG tolerance makes the basis worse here, not better. (The cold dual used to take 31.8 s on the Kaggle T4 host and 14-15 s
here before the LU speed-up; it is 11.8 s here now, which is why the PDHG advantage of 0.5 s approximate versus 12 s exact
does not translate into an exact-answer speedup.) PDHG alone remains the fast route to a near-optimal dfl001
(objective error 4e-5 at 1e-4, 5e-7 at 1e-6 against HiGHS), separate protocol line.

## packing_100000 (50,000 rows, 100,000 columns, 500k nonzeros; synthetic, gpu/gen_lp.py)
- PDHG tol 1e-6: 2.3 s, objective -68881.0220 against HiGHS (interior point) -68881.0219 in 2,331 s.
- Cold dual simplex: no answer in 120 s (34,656 iterations).
- Handoff: the warm basis was still inside its **initial LU factorization after 11.5 minutes (1.9 GB resident)** and was killed.
  The 50,000-row basis built from a PDHG point has a dense nucleus. The factorization has no time check, so the 120 s limit
  was not honoured: a known defect of the prototype path (a basis that cannot be factored in reasonable time should fall back
  to the slack basis; not done).
Conclusion: no crossover speed-up demonstrated on either case. This is the honest outcome; a usable crossover needs a
real basis-identification step (push/pull to a vertex), not a nearest-to-bound crash.
