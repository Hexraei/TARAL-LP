# Published application models: production and transportation

These two small linear programs complement the existing
[Williams refinery example](../refinery_williams.py). They demonstrate model
translation, solution output and independent checks, not industrial readiness
or real operational savings. They are textbook illustrations, not field data.

## Steel production planning

Source: Robert Fourer, David M. Gay and Brian W. Kernighan, *AMPL: A Modeling
Language for Mathematical Programming*, second edition (2003), chapter 1,
section 1.1, pp. 2-3.
[Publisher-hosted chapter](https://ampl.com/wp-content/uploads/Chapter-1-Production-Models-AMPL-Book.pdf).

Choose bands B and coils C (tons) to maximize 25 B + 30 C dollars of illustrative
weekly profit. Mill time: B/200 + C/140 <= 40 hours. Bounds: 0 <= B <= 6000,
0 <= C <= 4000. The chapter derives B=6000, C=1400 and profit 192000. The MPS
stores 1/140 as a double-precision coefficient; this is numerical reproduction,
not rational-arithmetic proof.

## Transportation / logistics

Source: George B. Dantzig, *Linear Programming and Extensions* (1963), chapter
3.3, as cited in the [GAMS transportation model](https://gams.com/latest/gamslib_ml/libhtml/gamslib_trnsport.html).
The model library also cites Richard E. Rosenthal's 1988 GAMS tutorial.

Choose nonnegative shipments from Seattle and San Diego to New York, Chicago
and Topeka. Plant capacities are 350 and 600 cases. Market demands are 325,
300 and 275 cases. Distances (thousands of miles) are [2.5,1.7,1.8] and
[2.5,1.8,1.4]. Freight is 90 dollars per case per thousand miles; the objective
uses 90*distance/1000, so its unit is thousands of dollars. Minimize shipment
cost, with each plant's total <= capacity and each market's total >= demand.

The independently measured optimum is 153.675 thousand dollars. One returned
shipment plan is Seattle [50,300,0], San Diego [275,0,275]. Other optimal plans
can exist; matching a particular point is not the acceptance rule.

## Reproduce and inspect

Build the current solver, install NumPy/HiGHS for the independent benchmark
checks, then run:

```sh
python examples/literature_lp/run.py --engine ./taral --out /tmp/literature-results
```

The runner retains MPS hashes, returned variable values, engine JSON and logs,
source/binary hashes, the separately parsed original-model checks and HiGHS
reference objectives. It requires optimal status, objective agreement and
row/bound checks at 1e-8. No general runtime advantage is claimed from tiny
millisecond cases. [Measured record](../../results/literature_lp_20261003/).

Application coverage remains partial: these are two LP case studies, not
mixed-integer scheduling, power dispatch or a full supply-chain implementation.
