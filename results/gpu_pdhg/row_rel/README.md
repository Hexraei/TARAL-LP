# Evidence for gpu/SCALING_RESULTS.md (Oct 3, RTX 5050 laptop)

- `default_ledgers/`: 94-file sweeps at 1e-4 and 1e-6 with the default binary (`--row-rel 0`), same flags as `out/gpu/run_harness.sh`; compare with `../netlib_94files_tol1e-{4,6}.csv` (`tools/cmp.py OLD NEW`).
- `row_rel100_ledger/`: the 1e-4 sweep with `--row-rel 100` (shows the bnl1, perold, pilot.we regressions). The 1e-6 sweep at `--row-rel 100` was halted and is not included.
- `large_cases/`: one JSON line per run on pds-100, fome21, ken-18, osa-60, rail2586, rail4284 (`res_base` = old binary, `res_new` = `--row-rel 100`, the rest are the scaling/stop-rule variants, built from env-switch experiment copies that are NOT in the repo); `highs_reference.json` = HiGHS IPM objectives. Rows come from `tools/drv.py` (numpy recheck of the written x).
- The six large MPS files are not stored: pds-100, fome21, rail* from plato.asu.edu/ftp/lptestset, ken-18 and osa-60 from netlib lp/data/kennington (emps-expanded).
