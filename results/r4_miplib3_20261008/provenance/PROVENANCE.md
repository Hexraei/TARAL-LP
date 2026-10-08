# R4 MIPLIB3 run provenance (Oct 8, 2026)

- Runner: m3run_vm_v3.py (protocol exactly as run on the VM; sha256 e4f9589c58d87f0031e8e81b8805f8d7ddc6b916bf10012e6e73163a83d5f896).
  NOTE: v3's classifier emits the v2 label "incumbent_below_ref"; the v3 rename to the neutral
  "incumbent_differs_from_ref" is label-only. The ledger carries the old label.
  Launch-mechanics repairs were applied by the engine on the VM side (self-containment, source
  archive pin of c4800e9, credential handling) per the Oct 8 00:17 review; protocol, classifier
  and verdicts are this file's.
- Runner v4: m3run_vm_v4.py (sha256 ce5bdb5992ff96c3e37886b947072943be8718c80ec863acc6127383a4467707):
  identical except mip_rel_gap=0 pinned on the HiGHS reference (the misc06 false-WRONG fix,
  engine-confirmed). Use v4 for any rerun.
- Build: c4800e9 (build_commit.txt), engine's verified source archive.
- Instance list: list.txt (65 MIPLIB 3 instances; index minus mps_format/references).
- source_manifest.sha256: sha256 of the c4800e9 engine SOURCE files as built on the VM (12 src/
  files, from the final archive). Per-MODEL sha256 for all 65 instances is the 'sha256' column of
  the merged ledger (ledger_c4800e9_merged.csv in the final archive). Engine separately recovered
  all 65 raw models and hash-matched them; air04 differs from the engine's older cache only in the
  MPS comment line (BEST SOLN 56137 vs 56138), parsed model field-identical.
- Protocol: taral --time-limit 300 --json --sol; HiGHS 1.15.1, 1 thread, 300 s. Resume-safe
  sharded runner (s0/s1). References: ref_highs_s0/s1.json in the final archive.
- Final archive (with parent): r4_results tar.gz sha256 a1d647958a85d7314e90e3fdedea781a6eab5d04f8bedd785eef58c60190e6cd
  (668748 bytes): merged + per-shard ledgers, raw .sol/.json for all 65, point-audit JSON,
  tight misc06 reference JSON, run logs.
