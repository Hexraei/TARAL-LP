Independent second-host run package (private). Docker build context = this directory. Build line and tests are the ones shipped here.
Contents: src/ (main 009fbf59, sanitized), tests/ + tools/adv (cloud/adversarial-tests 7d0af4d), gate/ (public sanitized 93-case gate, pinned netlib corpus + sha256 manifest), mit/ (rail2586, rail4284 solver-only runner), Dockerfile, run_all.sh.
Render: Background Worker, Docker, plan 4c-16g. Output = service logs; run ends with SECOND_HOST_DONE then idles. Delete the service after saving logs.
Suite exit 0 means only KNOWN failures (tests/known_failures.tsv), not that every case is correct.
