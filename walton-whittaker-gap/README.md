# Spectral gap [9/8, 23/20] for the clusters of the charged Walton–Whittaker tiling

* `FULL_PROOF.md` — the complete, self-contained proof (statements, lemmas, verbatim code, outputs, assessment). Start here.
* `PROOF.md` — shorter summary of the proof (two counting lemmas + a finite certificate checked in exact rational arithmetic).

```
python -m ww.verify              # proof: cells (closure under inflation), exact check of both certificates, ~35 s
python -m ww.verify --all-seeds  # the same, also for the chains that start next to holes of the tiling
python -m ww.counts              # exact cluster sizes, recurrence, gap labels
python -m tools.algebra          # exact algebra of the gap labels (factorization, irreducibility, Sturm, digits)
python -m tests.equivariance     # the window map commutes with rotations / species flip
python -m tests.crosscheck 6 5   # independent brute-force check on actual clusters (dense eigensolver is NOT part of the proof)
python -m tools.build_full_proof # rebuild FULL_PROOF.md from docs/FULL_PROOF.template.md, the sources and logs/
```

* `ww/`      tiling rules, windows, closure, cells, the two forms, exact verifier, counting
* `data/`    rational certificates (896 + 1028 positive definite pieces), gzip JSON
* `tools/`   how the certificates were produced (constraints → SDP with cvxpy/SCS → rounding), exact algebra, document builder; not needed to verify
* `tests/`   independent cross-checks (numpy brute-force inflation of the actual clusters), census, sizes, structure
* `logs/`    outputs of all programs, quoted verbatim in `FULL_PROOF.md`
* `docs/`    template of `FULL_PROOF.md`

Requirements: Python ≥ 3.9, numpy (verification); scipy + cvxpy + scs (regenerating certificates); sympy (optional, irreducibility check in `ww.counts`).
Run all commands from this directory.
