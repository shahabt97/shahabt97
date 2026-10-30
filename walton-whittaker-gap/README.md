# Spectral gap [9/8, 23/20] for the clusters of the charged Walton–Whittaker tiling

`PROOF.md` contains the proof (two counting lemmas + a finite certificate checked in exact rational arithmetic).

```
python -m ww.verify              # proof: cells (closure under inflation), exact check of both certificates, ~35 s
python -m ww.counts              # exact cluster sizes, recurrence, gap labels
python -m tests.equivariance     # the window map commutes with rotations / species flip
python -m tests.crosscheck 6 5   # independent brute-force check on actual clusters (dense eigensolver is NOT part of the proof)
```

* `ww/`      tiling rules, windows, closure, cells, the two forms, exact verifier, counting
* `data/`    rational certificates (896 + 1028 positive definite pieces), gzip JSON
* `tools/`   how the certificates were produced (constraints → SDP with cvxpy/SCS → rounding); not needed to verify
* `tests/`   independent cross-checks (numpy brute-force inflation of the actual clusters)

Requirements: Python ≥ 3.9, numpy (verification); cvxpy + scs (regenerating certificates); sympy (optional, irreducibility check).
Run all commands from this directory.
