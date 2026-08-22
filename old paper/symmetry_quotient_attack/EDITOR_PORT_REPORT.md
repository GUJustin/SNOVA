# Editor-port and final-audit report

The editor's Git history was diffed against the prior attack manuscript and
its applicable corrections were ported into the later referee revision.

## Validated editor corrections

- fixed outer maps are conditioned upon, not assigned random rank density;
- projective quotient-dual counting in the affine-rank union bound;
- exact public quotient/affine/coefficient-map preflights;
- coordinate-aligned reserved-line freshness;
- finite-cardinality homotopy start and separator proof;
- exact target rejection sampling;
- explicit `ell=2` quotient-rank preflight;
- adaptive-mixture table terminology;
- rank-deficient three-plus-one handling;
- bibliography placement and metadata.

## Additional defects found by the final audit

- The full `A^v` zero-offset repair theorem from the editor was accidentally
  left behind when the newer manuscript was assembled. It is restored.
- The later Just Guess section incorrectly asserted independence of Frobenius
  channel coefficients under byte-biased native sampling. It is replaced by
  an exact conditional-atom theorem and regenerated ledgers.
- A stale sentence still referred to normalization by a structural-preflight
  density. It now states that the ledger is conditional on the exact
  structural preflight.
- The checksum manifest was stale and included cache files. It is regenerated
  only after the final clean build.

See `FINAL_CORRECTNESS_AUDIT.md` for the complete audit and residual limits.
