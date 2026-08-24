# A Symmetry-Quotient Reduction for Odd-Characteristic SNOVA: Exact Structure and Conditional Recovery Ledgers

This directory contains the manuscript, source, and reproducibility artifact.
The paper specifies exact common-column reductions and conditional recovery
routes. It does not claim a production-size attack run or a complete-gate
security crossing.

## Requirements

- Python 3.10 or newer
- NumPy for the verifier-correspondence and reduced-parameter tests
- `latexmk`, BibTeX, and a TeX distribution with `natbib`, `cleveref`,
  `mathtools`, and `tabularx`

## Build and verify

```sh
make verify
make
```

`make verify` runs four separately identified classes of checks:

- the final 30-row PXL ledger, including the six `ell2_supported`, twelve
  `ell4_supported`, and twelve `ell4_H_struct` estimates across both matrix
  exponents;
- formula, target-sampling, dimension, and repair-floor consistency checks;
- exhaustive `F_19` scalar-netlist and quadratic-root checks, plus tower
  identity/recurrence checks;
- implementation-composition checks: an official Level-I KAT
  verifier/reduction harness and an unofficial reduced-parameter end-to-end
  forgery.

The KAT harness expands the public key from `pk`, validates the supplied KAT
signature and all 80 hash coordinates, checks the common-column decomposition,
and obtains a rank-50 residual system in 102 variables. It does not solve that
system and does not output a forgery.

The reduced test uses `q=19` and the Version 2.3 indexed public XOF, fixed ABQ,
message prehash, packing, public-map, signature-layout, and rejection
conventions at the unofficial shape `(v,o,ell,r)=(2,1,2,2)`. It imposes the
zero-offset relation `X_i=[u_i|0]`, interpolates the complete restricted map,
checks that its image equals the explicit rank-three quotient image, filters
one target-consistency coordinate before enumeration, and searches public
three-dimensional slices. The recomputed main case and eight fresh-key cases
pass positive and negative checks; a separate stress checker reruns 24 fresh
forgeries and a 100-key rank census. Verification uses a distinct literal path
in the shared KAT-anchored Python transcription, not the upstream C verifier.
This validates composition at reduced size, not independent verifier
correspondence or production feasibility.

## Regenerate deterministic ledgers

```sh
make regenerate
make verify
```

The final PXL artifact is `artifact/pxl_final_ledger.py` with pinned output in
`artifact/pxl_final_ledger.json`. Its checker recomputes all formulas, route
premises, trial-probability bounds, five PXL components, target-charge
convention, work estimates, and storage proxies. The `2^32` retained-target
charge is a stipulated model convention, not a measured or circuit-level
bound.

The formula checkers are internal consistency tools. They do not establish
the ordinary PXL regularity premise, the instancewise extraction premise, the
structured `H_struct` premise, the official frequency of `H_off`, the
idealized XOF-transcript model, or fixed-key spectrum events. The artifact
contains no production-size nonlinear solve or forgery.

Legacy homotopy and Just Guess files remain in `artifact/` as research
records. They are not part of the manuscript's current recovery headline.

The root `SHA256SUMS` file is the release manifest. It covers the files listed
in that manifest. It is not an inventory of every archived research file in
this directory. Verify the listed files with `sha256sum -c SHA256SUMS`.

`REVISION_RESPONSE_2026-07-31.md` is a revision record, not independent
evidence.
