# The Symmetry-Quotient Attack on Odd-Characteristic SNOVA

This package contains the manuscript source and deterministic numerical
artifact.

## Reproduce the ledger

Python 3.10 or newer is sufficient:

```text
python3 artifact/generate_ledger.py
python3 artifact/verify_final_results.py
```

The generator writes `artifact/all_nine_ledger.json`. The verifier imports none
of the generator code: it independently transcribes the paper formulas,
recomputes all nine direct and adaptive rows, and then compares the results.
It also checks the scalar multiplier on all 361 canonical pairs, the SHA-256
digest, field-tower identities and determinant, target sampling, separator
choices, the numerical consequences of the conditional coefficient bounds,
dimension inequalities, spectrum-tail factors, and repair floors.

The checker does not prove the affine reduction, the random-XOF idealization,
the cited symbolic-homotopy theorem, the correspondence with the official
SNOVA verifier, the fixed-layout coefficient-rank hypotheses, or any upper
bound on `kappa_hom`. The package also omits the pinned implementation needed
to regenerate the fixed outer maps. It therefore establishes conditional
formula-level results for the nine parameter shapes, not fixed-key coverage
of every official key.

## Build the paper

A TeX distribution with `latexmk`, `biblatex`, and `biber` is required:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error paper.tex
```

Verify the packaged publication inputs with:

```text
sha256sum -c SHA256SUMS
```

## Artifact scope

- `SHA256SUMS`: hashes of the publication inputs.
- `artifact/generate_ledger.py`: standalone formula generator.
- `artifact/verify_final_results.py`: independent deterministic checker.
- `artifact/all_nine_ledger.json`: generated direct/adaptive ledger.
- `artifact/f19_multiplier_netlist.json`: 150-gate scalar multiplier.
- `artifact/field_tower_circuits.json`: extension-field circuit record.

All numerical security comparisons are unit-leading AXN
pointwise-multiplication-schedule ledgers conditional on the exact public
preflights stated in the paper. They exclude the route-dependent multiplier
`kappa_hom`.
