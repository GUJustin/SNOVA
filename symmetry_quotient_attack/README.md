# The Symmetry-Quotient Attack on Odd-Characteristic SNOVA

This package contains the complete referee-revised paper and deterministic
reproducibility artifact.  It contains exactly one LaTeX manuscript source,
`paper.tex`.

## Build

Requirements: Python 3.10 or newer, `latexmk`, Biber, and a TeX distribution
including `biblatex`, `cleveref`, `mathtools`, and `tabularx`.

```sh
make
```

## Verify

```sh
make verify
```

The verification target:

- independently recomputes all nine corrected direct and adaptive homotopy
  ledgers, conditional on the exact public structural preflights;
- checks the exact official target-rejection sampling probability;
- exhaustively checks the 150-gate F19 multiplier on all 361 canonical pairs;
- checks the F19^2 and F19^4 tower identities and circuit recurrences;
- checks dimensions, spectrum tails, separator choices, adaptive mixtures,
  and repair floors;
- exhaustively checks the F19^2 quadratic-root routine; and
- recomputes the expected-tree and capped-tree Level-I Just Guess ledgers and
  their cross-channel second-moment bounds.

The checker does not prove the affine reduction, the random-XOF idealization,
the cited symbolic-homotopy theorem, the correspondence with the pinned
SNOVA verifier, the fixed-layout coefficient-rank hypotheses, or an upper
bound on `kappa_hom` or `kappa_JG`.

## Regenerate generated ledgers

```sh
make regenerate
make verify
```

`artifact/generate_ledger.py` is a self-contained generator for the corrected
all-nine ledger.  It assigns no random probability to SNOVA's fixed outer
maps: the l=4 figures are conditional on the exact public structural
preflight stated in the paper.
## Editor-correction port

- `EDITOR_PORT_REPORT.md` gives the semantic diff and explains every
  mathematical correction ported from the editor's manuscript.
- `EDITOR_PORT.patch` is the unified source diff from the preceding
  referee-revision `paper.tex` to this corrected version.

