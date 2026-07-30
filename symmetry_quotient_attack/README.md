# The Symmetry-Quotient Attack on Odd-Characteristic SNOVA

This package contains the final correctness-audited paper and deterministic
reproducibility artifact. It contains exactly one LaTeX manuscript source,
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

- independently recomputes all nine direct and adaptive homotopy ledgers,
  conditional on the exact public structural preflights;
- checks exact official target-rejection sampling;
- exhaustively checks the 150-gate F19 multiplier on all 361 canonical pairs;
- checks the F19^2 and F19^4 tower identities and circuit recurrences;
- checks dimensions, spectrum tails, separator choices, adaptive mixtures,
  and both repair-floor tables;
- exhaustively checks the F19^2 quadratic-root routine;
- exhaustively enumerates the official ell=2 conditional Frobenius-channel
  atoms; and
- recomputes the expected-operation and capped-tree Level-I Just Guess
  ledgers from the corrected conditional atom bound.

The checker does not prove the affine reduction, the random-XOF idealization,
the cited symbolic-homotopy theorem, correspondence with the pinned SNOVA
verifier, realized fixed-key coefficient-rank hypotheses, or upper bounds on
`kappa_hom`, `kappa_JG`, or `kappa_JG_cap`.

## Regenerate generated ledgers

```sh
make regenerate
make verify
```

`artifact/generate_ledger.py` is self-contained and assigns no random
probability to SNOVA's fixed outer maps. The ell=4 figures are conditional on
the exact public structural preflight stated in the paper.

`artifact/research/verify_l2_channel_conditional_atoms.py` regenerates the
exact weighted conditional-atom ledger used by the Just Guess second moment.

## Audit documents

- `FINAL_CORRECTNESS_AUDIT.md` records the adversarial re-derivation, the two
  additional defects found, and all residual limitations.
- `EDITOR_PORT_REPORT.md` records which editor corrections were validated and
  ported.
- `REFEREE_RESPONSE.md` explains how the paper addresses the Mode-A report.
- `RESEARCH_AUDIT.md` documents the focused solver research and revised Just
  Guess sensitivity.
