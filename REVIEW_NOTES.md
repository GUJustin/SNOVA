# Editorial and correctness review

## Scope

I reviewed both manuscripts line by line at the LaTeX-source level, followed the dependency chain of every theorem and proof, checked notation and hypotheses against later uses, recomputed the numerical claims that drive the concrete attack, audited citations and cross-references, rebuilt both papers from clean states, and visually inspected every rendered page.

This is a very strong manual and computational audit, not a machine-checked formal proof.

## Main corrections

### `alignment_trap`

- Removed a duplicated definition of the centralizer and normalizer.
- Made the simple-module hypothesis used by the Jacobson density theorem explicit.
- Expanded the normalizer-uniqueness proof to justify its strict inequality.
- Stated explicitly the transpose compatibility of the block embedding used in the Hermitian-factorization proof.
- Corrected `an w-dimensional` to `a w-dimensional`.
- Calibrated the application and conclusion so that the paper points to explicit forgery systems in the companion manuscript without overstating what the structural theorem alone proves.
- Improved paragraphing, line wrapping, and several transitions.

### `symmetry_quotient_attack`

- Corrected the reproduction command from the nonexistent `artifact/verify_all.py` to `artifact/verify_final_results.py`.
- Added `make verify` and `make regenerate` targets and documented them.
- Added an independent exact check of the official target-rejection sampling probability across all nine parameter rows. The minimum is `0.1524622985...`, so the stated lower bound `0.152462` is valid.
- Clarified the finite-field homotopy theorem: the separating form is sampled uniformly from the full coefficient space, each forbidden annihilator is a proper linear subspace, and its probability is at most `1/|k|`. Added the general success probability `1 - B^2/|k|`.
- Made the one-sided nature of the homotopy output explicit: substitution and a Jacobian-rank test remove spurious and singular outputs.
- Corrected the structural-preflight statement from a misleading simultaneous claim across six parameter shapes to a per-shape joint-event claim.
- Made the rank hypothesis in the zero-offset target split an explicit public preflight rather than an unconditional assertion.
- Relabeled the `l = 4` results table: the former `direct fallback` column actually reported the adaptive mixture that uses the fallback, not the fallback's own exponent.
- Distinguished unit-leading arithmetic ledgers and nominal gaps from complete NIST gate-count implementations. The excluded multiplier `kappa_hom` is now highlighted, especially for the narrow `0.30733`-bit Level-I nominal gap.
- Clarified what is and is not charged by the extension-field multiplication ledger: pointwise tower multiplications are charged; transforms, evaluation/interpolation overhead, inversion, factorization, control, and memory remain in `kappa_hom`.
- Tightened the three-plus-one theorem proof around consistent rank-deficient systems and incomplete trials.
- Moved the bibliography to the conventional end position after the appendices.
- Fixed grammar, punctuation, duplicated wording, and several dense or ambiguous transitions.

### Bibliography

- Corrected Eric Schost's first-name accent.
- Corrected the Safey El Din-Schost article DOI from `10.1016/j.jsc.2017.07.001` to `10.1016/j.jsc.2017.08.001` in both bibliographies.

## Mathematical and numerical checks

I found no remaining fatal contradiction in either manuscript's stated theorem chain.

For the concrete paper, the deterministic verifier passes all of the following after regeneration:

- the 150-gate `F_19` multiplier on all 361 canonical input pairs;
- the `F_(19^2)` and `F_(19^4)` tower identities and circuit recurrences;
- the exact official-target rejection-sampling acceptance;
- accepted-root moment bounds;
- conservative `B^2` separating-form extension optima;
- complete-square dimension inequalities for all nine parameter rows;
- projective `l = 4` Jacobian-preflight bounds;
- adaptive all-nine ledgers, the optional 55-to-16 orbit reduction, and repair floors.

I independently recomputed the key table values and repair percentages from the formulas in the paper. They agree with the artifact and the manuscript.

## Build and presentation checks

- Both papers build from a clean state with `latexmk`.
- Final logs contain no LaTeX warnings, unresolved references, unresolved citations, overfull boxes, or underfull boxes.
- Labels, references, citations, and environments were checked for missing or duplicated entries.
- All 13 pages of the structural paper and all 40 pages of the concrete paper were rendered and visually inspected.
- No clipping, overlap, malformed glyphs, broken tables, stranded headings, or visibly bad page breaks remain.
- All fonts are embedded, and text extraction contains no Unicode replacement characters.

## Residual caveats that should remain visible to readers

- The review is not a proof-assistant verification.
- The package does not contain a production-size end-to-end forger or a measured symbolic-homotopy implementation.
- The concrete exponents are unit-leading arithmetic ledgers, not complete wall-clock or NIST-gate estimates; `kappa_hom`, memory, transforms, factorization, and control overhead remain explicit exclusions.
- The Level-I direct `l = 4` nominal gap is only about `0.30733` bits before `kappa_hom`, so it should not be presented as an implementation margin.
