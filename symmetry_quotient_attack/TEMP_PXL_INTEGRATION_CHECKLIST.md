# Temporary PXL integration checklist

This checklist is anchored to the current canonical
`symmetry_quotient_attack/paper.tex`. The line numbers refer to the version
inspected on 12 August 2026. Recheck them immediately before integration.

The proposed prose is in
`TEMP_PXL_REPLACEMENT_DRAFT_DO_NOT_INCLUDE.tex`. Do not include that file as a
whole. Move each approved block into the canonical manuscript at the anchors
below.

## Invariants

- Do not change the statement or label of `thm:symmetry-quotient` near line
  690.
- Do not change the statement or label of `lem:canonical-residual` near line
  762.
- Do not change the statement or label of `prop:exact-affine` near line 793.
- Preserve `thm:collision-identity` and `thm:sq-accepted-root`. Replace their
  prose only with the versions that allow `h <= K` and say full column rank.
- Keep the final verifier check in every recovery theorem.
- For the Level I `ell=2` slice, require
  `dim_A Span_A(pi_v(L)) = h`. Base field injectivity of `pi_v|_L` is not
  enough. Require `pi_v(a) notin Span_A(pi_v(L))` for the affine coset.
- Sample a uniform affine coset first and discard a nontransverse coset. Do not
  sample from a conditioned coset distribution and then apply the uniform
  coset root theorem. For `h = v - 2`, subtract the exact failure probability
  `|A|^-2 = 19^-4` from the root lower bound.
- For the experiment certified by the current root theorem, sample a fresh
  retained target and affine coset in every trial. Put target filtering and
  every auxiliary cost inside the retry numerator. Do not claim that one fixed
  target can be reused across all cosets.
- Define `p_dir = max(p_reg - 19^-4, 0)` as a certified lower bound from a
  union bound. Do not call it the exact usable trial probability.
- Do not turn either `H_route` or `H_PXL` into a theorem. `H_route` is an
  instancewise premise about the actual restricted system. It does not assert
  that this system has the uniform MQ distribution.
- Keep the published PXL assumptions separate from the manuscript's added
  solver completeness premise. Charge root extraction and every cost not in
  the five published PXL terms through symbolic `C_aux`.
- Keep five logical categories separate: exact algebra, the ideal transcript
  spectrum statement, ordinary MQ heuristics, the instancewise route premise,
  and unit-leading asymptotic estimates.
- Do not call a conditional estimate an upper bound or a complete ledger.
- Keep every `ell=4` numerical row out of the draft until the structured family
  proof and the constrained ledger are both complete.

## Front matter

- Lines 69 to 97, abstract. Replace the homotopy exponents and the unmeasured
  `kappa_hom` sentence. State that the quotient and coset calculation are
  exact. Identify the ideal transcript spectrum statement, ordinary MQ
  assumptions, and unit-leading PXL estimate separately. Use only approved
  `ell=2` numerical rows.
- Lines 129 to 200, introduction and results. Replace the claim that every
  route ends in a square system. The exact reduction permits any affine slice.
  Present the Level I `48` equation and `46` variable route as the conditional
  concrete estimate.
- Lines 209 to 228, logical status table. Replace the homotopy row with separate
  rows for affine coefficient support, the ideal transcript spectrum condition,
  ordinary MQ `H_PXL`, instancewise `H_route`, the pending structured `ell=4`
  premise, the unit-leading estimate, the dense storage proxy, and the final
  exact check.
- Lines 230 to 283, attack overview. Replace the homotopy narrative with the
  target filter, transverse affine slice, full residual PXL solve, and final
  verifier check.

## Background and exact reduction

- Lines 464 to 591, related work. Add Furue and Kudo's PXL paper and state its
  affine semi-regularity, projected semi-regularity, and Macaulay rank
  assumptions. Position the present work as applying that standard model to an
  exact quotient system.
- Lines 593 to 899, exact quotient. Preserve the quotient theorem and all public
  rank checks. Generalize `cor:exact-square-restriction` to an `h` dimensional
  affine restriction with `h <= K`, or add a new corollary directly after it.
  Do not claim root existence in this section.
- Lines 901 to 1186, current preparation section. Replace the square-only
  preparation with `prop:pxl-level1-transverse` and
  `lem:pxl-uniform-coset`. Keep exact target rejection sampling and the fixed
  key rank checks.
- Port the proved `ell=2` affine coefficient support proposition with its
  complete proof. The temporary fragment now includes the incidence
  calculation, the alternating lift and its one dimensional kernel, and the
  descent to the base field with the official coefficient constraints.
- The `ell=4` proof is separate. Port the fixed-map channel surjectivity theorem
  only with its qualification on descent. Surjectivity descends only when the
  spaces, maps, and coefficient conditions arise by scalar extension from the
  base field.
- State what affine support proves. It proves surjectivity of the coefficient
  map. It does not prove uniform coefficients after byte reduction,
  semiregularity, or a Macaulay rank prediction.

## Root theorem and main result

- Lines 1189 to 1417, root analysis. Keep `thm:collision-identity` and
  `thm:sq-accepted-root`, but state them for `h <= K`. Replace nonsingular or
  full Jacobian language with full column rank when `h < K`. Insert
  `eq:pxl-alpha64`, `eq:pxl-eta-level1`, `eq:pxl-preg-level1`, and the
  certified transversality lower bound `eq:pxl-pdir-level1` after the general
  theorem.
- Lines 1419 to 1495, main recovery theorem. Replace the all-nine homotopy
  theorem with the conditional Level I `ell=2` PXL model. Keep a text
  placeholder for `ell=4` until the structured family proof and the
  constrained ledger are complete. Arithmetic alone does not
  establish the `ell=4` coefficient model.
- Lines 1497 to 1546, evidence table. Add rows for an exact slice checker, the
  fixed key spectrum certificate, the symbolic PXL ledger, and any small-size
  rank experiments. State that no production solve is included.

## Design section and conclusion

- Lines 1550 to 1753, design implications. Keep the exact output floors. Remove
  any conclusion that depends on homotopy. Add only consequences that follow
  from the exact quotient dimensions.
- Lines 1755 to 1776, conclusion. State the exact structural weakness first.
  State the Level I PXL result as a unit-leading conditional estimate. Do not
  present the exponent as a proved cost bound, a complete operation count, or
  a wall clock claim.

## Appendices

- Lines 1790 to 2608, homotopy and recovery appendices. Remove the homotopy
  adaptation and its cost formulas. Insert the named definitions
  `def:pxl-regularity` and `def:pxl-route`, the Hilbert series size `eq:pxl-A`,
  all five cost terms, the certified transverse retry lower bound, and the
  augmented auxiliary-cost expression.
- Lines 2738 to 2972, Just Guess appendix. Remove it from the main claim. If it
  is retained for comparison, mark it as a secondary sensitivity analysis and
  do not use it to support the headline Level I result.
- Lines 2974 to 3183, streamed XL appendix. Replace it with the complete PXL
  derivation. Charge `C1`, `C2`, `C3`, fixing, and second linearization.
  Explain why `19^k` appears in the last two terms but not in preprocessing.
- Lines 3185 to 3245, fixed key certificate. Keep this section. Update its
  cross references so it certifies the spectrum inequality used by
  `p_reg`.
- Lines 3246 to 3300, numerical appendix. Replace the homotopy rows with the
  independently checked PXL parameters. Record `h`, `K`, `k`, `D`, `A`,
  `omega`, each of the five PXL terms, `p_reg`, `p_dir`, field operation
  cost, gate exponent, target generation cost inside the retry numerator, the
  uncharged auxiliary costs, and final headroom. Optimize and report `k`, `D`,
  and `A` separately for each value of `omega`.
- Label `5 B^2 binom(k + D, D)` as a unit-leading dense proxy for symbolic
  preprocessing storage. State that it omits representation overhead,
  temporary arrays, and storage traffic. Do not label it as total dense state,
  total memory, or peak memory.

## Bibliography and artifact work after approval

- Do not edit `references.bib` during drafting. At integration time, add the
  final bibliographic entry for Furue and Kudo's PXL paper and use one stable
  citation key throughout.
- Add a checker that recomputes the coefficient series, degree prediction,
  matrix size, all five PXL terms, the transversality subtraction, and the root
  retry denominator.
- Add an independent checker that evaluates the same formulas without importing
  the main ledger module.
- Refresh `SHA256SUMS` only after the manuscript and artifact contents are
  final.
- Run the full artifact verification, LaTeX build, reference check, and a scan
  for stale `kappa_hom`, Just Guess headline, and homotopy headline language.

## Numerical release gate

Do not release a numerical row until two independent computations agree on:

- the chosen `k` and `D`;
- the coefficient `A`;
- `C1`, `C2`, `C3`, `C_fix`, and `C_lin`;
- `p_reg` from the accepted root theorem;
- `p_dir = max(p_reg - 19^-4, 0)` as a certified lower bound for the Level I
  `ell=2` route;
- the Boolean gate exponent after the factor `150`;
- the additive target generation cost;
- the fact that setup, extraction, control, storage traffic, and final
  verification remain outside the unit-leading estimate, unless a separate
  bound is supplied;
- the separate minimizing parameters and symbolic storage value for each
  value of `omega`;
- the final Level I headroom.
