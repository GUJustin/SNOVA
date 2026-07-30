# Final adversarial correctness audit

**Date:** 30 July 2026

## Bottom line

The editor's major corrections were mathematically sound and improved the
paper. They did not invalidate the symmetry quotient, exact output split,
accepted-root moment theorem, complete-square reductions, or conditional
homotopy ledgers.

The audit did find two defects in the subsequently combined manuscript:

1. an editor correction to the zero-offset dimension theorem was accidentally
   omitted; and
2. the newer Just Guess section used a false channel-independence assertion.

Both are corrected in this release. The first strengthens the repair barrier;
the second weakens the Just Guess numerical sensitivity by about `0.831316`
bits while making the proof valid under the native byte-biased XOF model.

## Re-derived editor changes

### Fixed outer map

The outer map, relation, offsets, and coordinate layout are deterministic.
The paper now conditions on exact coefficient-map and realized-rank
preflights. Random-XOF probabilities apply only to sampled native public-form
coefficients. No structural success probability is assigned to the fixed
outer map or divided out of the `ell=4` ledger.

### Structural union bounds

The affine-rank failure event is indexed by projective nonzero elements of the
quotient dual, giving `(q^(M-K)-1)/(q-1)` directions. The intersection event
is indexed by `A`-lines in `W+O_pub`. Maximum-atom pivot bounds and union
bounds require no independence between the two bad events.

### Finite-cardinality homotopy

The start system uses enough distinct field nodes for the Vandermonde
conditions. A uniformly sampled separating-form coefficient vector lies in
the annihilator of each fixed nonzero forbidden vector with probability at
most `1/|E|`; at most `B^2` vectors give success at least `1-B^2/|E|`.
Complete substitution and a Jacobian test make failure one-sided. The theorem
continues to expose all soft-O constants through `kappa_hom`.

### Exact targets and completion

The official target sampler is made exactly uniform by chunkwise rejection.
The fourth-block completion theorem now rejects inconsistent systems and
marks consistent rank-deficient systems incomplete unless separately priced.

## Newly corrected issues

### Full zero-offset domain

Zero offsets require no reserved line. The common column ranges over the full
`A^v` public-vinegar space, of base-field dimension `dv`. The necessary
condition for eliminating the chosen-message square by dimension is therefore
`rho_Q>dv`, not `rho_Q>d(v-1)`. The exact nine-row table is checked by an
independent script and ranges from `101.5625%` to `144.00%` more outputs.

### Conditional Frobenius-channel distribution

The native coefficients are independent before, but not after, the channel
transform. The audit derived the official transform explicitly and performed
weighted exhaustive enumeration over all `19^3` diagonal and `19^4`
off-diagonal native blocks. The common conditional atom bound is `49/829`.
The Just Guess second moment and every associated exponent are regenerated
from this value.

## Computational checks performed

- independent recomputation of every all-nine direct and adaptive ledger row;
- exact repair-floor tables;
- exhaustive `F_19` multiplier verification on all 361 inputs;
- exact `F_(19^2)` and `F_(19^4)` tower identities and circuit recurrences;
- Tonelli-Shanks verification on all 361 discriminants;
- exact conditional-channel atom enumeration;
- exact target-sampling acceptance;
- spectrum and structural coefficient union bounds;
- separator extension-degree optimization;
- regenerated Just Guess expected and capped ledgers;
- clean LaTeX build, log inspection, PDF rendering, and visual page review;
- clean checksum manifest and ZIP integrity test.

## Claims that remain conditional

The exact algebraic reduction is conditional on public preflights whose
realized official-key certificates are not distributed. The random-XOF
spectrum theorem is a model theorem. Homotopy costs retain `kappa_hom`.
Just Guess retains `kappa_JG` or `kappa_JG_cap`. The release does not claim a
production end-to-end forgery or an unconditional complete-gate break of all
nine rows.
