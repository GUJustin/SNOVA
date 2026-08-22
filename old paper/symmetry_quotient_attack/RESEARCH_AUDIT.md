# Final focused research and correctness audit

## Objective

Stress-test the editor's corrections and the later Level-I Just Guess route,
then retain only claims supported by exact algebra, a stated random model, or
an explicit sensitivity multiplier.

## Editor corrections

The following high-impact editor changes were re-derived and found correct:

- the outer map is deterministic and receives no random rank probability;
- the affine-rank union bound is over projective quotient-dual directions;
- the reserved line must be coordinate-aligned under biased byte sampling;
- the finite-cardinality homotopy separator succeeds with probability at least
  `1-B^2/|E|` when sampled from the full coefficient space;
- the target-rejection sampler is exact, with minimum acceptance
  `0.1524622985242944` over the official output lengths;
- rank-deficient fourth-block completion must be declared incomplete unless a
  separately priced enumeration procedure is used.

The audit found one missed editor port: the zero-offset repair theorem must
use the full `A^v` public-vinegar space. The corrected output-increase range is
`101.5625%--144.00%`.

## Just Guess correction

The earlier Just Guess draft treated the diagonal and cross Frobenius
channels as independent after conditioning on the diagonal solver transcript.
That is false for the native byte-biased coefficient distribution.

For `S=[[1,2],[2,15]]` over `F_19[u]/(u^2+1)`, eigenvectors
`(1,13+u)` and `(1,13-u)` give the exact block transforms printed in the
paper. Exhaustive weighted enumeration gives:

    diagonal max conditional atom = 49/829
    off-diagonal max conditional functional atom = 169246/2971565

Conditioning on all diagonal coefficients remains blockwise because the
transform is local to each independent native matrix block. For distinct
sign-canonical diagonal roots, one nonzero Hermitian functional per base form
then yields pair-collision probability at most `(49/829)^16`.

At `c=1/4`, the corrected second moment is

    Pr[cross match] >= c / (1 + c*(19*49/829)^16)
                    = 0.09613444074026915...

This adds `3.378802812` bits. The resulting sensitivities are:

    expected-operation route:
      135.425324835 + log2(kappa_JG)
      break-even log2(kappa_JG) < 7.574675165

    capped-tree route:
      141.092306187 + log2(kappa_JG_cap)
      break-even log2(kappa_JG_cap) < 1.907693813

The per-`z` `2^40` filtering amount is a declared reserve, not a proved bound
for every candidate population. Any excess belongs to `kappa_JG`.

## Remaining limitations

The package still does not contain:

1. official-v2.3 public keys and end-to-end equality tests against the pinned
   verifier;
2. realized structural-preflight certificates for every official fixed key;
3. a proof or production measurement of the Just Guess candidate-abundance
   multiplier;
4. a complete end-to-end Boolean implementation of symbolic homotopy.

These are stated limitations. They are not renamed as theorems or hidden in
unqualified "breaks all nine" language.
