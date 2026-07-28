# Frobenius-orbit-complete homotopy certificate

For h=K-2 there are 55 family-count patterns.  At a descent root, q-Frobenius cyclically rotates the four eigenblocks and carries every nonzero Jacobian minor to a nonzero conjugate minor.  One representative from each C4 orbit therefore suffices.  The exact orbit count is 16.

The finite-cardinality threshold uses arbitrary distinct start nodes and a uniformly random coefficient-vector separator: |k| >= max(max_j sum_i d_ij, 8 B^2, 8h).  All omitted equations and descent/rejection/verifier checks are charged as filters.

| Level | (h,K) | patterns -> orbits | compact AXN | compact margin | robust AXN | robust margin | gain over 55 | peak one-core output |
|:--:|:--:|:--:|--:|--:|--:|--:|--:|--:|
| I | (48,50) | 55 -> 16 | 137.74021 | 5.25979 | 138.32403 | 4.67597 | 1.75492 | 1.890 PiB |
| III | (68,70) | 55 -> 16 | 179.17452 | 27.82548 | 179.75835 | 27.24165 | 1.76217 | 2.076 ZiB |
| V | (88,90) | 55 -> 16 | 220.24857 | 51.75143 | 220.83240 | 51.16760 | 1.76629 | 2416.552 YiB |

## Random-XOF excess-spectrum tails

| Level | t by shape | log2 bad compact | log2 bad robust |
|:--:|:--:|:--:|:--:|
| I | 74/74 | -102.03/-102.03 | -105.35/-105.35 |
| III | 110/114 | -168.01/-184.78 | -171.32/-188.09 |
| V | 138/146 | -200.44/-233.98 | -203.76/-237.30 |

The packed-output values are straightforward dense parametrization ceilings, not peak-memory theorems.  Sequential processing needs only the largest individual parametrization plus solver workspace, but these ceilings are still very large.
