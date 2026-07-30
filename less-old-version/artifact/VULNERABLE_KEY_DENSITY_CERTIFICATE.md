# Random-XOF vulnerable-key density certificate

This certificate combines three exact ingredients for each official `ell=4` shape:

1. the disjoint-coordinate public preflight bound (full quotient/source rank, outer-map ranks, and projection injectivity);
2. the blockwise anti-concentration bound for a nonzero core-Jacobian witness polynomial; and
3. the robust seeded-spectrum tail `Pr[bar_eta_L > 1/2]`.

The cross-preflight coordinates are disjoint from the internal `V' x V'` coordinates, so those probabilities multiply.  The Jacobian and spectrum events share the internal coordinates, so the exact certified lower bound is

```text
p_vuln >= p_cross * (p_Jac - eps_spectrum).
```

All probability arithmetic in the JSON is exact rational arithmetic.  Logarithms below are presentation values.

| Shape | p_cross | p_Jac | eps_spectrum | p_vuln lower | penalty (bits) | per-key robust | normalized robust | margin | normalized + 2^8 | margin |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| I-a | 0.942322268 | 0.136912282 | 3.262e-13 | 0.129015492 | 2.954 | 134.077 | 137.031 | 5.969 | 145.031 | -2.031 |
| I-b | 0.942322268 | 0.136912282 | 3.262e-13 | 0.129015492 | 2.954 | 134.077 | 137.031 | 5.969 | 145.031 | -2.031 |
| III-a | 0.942322268 | 0.123579079 | 1.131e-33 | 0.116451318 | 3.102 | 196.246 | 199.348 | 7.652 | 207.348 | -0.348 |
| III-b | 0.942322268 | 0.123579079 | 1.131e-33 | 0.116451318 | 3.102 | 196.246 | 199.348 | 7.652 | 207.348 | -0.348 |
| V-a | 0.942322268 | 0.212366176 | 5.110e-49 | 0.200117376 | 2.321 | 247.432 | 249.753 | 22.247 | 257.753 | 14.247 |
| V-b | 0.942322268 | 0.212366176 | 5.110e-49 | 0.200117376 | 2.321 | 247.432 | 249.753 | 22.247 | 257.753 | 14.247 |

## Interpretation

The per-key column is the attack ledger conditional on the public certificate.  The normalized column additionally pays the inverse certified key density, equivalently measuring work per successful forgery over a random generated public key in the idealized distribution.  It is not a claim that an adversary can force a target installation to regenerate its key.

The `+ 2^8` column exposes a factor-256 implementation sensitivity on the leading homotopy term.  It is a sensitivity point, not a theorem that all hidden constants are at most 256.
