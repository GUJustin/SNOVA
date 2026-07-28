# Direct core-Jacobian probability certificate

For the recommended original-family cores, a nonzero left-kernel combination supported on one public base form forces at least `s-1` independent A-valued equations on disjoint off-diagonal coordinate-pair blocks.  Each A-valued equation has F_19-rank four.  With `rho=14/256`,

```text
delta_core <= product_i (1 + (|A|^t_i-1) rho^(4(s-1))) - 1.
```

Here `t_i` is the exact number of selected core rows drawn from public form `i`.  The row-template map is injective: loop multipliers appear in distinct diagonal blocks, adjacent-edge multipliers appear in an exposed first/last block, and the opposite-edge multiplier appears in its own H component.  Thus every active form contributes a nonzero coefficient map.

| Level | s | rows by public form | log2 core failure | core success |
|:--:|--:|:--|--:|--:|
| I | 10 | `[6, 6, 6, 6, 6]` | -46.66 | 0.999999999999991 |
| III | 11 | `[5, 5, 5, 5, 5, 4, 4]` | -80.43 | 1.000000000000000 |
| V | 13 | `[3, 3, 3, 3, 3, 3, 3, 3, 2]` | -147.27 | 1.000000000000000 |

Combining this internal-core event with the exact aggregate-spectrum tail by a union bound, and multiplying by the disjoint structural-preflight probability, gives:

| Shape | robust per certified key | vulnerable-key density | normalized robust | margin | normalized robust + 2^8 | margin |
|:--:|--:|--:|--:|--:|--:|--:|
| I-a | 134.077 | 0.942322267954 | 134.163 | 8.837 | 142.163 | 0.837 |
| I-b | 134.077 | 0.942322267954 | 134.163 | 8.837 | 142.163 | 0.837 |
| III-a | 196.246 | 0.942322267955 | 196.332 | 10.668 | 204.332 | 2.668 |
| III-b | 196.246 | 0.942322267955 | 196.332 | 10.668 | 204.332 | 2.668 |
| V-a | 247.432 | 0.942322267955 | 247.517 | 24.483 | 255.517 | 16.483 |
| V-b | 247.432 | 0.942322267955 | 247.517 | 24.483 | 255.517 | 16.483 |

The normalized ledger is a time-success ratio over a randomly generated public key in the pinned random-XOF idealization.  It does not mean that an adversary can resample a fixed target key.  The `+2^8` column is an exposed implementation-sensitivity point, not a proof that the suppressed homotopy factor is at most 256.
