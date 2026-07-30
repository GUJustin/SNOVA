# Structural public-preflight probability certificate

This certificate replaces the earlier single-minor anti-concentration bound by direct kernel-union bounds.  The source, projection, outer-map, and internal-core byte regions are disjoint in the pinned expansion.

For the affine source, the raw `(j,k)=(1,0)` labels quotient by the symmetric-square image to one `wedge^2(A)` per public base form.  Every nonzero dual alternating tuple imposes `4(v-1)` independent affine constraints on the fresh `W x V'` symbols.  Hence

```text
delta_src <= ((19^(6m1)-1)/18) * (14/256)^(4(v-1)).
```

For projection injectivity, the bad `A`-line `W` forces all `10m1` symmetric `W x W` symbols to vanish.  Any other `A`-line in `J=W+O` forces `16m1` independent `W x O` constraints.  Conditional on these source facts, one full-rank outer-map event implies both `rank Q=K` and `rank[Q Lambda]=M`.

| Shape | affine rank proved / needed | log2 delta_src | log2 delta_proj | outer success | complete cross success | penalty bits |
|:--:|:--:|--:|--:|--:|--:|--:|
| I-a | 30 / 30 | -329.54 | -209.63 | 0.942322268 | 0.942322268 | 0.0857 |
| I-b | 30 / 30 | -329.54 | -209.63 | 0.942322268 | 0.942322268 | 0.0857 |
| III-a | 42 / 42 | -479.81 | -293.49 | 0.942322268 | 0.942322268 | 0.0857 |
| III-b | 42 / 30 | -446.27 | -293.49 | 0.942322268 | 0.942322268 | 0.0857 |
| V-a | 54 / 54 | -596.54 | -377.34 | 0.942322268 | 0.942322268 | 0.0857 |
| V-b | 54 / 54 | -630.08 | -377.34 | 0.942322268 | 0.942322268 | 0.0857 |

All displayed probabilities are exact rationals in the JSON.  The decimal cross-preflight values differ only far beyond the displayed digits because the source and projection failures are already negligible; the visible `0.942322...` factor is the conservative biased-product full-rank bound for the outer map.
