# Nine-shape quotient-rank certificate (Task 2)

**Claim under test (paper §Implementation and Algebraic Validation, "All-parameter
quotient-rank certificates"):** for the simple public relation `rho = e_0`, the
symmetry-reduced feature matrix `E_R` attains full column rank `K = m1*C(l+1,2)`
for all nine `q=19` Version 2.3 shapes.

**Verdict: CONFIRMED. All nine ranks equal K exactly. No mismatch.**

## Method

- Constants: the fixed public `A,B,q1,q2,S` reconstructed deterministically by
  `reconstruct_abq` (`shake_256(b'SNOVA_ABQ')` + the documented `improve`/coefficient
  repair). These are the official Version 2.3 public constants; this certificate is
  KAT-independent because `E_R` depends only on `(A,B,q1,q2,S)` and the relation
  `rho`, not on the per-key `P22` block (which needs the KAT).
- Relation: `rho = e_0` (`rho[0]=1`, else 0), matching the paper's "simple public
  relation".
- Rank: exact over F_19 (no floating point). Computed **twice, independently**:
  1. `repro/symmetry_attack_validation.py --rank-only` (its own RREF-mod-19), and
  2. the designated tool `rank_cert.py` run on the dumped `E_R` matrices (Fermat-inverse
     RREF mod 19).
  Both agree on every shape.

## Certificate

| Level | (v,o,l,r)      | matrix size | rank | K = m1·C(l+1,2) | match |
|-------|----------------|-------------|------|-----------------|-------|
| I     | (28, 5, 4, 4)  | 80 × 50     | 50   | 50              | OK    |
| I     | (48,16, 2, 2)  | 64 × 48     | 48   | 48              | OK    |
| I     | (28, 4, 4, 5)  | 80 × 50     | 50   | 50              | OK    |
| III   | (40, 7, 4, 4)  | 112 × 70    | 70   | 70              | OK    |
| III   | (72,24, 2, 2)  | 96 × 72     | 72   | 72              | OK    |
| III   | (38, 5, 4, 5)  | 100 × 70    | 70   | 70              | OK    |
| V     | (50, 9, 4, 4)  | 144 × 90    | 90   | 90              | OK    |
| V     | (96,32, 2, 2)  | 128 × 96    | 96   | 96              | OK    |
| V     | (52, 6, 4, 6)  | 144 × 90    | 90   | 90              | OK    |

These (size, rank) pairs reproduce the table in the paper's validation section
exactly.

## Reproduce

```bash
cd SNOVA/files/SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19
../.venv/bin/python repro/symmetry_attack_validation.py --rank-only \
    --out results/validation_rank_only.json
# independent cross-check with the designated tool:
../.venv/bin/python - <<'PY'
import sys, numpy as np; sys.path.insert(0,'repro')
import symmetry_attack_validation as V
for p in V.PARAMS:
    abq=V.reconstruct_abq(p); rho=np.zeros(p.r,dtype=np.int64); rho[0]=1
    E,_=V.build_E(p,abq,rho); np.save(f"/tmp/{p.name}.npy", np.asarray(E,dtype=np.int64))
PY
for f in /tmp/*-l*.npy; do ../.venv/bin/python ../rank_cert.py "$f"; done
```

*Environment note:* the system `python3` has a broken numpy (x86_64 build on arm64);
use the arm64 venv at `SNOVA/files/.venv`.
