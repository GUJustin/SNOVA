# Low-degree sanity check (Concern #3 / paper §Validation, lines 766–784)

**Provenance flag (important).** The paper's "Low-degree sanity check" paragraph
reports Jacobian / strip / polar / Hessian statistics, but **no script in the
shipped reproducibility package computes any of them** (verified by exhaustive
search of `SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19/`: the only regex
hit for jacobian|hessian|polar|strip|Macaulay was an unrelated `.split` in the
KAT parser). Those numbers were therefore **unverified** by shipped code.

This report reconstructs them here from the OFFICIAL key.
Tool: `experiments/task3_lowdegree.py`, raw log `experiments/logs/task3_full.log`,
machine-readable `results/task3_lowdegree.json`. All arithmetic is exact over
F_19. The 50 residual quadratic forms are taken from the same verified builder
used for the solving-degree core: emitted polynomials are checked against the
reference verifier `direct_output` on random points before any statistic is
computed, and each of the 20 additional keys re-runs that same verifier net.

## Results — official Level-I 50-in-102 residual (from the KAT key)

| Quantity | Paper | This run | Match |
|----------|-------|----------|-------|
| Planted-root Jacobian rank | 50 | **50** | yes |
| Random coord 50-strips invertible / 10,000 | 9,455 | **9,466** | yes (both ≈ rmt) |
| Random-matrix probability ∏_{k=1}^{50}(1−19^{−k}) | 0.9445987429… | **0.9445987429** | exact |
| Polar-matrix output rank (102 coord + root + 10,000 random dirs) | 50 (all) | **50** (all 10,103) | yes |
| Hessian min rank over 22,100 wt≤2 + 20,000 random combos | 100 | **100** (42,100 tested) | yes |

The empirical strip fraction 0.9466 differs from the paper's 0.9455 only by
Monte-Carlo sampling (different RNG seed); both are consistent estimates of the
theoretical invertibility probability 0.9445987429 (σ for 10,000 draws ≈ 0.0023).

## Results — 20 independently generated official-shaped keys

Each key built from a fresh 48-byte seed through the identical key-gen + residual
path (`build_from_seed`), verified against the reference verifier.

- **Every** planted root had Jacobian rank **50** (20/20).
- Coordinate-strip full-rank counts per key ranged 9,412–9,491; aggregate
  **188,944 / 200,000 = 0.94472**, sitting on the theoretical 0.9445987429.
- The paper's single figure "9,418 / 10,000" falls inside the observed per-key
  range and is ~1.2σ below the mean of 9,447 — statistically consistent. (The
  paper's phrasing is ambiguous between a per-key draw and an aggregate; either
  reading is consistent with this run.)

## Not yet reproduced

- **Degree-one Macaulay full column rank (5,150 columns), "no unexpected
  degree-one syzygy."** Not recomputed by an explicit Macaulay-matrix rank here.
  Partial corroboration: the msolve F4 run on the official core shows **0 zero
  reductions at degree 3** (the degree-one multiplication level), consistent with
  the absence of unexpected degree-one syzygies. A direct 187k×5,150 sparse rank
  over F_19 on the 50-in-102 residual is the clean check and is left as a
  follow-up.

## Verdict

The paper's low-degree sanity-check numbers were **not backed by shipped code**
but are **reproducible and correct**: reconstructed independently from the
official key (and 20 further keys), every named statistic matches within
Monte-Carlo tolerance, and the exact quantities (Jacobian rank 50, polar rank 50,
Hessian min rank 100, rmt probability) match exactly. The paragraph's claims
stand; only the reproducibility package needs `task3_lowdegree.py` added so the
numbers are no longer unverified.
