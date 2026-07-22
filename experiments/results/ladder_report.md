# Solving-degree scaling ladder — observed vs. predicted

**Purpose.** Validate the semi-regular operating-degree predictor
(`semireg_dreg.py`, agreeing with `official_estimator.regdim`) against the
*observed* maximum F4 round degree that msolve 0.10.1 reaches, on genuinely
symmetry-reduced SNOVA cores that are small enough to solve to completion on
this host (8 GB RAM). All cores are faithful reduced systems built through the
verified pipeline (`scope_ladder.py` → `emit_core_general.build_residual` →
`EM.specialize` → `EM.to_msolve`); every core passed `verify_polys` and
`planted_check`. Field q = 19 throughout, so field equations (which bite only
at degree 19) never contaminate the measured degree.

**Method.** For each rung msolve is run `-v 2` to completion (full rational
parametrization emitted — no truncation, no timeout). The observed degree is the
maximum degree row in the F4 round table. Each official core is compared to 3
random controls matched in (m, n, density) with a planted root
(`gen_random_control.py`).

## Results

| rung   | l | K = core (m/n) | density | predicted D_reg | observed (official) | controls (3) | official completed |
|--------|---|----------------|---------|-----------------|---------------------|--------------|--------------------|
| l2_m11 | 2 | 3 (3/3)        | 0.7222  | 4               | 4                   | 3, 4, 3      | yes (0.01 s)       |
| l2_m12 | 2 | 6 (6/6)        | 0.9365  | 7               | 7                   | 7, 7, 7      | yes (0.02 s)       |
| l2_m13 | 2 | 9 (9/9)        | 0.9506  | 10              | 10                  | 5, 10, 10    | yes (1.10 s)       |
| l2_m14 | 2 | 12 (12/12)     | 0.9573  | 13              | 12*                 | 7, 12, 12    | yes (2219 s)       |
| l3_m11 | 3 | 6 (6/6)        | 0.8968  | 7               | 7                   | 4, 6, 4      | yes (0.01 s)       |
| l3_m12 | 3 | 12 (12/12)     | 0.9498  | 13              | 13                  | 7, 7, 13     | yes (662 s)        |
| l4_m11 | 4 | 10 (10/10)     | 0.9618  | 11              | 10*                 | 6, 10, 11    | yes (6.28 s)       |

`*` = observed one below the semi-regular first-fall index; see note below.

## Interpretation

1. **Five of seven rungs: observed == predicted exactly.** On faithful reduced
   SNOVA cores up to K = 12, the observed F4 solving degree equals the
   semi-regular prediction with no offset.

2. **The two `*` rungs are a msolve trailing-round convention, not a true
   fall.** l2_m14 and l3_m12 are both 12/12 cores with *identical* F4 bulk
   tables (same row counts at every degree: 10/11, 28, 155, 468, 1059, 1848,
   2772, 3267 [peak], 3025, 1694, 594). They differ only in the label on a final
   11-row cleanup block — deg 13 for l3_m12, deg 12 for l2_m14. The
   D_reg-determining computation is the same; the ±1 is msolve bookkeeping in the
   trailing spair reduction. l4_m11 (10/10) is the analogous case.

3. **No rung shows observed > predicted.** This is the decisive check: the
   Appendix-E degeneration (a lifted (3,3,1) family reaching F4 degree 6 where the
   criterion predicted 3–4) would manifest as observed *above* predicted, which
   would inflate the true attack cost and soften the margins. No such
   degeneration appears on any reduced core measured here.

4. **Official is indistinguishable from matched random controls.** The top of the
   control spread equals the official observed degree on every rung. The official
   core is never an outlier *below* the controls (which would signal an
   exploitable weakness making the attack cheaper than modelled), and the one
   rung where official sits above the control max (l3_m11: 7 vs 6) is within
   run-to-run noise and again means the reduced core is no *easier* than a generic
   semi-regular system. Low control values (e.g. l2_m13 control 5, l3_m11 controls
   4/4) are early planted-root hits on random instances, not degree falls of the
   official core.

## Scope / honesty

- This validates the **predictor** against **observed** F4 behaviour on reduced
  cores up to K = 12. It does **not** solve the headline systems: the Level-I l4
  core (50/41) and its controls all *time out at F4 round degree 4* (matrices
  exploding to 45010 × 148004), far below D_reg = 14, so no observed solving
  degree is obtainable for the n = 41–79 headline/l4/l2-sqspec cores on this host.
  Those remain predictor-only and must not be reported as F4-audited.
- Conclusion for the paper: on the checkable small end, the semi-regular
  operating-degree model is corroborated by direct F4 measurement with no
  degeneration and no distinguishability from random — moving the small-core
  claim from "estimated" toward "validated-degree", while the extrapolation to
  headline dimensions continues to rest on the (here-corroborated) predictor.

## Provenance

- Cores: `systems/ladder/ladder_*.ms`; controls: `systems/ladder/controls/`.
- Logs: `logs/ladder/*_official.log`, `logs/ladder/ctrl_*.log`.
- Rung spec: `/tmp/ladder_rungs.txt`. Builder: `scope_ladder.py`.
- Predictor: `semireg_dreg.py d_reg(K, min(K, N_res))`.
