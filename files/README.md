# SNOVA cryptanalysis experiment tooling

Companion utilities for the symmetry-quotient forgery paper. These support the
two reviewer concerns that need computation rather than prose.

## What runs where

| Tool | Concern | Runs in sandbox? | Needs from you |
|------|---------|------------------|----------------|
| `rank_cert.py` | #3 generic-rank certificates for III/V/v2.4 | YES | the reduced matrices (E_R-bar, W_R L_R) as .npy/.csv |
| `semireg_dreg.py` | #1 predicted operating degree (determined systems) | YES | just (m,n) |
| `solving_degree.magma` | #1 OBSERVED solving degree | NO (needs Magma) | residual system in `system.mag` |
| `run_msolve.sh` | #1 OBSERVED solving degree | NO (needs msolve) | residual system in `system.ms` |

## The gating experiment (#1)

Appendix E already shows a (3,3,1) lifted family hitting F4 degree 6 where the
formal criterion predicted ~3-4. To know whether the ONE-BASE headline systems
behave that way, run `solving_degree.magma` (or `run_msolve.sh`) on the actual
50-in-102 (Level I) and 70-in-146 (Level III) residual systems your reduction
emits, and compare the observed max step degree to `semireg_dreg.py`'s
prediction for the corresponding determined subsystem. If observed > predicted,
the true cost is above the table and the margins soften; if it matches, the
"below target" claim firms up from "estimated" toward "validated-degree".

## Honesty scope

`semireg_dreg.py` is the standard semi-regular predictor for DETERMINED systems.
It does NOT run the Hashimoto (a,k) optimization and does NOT reproduce the
paper's gate costs — that is your reproducibility package. Per-set operating
degrees for the underdetermined headline systems require that pipeline first.
