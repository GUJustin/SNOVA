# Correspondence tests

These tests are deliberately separate from the formula ledgers.

- `official_key_reduction_harness.py --check` parses one pinned official
  `(28,5,4,4)` KAT from public fields, verifies the supplied signature in a
  Python Version 2.3 transcription, checks the common-column decomposition,
  and confirms quotient/affine ranks 50/30. It leaves the 50-in-102 residual
  system unsolved and emits no forgery.
- `reduced_parameter_end_to_end_forgery.py --check` uses the deliberately
  unofficial shape `(2,1,2,2)`. It interpolates the complete zero-offset
  restriction `X_i=[u_i|0]`, checks equality of the interpolated and explicit
  rank-three quotient images, fixes each SHAKE target before the public slice
  search, and reconstructs a non-planted serialized signature. `--check`
  recomputes the main trial and eight fresh-key regressions before comparing the
  deterministic JSON.
- `verify_reduced_parameter_stress.py` reruns 24 further fresh-key forgeries,
  reconstructs their targets through the literal evaluator, checks changed-key
  rejection, obtains rank three on 100 further reduced keys, and round-trips
  base-19 vectors of lengths 1 through 100.

The reduced forger and literal verifier path share the KAT-anchored Python
scheme helpers, so these are composition regressions, not independent verifier
correspondence or production-size cost evidence. The original supplied
standalone demo required two convention corrections before its attack logic
could be integrated; `REDUCED_ATTACK_REVIEW.md` records them.

The two transcript scripts accept `--write` to regenerate their committed JSON.
With `--check` they do not modify files. The shared `snova_v23_reference.py` is a
local transcription pinned to the commit named in `../official/COMMIT`; it is
not the unmodified official C verifier.
