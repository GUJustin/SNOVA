# SNOVA Level-I rank-24 / KAT-census evidence

This supplement contains the exact finite-field data used by the strengthened
Level-I Just-Guess section.

* `rank24/rank24_certificate_summary.json`: statement and per-chart totals.
* `rank24/full24_prep_e*.npz`: chart kernels and transfer blocks over
  A = F_19[u]/(u^2+1), encoded as `a+19*b`.
* `rank24/verify_all_rank24.py`: reruns all 120 single-block and 560 pair-block
  rank checks.  Every rank decision is exact over F_361.  The dense Schur
  multiplier uses the exact integer embedding `a+ib` only for accumulation;
  all integer sums are far below 2^53 and are rounded then reduced mod 19.
* `census/m15_kat100_census_summary.json`: corrected 100/100 deterministic
  (15,4,7) staircase census.  This is finite-corpus evidence, not a key-
  distribution theorem.
* `census/PQCsignKAT_SNOVA_48_16_19_2.rsp`: pinned KAT corpus used by the census.
* `census/census2_*.txt`: per-key run logs.

The rank-24 checker can be run with:

```
python rank24/verify_all_rank24.py
```

The KAT file SHA-256 is
`cb4951074a1a28366617cf3cc2f3d64572dff2daab89045f5bf0f9606f3c5627`.
The count-0 serialized public-key SHA-256 is
`ce99d967cfdc849441981be59fe50cbb92c4faacf2b523e0526f62c470ed9d7a`.


## Final conditional joint-pencil strengthening

The `conditional_joint/` directory supports the strengthened Level-I success
analysis. `verify_conditional_joint_bound.py` exactly enumerates the native
conditional block distributions, verifies

* diagonal atom `49/829`;
* off-diagonal scalar-functional atom `169246/2971565`;
* off-diagonal full two-coordinate atom `38416/11886083 < (49/829)^2`;
* the exact projective union bound for the rank-126 cross-tensor argument,
  whose base-2 logarithm is below `-114.684`; and
* the revised clean/literal success-adjusted Boolean exponents `132.447` and
  `132.340`.

`joint_pencil_basis.npz` reconstructs the actual count-0 joint D/H polar
pencil. `reconnaissance_joint_rank_fast.py` is explicitly reconnaissance only;
its random tests are not used in any theorem.

The `ledger/` directory contains the audited `(15,4,7)` Boolean ledger used by
the final paper.
