# Napkin Runbook

## Curation Rules

- Read only at session start and write only at session end.
- Keep recurring, high-value notes only.
- Re-prioritize on every write and keep at most 10 items per category.
- Give every item a concrete `Do instead:` action.

## Complexity Claims (Highest Priority)

1. **[2026-08-12] Separate necessary floors from constructive upper bounds**
   Do instead: Rule out a route only with operations every implementation of
   that named kernel must perform. Use a complete constructive count for a
   positive security claim. A failing conservative upper envelope is not a
   lower-bound exclusion.

2. **[2026-08-12] Anchor homotopy refinements to the ledger's exact baseline**
   Do instead: Start from `exact_homotopy_factor`, whose numerator charges
   `(2e-1) * g_A` and whose transforms remain in `kappa_hom`. State explicitly
   which baseline factor each new circuit or operation ratio replaces.

3. **[2026-08-12] Charge the complete PXL model and slice retries**
   Do instead: Use all three terms in PXL equation 4.8, optimize the guessed
   variables and predicted degree, multiply by the accepted root theorem's
   reciprocal probability, and only then apply the field operation gate cost.
   Report the PXL semi regularity assumptions separately from proved
   coefficient support.

## Execution and Release Hygiene

1. **[2026-08-13] Solver-headline revisions require matching verifier updates**
   Do instead: When replacing the homotopy tables with PXL tables, update
   `artifact/verify_ledgers.py` in the same release slice. Until then, expect
   `make verify` to fail because it still searches `paper.tex` for the removed
   homotopy rows.

2. **[2026-08-12] Full artifact runs create untracked Python 3.14 bytecode**
   Do instead: After `make regenerate` or `make verify`, remove only the
   untracked `*.cpython-314.pyc` files. Preserve the tracked Python 3.13
   bytecode already present in the repository.

3. **[2026-08-12] Verify the root checksum manifest after artifact changes**
   Do instead: Run `sha256sum -c SHA256SUMS` from
   `symmetry_quotient_attack/` and refresh every stale distributed-file entry,
   including existing entries affected by earlier changes.
