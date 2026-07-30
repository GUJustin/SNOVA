# Response to the Mode-A referee report

The revision accepts the report's central criticism: the previous
unit-leading homotopy ledger was not a complete NIST classical-gate count.
The abstract, introduction, main theorem, numerical tables, and conclusion
now use conditional break-even language and display the exact allowable
`kappa_hom` thresholds rather than claiming unconditional all-nine gate
breaks.

## Mathematical and implementation-scope changes

1. The exact symmetry quotient, canonical output split, accepted-root theorem,
   complete-square reductions, and repair floors are retained.
2. The arbitrary reserved-line decomposition is replaced by a coordinate raw
   A-line and a complementary raw-coordinate summand.
3. For the three ell=2 rows, `rank(Q_R)=K` is an explicit fixed-key public
   preflight; no unproved random-key density is charged as a theorem.
4. The unreproducible numerical 16-orbit branch was removed.
5. The unreproducible low-state XL numerical table was removed; only the
   conditional method remains.
6. The square-row rejection statement is now a conservative lower bound, not
   an exact density claim.
7. The package description and checker description now match the distributed
   files and distinguish formula checking from official-verifier testing.
8. Prior-work discussion separates Ran's wedge attack and expands the
   comparison with Jin et al.


## Corrections ported from the independent editor's manuscript

A subsequent line-by-line editor pass identified a serious modeling error in
the earlier structural-preflight ledger.  SNOVA's outer map is fixed, not a
random object to which a full-rank probability may be assigned.  The revised
paper therefore:

- conditions on exact public checks `rank(Q_R)=K` and on fixed-layout
  coefficient-map ranks;
- applies random-XOF probability only to genuinely sampled native public-form
  coefficients;
- uses `(q^(M-K)-1)/(q-1)` projective quotient directions in the affine-rank
  union bound; and
- removes the former `p_pre` normalization from all ell=4 ledgers.

This lowers each conditional ell=4 exponent by 0.085708 bits while making the
claim weaker in the correct logical sense: the ledger is conditional on an
exact public structural preflight, not a random-key theorem for the fixed
outer map.  The all-nine adaptive maxima are unchanged because the ell=2 rows
remain levelwise dominant.

The editor's strengthened finite-cardinality homotopy proof, exact target
sampling calculation, coordinate-aligned freshness condition, corrected
three-plus-one completion proof, table terminology, bibliography placement,
and metadata corrections have also been ported.

## New focused research result

The weakest numerical row, Level-I ell=2, now has a separate Just Guess route.
The published expected field-operation schedule is instantiated with explicit
Boolean circuits for F19^2 arithmetic and quadratic roots.  One trial has
ledger exponent 132.046522.  Conditional on a transcript-checkable candidate
regularity event, a new cross-channel second-moment lemma gives success
probability at least 0.171052 and total exponent

    134.594009 + log2(kappa_JG).

Thus the Level-I break-even condition is

    log2(kappa_JG) < 8.405991.

A separately generated capped schedule explores at most 63 quadratic nodes
and 64 leaves per guess, aborting nonunique residual linear branches.  It has
trial exponent 137.713503 and success-adjusted exponent

    140.260990 + log2(kappa_JG_cap),

with break-even threshold 2.739010 bits.  This is a bounded-work sensitivity,
not a proof of the capped candidate event.

The revision does not rename either regularity multiplier into a theorem.
A production fixed-key transcript remains the cleanest way to measure or
certify it.
