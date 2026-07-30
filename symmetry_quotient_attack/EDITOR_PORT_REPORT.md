# Editor correction diff and port report

## Inputs compared

This port used two independent comparisons:

1. the editor repository's principal review commit
   `40f1a82` (`Thorough review pass SDK-2`) against its parent `f014001`;
2. the editor-corrected concrete manuscript against the later
   referee-revision manuscript containing the new Just Guess section.

The resulting source-level port changes `paper.tex` by 329 insertions and 205
deletions, and `references.bib` by 30 insertions and 19 deletions.  Purely
stylistic rewrites were not copied wholesale where they would overwrite newer
research; the mathematical corrections were ported semantically.

## Critical correction ported

### The fixed outer map cannot be assigned a random rank probability

The earlier referee-revision source folded a factor approximately
`0.942322...` into the structural preflight.  That factor treated SNOVA's
fixed outer map as though it were sampled independently in the random-XOF
experiment.  This was not justified.

The corrected theorem now fixes the deterministic outer map `E`, relation
`R`, complete offset vector `Gamma`, and native layout.  It requires exact
public or symbolic preflights:

- `rank(Q_R)=K`;
- `t_src >= d(v-1)`;
- `t_W >= m_1 binom(d+1,2)`;
- `t_oth >= m_1 d^2`.

Only the genuinely sampled native public-form coefficients receive a
random-XOF probability.  The affine-rank union bound ranges over the
projective quotient dual and therefore has factor

`(q^(M-K)-1)/(q-1)`,

not the former source-feature count.  The six ell=4 ledgers are consequently
conditional on passing the exact public structural preflight.  No
probability is assigned to the fixed outer map.

Removing the invalid normalization lowers each conditional ell=4 exponent by
`0.08570756` bits.  Representative changes are:

| Row | Previous direct | Corrected direct | Previous adaptive | Corrected adaptive |
|---|---:|---:|---:|---:|
| `(28,5,4,4)` | 142.692668 | 142.606960 | 128.331985 | 128.246278 |
| `(40,7,4,4)` | 185.536587 | 185.450879 | 171.702386 | 171.616678 |
| `(50,9,4,4)` | 227.646114 | 227.560406 | 214.644029 | 214.558321 |

The all-nine adaptive maxima remain
`132.188072 / 183.794625 / 233.843987`, because the ell=2 rows dominate.
The direct all-nine maxima become
`142.606960 / 185.450879 / 233.843987`.

## Other mathematical corrections ported

### Finite-cardinality homotopy

The proof now states that the separating form is sampled uniformly from the
full coefficient space.  Each forbidden annihilator is a proper linear
subspace and has probability at most `1/|k|`; a union bound gives success at
least `1-B^2/|k|`.  The start system, its exactly `B` simple roots, and the
one-sided substitution plus Jacobian-rank filtering are explicit.

### Coordinate-aligned freshness

The reserved line is now a native raw-coordinate line `W=A e_j`, with a
complementary raw-coordinate summand.  No arbitrary A-linear basis change is
used before invoking independence of biased byte-reduced symbols.  This is
the condition under which the disjoint byte-region argument is valid.

### Exact official-target sampling

The target XOF chunks are modeled separately from key expansion.  Rejection
sampling is written explicitly for each bit chunk, and exact evaluation over
all nine official target lengths gives minimum acceptance
`0.1524622985...`, hence the stated lower bound `0.152462`.

### Rank and completion logic

- `rank(Q_R)=K` is an explicit public preflight for the zero-offset ell=2
  route.
- The three-plus-one proof now declares a trial incomplete whenever a
  consistent completion system is rank deficient; the uniqueness conclusion
  is no longer circular.
- Generic quotient rank is denoted `rho_Q`, avoiding collision with the byte
  atom `rho_byte`.

## Presentation and accounting corrections ported

- The ell=4 column formerly called `direct fallback` is now correctly called
  `adaptive mixture`.
- Extension-field ledgers explicitly charge pointwise tower multiplications;
  transforms, evaluation/interpolation, inversion, factorization, control,
  and memory remain in `kappa_hom`.
- The bibliography is after the appendices.
- The notation table distinguishes public parameters, slice dimensions,
  extension degrees, quotient ranks, and byte atoms.
- The Safey El Din--Schost DOI and author accent, Thinh Hung Dang's name, and
  several publication records were corrected.
- The Makefile and README now expose honest `verify` and `regenerate` targets,
  and artifact descriptions distinguish generated checks from external
  official-verifier tests that are not included.

## Newer material deliberately retained

The editor's source predates the referee-driven Just Guess research.  The
following newer material was preserved and rechecked rather than replaced:

- the Level-I Just Guess plus quotient-filter route;
- its expected-tree and capped-tree ledgers;
- the explicit `kappa_JG` and `kappa_JG_cap` regularity caveats;
- the referee-response and focused research-audit documents;
- removal of the unreproducible numerical orbit and streamed-XL claims.

## Validation

From a clean source directory:

- `make regenerate` passes;
- `make verify` passes;
- `latexmk` produces a 45-page PDF;
- the log has no warnings, unresolved references or citations, overfull
  boxes, or underfull boxes;
- all 45 pages were rendered and visually inspected;
- all-nine and Just Guess ledgers reproduce from the distributed generators.
