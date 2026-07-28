# Third Full Correctness and MQ-Accessibility Audit

## Bottom line

I performed another theorem-by-theorem, proof-dependency, numerical, literature,
cost-model, source, and presentation audit, using the previous 32-page
readability revision as the baseline.

I did not find a fatal flaw in the central cryptanalysis, and none of the
headline exponents changed. The symmetric-square quotient, exact affine
reduction, stable-kernel construction, fixed-skew identity, rejection-filter
bounds, collision and retry formulas, and Las Vegas verification wrapper all
survived independent rederivation.

The revision is now 37 pages. The added pages are not padding: they replace
specialist shorthand with a running example, a self-contained MQ primer,
operational explanations of the linear algebra, and explicit separation of
exact statements from rank assumptions, solver assumptions, and circuit-cost
conversions.

The principal remaining qualifications are substantive rather than editorial:

1. the production global polar and cross-polar rank spectra are supported by
   large finite campaigns, not exhaustively certified over every projective
   direction;
2. production-size special-XL and generic-MQ behavior is modeled rather than
   established by full-size solver runs; and
3. the AXN/NAND conversions price the arithmetic schedule charged by the cited
   estimators, not a complete implementation including pivoting, control,
   memory, indexing, and data movement.

## New correctness and precision fixes

### 1. The symmetry theorem gives a cap; full quotient rank is a separate check

The symmetric-square theorem proves that the ordered feature map has rank at
most

`K = m_1 binom(d+1,2)`.

It does not by itself prove that the concrete relation-dependent matrix has
rank exactly `K`. The paper now says this at the first informal statement of
the result and consistently treats `rank Q_R = K` as an exact public preflight
condition. This avoids turning an upper bound into a generic-equality claim.

### 2. The affine base point is target-dependent

The solution of the affine constraints is now written `a(t)` rather than `a`.
This makes clear that conditioning on the affine coordinates of a uniform
target changes the translate but not the residual quadratic polar maps. The
probability argument no longer relies on the reader silently carrying this
dependence.

### 3. Fixed skew ensures acceptance, not root existence

The fixed-skew identity proves that every block is nonsymmetric for every point
satisfying the column relation. Together with full affine rank, every target
passes the *linear consistency* stage and leaves a residual MQ system. It does
not alone imply that the residual system has a root. The proposition and its
operational explanation now state this distinction explicitly; root existence
is supplied separately by the polar-rank analysis.

### 4. Three different success events are kept separate everywhere

The paper now distinguishes, from the MQ primer onward:

1. at least one base-field root exists;
2. exactly one base-field root exists on the selected affine slice; and
3. the chosen dehomogenized Macaulay matrix has a one-dimensional right kernel
   from which the root can be read.

The moment and rank-spectrum theorems address the first two. The third remains
part of the production special-XL model. Table headings now say
"all-target root" rather than language that could be read as "one solver call
succeeds."

### 5. The ordered-label comparison is described as an estimator baseline

The verifier really does compute an ordered-label feature vector, so the
baseline is operationally meaningful. The counterfactual is treating the two
reversed labels as independent generic features in the solver estimator after
the public symmetry identity is known. The appendix now states exactly this,
rather than suggesting that the ordered coordinate system itself is fictitious.

### 6. The Macaulay-kernel caveat is more precise

The old wording could suggest that "syzygies" automatically create right-kernel
vectors. The revision instead identifies the actual unresolved mechanisms:
additional extension-field roots or failure of the selected Macaulay rows to
span all expected relations at the modeled multidegree. The recovery lemma,
right-versus-left-kernel convention, and Las Vegas checks are now explicit.

### 7. The `q^tau` retry heuristic is scoped to the selected positive deficits

For the cost-driving choices `tau=10` and `tau=14`, the exact full-rank retry
factor is essentially `q^tau`. The unused edge case `tau=0` behaves
differently. The explanatory text now restricts the heuristic to the positive
values actually used.

### 8. Notation collisions were removed

Several overloaded symbols were individually standard but collectively costly
for nonexperts. The revision now uses:

- bold `v` for the offset tuple, rather than the same letter used for a vector
  space or the SNOVA vinegar count;
- `mathcal F_{R,v}` for the feature-linear matrix, rather than `L`, which is
  reserved for an affine-slice direction space;
- `D_j` for structured offset multipliers, rather than `T`, which is reserved
  for directional derivatives;
- `Phi_{b,L}` for the cross-polar map, rather than `A_b`, which collided with
  the public algebra `A`; and
- `t_rej` for the rejection threshold, rather than `t`, which is the target.

The notation table and every theorem were updated consistently.

### 9. A real bibliography error was corrected

The EUROCRYPT 2026 exterior-algebra paper is by **Lars Ran**, not Kexin Ran.
The entry is corrected. The Hashimoto title correction from the preceding pass
is retained, and all 19 bibliography records now pass Biber's datamodel
validation.

## Readability reconstruction for readers new to MQ

### The abstract now gives the attack in ordinary language

The abstract first says what an MQ forgery is, then explains the duplicated
ordered features, gives the 16-to-10 and 4-to-3 collapses, and only then states
the construction and costs. It closes by naming the exact parts and the
remaining modeling assumptions, rather than asking a reader to infer the scope
from technical qualifiers later in the paper.

### The introduction has a running attack narrative

The opening pages now contain:

- the public map `P: F_q^N -> F_q^M` and the equation a forger must solve;
- a four-label toy example in which `(1,S)` and `(S,1)` are the same polynomial;
- a boxed statement of the core symmetry vulnerability;
- an attack data-flow table from public key to accepted forgery;
- four clearly named jobs: quotient, affine completion, solver-compatible
  slice, and verifier acceptance;
- a complete Level-I dimension example, tracking the raw feature space,
  quotient, affine constraints, stable kernel, selected slice, and retry
  deficit; and
- a claim-status table separating exact algebra, public preflight, finite
  evidence, and modeled costs.

This lets an MQ newcomer know what the paper is trying to accomplish before
encountering tensor products, invariant kernels, or Hilbert series.

### A self-contained MQ primer now precedes the attack machinery

The primer explains:

- finite-field MQ maps, roots, fibers, and affine slices;
- why "more variables than equations" does not by itself give an efficient
  algorithm;
- homogeneous quadratic parts and polar forms;
- how polar rank controls deviations from random-map fiber behavior;
- projective coefficient directions and why scalar multiples need not be
  retested;
- XL as multiplying equations by monomials and linearizing the resulting
  monomial columns;
- the construction and orientation of a Macaulay matrix;
- why a root gives a right-kernel evaluation vector;
- why multidegrees and separate homogenizers reduce the special-XL matrix; and
- why a unique base-field root is not the same statement as a one-dimensional
  Macaulay kernel.

A notation table immediately before the technical development defines every
recurring space and parameter.

### The affine-column framework is explained operationally

The revision now tells the reader what each matrix does:

- the column relation replaces all signature columns by affine functions of one
  common vector;
- substitution creates quadratic, linear, and constant feature parts;
- `Q_R` records the surviving quadratic output directions;
- a left-kernel basis kills those directions and exposes ordinary affine
  constraints; and
- solving those constraints leaves a residual MQ system on an affine translate.

A two-output toy elimination example illustrates the quotient-space argument
behind the exact rank formula.

### The stable lift is presented as a design problem, not a magic formula

The offsets have two competing jobs: fill the missing affine output directions
and preserve an `A`-linear direction space on which the offset-linear feature
terms are constant. The text now explains the four-step path

`structured offsets -> A-stable row space -> A-linear kernel -> four eigenblocks`

before stating the formal theorem. It also explains the splitting field,
primitive idempotents, and block multidegrees at the point where they are used.

### The probability section begins with the random-map baseline

Before the exact moment theorem, the reader is shown that a random map on an
`h`-dimensional slice targeting `K` coordinates should have mean fiber size
`mu=q^(h-K)`. The real issue is then framed as excess collisions, not merely
expected root count. The factorial second moment is described as an exact count
of ordered colliding pairs, and each corollary is labeled by the operational
question it answers: nonempty fiber, unique root, planted-versus-conditioned
sampling, or all-target availability.

### The cost section separates four layers

A new interpretation table distinguishes:

1. the exact reduced MQ instance;
2. the cited solver's field-operation model;
3. conversion of the charged arithmetic to AXN or NAND gates; and
4. costs omitted from a complete solver implementation.

The paper explains base-2 exponents, arithmetic overhead budgets, and the role
of the AES calibration without calling the converted arithmetic a complete
attack circuit.

### The evidence section decodes experimental jargon

The evidence hierarchy now distinguishes public preflight, finite rank
campaigns, tractable Macaulay experiments, and circuit-primitive validation.
Terms such as coefficient-basis direction, weight-two direction, dense
direction, polar radical, homogeneous Jacobian, planted-solvable matrix, and
projective distinctness are explained when first used.

## Mathematical audit performed

The following were rederived independently from the definitions in the paper:

- `rank[Q_R F_R,v] = rank Q_R + rank(W_R F_R,v)`;
- the symmetric-square factorization under transpose symmetry;
- closure and dimension of the structured offset row space;
- `A`-linearity and dimension of the stable kernel;
- preservation of block multidegrees under translation and blockwise
  homogenization;
- the fixed-skew identity for the official `ell=2` matrix;
- the accepted-target union and double-counting bounds;
- the first and factorial second moments of the number of roots in a random
  target-coset pair;
- the rank-spectrum identity obtained by double-counting `(b,y)` annihilator
  pairs;
- the variance, Paley-Zygmund, singleton, conditional multiple-root, and
  planted-distribution inequalities;
- equivalence of full cross-polar rank and surjectivity of every nonzero
  directional derivative;
- the full-space character-sum availability bound; and
- exact recovery from a one-dimensional dehomogenized Macaulay right kernel,
  followed by the public-verifier Las Vegas wrapper.

No theorem statement had to be weakened beyond the precision changes described
above.

## Numerical checks rerun

Independent calculations reproduce:

- special-XL model exponents:
  `128.81713295`, `171.68561994`, `214.12301141`;
- special-XL AXN exponents:
  `132.17186724`, `175.04035423`, `217.47774570`;
- special-XL NAND exponents:
  `133.04030634`, `175.90879334`, `218.34618481`;
- fixed-skew generic-MQ model exponents:
  `129.44029908`, `183.43802955`, `237.84309628`;
- Level-V filtered generic-MQ model exponent:
  `235.19473743`;
- corresponding AXN exponents:
  `133.74642790`, `187.74415790`, `242.14922790`, and
  `239.50086790` for the filter;
- fixed-skew arithmetic-budget rank thresholds `46`, `68`, and `89`;
- fixed-skew all-target-root thresholds `96`, `144`, and `192`;
- filtered arithmetic-budget and all-target-root thresholds `177` and `240`;
- accepted-fiber fraction `0.9999999990155023` at filter rank `254`;
- conditional multiple-root logarithms below `-43.479` bits for `tau=10` and
  below `-60.471` bits for `tau=14`; and
- planted-distribution total-variation logarithms below `-42.479` and
  `-59.471` bits in those cases.

The independent bounded Hilbert-series enumeration checked `3,690`, `11,217`,
and `27,449` candidate multidegrees for `m_1=5,7,9`, respectively. It found no
zero trigger coefficients and recovered the same unique optima up to block
permutation:

- `(s,d)=(10,(3,3,3,4))`, coefficient `-14,677,784`;
- `(s,d)=(15,(4,4,6,6))`, coefficient `-124,983,335,190`; and
- `(s,d)=(19,(5,6,6,6))`, coefficient
  `-1,164,044,975,088,360`.

## Literature and bibliography audit

The references were checked against primary publisher, author, NIST, and IACR
records where available. Besides correcting Lars Ran's name, the revision adds
or updates the current third-round NIST context, the recent `Just Guess`
solver, generalized variable partitioning, and the Sakata-Takagi
underdetermined-MQ result. The text does not claim to implement generalized
variable partitioning; it presents the chosen generic-MQ estimates as pinned,
reproducible baselines.

## Source and PDF quality checks

- 37-page US Letter PDF.
- Clean normal build: no unresolved citations or references, LaTeX/package
  warnings, overfull boxes, or underfull boxes.
- Optional anonymous and coordinated-disclosure branches also compile cleanly.
- 108 labels, all unique, with no missing references.
- 19 bibliography entries, all cited, with no missing or duplicate keys.
- Biber datamodel validation passes without warnings.
- All fonts are embedded and subset.
- Text extraction contains no replacement characters or NUL bytes.
- All 37 pages were rendered at 180 dpi with both PDFium and Poppler and
  visually inspected.
- The two renderers agree on every page's nonwhite-content bounding box to
  within one pixel.
- No clipped text, overlap, malformed glyph, broken table, bad float placement,
  or anomalous blank page was found.

## External limitation

The uploaded materials do not contain the complete evidence directory or the
released arithmetic netlists described in the artifact appendix. I could audit
the mathematics, source, bibliography, dimensions, cost formulas, independent
numeric calculations, and PDF, but I could not independently rerun the full
key-generation campaign, all production rank records, the complete verifier
regression suite, or the netlist-level validation from the stated release.
