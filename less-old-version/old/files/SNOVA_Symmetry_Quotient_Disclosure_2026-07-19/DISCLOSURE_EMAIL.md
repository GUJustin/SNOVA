Subject: Coordinated disclosure: symmetry-quotient forgery attack on odd-characteristic SNOVA

To: SNOVA submission team  
Cc: NIST PQC technical contacts

Dear SNOVA team and NIST colleagues,

I am writing to disclose a structural forgery reduction affecting the nine q=19 parameter shapes implemented in the public SNOVA Version 2.3 draft, and conditionally affecting the six Version 2.4-preview shapes if their eventual implementation retains symmetric scalar-expanded public matrices.

Version 2.3 analyzes the Beullens affine-column attack using one homogeneous feature coordinate for each ordered power pair `(a,b)`, for a total of `m1*l^2` coordinates. For every symmetric public base matrix `P_i` and symmetric power matrix `S`, however,

```text
u^T (I_n tensor S^a) P_i (I_n tensor S^b) u
  = u^T (I_n tensor S^b) P_i (I_n tensor S^a) u.
```

Equivalently, with `A = F_q[S]` viewed as an `l`-dimensional `F_q`-space, the homogeneous feature map factors through `m1` copies of `Sym^2_{F_q}(A)`. Its rank is therefore at most

```text
K = m1 * binom(l+1,2),
```

not `m1*l^2`. This is a universal identity for every key generated under the symmetric-matrix rule, not a weak-key event or a search for an exceptional low-rank relation.

After eliminating the complementary affine-linear output combinations, the Version 2.3 shapes reduce to:

- Level I: 50 quadratics in 102 variables; 48 in 112; 50 in 98.
- Level III: 70 in 146; 72 in 168; 70 in 142.
- Level V: 90 in 182; 96 in 224; 90 in 178.

Using the stricter `p1=1` convention in the current public analysis code and the same Hashimoto/Wiedemann semi-regular methodology as the submission, the estimated gate costs are:

- Level I: 2^138.94, 2^130.30, and 2^138.94;
- Level III: 2^184.29 for all three shapes;
- Level V: 2^227.95, 2^238.77, and 2^227.95.

Thus all nine q=19 Version 2.3 shapes are estimated below their nominal categories, with shortfalls from 4.06 to 44.05 bits. These are heuristic algebraic-complexity estimates in the same sense as the submission's estimates; I have not carried out a production-size forgery.

The attack also handles the verifier's format check. For square signature blocks, public affine offsets can force a fixed nonzero skew entry in every block for every residual assignment, without changing the homogeneous reduction. Rectangular blocks are not subject to a symmetry test.

The attached artifact records two exact Level-I normal forms. First, the source-level KAT audit reconstructs the public key byte-for-byte and verifies exact equivalence between all 80 verifier equations and 50 quadratics in 102 variables after 30 affine eliminations. Second, the official feature map has six invertible 80-by-80 cross-column blocks; an explicit two-column restriction leaves 50 self-quadrics on a 52-dimensional affine space.

I would appreciate your review of four points:

1. Do you reproduce the symmetric-square quotient and the residual dimensions?
2. Is this quotient already incorporated into a newer internal analysis not reflected in Version 2.3 or the public Version 2.4-preview scripts?
3. Do you see any verifier or signature-validity condition that invalidates the fixed-skew affine family?
4. Which semi-regular convention and concrete MQ cost model do you intend to use for the updated Round-3 package?

I propose an initial 14-day confidential review period, followed by a coordinated public-disclosure date. I am happy to provide the official-KAT harness directly and compare transcripts with your implementation.

Attached:

- confidential disclosure manuscript;
- public and anonymous manuscript builds;
- LaTeX source archive;
- reproducibility archive and recorded results;
- concise claim-and-scope statement.

Best,

Justin Thaler  
a16z crypto research  
Georgetown University (on leave)
