# Claim and scope statement

## Exact algebraic claims

For every odd-characteristic SNOVA public key whose scalar-expanded public base matrices are symmetric, the ordered power-pair feature map in the Beullens affine-column attack factors through

\[
\bigoplus_{i=1}^{m_1}\operatorname{Sym}^2_{\mathbb F_q}(\mathbb F_q[S]).
\]

Its homogeneous rank is therefore at most

\[
K=m_1\binom{\ell+1}{2},
\]

independently of the emulsifier choices and independently of any weak-key event. Left-kernel elimination yields the residual MQ systems stated in the paper. Fixed affine offsets deterministically bypass the symmetric-block rejection without changing the homogeneous quotient.

## Exact computational certificates

The recorded Level-I audit reconstructs the official Version 2.3 KAT public key byte-for-byte, verifies the 80-by-50 quotient rank, obtains 30 independent affine constraints, and checks exact equivalence with 50 quadratics in 102 variables. The official feature-map certificate has four rank-50 self maps, six invertible 80-by-80 cross maps, and full 80-by-680 rank. An explicit two-column restriction leaves 50 self-quadrics on a 52-dimensional affine space.

## Heuristic complexity claims

The displayed gate exponents use the same finite-field semi-regular, Hashimoto/Wiedemann, sparse-linear-algebra, and field-to-gate methodology as the public SNOVA analysis. They are comparative parameter estimates, not proven fixed-parameter running-time upper bounds and not measured attacks.

All nine implemented Version 2.3 q=19 shapes are estimated below their category targets under the stricter current-code convention. The six Version 2.4-preview rows are conditional on a corresponding implementation retaining symmetric public matrices; they are not claimed to be final Round-3 parameters.

## Not claimed

The work does not claim a practical production-size forgery, key recovery, a break of the unsymmetrized q=16 parameters by this quotient, an attack on a definitive NIST-hosted Round-3 package, or any stronger solver extrapolation omitted from the public paper.
