# IACR ePrint submission metadata

## Title

Symmetry-Quotient Forgery Attacks on Odd-Characteristic SNOVA

## Author

Justin Thaler

## Affiliations

a16z crypto research; Georgetown University (on leave)

## Suggested category and topics

- Attacks and cryptanalysis
- Post-quantum cryptography
- Multivariate cryptography

## Abstract

The public odd-characteristic variants of SNOVA use symmetric scalar-expanded public quadratic matrices. We show that this symmetry creates a universal quotient inside the affine-column forgery framework of Beullens. The existing analysis assigns one homogeneous feature coordinate to every ordered power pair (a,b). For symmetric S and P_i, the (a,b) and (b,a) quadratic evaluations coincide. Thus, with A=F_q[S] viewed as an l-dimensional F_q-space, the homogeneous map factors through Sym^2_{F_q}(A) and has rank at most m_1*binom(l+1,2), rather than m_1*l^2.

For the nine q=19 parameter sets in the public Version 2.3 draft, the quotient gives residual systems ranging from 48 quadratics in 112 variables to 96 quadratics in 224 variables. Under the same semi-regular Hashimoto/Wiedemann methodology used in the SNOVA analysis, with the stricter convention in the public analysis code, the estimated costs range from 2^130.30 to 2^238.77 gates. Every set is below its nominal NIST category, with estimated shortfalls from 4.06 to 44.05 bits. Fixed affine offsets bypass the verifier's rejection of symmetric signature blocks without changing the homogeneous reduction. A source-level reconstruction of the official Level-I KAT gives exact equivalence between the 80 verification equations and 50 quadratics in 102 variables; an independent two-column restriction gives 50 self-quadrics on a 52-dimensional affine space.

The algebraic reductions and public certificates are exact. The displayed work factors remain heuristic, as are the corresponding estimates in the SNOVA materials. We therefore claim a security-category or parameter break of the public odd-characteristic drafts under the candidate's accepted methodology, not a practical full-size forgery. This quotient does not apply to the unsymmetrized q=16 public matrices.

## Keywords

SNOVA; multivariate signatures; post-quantum cryptography; forgery attack; symmetric square; underdetermined MQ; NIST PQC

## Suggested comments

20 pages, 4 tables, 1 figure. Includes a reproducibility artifact with exact finite-field rank certificates, the official Level-I KAT verifier, and cost-table code. Version 2.4-preview estimates are conditional on retaining symmetric public matrices and are not claimed to describe a definitive Round-3 package.

## Files

- PDF: `snova_symmetry_quotient_eprint.pdf`
- Source archive: `SNOVA_Symmetry_Quotient_Eprint_Source_2026-07-19.zip`
- Artifact archive: `SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19.zip`
