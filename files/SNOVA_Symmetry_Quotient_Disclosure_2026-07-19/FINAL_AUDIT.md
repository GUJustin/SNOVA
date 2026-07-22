# Final SNOVA ePrint audit

## Editorial disposition

The 25-page newer manuscript was used as the authoritative source. The public paper was reduced to 20 pages and now contains only:

- the universal symmetric-square quotient;
- the residual-MQ reduction and candidate-relative cost estimates;
- the deterministic verifier-format bypass;
- the official Level-I 50-in-102 KAT certificate;
- the official Level-I 50-in-52 two-column certificate;
- the explicitly conditional Version 2.4-preview table;
- reproducibility, novelty, countermeasures, and a compact claim checklist.

The following material was removed from the public paper rather than relegated to an appendix:

- the failed extension-field-linear MinRank shortcut;
- preliminary multibase solving-degree extrapolations;
- negative radical, adjoint-algebra, and low-rank searches;
- QR-UOV, MAYO, and TSUOV implications;
- other research-program dead ends.

## Correctness and source fixes

- Corrected the quoted Hashimoto formula from an undefined `m` to the paper's `m_1`.
- Tightened the symmetric-square theorem so that `wedge^2(A)` is identified as the kernel of the canonical ordered-to-symmetric quotient, not asserted to be the entire kernel of every induced public map.
- Clarified that the 50-in-52 construction is an exact verifier-preimage normal form and is not used for the rejection-safe headline work factor.
- Corrected the PKC 2026 symmetric-algebra authors, part number, DOI, and publication metadata.
- Completed the EUROCRYPT 2026 publication metadata for Ran's wedge paper.
- Replaced the wedge-attack citation with its Cryptology ePrint record.
- Reduced the bibliography to the sixteen entries actually cited.
- Rephrased the Round-3 scope to distinguish NIST's candidate listing from a NIST-hosted definitive specification package.

## Build and PDF audit

- Public PDF: 20 pages, named author metadata.
- Anonymous PDF: 20 pages, anonymous author metadata.
- Disclosure PDF: 20 pages, dated confidential notice.
- No undefined references or citations.
- No overfull or underfull boxes.
- Only a benign `amsmath` accent warning remains.
- All fonts are embedded and subset.
- PDF preflight passed for all three builds.
- All 20 public pages were rendered and visually inspected; the first pages of both variants were separately inspected.

## Reproducibility status

The KAT-independent rank and estimator scripts were rerun successfully in the release environment. The official Level-I KAT response file is not redistributed in the artifact; the package pins the audited repository commit and expected filename, includes the recorded exact KAT output, and includes executable scripts that rerun the byte-level and direct-verifier checks once that public input is supplied.
