# Final release audit

## Frozen headline

The fixed-core-free, selected-degree-free branch uses 16 Frobenius-orbit representatives rather than 55 separate family-count patterns. Its robust per-certified-key AXN exponents are 138.32403, 179.75835, and 220.83240. Charging the exact random-XOF structural-preflight lower bound greater than 0.9423222679 gives random-key-normalized exponents 138.40973, 179.84405, and 220.91810.

## Consistency checks

- All obsolete 139.50 / 187.84 / 222.06 headline values were removed from the source and all PDF variants.
- Abstract, introduction, main theorem, conclusion, and artifact appendix now use the same 16-orbit theorem and normalized costs.
- The separately conditioned low-output eigenblock-core figures remain 134.077, 196.246, and 247.432 per certified key.
- The discarded core-Jacobian probability calculation remains isolated under `artifact/rejected/` and is not cited as a theorem.

## Build checks

- Public, anonymous, and disclosure variants compile through Biber and repeated LaTeX passes.
- Final logs contain no unresolved citations or references, package/LaTeX warnings, overfull boxes, or underfull boxes.
- Public and anonymous PDFs contain 77 pages; disclosure adds only the intended confidentiality banner when page flow requires it.
- The public PDF was rendered page-by-page at 120 dpi and reviewed through eight contact sheets, including the new theorem and tables.

## Artifact checks

The following generators were rerun successfully in the frozen directory:

- structural-preflight probability;
- Frobenius-orbit sweep;
- normalized orbit-attack density;
- exhaustive eigenblock-core frontier; and
- eigenblock witness certificate.

## Remaining scientific qualification

The orbit-complete branch is theorem-backed in operation count but has extremely large dense parametrization ceilings; it is not presented as a low-memory implementation. The low-output eigenblock-core branch still conditions on an exact public nonzero-Jacobian preflight. The streamed-XL branch remains the low-state alternative and retains its explicitly stated production selected-degree qualification.
