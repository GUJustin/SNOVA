# Fixed-core-free orbit-attack density certificate

The structured source and projection preflights succeed with the exact lower bound in `STRUCTURAL_PREFLIGHT_PROBABILITY_CERTIFICATE.md`.  The orbit sweep needs no fixed core-Jacobian event: at every descent root where the full restricted Jacobian has rank `h`, one of the 16 Frobenius-orbit representatives contains a nonsingular square core.  The only internal-coefficient failure event is therefore the seeded aggregate-spectrum tail.

| Shape | `(h,K)` | 55 -> 16 | compact per key | compact normalized | margin | robust per key | robust normalized | margin | peak one-core output |
|:--:|:--:|:--:|--:|--:|--:|--:|--:|--:|--:|
| I-a | (48,50) | 55 -> 16 | 137.740 | 137.826 | 5.174 | 138.324 | 138.410 | 4.590 | 1.890 PiB |
| I-b | (48,50) | 55 -> 16 | 137.740 | 137.826 | 5.174 | 138.324 | 138.410 | 4.590 | 1.890 PiB |
| III-a | (68,70) | 55 -> 16 | 179.175 | 179.260 | 27.740 | 179.758 | 179.844 | 27.156 | 2.076 ZiB |
| III-b | (68,70) | 55 -> 16 | 179.175 | 179.260 | 27.740 | 179.758 | 179.844 | 27.156 | 2.076 ZiB |
| V-a | (88,90) | 55 -> 16 | 220.249 | 220.334 | 51.666 | 220.832 | 220.918 | 51.082 | 2416.552 YiB |
| V-b | (88,90) | 55 -> 16 | 220.249 | 220.334 | 51.666 | 220.832 | 220.918 | 51.082 | 2416.552 YiB |

The normalized exponent pays the inverse certified key density and is therefore a time-success-ratio ledger over a random generated public key.  The dense parametrization ceilings are not low-memory claims; this route is included to remove fixed-core and selected-Macaulay-degree assumptions, while the eigenblock-core and streamed-XL branches supply lower-output or lower-state points.
