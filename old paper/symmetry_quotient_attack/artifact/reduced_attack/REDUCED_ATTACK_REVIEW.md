# Integration note for the supplied reduced-attack demonstration

This note records why `reduced_attack.zip` was not copied into the release
verbatim. It is a revision record, not independent evidence. The reviewed ZIP
has SHA-256
`b63a571218d061e2d8662615babaa8042731d3a0bfa060511269a048bbaa3de1`.

## Defects in the supplied files

| Supplied location | Finding | Concrete consequence | Release correction |
|---|---|---|---|
| `reduced_attack/attack.py:71-72,183-205` | The public seed is expanded by a direct SHAKE256 call. Line 184 also allocates 20 symbols although lines 190-204 consume 18. | The committed public matrices differ from the pinned Version 2.3 matrices in 32 of 36 scalars. Under the KAT-anchored transcription, the committed signature evaluates to `[18,10,8,3]`, not its recorded target `[1,6,11,3]`. | The integrated test uses `snova_v23_reference.py:202-209,320-358`: indexed 168-byte SHAKE128 blocks over `seed || LE64(block)` and the exact 18-symbol count. |
| `reduced_attack/attack.py:379-382,627-635` | The target helper accepts a digest, but the CLI passes a 40-byte text string directly. | The committed run is a digest-level toy with a nonstandard digest length, not the message-level Version 2.3 interface. | The integrated test calls `snova_v23_reference.py:647-654`, which computes the 64-byte SHAKE256 message digest before the target hash. |
| `reduced_attack/verify_transcript.py:47-50,61` | Parsed JSON equality is described as “bit-for-bit” reproduction. | The committed CRLF file and regenerated LF file are semantically equal but bytewise different. | The integrated checker emits and compares one deterministic canonical JSON representation. |
| `reduced_attack/verify_transcript.py:52-54` | The eight-row batch ledger is trusted rather than recomputed. | A forged `success: true` field would pass this part of the checker. | `reduced_parameter_end_to_end_forgery.py --check` reruns the main trial and all eight fresh-key trials before comparing the committed output. Any failed trial exits nonzero. |
| `reduced_attack/attack.py:617-620` | Several success and input-boundary fields are literal constants. | Those fields are descriptions, not checks. | The integrated output is produced only after executed positive, negative, rank, image, interpolation, and reconstruction checks; the public attack boundary is also visible in the `forge(public_key, message, rng, ...)` signature. |
| `reduced_attack/README.md:38-40` | The text says the new script “closes” the reduced composition gap. | The preceding release already contained a non-planted end-to-end test at the same unofficial shape. | The manuscript describes the incremental value accurately: a zero-offset realization with complete restricted-map interpolation, explicit rank-three quotient/image agreement, one-coordinate target filtering, and public three-dimensional slice search. |

The supplied Python bytecode cache is omitted.

## What survived adversarial testing

After changing only the public XOF/count and the message prehash, the supplied
zero-offset attack logic continued to forge. The integrated version therefore
ports that attack logic onto the package's KAT-anchored scheme helpers rather
than retaining a second 650-line scheme transcription. It checks:

- equality of staged and literal public-map evaluators on unrestricted inputs;
- exact interpolation of the complete zero-offset restricted verifier;
- rank three for both the interpolated coefficient image and the explicit
  `rho=(1,0)` symmetric quotient, plus equality of those output images;
- one target-consistency coordinate before any root search;
- public three-dimensional slice enumeration, reconstruction, canonical
  base-19 serialization, and the `ell=2` rejection rule;
- acceptance of the main forgery and eight fresh-key forgeries, and rejection
  after changing the message, public key, signature, or salt.

The separate executable stress checker additionally passed 24/24 fresh-key
forgeries and obtained restricted-map rank three on 100/100 further reduced
keys; it also round-tripped base-19 vectors of every length from 1 through 100.
The committed corrected main trial has public key
`e74b264b7fadf91c87a567a16396f6daa904`, signature
`c4b76522ef5e0004000000000000000000000000000000`, salt counter 4, and
verifier target/output `[2,15,17,2]`.

This repaired result is genuine toy-size composition evidence at
`(v,o,q,ell,r)=(2,1,19,2,2)`. It remains neither a published-parameter forgery
nor evidence about production solver feasibility, `H_hom`, `kappa_hom`,
`kappa_JG`, memory, or a complete gate count. The verifier path is a distinct
literal evaluator in the shared KAT-anchored Python transcription; it is not an
invocation of the upstream C verifier.
