# Response to the 31 July 2026 adversarial reviews

This response records revision decisions. It is not evidence for the paper's
claims.

The criticism assessment is right that one adversarial report overread the
manuscript: the paper specified a conditional algorithm and expressly
disclaimed a completed production-size run. The narrower artifact criticism
was nevertheless correct: the prior release checked many constituent formulas
but did not instantiate the reduction-and-forgery composition.

The revision makes four substantive changes:

1. The symmetric-square quotient is now stated only after the attacker-imposed
   common-column relation. For the unrestricted verifier, swapping power labels
   also transposes the signature-column indices.
2. The structural companion is corrected accordingly. Its unrestricted
   per-base cap is `binom(d*r+1, 2)`, nonbinding for all nine published shapes;
   it no longer claims an unrestricted SNOVA no-go theorem or concrete break.
3. Claim labels distinguish exact reductions, idealized-transcript theorems,
   the unverified inherited homotopy hypothesis `H_hom`, conditional recovery
   algorithms, cost sensitivities, and supplied implementation evidence. The
   evidence table now reports availability.
4. The artifact adds two composition tests: a public-key-only Level-I KAT
   verifier/reduction harness and a genuine non-planted end-to-end forgery at
   an unofficial reduced shape. The reduced forger and literal evaluator use
   distinct paths in one shared Python transcription, so this is composition
   evidence rather than independent verifier correspondence. Neither test is
   presented as production-size feasibility evidence.

The official-key harness also found and corrected a flaw in an archived KAT
script: the official `be_invertible_by_add_aS` loop updates its matrix
cumulatively, whereas the archived transcription reset to `M + aS` on each
iteration. The archived evaluator matched only 65 of 80 KAT target
coordinates. It is not reused as evidence.

An additional supplied zero-offset demo was also reviewed before integration.
Its committed transcript was not Version 2.3-correct because it used direct
SHAKE256 rather than the indexed SHAKE128 public expansion and passed raw text
instead of the message interface's 64-byte digest. The release does not retain
that duplicate transcription. It ports the corrected zero-offset attack logic
onto the KAT-anchored helpers, recomputes the main and eight batch trials, and
adds 24-forgery/100-rank stress checks. This strengthens only the existing
toy-size composition evidence; it does not change the security claim.

The remaining limitations are explicit: no production solver or forgery, no
proof or test of `H_hom` for a published residual ideal, no complete Boolean
homotopy implementation, no bound on `kappa_hom`, no proved Just Guess
candidate abundance, no all-shape official-key preflight census, no fixed-key
spectrum certificate, and no official-C-verifier forgery transcript.
