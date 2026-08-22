# Reduced-parameter end-to-end forgery report

## Result

The artifact executes the symmetry-quotient forgery pipeline on the reduced
SNOVA row

```
(v,o,q,l,r) = (2,1,19,2,2).
```

It preserves the pinned odd-prime symmetric reference formulas and byte
formats, changing only the dimensions. The attack receives a serialized public
key and a chosen digest; it receives no secret-key material.

## Executed pipeline

1. Parse and expand the serialized public key.
2. Restrict signatures to the zero-offset common-column relation.
3. Interpolate the complete restricted public verifier.
4. Compute the rank-3 symmetry quotient and its one consistency coordinate.
5. Search 16-byte salts until the target passes consistency.
6. Enumerate a public three-dimensional slice over F19.
7. Reconstruct and serialize the signature.
8. Verify with a separately organized direct verifier.
9. Reject a signature mutation and a salt mutation as negative controls.

## Committed transcript

- Public-key size: 18 bytes.
- Signature size: 23 bytes.
- Quotient rank: 3.
- Target-consistency dimension: 1.
- Salt counter: 4.
- Slices tested: 1.
- Independent direct verifier: accepts.
- Canonical serialization: yes.
- Secret material passed to attack: no.

Eight additional deterministic toy keys were forged successfully. The checker
also reproduces the principal transcript bit for bit.

## Scope

This is an end-to-end reduced-parameter instantiation. It validates the
composition from public-key bytes through quotient extraction, solving,
reconstruction, serialization, rejection, and verification. It is not a
production-size attack, does not measure the homotopy or Just Guess multipliers,
and is not a compiled invocation of the upstream C code. The implementation is
a clean Python transcription of the pinned formulas with independently
organized attack and verifier evaluators.

## Reproduction

From the package root:

```
python3 artifact/reduced_attack/verify_transcript.py
```

To regenerate the main and batch transcripts:

```
python3 artifact/reduced_attack/attack.py --out artifact/reduced_attack/transcript.json
python3 artifact/reduced_attack/run_batch.py
```
