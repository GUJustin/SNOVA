# Reduced-parameter end-to-end SNOVA forgery

This directory contains an executable demonstration of the symmetry-quotient
attack on the reduced row

```
(v,o,q,l,r) = (2,1,19,2,2).
```

It is a clean Python transcription of the odd-prime, symmetric, `l=r=2` path
of the SNOVA reference implementation pinned at commit
`9da14981336ede257c41ef53cc069989051e8181`.  It preserves the reference
implementation's q=19 field, public S matrix, fixed-ABQ expansion, hidden-UOV
public-key generation, base-19 public-key/signature serialization, verifier
map, target hash, salt, and rejection rule.  Only the dimensions are reduced.

The attack receives only the serialized public key and a chosen digest.  It
imposes the zero-offset common-column relation, interpolates the complete
restricted verifier, finds its rank-3 quotient and one target-consistency
coordinate, searches salts, enumerates public 3-dimensional slices, serializes
the result, and checks it using a separately organized direct verifier.

Run:

```bash
python3 attack.py --out transcript.json
python3 verify_transcript.py
python3 run_batch.py
```

The committed transcript records an 18-byte public key and a 23-byte forged
signature accepted by the direct verifier.  The verifier also rejects a
mutated signature and a mutated salt.  `batch_results.json` records successful
forgeries on eight additional deterministic keys.

## Scope

This closes the reduced-parameter composition gap: public-key expansion,
quotient extraction, target consistency, solving, signature reconstruction,
serialization, rejection, and verification all execute together.  It is not a
production-size attack, does not measure `kappa_hom` or `kappa_JG`, and is not a
compiled invocation of the upstream C implementation.  The verifier and
key-generation logic are independent Python transcriptions of the pinned C
formulas, with two differently organized public-map evaluators cross-checked
on unrestricted and restricted inputs.
