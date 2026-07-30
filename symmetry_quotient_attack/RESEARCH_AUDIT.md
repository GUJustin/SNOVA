# Focused strengthening audit after the Mode-A report

## Objective

Replace or sharply narrow the unknown homotopy multiplier on the Level-I
bottleneck, while preserving one-sided forgery soundness and honest resource
accounting.

## Routes investigated

- Random projected square cores: valid but no exponent improvement.
- Direct Just Guess on the full 48-equation F19 quotient: too expensive.
- Mixed inclusion of cross-channel equations: no favorable classical
  tradeoff in the scanned feasible parameter range.
- Worst-case binary-tree Just Guess: finite and reproducible, but leaves only
  2.739 bits for the capped candidate-regularity multiplier.
- Expected-tree Just Guess on the 16 A-valued diagonal equations, streaming
  an A^8 family and filtering through the cross channel: best result.

## Best retained result

For A=F_{19^2}, `(n,m,k,p)=(64,16,5,6)` satisfies the published Just Guess
structural inequalities.  With the verified 692-gate A multiplier, 84-gate A
adder, and a 16,886-gate all-roots routine for monic quadratics over A, the
published expected schedule gives 3,179,584 AXN gates per guess.  Streaming
A^8 assignments and all A^5 guesses gives exponent 132.046522 before success
repetitions.

For any set of `c*19^16` sign-canonical diagonal roots, a fresh cross channel
and uniform cross target satisfy

    Pr[match] >= c / (1 + c*(19*14/256)^16).

At c=1/4 this is at least 0.1710527744, adding 2.547487 bits.  Hence the
success-adjusted expression is

    134.594009 + log2(kappa_JG).

## Remaining gap

The transformation-selected search family depends on the diagonal system.
The current proof does not establish that one trial produces at least
`19^16/4` sign-canonical roots with a bounded expected number of
transformation/linearization retries.  This is isolated as `kappa_JG`.
Attempts to remove it by a naive diagonal second moment fail because the
search family is transformation-dependent.  A proof would need either:

1. a fresh-coefficient decoupling theorem for the Just Guess transformation;
2. a distributional invariance theorem for the transformed random system; or
3. a fixed-key production transcript measuring candidate abundance and
   restarts.

No stronger unconditional claim is made.
