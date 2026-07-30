#!/usr/bin/env python3
"""Independent deterministic checks for the final SNOVA paper artifact."""
from __future__ import annotations
import hashlib, json, math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = 19
RHO = Fraction(14, 256)
A4 = Q**4
REF = {'I': 143.0, 'III': 207.0, 'V': 272.0}
EXPECTED_SIMPLE = {
    'I': 142.6926676651965,
    'III': 185.5365868407971,
    'V': 233.84398678063826,
}
EXPECTED_OPT = {
    'I': 132.18807210764757,
    'III': 183.79462510032272,
    'V': 233.84398678063826,
}
EXPECTED_COMPACT = {
    'I': 132.18807210764757,
    'III': 185.4345363337288,
    'V': 237.78324097428373,
}

PARAMETER_ROWS = (
    (28, 5, 4, 4), (48, 16, 2, 2), (28, 4, 4, 5),
    (40, 7, 4, 4), (72, 24, 2, 2), (38, 5, 4, 5),
    (50, 9, 4, 4), (96, 32, 2, 2), (52, 6, 4, 6),
)


def eval_netlist(net, a, b):
    vals = [bool((a >> i) & 1) for i in range(5)] + [bool((b >> i) & 1) for i in range(5)]
    for op, x, y in net['gates']:
        xv = vals[int(x[1])] if x[0] == 'w' else bool(x[1])
        yv = vals[int(y[1])] if y[0] == 'w' else bool(y[1])
        vals.append(xv and yv if op == 'AND' else xv ^ yv if op == 'XOR' else not (xv ^ yv))
    return sum((1 << i) for i, j in enumerate(net['outputs']) if vals[j])


def circuit_checks():
    net = json.loads((HERE / 'f19_multiplier_netlist.json').read_text())
    payload = json.dumps(
        {'ninputs': net['ninputs'], 'gates': net['gates'], 'outputs': net['outputs']},
        separators=(',', ':'), sort_keys=True,
    ).encode()
    assert hashlib.sha256(payload).hexdigest() == net['sha256']
    assert net['gate_count'] == 150
    for a in range(Q):
        for b in range(Q):
            assert eval_netlist(net, a, b) == (a * b) % Q
    d = json.loads((HERE / 'field_tower_circuits.json').read_text())
    assert d['F19_2']['multiplication'] == 692
    assert d['F19_4']['multiplication'] == 2628


def target_sampling_checks():
    """Verify the exact rejection-sampling acceptance for official targets.

    The pinned q=19 expansion emits 15 field elements from each full
    eight-byte chunk.  A final partial chunk contains the remaining field
    elements in the minimum whole number of bytes.
    """
    probabilities = []
    for _, o, ell, r in PARAMETER_ROWS:
        digits = o * ell * r
        byte_length = math.ceil(8 * digits / 15)
        full_chunks, remaining_digits = divmod(digits, 15)
        probability = Fraction(1, 1)
        full_bound = (2**64 // Q**15) * Q**15
        probability *= Fraction(full_bound, 2**64) ** full_chunks
        if remaining_digits:
            remaining_bytes = byte_length - 8 * full_chunks
            bits = 8 * remaining_bytes
            bound = (2**bits // Q**remaining_digits) * Q**remaining_digits
            probability *= Fraction(bound, 2**bits)
        probabilities.append(probability)
    minimum = min(probabilities)
    assert abs(float(minimum) - 0.1524622985242944) < 1e-16
    assert minimum > Fraction(152462, 10**6)


def separator_optimum(info, base_size, B, D, gate):
    """Recompute the conservative finite-cardinality B^2 bad-vector ledger."""
    H = B * B
    candidates = []
    for r in range(1, 300):
        E = base_size**r
        if E < D or E <= H or 2 * r - 1 > base_size:
            continue
        succ = 1 - Fraction(H, E)
        candidates.append((Fraction((2 * r - 1) * gate, 1) / succ, r, E, succ))
    z = min(candidates)
    assert z[1] == info['extension_degree']
    assert z[0] == Fraction(*info['gate_factor'])
    assert H == info['separator_bad_hyperplanes']


def spectrum_checks(rep):
    for row in rep['rows'].values():
        for key in ('diagonal', 'complete', 'fast', 'complete_square'):
            if key not in row:
                continue
            z = row[key]
            eta = Fraction(*z['eta'])
            h, K = z['h'], z['K']
            assert eta - Fraction(Q**h - 1, Q**K) == Fraction(1, 2**128)


def dimension_checks(rep):
    """Check that the direct complete-square slices really fit by dimensions."""
    for row in rep['rows'].values():
        v, o, d, r = row['parameters']
        if d == 4:
            M = o * r * d
            m1 = (o * r + d - 1) // d
            K = m1 * d * (d + 1) // 2
            lower = d * (v - 1) - (M - K)
            assert lower >= K
            assert row['complete_square']['h'] == K
        else:
            m = o
            K = 3 * m
            assert m % 2 == 0
            assert 2 * row['complete']['s'] == K
            assert row['complete']['h'] == K


def route_checks(rep):
    for row in rep['rows'].values():
        if 'complete' in row:
            m = row['parameters'][1]
            z = row['diagonal']
            separator_optimum(z['homotopy'], Q**2, 2**z['s'], 2 * z['s'], 692)
            z = row['complete']
            B = 2**(2 * m) * math.comb(m, m // 2)
            separator_optimum(z['homotopy'], Q**2, B, 2 * z['s'], 692)
            assert row['compact_output']['fast_branch_probability_lower'] > 0.997
        else:
            z = row['fast']
            s, a, b = z['profile']
            separator_optimum(z['homotopy'], Q**4, 2**a * 20**b, 2 * a + 20 * b, 2628)
            z = row['complete_square']
            separator_optimum(z['homotopy'], Q**2, 2**z['K'], 2 * z['K'], 692)
            # Conservative projective right-kernel union bound.
            lines = Fraction(A4**s - 1, A4 - 1)
            fail = lines * RHO**(4 * s)
            assert abs(float(fail) - row['fast']['jacobian_failure_upper']) < 1e-16


def adaptive_check(rep):
    for row in rep['rows'].values():
        if 'complete' in row:
            f, F, ad = row['diagonal'], row['complete'], row['optimized']
            if row['level'] != 'I':
                assert abs(ad['normalized_log2_AXN'] - F['total_log2_AXN']) < 1e-12
                continue
            eps = min(.5, 2**f['spectrum_failure_log2'] + 2**F['spectrum_failure_log2'])
            q = min(1., f['jacobian_failure_upper'] / (1 - eps))
            W = 2**f['per_good_key_log2_AXN'] + q * max(
                0., 2**F['per_good_key_log2_AXN'] - 2**f['per_good_key_log2_AXN'])
            solve = math.log2(W / (1 - eps))
            t = f['target_filter_log2_AXN']
            tot = max(solve, t) + math.log2(1 + 2**(min(solve, t) - max(solve, t)))
            assert abs(tot - ad['normalized_log2_AXN']) < 2e-9
        else:
            f = row['fast']
            for fallback_key, adaptive_key in [('complete_square', 'simple_adaptive'), ('orbit_fallback', 'optimized')]:
                F, ad = row[fallback_key], row[adaptive_key]
                eps = min(.5, 2**f['spectrum_failure_log2'] + 2**F['spectrum_failure_log2'])
                q = min(1., f['jacobian_failure_upper'] / (1 - eps))
                W = 2**f['per_good_key_log2_AXN'] + q * max(
                    0., 2**F['per_good_key_log2_AXN'] - 2**f['per_good_key_log2_AXN'])
                norm = W / (f['structural_probability'] * (1 - eps))
                assert abs(math.log2(norm) - ad['normalized_log2_AXN']) < 3e-9


def headline_checks(rep):
    for title, exp in [
        ('all_nine_simple_complete_square', EXPECTED_SIMPLE),
        ('all_nine_optimized', EXPECTED_OPT),
        ('all_nine_compact_output', EXPECTED_COMPACT),
    ]:
        for level, x in exp.items():
            z = rep[title][level]
            assert abs(z['exponent'] - x) < 1e-10
            assert abs(z['headroom'] - (REF[level] - x)) < 1e-10
    assert all(EXPECTED_SIMPLE[L] < REF[L] for L in REF)


def orbit_count_check():
    fam = [(a, b) for a in range(4) for b in range(a, 4)]
    def pair(a, b):
        a %= 4; b %= 4
        return (a, b) if a <= b else (b, a)
    perm = [fam.index(pair(a + 1, b + 1)) for a, b in fam]
    vec = []
    for i in range(10):
        z = [0] * 10; z[i] = 2; vec.append(tuple(z))
    for i in range(10):
        for j in range(i + 1, 10):
            z = [0] * 10; z[i] = z[j] = 1; vec.append(tuple(z))
    unseen = set(vec); n = 0
    while unseen:
        x = next(iter(unseen)); orbit = set(); y = x
        for _ in range(4):
            orbit.add(y); z = [0] * 10
            for i, a in enumerate(y): z[perm[i]] = a
            y = tuple(z)
        unseen -= orbit; n += 1
    assert len(vec) == 55 and n == 16


def repair_floor_checks():
    rows = [
        (28,5,4,4),(48,16,2,2),(28,4,4,5),
        (40,7,4,4),(72,24,2,2),(38,5,4,5),
        (50,9,4,4),(96,32,2,2),(52,6,4,6),
    ]
    combined, zero = [], []
    for v, o, d, r in rows:
        x = o * r; M = d * x; m1 = (x + d - 1) // d
        Kord = min(M, m1 * d * d); C = d * (d + 1) // 2
        xc = max(v, (Kord + d - 1) // d, d * (((Kord + C - 1) // C) - 1) + 1)
        combined.append(100 * (d * xc - M) / M)
        msq = d * (v - 1) // C + 1
        xz = max(v, d * (msq - 1) + 1, (Kord + d - 1) // d,
                 d * (((Kord + C - 1) // C) - 1) + 1)
        zero.append(100 * (d * xz - M) / M)
    assert min(combined) == 45
    assert abs(max(combined) - 60.714285714285715) < 1e-12
    assert min(zero) == 96.875 and max(zero) == 128


def main():
    circuit_checks()
    target_sampling_checks()
    rep = json.loads((HERE / 'primary_ledger.json').read_text())
    spectrum_checks(rep)
    dimension_checks(rep)
    route_checks(rep)
    adaptive_check(rep)
    headline_checks(rep)
    orbit_count_check()
    repair_floor_checks()
    print('Independent final-paper checks passed')
    print('- 150-gate F19 multiplier on all 361 canonical pairs')
    print('- exact official-target rejection-sampling acceptance')
    print('- sharp accepted-root denominator a+eta')
    print('- conservative B^2 finite-cardinality separator optima')
    print('- direct complete-square dimension inequalities for all nine rows')
    print('- conservative projective l=4 Jacobian-preflight bounds')
    print('- adaptive all-nine ledgers, optional 55-to-16 orbit count, and repair floors')


if __name__ == '__main__':
    main()
