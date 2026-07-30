#!/usr/bin/env python3
"""Constructive AXN upper bounds for F19, F19^2, and F19^4 arithmetic.

Every base-field circuit in this file is an explicit network over two-input
AND/XOR/XNOR gates with free constants, wires, fan-out, and structural sharing.
The F19 addition, subtraction, negation, and multiplication circuits are
exhaustively checked on every canonical input.  The extension-field bounds are
transparent Karatsuba compositions.

The dedicated F19 multiplier is substantially smaller than the earlier
multiply-accumulate primitive because symbolic homotopy needs a plain
pointwise multiplication.  It computes the exact nine-bit product p=l+32h,
uses the identity 32=13 (mod 19), evaluates 13h mod 19 by a five-output ANF
valid on h=0,...,10, and performs a final bounded reduction.

The F19^4 tower is connected to the paper's sparse basis
F19[t]/(t^4-t-1) by exact identities checked below.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

Q = 19


@dataclass(frozen=True)
class Wire:
    idx: int


class Circuit:
    def __init__(self, ninputs: int):
        self.ninputs = ninputs
        self.gates: list[tuple[str, tuple[str, int | bool], tuple[str, int | bool]]] = []

    def _gate(self, op: str, a: Wire | bool, b: Wire | bool) -> Wire:
        ai = ("w", a.idx) if isinstance(a, Wire) else ("c", bool(a))
        bi = ("w", b.idx) if isinstance(b, Wire) else ("c", bool(b))
        self.gates.append((op, ai, bi))
        return Wire(self.ninputs + len(self.gates) - 1)

    def xor(self, a: Wire | bool, b: Wire | bool) -> Wire:
        return self._gate("XOR", a, b)

    def xnor(self, a: Wire | bool, b: Wire | bool) -> Wire:
        return self._gate("XNOR", a, b)

    def and_(self, a: Wire | bool, b: Wire | bool) -> Wire:
        return self._gate("AND", a, b)

    def not_(self, a: Wire | bool) -> Wire:
        return self.xnor(a, False)

    def eval(self, bits: Iterable[int | bool]) -> list[bool]:
        vals = [bool(x) for x in bits]
        if len(vals) != self.ninputs:
            raise ValueError("wrong input length")
        for op, a, b in self.gates:
            av = vals[int(a[1])] if a[0] == "w" else bool(a[1])
            bv = vals[int(b[1])] if b[0] == "w" else bool(b[1])
            if op == "XOR":
                vals.append(av ^ bv)
            elif op == "XNOR":
                vals.append(not (av ^ bv))
            elif op == "AND":
                vals.append(av and bv)
            else:
                raise AssertionError(op)
        return vals

    @property
    def count(self) -> int:
        return len(self.gates)

    def digest(self, outputs: list[Wire]) -> str:
        payload = json.dumps(
            {
                "ninputs": self.ninputs,
                "gates": self.gates,
                "outputs": [w.idx for w in outputs],
            },
            separators=(",", ":"),
            sort_keys=True,
        ).encode()
        return hashlib.sha256(payload).hexdigest()


def bits5(x: int) -> list[int]:
    return [(x >> i) & 1 for i in range(5)]


def value(bits: list[bool], wires: list[Wire]) -> int:
    return sum((1 << i) for i, w in enumerate(wires) if bits[w.idx])


def half_adder(c: Circuit, a: Wire, b: Wire) -> tuple[Wire, Wire]:
    return c.xor(a, b), c.and_(a, b)


def full_adder(c: Circuit, a: Wire, b: Wire, carry: Wire | bool) -> tuple[Wire, Wire]:
    t = c.xor(a, b)
    s = c.xor(t, carry)
    u = c.and_(a, b)
    v = c.and_(carry, t)
    cout = c.xor(u, v)  # the two terms are disjoint
    return s, cout


def add5_to_6(c: Circuit, a: list[Wire], b: list[Wire]) -> list[Wire]:
    out = [c.xor(a[0], b[0])]
    carry = c.and_(a[0], b[0])
    for i in range(1, 5):
        s, carry = full_adder(c, a[i], b[i], carry)
        out.append(s)
    out.append(carry)
    return out


def or2(c: Circuit, a: Wire, b: Wire) -> Wire:
    # a OR b = (a XOR b) XOR (a AND b).
    t = c.xor(a, b)
    u = c.and_(a, b)
    return c.xor(t, u)


def add_conditional_13_mod32(c: Circuit, x: list[Wire], q: Wire) -> list[Wire]:
    """Return x + 13*q modulo 32."""
    y: list[Wire] = []
    y.append(c.xor(x[0], q))
    carry = c.and_(x[0], q)
    y.append(c.xor(x[1], carry))
    carry = c.and_(x[1], carry)
    for i in (2, 3):
        s, carry = full_adder(c, x[i], q, carry)
        y.append(s)
    y.append(c.xor(x[4], carry))
    return y


def build_add19() -> tuple[Circuit, list[Wire]]:
    c = Circuit(10)
    inputs = [Wire(i) for i in range(10)]
    a, b = inputs[:5], inputs[5:]
    z = add5_to_6(c, a, b)

    # For valid inputs z<=36.  Hence z_5 and z_4 cannot both be one, and
    # z>=19 iff z_5 XOR [z_4 AND (z_3 OR z_2 OR (z_1 AND z_0))].
    z10 = c.and_(z[1], z[0])
    z32 = or2(c, z[3], z[2])
    low_ge_3 = or2(c, z32, z10)
    middle = c.and_(z[4], low_ge_3)
    reduce_bit = c.xor(z[5], middle)
    return c, add_conditional_13_mod32(c, z[:5], reduce_bit)


def build_neg19() -> tuple[Circuit, list[Wire]]:
    """Return -x modulo 19 for canonical x in [0,18]."""
    c = Circuit(5)
    x = [Wire(i) for i in range(5)]

    # Compute 19-x modulo 32 as (~x)+20, then repair x=0.
    nx0 = c.not_(x[0])
    nx1 = c.not_(x[1])
    nx2 = c.not_(x[2])
    nx3 = c.not_(x[3])
    nx4 = c.not_(x[4])
    y0 = nx0
    y1 = nx1
    y2 = x[2]
    carry3 = nx2
    y3 = c.xor(nx3, carry3)
    carry4 = c.and_(nx3, carry3)
    t4 = c.xor(nx4, carry4)
    y4 = c.not_(t4)
    y = [y0, y1, y2, y3, y4]

    z01 = c.and_(nx0, nx1)
    z23 = c.and_(nx2, nx3)
    z0123 = c.and_(z01, z23)
    is_zero = c.and_(z0123, nx4)
    y[0] = c.xor(y[0], is_zero)
    y[1] = c.xor(y[1], is_zero)
    y[4] = c.xor(y[4], is_zero)
    return c, y


def _sum_columns_truncated(c: Circuit, columns: list[list[Wire]], width: int) -> list[Wire]:
    """Sum weighted bit columns modulo 2^width, dropping the final carry."""
    cols = [list(col) for col in columns[:width]]
    for i in range(width - 1):
        while len(cols[i]) >= 3:
            a, b, d = cols[i].pop(), cols[i].pop(), cols[i].pop()
            s, carry = full_adder(c, a, b, d)
            cols[i].append(s)
            cols[i + 1].append(carry)

    # Only parity matters in the final column because the 2^width carry is
    # deliberately discarded.
    while len(cols[width - 1]) > 1:
        a, b = cols[width - 1].pop(), cols[width - 1].pop()
        cols[width - 1].append(c.xor(a, b))

    out: list[Wire] = []
    carry: Wire | None = None
    for i in range(width):
        xs = cols[i]
        final = i == width - 1
        if carry is None:
            if not xs:
                raise AssertionError("unexpected empty product column")
            if len(xs) == 1:
                out.append(xs[0])
            elif final:
                out.append(c.xor(xs[0], xs[1]))
            else:
                s, carry = half_adder(c, xs[0], xs[1])
                out.append(s)
        else:
            if not xs:
                out.append(carry)
                carry = None
            elif len(xs) == 1:
                if final:
                    out.append(c.xor(xs[0], carry))
                    carry = None
                else:
                    s, carry = half_adder(c, xs[0], carry)
                    out.append(s)
            elif final:
                t = c.xor(xs[0], xs[1])
                out.append(c.xor(t, carry))
                carry = None
            else:
                s, carry = full_adder(c, xs[0], xs[1], carry)
                out.append(s)
    return out


def _reduce_0_31(c: Circuit, x: list[Wire]) -> list[Wire]:
    """Return x mod 19 for a five-bit x in [0,31]."""
    x10 = c.and_(x[1], x[0])
    x32 = or2(c, x[3], x[2])
    low_ge_3 = or2(c, x32, x10)
    reduce_bit = c.and_(x[4], low_ge_3)
    return add_conditional_13_mod32(c, x, reduce_bit)


def _map_13h(c: Circuit, h: list[Wire]) -> list[Wire]:
    """Return 13*h mod 19 for the only reachable values h=0,...,10.

    The following algebraic-normal-form circuit uses five shared quadratic
    monomials.  Values on h=11,...,15 are don't-cares.
    """
    h0, h1, h2, h3 = h
    m01 = c.and_(h0, h1)
    m02 = c.and_(h0, h2)
    m12 = c.and_(h1, h2)
    m03 = c.and_(h0, h3)
    m13 = c.and_(h1, h3)

    def xor_many(xs: list[Wire]) -> Wire:
        y = xs[0]
        for z in xs[1:]:
            y = c.xor(y, z)
        return y

    return [
        xor_many([h0, h1, m01, m02, m12, h3, m03]),
        xor_many([h1, m01, h2, m02, m12, m03, m13]),
        xor_many([h0, h1, h2, m03, m13]),
        xor_many([h0, m01, h2, m02, m12, h3, m13]),
        m13,
    ]


def build_mul19() -> tuple[Circuit, list[Wire], dict[str, int]]:
    """Return a 194-gate circuit for canonical F19 multiplication."""
    c = Circuit(10)
    inputs = [Wire(i) for i in range(10)]
    a, b = inputs[:5], inputs[5:]

    columns: list[list[Wire]] = [[] for _ in range(9)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            columns[i + j].append(c.and_(ai, bj))
    product = _sum_columns_truncated(c, columns, 9)
    after_product = c.count

    # For canonical inputs p<=18^2=324<2^9, so the truncated nine-bit product
    # is exact.  Write p=l+32h, with h in {0,...,10}.
    low = product[:5]
    high_correction = _map_13h(c, product[5:9])
    after_map = c.count

    # low + correction <= 31+16=47.  Fold its possible 32 once using
    # 32=13 mod 19, leaving [0,31], and reduce once more.
    z = add5_to_6(c, low, high_correction)
    folded = add_conditional_13_mod32(c, z[:5], z[5])
    out = _reduce_0_31(c, folded)
    stages = {
        "nine_bit_product": after_product,
        "after_13h_map": after_map,
        "complete_multiplier": c.count,
    }
    return c, out, stages


# Polynomial arithmetic in F19[t]/(t^4-t-1).
def poly_add(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((x + y) % Q for x, y in zip(a, b))


def poly_mul(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    raw = [0] * 7
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            raw[i + j] = (raw[i + j] + x * y) % Q
    for k in range(6, 3, -1):
        x = raw[k] % Q
        raw[k] = 0
        raw[k - 3] = (raw[k - 3] + x) % Q
        raw[k - 4] = (raw[k - 4] + x) % Q
    return tuple(x % Q for x in raw[:4])


def det_mod(matrix: list[list[int]]) -> int:
    a = [row[:] for row in matrix]
    n = len(a)
    det = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col] % Q), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        p = a[col][col] % Q
        det = det * p % Q
        inv = pow(p, Q - 2, Q)
        for r in range(col + 1, n):
            f = a[r][col] * inv % Q
            for j in range(col, n):
                a[r][j] = (a[r][j] - f * a[col][j]) % Q
    return det % Q


def main() -> None:
    add, add_out = build_add19()
    neg, neg_out = build_neg19()
    mul, mul_out, mul_stages = build_mul19()

    # Compose subtraction as a + (-b); this remains a transparent upper bound.
    sub19 = add.count + neg.count
    for a in range(Q):
        for b in range(Q):
            add_value = value(add.eval(bits5(a) + bits5(b)), add_out)
            mul_value = value(mul.eval(bits5(a) + bits5(b)), mul_out)
            neg_b = value(neg.eval(bits5(b)), neg_out)
            sub_value = value(add.eval(bits5(a) + bits5(neg_b)), add_out)
            assert add_value == (a + b) % Q
            assert mul_value == (a * b) % Q
            assert sub_value == (a - b) % Q
    for x in range(Q):
        assert value(neg.eval(bits5(x)), neg_out) == (-x) % Q

    add19 = add.count
    neg19 = neg.count
    mul19 = mul.count

    add2 = 2 * add19
    sub2 = 2 * sub19
    # Complex Karatsuba over u^2=-1.  The two input sums cost 2 additions;
    # reconstruction uses two shared negations and three additions.
    neg2_outputs = 2 * neg19
    add2_outputs = 3 * add19
    mul2 = 2 * add19 + 3 * mul19 + neg2_outputs + add2_outputs
    mac2_plus = mul2 + add2
    mac2_minus = mul2 + sub2

    # Exact tower identities in the sparse t-basis used by the paper.
    one = (1, 0, 0, 0)
    minus_one = (Q - 1, 0, 0, 0)
    u = (1, 2, 15, 5)
    delta = poly_add(one, u)
    v = (4, 13, 13, 1)
    uv = poly_mul(u, v)
    assert poly_mul(u, u) == minus_one
    assert poly_mul(v, v) == delta
    assert poly_mul(uv, uv) == tuple((-x) % Q for x in delta)
    basis = [one, u, v, uv]
    basis_matrix = [[basis[col][row] for col in range(4)] for row in range(4)]
    basis_det = det_mod(basis_matrix)
    assert basis_det != 0

    # Norm(1+u)=(1+u)(1-u)=2 is nonsquare in F19, so v^2-(1+u)
    # is irreducible over F19^2.
    norm_delta = 2
    assert pow(norm_delta, (Q - 1) // 2, Q) == Q - 1

    delta_mul2 = sub19 + add19
    add4 = 2 * add2
    sub4 = 2 * sub2
    mul4 = 2 * add2 + 3 * mul2 + delta_mul2 + add2 + 2 * sub2
    mac4_plus = mul4 + add4
    mac4_minus = mul4 + sub4

    report = {
        "model": "unit-cost two-input AND/XOR/XNOR; free wires/constants/fan-out",
        "scalar": {
            "add_mod_19": add19,
            "neg_mod_19": neg19,
            "sub_mod_19": sub19,
            "multiplication": mul19,
            "multiplication_stages": mul_stages,
            "multiplication_netlist_sha256": mul.digest(mul_out),
            "validated_add_pairs": Q * Q,
            "validated_mul_pairs": Q * Q,
            "validated_sub_pairs": Q * Q,
            "validated_negations": Q,
        },
        "F19_2": {
            "representation": "F19[u]/(u^2+1)",
            "addition": add2,
            "subtraction": sub2,
            "multiplication": mul2,
            "multiplication_output_negations": neg2_outputs,
            "multiplication_output_additions": add2_outputs,
            "mac_plus": mac2_plus,
            "mac_minus": mac2_minus,
            "conservative_signed_mac": max(mac2_plus, mac2_minus),
        },
        "F19_4": {
            "paper_basis": "F19[t]/(t^4-t-1)",
            "tower": "F19^2[v]/(v^2-(1+u))",
            "u_t_coordinates": list(u),
            "v_t_coordinates": list(v),
            "uv_t_coordinates": list(uv),
            "tower_basis_determinant_mod_19": basis_det,
            "norm_1_plus_u": norm_delta,
            "delta_multiplication": delta_mul2,
            "addition": add4,
            "subtraction": sub4,
            "multiplication": mul4,
            "mac_plus": mac4_plus,
            "mac_minus": mac4_minus,
            "conservative_signed_mac": max(mac4_plus, mac4_minus),
            "previous_bound": 6081,
        },
    }
    netlist = {
        "model": report["model"],
        "field": "F19",
        "canonical_input_encoding": "two little-endian five-bit integers in 0,...,18",
        "canonical_output_encoding": "one little-endian five-bit integer in 0,...,18",
        "ninputs": mul.ninputs,
        "gates": mul.gates,
        "outputs": [w.idx for w in mul_out],
        "gate_count": mul.count,
        "sha256": mul.digest(mul_out),
    }
    netout = Path(__file__).with_name("snova_f19_mul194_netlist.json")
    netout.write_text(json.dumps(netlist, separators=(",", ":"), sort_keys=True) + "\n")
    report["scalar"]["multiplication_netlist_file"] = netout.name
    out = Path(__file__).with_name("snova_field_tower_circuits_v5.json")
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    print(f"Wrote {out}")
    print(f"Wrote {netout}")


if __name__ == "__main__":
    main()
