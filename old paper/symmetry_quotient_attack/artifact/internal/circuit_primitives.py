#!/usr/bin/env python3
"""Shared primitives for the committed F19 AXN circuit generator.

The current generator imports the circuit container, bit adders, canonical
negation circuit, and polynomial-basis helpers from this module.  This module
does not emit artifacts on its own.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
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
