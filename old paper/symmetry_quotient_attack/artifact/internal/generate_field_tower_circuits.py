#!/usr/bin/env python3
"""Generate the committed AXN circuit artifacts for F19, F19^2, and F19^4.

The F19 multiplier first maps each canonical residue to balanced signed
magnitude in {0,...,9}, multiplies two four-bit magnitudes exactly, reduces the
seven-bit product modulo 19, then conditionally restores the sign.  All
canonical input pairs are exhaustively verified.
"""
from __future__ import annotations

import json
from pathlib import Path

import circuit_primitives as core

HERE = Path(__file__).resolve().parent
ARTIFACT_DIR = HERE.parent
Circuit = core.Circuit
Wire = core.Wire
Q = core.Q
bits5 = core.bits5
value = core.value



def build_add19() -> tuple[Circuit, list[Wire]]:
    """A 42-gate canonical F19 adder with shared compare/carry logic."""
    c=Circuit(10)
    inp=[Wire(i) for i in range(10)]
    z=core.add5_to_6(c,inp[:5],inp[5:])
    # z<=36.  Let low indicate that the low four bits encode at least 3.
    # The same nodes also give the carries when 13 is added conditionally.
    z10=c.and_(z[1],z[0])
    a2=core.or2(c,z[2],z10)
    low=core.or2(c,z[3],a2)
    middle=c.and_(z[4],low)
    reduce_bit=c.xor(z[5],middle)
    c1=c.and_(reduce_bit,z[0])
    y0=c.xor(z[0],reduce_bit)
    y1=c.xor(z[1],c1)
    c2=c.and_(reduce_bit,z10)
    t2=c.xor(z[2],reduce_bit)
    y2=c.xor(t2,c2)
    c3=c.and_(reduce_bit,a2)
    t3=c.xor(z[3],reduce_bit)
    y3=c.xor(t3,c3)
    c4=c.and_(reduce_bit,low)
    y4=c.xor(z[4],c4)
    return c,[y0,y1,y2,y3,y4]

def balanced_decode(c: Circuit, x: list[Wire]) -> tuple[Wire,list[Wire]]:
    """For x in 0,...,18 return sign and magnitude min(x,19-x).

    The identities exploit don't-care inputs 19,...,31 and use 9 gates.
    """
    x0,x1,x2,x3,x4=x
    a=c.and_(x1,x3)
    b=c.and_(x2,x3)
    ab2=c.and_(a,x2)
    s=c.xor(a,b)
    s=c.xor(s,ab2)
    s=c.xor(s,x4)
    m0=c.xor(x0,s)
    m1=c.xor(x1,s)
    m2=x2
    m3=c.xor(x3,b)
    return s,[m0,m1,m2,m3]


def mul4_exact(c: Circuit,a:list[Wire],b:list[Wire])->list[Wire]:
    columns=[[] for _ in range(7)]
    for i,ai in enumerate(a):
        for j,bj in enumerate(b):
            columns[i+j].append(c.and_(ai,bj))
    return core._sum_columns_truncated(c,columns,7)



def add5_special_last_zero(c: Circuit, a: list[Wire], b: list[Wire | bool]) -> list[Wire]:
    """Add two five-bit words when b[4]=0, returning six bits."""
    if b[4] is not False:
        raise ValueError("special adder requires a zero top bit")
    out=[c.xor(a[0],b[0])]
    carry=c.and_(a[0],b[0])
    for i in range(1,4):
        si,carry=core.full_adder(c,a[i],b[i],carry)
        out.append(si)
    out.append(c.xor(a[4],carry))
    out.append(c.and_(a[4],carry))
    return out

def map_32_high(c: Circuit, high:list[Wire])->list[Wire|bool]:
    """For high h in {0,1,2}, return 13*h mod19 as five little-endian bits.

    Since 32 == 13 mod19, the values are 0,13,7.  State 3 is unreachable.
    """
    h0,h1=high
    u=c.xor(h0,h1)
    return [u,h1,u,h0,False]



def reduce_folded_product(c: Circuit, x: list[Wire]) -> list[Wire]:
    """Reduce the folded balanced-product values modulo 19 in 18 gates.

    The only reachable inputs are
      {0,...,10,12,...,18,20,21,23,...,30}.
    Therefore the reduction predicate is x4 AND (x3 OR x2); the omitted
    values 19,22,31 are don't-cares fixed by the preceding product circuit.
    Shared carry logic then adds 13 modulo 32 when reduction is required.
    """
    x0,x1,x2,x3,x4=x
    high=core.or2(c,x3,x2)
    t=c.and_(x4,high)
    x10=c.and_(x1,x0)
    y0=c.xor(x0,t)
    c1=c.and_(x0,t)
    y1=c.xor(x1,c1)
    c2=c.and_(t,x10)
    q2=c.xor(x2,t)
    y2=c.xor(q2,c2)
    a2=core.or2(c,x2,x10)
    c3=c.and_(t,a2)
    q3=c.xor(x3,t)
    y3=c.xor(q3,c3)
    y4=c.xor(x4,t)
    return [y0,y1,y2,y3,y4]


def reduce_direct_z(c: Circuit, z: list[Wire]) -> list[Wire]:
    """Reduce the exact six-bit pre-reduction word directly modulo 19.

    The preceding balanced-product circuit reaches only
      {0,...,10,12,...,18,20,21,23,...,30,35,37,44}.
    First reduce the low five bits with the 18-gate folded reducer.  When the
    sixth bit is set, the only low words are 3,5,12 and the desired residues
    are 16,18,6.  A ten-gate override implements precisely those three cases.
    All other six-bit inputs are don't-cares.
    """
    if len(z) != 6:
        raise ValueError("direct reducer expects six bits")
    r=reduce_folded_product(c,z[:5])
    h=z[5]
    x0,x1,x2,x3,_x4=z[:5]
    nx3=c.not_(x3)
    y0=c.xor(r[0],c.and_(h,x0))
    y1=c.xor(r[1],h)
    t=c.and_(h,x2)
    flip2=c.and_(t,nx3)
    y2=c.xor(r[2],flip2)
    y3=c.xor(r[3],c.and_(h,x3))
    # Multiplex r_4 (when h=0) with not x_3 (when h=1).
    d=c.xor(r[4],nx3)
    y4=c.xor(r[4],c.and_(h,d))
    return [y0,y1,y2,y3,y4]

def build_mul19_balanced()->tuple[Circuit,list[Wire],dict[str,int]]:
    c=Circuit(10)
    inp=[Wire(i) for i in range(10)]
    sa,a=balanced_decode(c,inp[:5]); sb,b=balanced_decode(c,inp[5:])
    after_decode=c.count
    sign=c.xor(sa,sb)
    product=mul4_exact(c,a,b)
    after_product=c.count

    # p = low + 32*high, p<=81 and high in {0,1,2}.  Replace 32*high by
    # (13*high mod19), then reduce the exact reachable six-bit word directly.
    corr=map_32_high(c,product[5:7])
    z=add5_special_last_zero(c,product[:5],corr)
    residue=reduce_direct_z(c,z)
    after_reduce=c.count

    # For canonical r in 0,...,18, the bitwise difference r XOR (-r mod 19)
    # has a particularly small form:
    #   bits 0,1: nonzero(r); bit 2: 0; bit 3: r_2;
    #   bit 4: 1 exactly on {1,2,3,16,17,18}.
    # The last predicate equals (r4 OR r1 OR r0) AND NOT(r3 OR r2).
    r0,r1,r2,r3,r4=residue
    low01=core.or2(c,r1,r0)
    low=core.or2(c,r4,low01)
    high=core.or2(c,r3,r2)
    high_zero=c.not_(high)
    d4=c.and_(low,high_zero)
    nonzero=c.xor(d4,high)  # disjoint on canonical residues
    out=[]
    for ri,di in ((r0,nonzero),(r1,nonzero)):
        out.append(c.xor(ri,c.and_(sign,di)))
    out.append(r2)
    out.append(c.xor(r3,c.and_(sign,r2)))
    out.append(c.xor(r4,c.and_(sign,d4)))
    after_neg=c.count
    stages={
        'balanced_decode_both':after_decode,
        'sign_and_exact_magnitude_product':after_product,
        'after_modular_reduction':after_reduce,
        'after_signed_recode':after_neg,
        'complete_multiplier':c.count,
    }
    return c,out,stages


def main()->None:
    add,add_out=build_add19(); neg,neg_out=core.build_neg19(); mul,mul_out,stages=build_mul19_balanced()
    for x in range(Q):
        assert value(neg.eval(bits5(x)),neg_out)==(-x)%Q
        for y in range(Q):
            vals=mul.eval(bits5(x)+bits5(y))
            got=value(vals,mul_out)
            assert got==(x*y)%Q,(x,y,got,(x*y)%Q)
            assert value(add.eval(bits5(x)+bits5(y)),add_out)==(x+y)%Q
    add19=add.count; neg19=neg.count; sub19=add19+neg19; mul19=mul.count
    add2=2*add19; sub2=2*sub19
    mul2=2*add19+3*mul19+2*neg19+3*add19
    add4=2*add2; sub4=2*sub2
    delta_mul2=sub19+add19
    # For the high coefficient, form p0+p1 once and subtract it from p2;
    # this saves one F19^2 subtraction versus two serial subtractions.
    mul4=2*add2+3*mul2+delta_mul2+add2+(add2+sub2)

    # Recheck the exact tower isomorphism used by the recurrence ledger.
    one=(1,0,0,0); minus_one=(Q-1,0,0,0); u=(1,2,15,5); delta=core.poly_add(one,u); vv=(4,13,13,1); uv=core.poly_mul(u,vv)
    assert core.poly_mul(u,u)==minus_one
    assert core.poly_mul(vv,vv)==delta
    assert core.poly_mul(uv,uv)==tuple((-x)%Q for x in delta)
    basis=[one,u,vv,uv]; mat=[[basis[col][row] for col in range(4)] for row in range(4)]
    det=core.det_mod(mat); assert det!=0

    report={
      'model':'unit-cost two-input AND/XOR/XNOR; free wires/constants/fan-out',
      'scalar':{
        'addition':add19,'negation':neg19,'subtraction':sub19,'multiplication':mul19,
        'multiplication_stages':stages,'validated_mul_pairs':Q*Q,
        'multiplication_payload_sha256':mul.digest(mul_out),
      },
      'F19_2':{
        'representation':'F19[u]/(u^2+1)','addition':add2,'subtraction':sub2,'multiplication':mul2,
      },
      'F19_4':{
        'paper_basis':'F19[t]/(t^4-t-1)','tower':'F19^2[v]/(v^2-(1+u))',
        'u_t_coordinates':list(u),'v_t_coordinates':list(vv),'uv_t_coordinates':list(uv),
        'tower_basis_determinant_mod_19':det,'addition':add4,'subtraction':sub4,'multiplication':mul4,
      }
    }
    net={
      'model':report['model'],'field':'F19','canonical_input_encoding':'two little-endian five-bit integers in 0,...,18',
      'canonical_output_encoding':'one little-endian five-bit integer in 0,...,18','ninputs':mul.ninputs,
      'gates':mul.gates,'outputs':[w.idx for w in mul_out],'gate_count':mul.count,'payload_sha256':mul.digest(mul_out),
    }
    (ARTIFACT_DIR/'f19_multiplier_netlist.json').write_text(json.dumps(net,separators=(',',':'),sort_keys=True)+'\n')
    (ARTIFACT_DIR/'field_tower_circuits.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__': main()
