#!/usr/bin/env python3
"""Scope the solving-degree scaling ladder (#1).

Builds FAITHFUL small l=2 SNOVA residual cores via the genuine symmetry-quotient
+ affine-column reduction (no KAT: a fresh key from a random seed, verified
self-consistently by the reference verifier net), emits the determined K-in-K
core, and prints the predicted semi-regular operating degree so we can compare
against an actual msolve solve.

This reuses emit_core_general.build_residual verbatim by synthesizing a
self-consistent KAT (sk = random seed, pk = the key that seed generates), so the
byte-for-byte check passes trivially and every downstream verified code path runs
unchanged.
"""
import sys, os, argparse, tempfile
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "files"))
sys.path.insert(0, str(REPO / "files" / "SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19" / "repro"))
sys.path.insert(0, str(HERE))

import symmetry_attack_validation as V
import emit_level1_core as EM
import emit_core_general as G
sys.path.insert(0, str(REPO / "files"))
from semireg_dreg import d_reg


def synth_kat(p, seed48: bytes, tmpdir: Path) -> Path:
    pk, _ = G.kat_public_key_general(seed48, p)
    path = tmpdir / f"synth_{p.name}.kat"
    path.write_text(f"sk = {seed48.hex()}\npk = {pk.hex()}\n")
    return path


def build_rung(v, o, l, r, seed=20260721, outdir=None):
    p = V.Params(f"ladder-{v}-{o}-{l}-{r}", v, o, l, r, 0)
    K = p.unordered
    rng = np.random.default_rng(seed)
    with tempfile.TemporaryDirectory() as td:
        seed48 = bytes(rng.integers(0, 256, size=48, dtype=np.uint8).tolist())
        kat = synth_kat(p, seed48, Path(td))
        data = G.build_residual(p, kat, seed)
    ok_v = EM.verify_polys(data)
    ok_p = EM.planted_check(data)
    nres = data["residual_vars"]
    core_n = min(K, nres)  # square-specialize down to determined (K-in-K) if possible
    newpolys, keep, zc = EM.specialize(data, core_n, seed=seed)
    dreg = d_reg(K, core_n)
    dens = G.density(newpolys, core_n)
    print(f"[rung] (v,o,l,r)=({v},{o},{l},{r})  m1={p.m1}  K={K}  N_res={nres}  "
          f"-> core {K}/{core_n}  verify={ok_v} planted={ok_p}  density={dens:.4f}  predicted D_reg={dreg}")
    if outdir is not None:
        outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
        tag = f"ladder_l{l}_m1{p.m1}"
        msp = outdir / f"{tag}_core_{K}in{core_n}.ms"
        EM.to_msolve(newpolys, core_n, msp)
        print(f"        emitted {msp}")
    return dict(K=K, core_n=core_n, nres=nres, dreg=dreg, verify=ok_v, planted=ok_p, m1=p.m1, l=l, density=dens)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rungs", default="3,1,2,2;4,2,2,2",
                    help="semicolon-separated v,o,l,r tuples")
    ap.add_argument("--emit", action="store_true", help="write .ms cores")
    ap.add_argument("--outdir", default=str(HERE / "systems" / "ladder"))
    a = ap.parse_args()
    for spec in a.rungs.split(";"):
        v, o, l, r = (int(x) for x in spec.split(","))
        build_rung(v, o, l, r, outdir=a.outdir if a.emit else None)
