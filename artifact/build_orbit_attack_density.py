#!/usr/bin/env python3
"""Combine structural preflight, seeded-spectrum, and orbit-complete homotopy.

The resulting table gives both per-certified-key and random-key-normalized
operation exponents for the fixed-core-free Frobenius-orbit sweep.
"""
from __future__ import annotations
import json, math, sys
from fractions import Fraction
from pathlib import Path
sys.set_int_max_str_digits(1000000)
ROOT=Path(__file__).resolve().parent

def load(n):
 return json.loads((ROOT/n).read_text())
def log2f(x): return math.log2(x.numerator)-math.log2(x.denominator)

pre=load('structural_preflight_probability.json')
orb=load('frobenius_orbit_sweep_certificate.json')
levels={x['level']:x for x in orb['levels']}
rows=[]
for sh in pre['results']:
 z=levels[sh['level']]
 cross=Fraction(sh['complete_cross_preflight_success']['numerator'],sh['complete_cross_preflight_success']['denominator'])
 # Order in the orbit certificate is the order of the two official shapes.
 shape_index=0 if sh['key'].endswith('-a') else 1
 tail=z['spectrum_tails'][shape_index]
 modes={}
 for mode in ('compact','robust'):
  eps=Fraction(tail[mode]['numerator'],tail[mode]['denominator'])
  density=cross*(1-eps)
  bits=float(z['costs'][mode]['bits'])
  penalty=-log2f(density)
  normalized=bits+penalty
  modes[mode]={
   'eta_numerator':z['costs'][mode]['eta_numerator'],
   'eta_denominator':z['costs'][mode]['eta_denominator'],
   'spectrum_failure':{'numerator':eps.numerator,'denominator':eps.denominator,'decimal':float(eps),'log2':log2f(eps)},
   'vulnerable_key_density_lower':{'numerator':density.numerator,'denominator':density.denominator,'decimal':float(density),'log2':log2f(density),'penalty_bits':penalty},
   'per_certified_key_bits':bits,
   'random_key_normalized_bits':normalized,
   'reference_bits':z['reference'],
   'normalized_margin_bits':z['reference']-normalized,
  }
 rows.append({'shape':sh['key'],'level':sh['level'],'h':z['h'],'K':z['K'],'patterns':z['all_pattern_count'],'orbits':z['frobenius_orbit_count'],'peak_output':z['peak_individual_output_human'],'modes':modes})

data={'model':'random-XOF idealization with coordinatewise modulo-19 bytes','algorithm':'fixed-core-free Frobenius-orbit-complete homotopy sweep','results':rows}
(ROOT/'orbit_attack_density.json').write_text(json.dumps(data,indent=2)+'\n')
md=['# Fixed-core-free orbit-attack density certificate','',
'The structured source and projection preflights succeed with the exact lower bound in `STRUCTURAL_PREFLIGHT_PROBABILITY_CERTIFICATE.md`.  The orbit sweep needs no fixed core-Jacobian event: at every descent root where the full restricted Jacobian has rank `h`, one of the 16 Frobenius-orbit representatives contains a nonsingular square core.  The only internal-coefficient failure event is therefore the seeded aggregate-spectrum tail.','',
'| Shape | `(h,K)` | 55 -> 16 | compact per key | compact normalized | margin | robust per key | robust normalized | margin | peak one-core output |',
'|:--:|:--:|:--:|--:|--:|--:|--:|--:|--:|--:|']
for r in rows:
 c=r['modes']['compact'];b=r['modes']['robust']
 md.append(f"| {r['shape']} | ({r['h']},{r['K']}) | {r['patterns']} -> {r['orbits']} | {c['per_certified_key_bits']:.3f} | {c['random_key_normalized_bits']:.3f} | {c['normalized_margin_bits']:.3f} | {b['per_certified_key_bits']:.3f} | {b['random_key_normalized_bits']:.3f} | {b['normalized_margin_bits']:.3f} | {r['peak_output']} |")
md += ['', 'The normalized exponent pays the inverse certified key density and is therefore a time-success-ratio ledger over a random generated public key.  The dense parametrization ceilings are not low-memory claims; this route is included to remove fixed-core and selected-Macaulay-degree assumptions, while the eigenblock-core and streamed-XL branches supply lower-output or lower-state points.']
(ROOT/'ORBIT_ATTACK_DENSITY_CERTIFICATE.md').write_text('\n'.join(md)+'\n')
print('wrote orbit_attack_density.json and ORBIT_ATTACK_DENSITY_CERTIFICATE.md')
for r in rows:
 print(r['shape'], 'compact',f"{r['modes']['compact']['random_key_normalized_bits']:.6f}", 'robust',f"{r['modes']['robust']['random_key_normalized_bits']:.6f}")
