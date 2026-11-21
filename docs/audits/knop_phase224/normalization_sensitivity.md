# Phase 2.24 normalization sensitivity (ray 82P, deck r0000, production
# trajectory committed states passes 13-16 via replay maxp runs)

Population: 4*56377 = 225508 cells. Committed delta x(n+1)-x(n), pass 15->16.

normalization         max       med       p90       p99     >2e-3
prod |d|/max(vo,1e-4vf) 3.55E-2  3.04E-5  3.11E-4  4.31E-3   5593
|d|/max(|vn|,1e-4vf)    3.68E-2  3.04E-5  3.11E-4  4.31E-3   5610
|d|/vmaxf(col peak)     3.60E-2  2.62E-5  3.06E-4  4.21E-3   5232
absolute |d|            1.87E-8  2.03E-15 1.25E-12 3.39E-11  0

Bad-cell row census (prod norm): j1=2030 j2=2764 j3=798 j4=0.
Denominator-floored cells (VOLD<=1e-4*vmaxf) among bad&sig: 0.
The denominators are the cells' own values (~1e-20..1e-26), not the
floor - the defect is relative-to-self on absolutely dark cells.
|d|>1e-9: 4 cells; |d|>1e-12: 5016 cells.

Interpretation: relative-normalization variants all report ~5.5K cells
over tolerance; every cell is dark (x << any flux-relevant scale).
The production test cannot separate a 3.5% wobble on a 1e-24 field
from a real transfer error.
