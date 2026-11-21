import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog

root = Path.cwd()
folder = root / 'research/pde_ledger_v3/_measurements'
archive = json.loads((folder/'mouth_bend_depth_participation.json').read_text())['runs']
sys.path.insert(0, str(root/'software/stage1_solver/tools'))
from mouth_bend_depth_participation import spectral_readout

def f(x):
    return np.sin(np.pi*x)**2

def pressure(m):
    lines = archive[f'bend_depth_pressure{m}']['stdout_lines']
    slices = [[float(s.split('=')[1]) for s in text.split()[1:]]
              for text in lines if text.startswith('attached ')]
    bounds = [0.]
    for midpoint,share in slices:
        bounds.append(2*midpoint-bounds[-1])
    assert abs(bounds[-1]-1)<1e-7
    bounds[-1]=1.
    def moment(lo,hi,power):
        total=0.
        for (_,share), x,y in zip(slices,bounds,bounds[1:]):
            l,r=max(x,lo),min(y,hi)
            if r>l:
                total+=share/(y-x)*quad(lambda w:f(w)**power,l,r)[0]
        return total
    return moment

def central(lo,hi,power):
    return quad(lambda w:(.02+f(w)**4)*f(w)**power,lo,hi)[0]/(.02+math.comb(8,4)/256)

radial=59.3830741481
def review(label,moment):
    e,c=moment(0,.25,1),moment(.375,.625,1)
    mf=moment(0,1,2)
    outer_max=f(.375)
    mass_c=min(c/outer_max,1-2*e)
    lower=e**2/(1-mass_c)+c**2/mass_c
    upper=.5*e+c+outer_max**2*(1-2*e-c)
    print(f'{label} mean_f_times_mu_E={4*e:.12g} mean_f_times_mu_C={4*c:.12g} '
          f'P_E/bdot={radial*e/6:.12g} P_C/bdot={radial*c/6:.12g} '
          f'C/E={c/e:.12g} M_modal={radial*mf/3:.12g} '
          f'M_lower={radial*lower/3:.12g} M_upper={radial*upper/3:.12g}')
    # Independent optimization over a mesh of available point masses, entirely
    # inside each region. The closure certificate above includes boundary atoms.
    for n in (1024,16384):
        x=(np.arange(n)+.5)/n
        fx=f(x)
        A=np.array([np.ones(n),fx*(x<.25),fx*((x>.375)&(x<.625))])
        low=linprog(fx**2,A_eq=A,b_eq=[1.,e,c],bounds=(0,None),method='highs')
        high=linprog(-fx**2,A_eq=A,b_eq=[1.,e,c],bounds=(0,None),method='highs')
        assert low.success and high.success
        print(f'{label} n={n} LP_lower={radial*low.fun/3:.12g} '
              f'LP_upper={-radial*high.fun/3:.12g}')
        assert abs(low.fun-lower)<3/n and abs(-high.fun-upper)<4/n
    for wall in (12.,14.):
        if low.fun < 3*wall/radial < -high.fun:
            mix=(3*wall/radial-low.fun)/(-high.fun-low.fun)
            weights=(1-mix)*low.x+mix*high.x
            assert max(abs(A@weights-np.array([1,e,c])))<1e-10
            print(f'{label} same_bands_mass={radial*(fx**2@weights)/3:.12g} '
                  f'mixture_fraction={mix:.10g}')

for name in ('030','015'):
    review('pressure'+name,pressure(name))
review('central',central)
for wall in (12.,14.):
    spectral_readout(wall,'witness'+str(int(wall)))
