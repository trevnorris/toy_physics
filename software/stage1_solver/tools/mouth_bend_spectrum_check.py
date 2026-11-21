"""Solve the three-by-three matrix on a frequency grid with a linear table interpolant."""
import json
from pathlib import Path
import re
import sys

import numpy as np
from numpy.polynomial import Polynomial

sys.path.insert(0, str(Path('software/stage1_solver/tools').resolve()))
from mouth_carrier_bend_budget import settings

point,energy,bend,_,arc = settings()
ws,wd = (float(point[k]) for k in ('omega_s','omega_d'))
g = float(bend['overlap'])*np.sqrt(2*arc['parameters']['I_wave']/ws)
gap = abs(float(energy['U_over_action']))
base = Path('_scratch/s11c')
def parse(rows,start):
    return [dict(re.findall(r'(\w+)=([^ ]+)',s)) for s in rows if s.startswith(start)]
def load(name):
    return (base/name/'stdout').read_text().splitlines()
geo, = parse(load('bend_projection_wave015'),'geometry ')
projection, = parse(load('bend_projection_wave015'),'support ')
raw = json.loads(Path('research/pde_ledger_v3/_measurements/mouth_open_reservoirs.json').read_text())
rows = raw['runs']['final_open_lines64']['stdout_lines']+load('bend_projection_rad64')
direct = float(projection['direct_K'])
k0 = direct+float(geo['tensile_both'])
x = np.linspace(-.016,.08,12001)
z = x+wd
def interp(rows,field):
    pairs = sorted((float(r['offset']),float(r[field])) for r in rows)
    return np.interp(x,*np.array(pairs).T)
forces = []
for tone in ('difference','sum'):
    samples = parse(rows,'tone='+tone+' ')
    for s in samples:
        s['imag'] = 2*float(s['power_per_b2'])/float(s['omega'])
    forces.append(interp(samples,'induced_real')+1j*interp(samples,'imag'))
wave = parse(load('bend_projection_wresp03'),'wave ')
common, lower, upper, emitted = (interp(wave,k) for k in ('central','lower','upper_real','upper_imag'))
for key in ['uniform_wall_mass','central_wall_mass']:
    mass = float(geo[key])
    kernel = np.zeros((len(x),3,3),complex)
    kernel[:,0,0] = wd**2-z**2
    kernel[:,0,1:] = g/2
    kernel[:,1:,0] = g/2
    kernel[:,1,1] = k0-mass*(z-ws)**2-forces[0]-common-lower
    kernel[:,2,2] = k0-mass*(z+ws)**2-forces[1]-common-upper-1j*emitted
    kernel[:,1,2] = kernel[:,2,1] = direct/2-common
    profile = np.linalg.inv(kernel)[:,0,0].imag
    maxima = np.flatnonzero((profile[1:-1]>profile[:-2])&(profile[1:-1]>profile[2:]))+1
    idx = maxima[np.argmax(profile[maxima])]
    half = profile[idx]/2
    left = np.interp(half,profile[:idx+1],x[:idx+1])
    right = np.interp(half,profile[:idx-1:-1],x[:idx-1:-1])
    # derivative numerator of the inherited equal-area Lorentzian readout, x in gap units
    a = (right-left)/(2*gap)
    u = Polynomial([0,1])
    derivative = u*((u+1)**2+9*a*a)**2 + 3*(u+1)*(u**2+a*a)**2
    roots = derivative.roots()
    count = sum(abs(q.imag)<1e-9 and derivative.deriv()(q.real)>0 for q in roots)
    print(f'sheet={key} raw_band_maxima={len(maxima)} delta={x[idx]:.10g} width={right-left:.12g} exact_readout_maxima={count} left_half={left:.12g} right_half={right:.12g}')
