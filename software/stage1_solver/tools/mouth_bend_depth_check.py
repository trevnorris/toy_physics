"""Independent checks of the fixed-wall phase response and sub-branch spectrum."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.interpolate import PchipInterpolator
from scipy.optimize import linprog
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path('software/stage1_solver/tools').resolve()))
import mouth_open_reservoirs as reservoir
from mouth_carrier_bend_budget import settings
from mouth_projected_bend_spectrum import fields, spectrum

point, state, bend, flow, archive = settings()
par = archive['parameters']
rad = archive['phase_A']['q_star'][0]
height = float(point['L'])/2
rec = Path('research/pde_ledger_v3/_measurements')
prior = json.loads((rec/'mouth_open_reservoirs.json').read_text())['runs']
old = json.loads((rec/'mouth_bend_coefficients.json').read_text())['runs']
geo, support = map(fields,old['bend_projection_wave015']['stdout_lines'])
mass = float(geo['breathing'])

# Independent grid linear program over descending depth cell masses.
n = 2048
x = np.arange(1,n+1)/n
weights = np.sin(np.pi*(x-0.5/n))**4
descending = sparse.diags([-np.ones(n-1),np.ones(n-1)],[0,1],shape=(n-1,n))
sol = linprog(-weights,A_ub=descending,b_ub=np.zeros(n-1),
              A_eq=np.ones((1,n)),b_eq=[1],bounds=(0,None),method='highs')
assert sol.success
print(f'LP cells={n} max_f2={-sol.fun:.12g} mass={-mass*sol.fun/3:.12g} '
      f'last_positive_depth={x[sol.x>1e-9][-1]:.12g}',flush=True)

wd, ws = float(point['omega_d']),float(point['omega_s'])
gap = abs(float(state['U_over_action']))
g = float(bend['overlap'])*np.sqrt(2*par['I_wave']/ws)
direct = float(support['direct_K'])
for cells, radiation, kernel, response in [(48,'open_reservoir_lines48','rad48','wresp05'),
                                           (64,'final_open_lines64','rad64','wresp03')]:
    lines = prior[radiation]['stdout_lines'] + old['bend_projection_'+kernel]['stdout_lines']
    rows = [fields(s) for s in old['bend_projection_'+response]['stdout_lines'] if s.startswith('wave ')]
    disp = [PchipInterpolator([float(row['offset']) for row in rows],
                              [float(row[k]) for row in rows],extrapolate=False)
            for k in ('central','lower','upper_real','upper_imag')]
    for m in (0., .01, .1, .3, .5, 1., 2., 3., 4., 5., 6., 7., -mass*sol.fun/3):
        try:
            shift,width,asym,count=spectrum(lines,ws,wd,g,gap,m,
                                float(geo['tensile_both'])+direct,(direct,*disp,1))
            print(f'spectrum cells={cells} M={m:g} shift={shift:.9g} FWHM/gap={width/gap:.9g} count={count}',flush=True)
        except (AssertionError, ValueError) as exc:
            print(f'spectrum cells={cells} M={m:g} error={exc}',flush=True)

# Capture both complex quadratures of the existing physical solve (source over
# two) in a finite real box. Differentiate the phase with respect to wall
# speed directly at real frequency. Check mass with the wall susceptibility.
cells = 32
d = reservoir.ReservoirDomain(rad,height,rad,cells,2,2,3,par)
d.set_current(float(flow['reference_inward_current'])/par['m'])
d.ac = d.metric(d.re,d.ze)
solutions = []
def capture(a,b):
    u = spsolve(a,b)
    solutions.append(u)
    return u
reservoir.spsolve = capture
static = d.response(0)[0].real
omega = .0001
dynamic = d.response(omega)[0].real
print(f'independent_force_mass={(dynamic-static)/omega**2:.12g}',flush=True)
a,b = np.split(solutions[-1],2)
phase_per_speed = ((d.psi.real*b-d.psi.imag*a)/abs(d.psi)**2/(-1j*omega)).real
def at(index):
    in_neck = (d.z_index>=d.nfeed)&(d.z_index<d.neck_end)
    return phase_per_speed[in_neck&(d.r_index==index)]
velocity = par['hbar']/par['m']/rad*(3*at(cells-1)-at(cells-2))/2
for k, q in enumerate(np.array_split(velocity,4),1):
    print(f'independent_quarter {k} speed={q.mean():.12g}',flush=True)
