"""Estimate axial pressure shares by independently dilating each quarter of the well."""
import sys
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.special import expit
sys.path.insert(0,str(Path('software/stage1_solver/tools').resolve()))
from mouth_carrier_bend_budget import settings
from mouth_smooth_dipole_relaxation import AxialModes

point, _, _, _, archived = settings()
p = archived['parameters']
radius, _ = archived['phase_A']['q_star']
half = float(point['L'])/2
h=.03
mesh=AxialModes(radius,p,h,7.2)
w=mesh.w
r=(np.arange(len(w)-1)+1)*h
def weight(lo, hi):
    return np.clip(np.minimum(hi,w+h/2)-np.maximum(lo,w-h/2),0,h)/h

def scalar_energy(profile,delta):
    transition = expit((2*half-2*np.hypot(w,p['eps_w']))/p['delta_parallel'])
    strength = p['mu_out']**2-p['mu_in']**2
    confine = p['mu_out']**2 - strength*transition*expit(
          2*(radius-r[:,None]+delta*profile[None,:])/p['delta_r'])
    A=mesh.kinetic+sparse.diags(confine.ravel())
    eig,=eigsh(A,sigma=0,which='LM',k=1,return_eigenvectors=False,v0=mesh.seed,tol=2e-12)
    return p['I_wave']*np.sqrt(eig)
def load(profile):
    small=.00025
    return (scalar_energy(profile,-small)-scalar_energy(profile,small))/(2*small)
total=load(weight(0,half))
full=load(np.ones_like(w))
quarters=np.array([load(weight(i*half/4,(i+1)*half/4)) for i in range(4)])
moment=load(weight(0,half)*np.sin(np.pi*w/half)**4)/total
radial_mass=p['rho0']*4*np.pi/3*radius**3*2*half
print(f'directed_dilations step={h} physical_side_load={total:.12g} all_w_load={full:.12g} '
      f'fraction_on_physical_side={total/full:.12g} quarters={quarters/total} sum={quarters.sum()/total:.12g} '
      f'moment={moment:.12g} modal_mass={radial_mass*moment/3:.12g}',flush=True)
assert abs(quarters.sum()/total-1)<1e-5
