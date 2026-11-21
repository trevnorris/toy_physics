"""Independent fixed-shape pressure projection, translated well, and force inversion."""
import sys
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import null_space
from scipy.sparse.linalg import eigsh, spsolve
from scipy.special import expit

sys.path.insert(0, str(Path('software/stage1_solver/tools').resolve()))
from mouth_carrier_bend_budget import settings
from mouth_open_reservoirs import ReservoirDomain
from mouth_smooth_dipole_relaxation import AxialModes

source, _, _, run, model = settings()
p = model['parameters']
R = model['phase_A']['q_star'][0]
H = float(source['L'])/2
nu = np.array([0., .01, .02])

# Translate the actual well and differentiate its two angular matrix elements
# by quadrature, without derivatives of its confining profile. Evaluate the
# even-frequency combination using its convergent Taylor series in omega^2.
h = .03
mesh = AxialModes(R, p, h, 7.2)
r = h*np.arange(1,len(mesh.w))[:,None]
w = mesh.w[None,:]
f = np.where(w < H, np.sin(np.pi*w/H)**2, 0)
depth = (p['mu_out']**2-p['mu_in']**2)*expit((2*H-2*np.sqrt(w**2+p['eps_w']**2))/p['delta_parallel'])
def potential(radius):
    return p['mu_out']**2-depth*expit(2*(R-radius)/p['delta_r'])
base = potential(r)
hs = mesh.kinetic+sparse.diags(base.ravel())
es,v = eigsh(hs,k=1,which='LM',sigma=1.,v0=mesh.seed,tol=1e-11)
ws = np.sqrt(es[0])
s = v[:,0]
hd = hs + sparse.diags(2*mesh.centrifugal)
from scipy.sparse.linalg import splu
G = splu((hd-sparse.eye(hd.shape[0])*es[0]).tocsc()).solve
beta = .0003
vpp,vp = np.zeros_like(base),np.zeros_like(base)
for angle,weight in zip(*np.polynomial.legendre.leggauss(20)):
    translated = potential(np.sqrt(r*r+(beta*f)**2-2*r*beta*f*angle))
    vpp += (translated-base)*weight/beta**2
    vp += (translated-base)*weight*np.sqrt(3)*angle/(2*beta)
q = vp.ravel()*s
moments = []
u = q.copy()
for n in range(10):
    u = G(u)
    moments.append(q@u)
def wave(omega):
    pos,neg = omega**2+2*ws*omega,omega**2-2*ws*omega
    chi = sum(mu*(pos**n+neg**n) for n,mu in enumerate(moments))
    return p['I_wave']/(2*ws)*(s@(vpp.ravel()*s)-chi)
for x in nu:
    print(f'translated_well h={h} omega={x:g} holding_per_b={wave(x):.12g}',flush=True)

# Independent generalized force from local wall pressure, tested against
# the complex-susceptibility sign and refined grid. Broaden the far layer.
cells = 32
d = ReservoirDomain(R,H,R,cells,2,4,3,p)
def expand(x):
    return x+40*R*(np.maximum(x/R-2,0)/4)**3
d.re = expand(d.re)
d.ze = H+expand(d.ze-H)
d.volume,d.conductance,d.inlet,d.outlet,ang,d.operator = d.metric(d.re,d.ze)
d.set_current(float(run['reference_inward_current'])/p['m'])
def stretch(x):
    return x+0.7j*np.maximum(x-10*R,0)**3/(36*R)**2
side = (d.r_index == cells-1) & (d.w_index < cells)
inner = (d.r_index == cells-2) & (d.w_index < cells)
def trace(v):
    return (9*v[side]-v[inner])/8
b = trace(d.psi)
profile = lambda x: np.sin(np.pi*x/H)**2
fmid = profile((np.arange(cells)+.5)*d.dw)
fedge = profile(np.arange(1,cells)*d.dw)
for x in nu:
    vol,_,_,_,ang,lap = d.metric(stretch(d.re),H+stretch(d.ze-H)) if x else d.metric(d.re,d.ze)
    wgt = vol/3
    stiffness = lap/3+sparse.diags(d.kinetic*ang+1.25*d.eos_k*(abs(d.psi)**8-d.bath_density**4)*wgt)
    g = 10*d.eos_k*abs(d.psi)**6*wgt
    c = g*d.psi.real*d.psi.imag
    sym = 1j*x*d.hbar*wgt
    mat = sparse.bmat([[stiffness+sparse.diags(g*d.psi.real**2),sparse.diags(c-sym)],
                       [sparse.diags(c+sym),stiffness+sparse.diags(g*d.psi.imag**2)]],format='csc')
    u = spsolve(mat,(d.spatial+1j*x*d.temporal)/2)
    a,v = map(trace,np.split(u,2))
    pressure = (1j*x*d.hbar*(b.real*v-b.imag*a)
                +2.5*d.eos_k*(d.bath_density**4-abs(b)**8)*(b.real*a+b.imag*v))
    grad = np.diff(b.real)*np.diff(a)+np.diff(b.imag)*np.diff(v)
    load = -4*np.pi*R**2/3*(d.dw*(fmid@pressure)-2*d.kinetic*(fedge@grad)/d.dw)
    power = sum(d.cuts(x,u,round(1.25*cells),round(.75*cells)))
    print(f'pressure n={cells} outer_R=46 omega={x:g} holding_per_b={load.real:.12g} '
          f'holding_imag_per_b={load.imag:.12g} emitted_power={power:.12g}',flush=True)
    assert abs(power+.5*x*load.imag)<1e-17

# Independent held-out residual from the rounded forward-prediction table.
def table(dir,kind):
    rows = [dict(i.split('=',1) for i in line.split()[1:])
            for line in (Path('_scratch/s11c')/dir/'stdout').read_text().splitlines()
            if line.startswith(kind+' ')]
    return rows
field = np.array([float(row['holding_per_b']) for row in table('bend_hold_fluid48','carrier')])
field += np.array([float(row['holding_per_b']) for row in table('bend_hold_wave015','wave')])
for row in table('bend_hold_readout','conditional'):
    reduced = np.array([float(row[f'F{i}']) for i in range(3)])/1e-5-field
    m1,m2 = (reduced[0]-reduced[1:])/nu[1:]**2
    print(f"fitted_printed_totals sheet={row['sheet']} K={reduced[0]:.12g} "
          f'M01={m1:.12g} M02={m2:.12g} high_holdout={(reduced@[3.,-4.,1.])*1e-5:.12g}',flush=True)
    assert abs(m1-m2)<1e-7
# A fourth-order wall memory term is detectable in the third reading, and a
# frequency dependence zero at all three probes is invisible to these data.
dispersive = nu**4
print(f'quartic_null_test={(dispersive@[3.,-4.,1.]):.12g}')
assert np.isclose(dispersive@[3.,-4.,1.],1.2e-7)
