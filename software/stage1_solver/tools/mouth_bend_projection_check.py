"""Direct angular quadrature of the translated well and fixed-box force checks."""
from pathlib import Path
import sys

import numpy as np
from scipy import sparse
from scipy.integrate import quad
from scipy.sparse.linalg import eigsh, spsolve
from scipy.special import expit

sys.path.insert(0, str(Path("software/stage1_solver/tools").resolve()))
from mouth_carrier_bend_budget import settings
from mouth_open_reservoirs import ReservoirDomain
from mouth_smooth_dipole_relaxation import AxialModes


ref, _, _, flow, cfg = settings()
p = cfg['parameters']
radius = cfg['phase_A']['q_star'][0]
length = float(ref['L'])
half = length/2
h = .06
extent = 7.2
mesh = AxialModes(radius, p, h, extent)
w = mesh.w[None,:]
r = h*np.arange(1,len(w[0]))[:,None]
shape = np.where(w<half,np.sin(np.pi*w/half)**2,0)
axial = expit((length-2*np.sqrt(w*w+p['eps_w']**2))/p['delta_parallel'])
cx,cw = np.polynomial.legendre.leggauss(16)

def blocks(beta):
    a = [np.zeros((len(r),len(w[0]))) for _ in range(3)]
    for c, wt in zip(cx,cw):
        q = np.sqrt(r*r+(beta*shape)**2-2*r*beta*shape*c)
        v = p['mu_out']**2 - mesh.depth*axial*expit(2*(radius-q)/p['delta_r'])
        for k, f in enumerate([1,np.sqrt(3)*c,3*c*c]):
            a[k] += v*wt*f/2
    return [sparse.diags(v.ravel()) for v in a]

v00, _, v11 = blocks(0)
hs = mesh.kinetic+v00
hd = mesh.kinetic+v11+sparse.diags(2*mesh.centrifugal)
def ground(mat, seed):
    val, vec = eigsh(mat,k=1,sigma=1.,which="LM",v0=seed,tol=1e-12)
    return float(val[0]), vec[:,0]
ls, s = ground(hs,mesh.seed)
ld, d = ground(hd,mesh.seed)
ws = np.sqrt(ls)
wd = np.sqrt(ld)
free = p['I_wave']/(2*ws)
print(f'grid h={h} radius={radius:.12g} L={length:.12g} ws={ws:.12g} wd={wd:.12g}')
for b in [.0016,.0008]:
    a, c, e = blocks(b)
    bent = sparse.bmat([[mesh.kinetic+a,c],[c,mesh.kinetic+sparse.diags(2*mesh.centrifugal)+e]],format="csc")
    energy = np.sqrt(ground(bent,np.r_[s,np.zeros_like(s)])[0])
    direct = free*2*np.dot(s,(a-v00)@s)/b**2
    src = c@s/b
    overlap = (d@src)**2
    residual_src = src-(src@d)*d
    def rat(source,z):
        return source@spsolve(hd-sparse.eye(len(s),format='csc')*z**2,source)
    selected = 2*free*overlap/(ld-ls)
    relaxed = direct-2*free*rat(src,ws)
    def ir(nu):
        return free*(rat(src,ws+nu)+rat(src,ws-nu))
    m_taylor = (ir(.001)-ir(0))/.001**2
    resolv = spsolve(hd-sparse.eye(len(s),format='csc')*ls,src)
    cubic = resolv@spsolve(hd-sparse.eye(len(s),format='csc')*ls,resolv)
    m_exact = 2*free*(resolv@resolv+4*ls*cubic)
    # A one-part-in-10^7 detuning lets us evaluate the regular spectator resolvent.
    rest = rat(residual_src,wd+1e-7)
    cross = direct/2-free*rest
    print(f'b={b} eigenenergy_K={2*p["I_wave"]*(energy-ws)/b**2:.12g} direct_K={direct:.12g} response_K={relaxed:.12g} selected={selected:.12g} Mslow_secant={m_taylor:.12g} Mslow_exact={m_exact:.12g} Kcross_at_wd={cross:.12g}')
theta = lambda w: np.sin(np.pi*w/half)**2
deriv = lambda w: np.pi/half*np.sin(2*np.pi*w/half)
area = 4*np.pi*radius**2
radial_mass = p['rho0']*(4*np.pi/3)*radius**3*length
def integrate(fun):
    return quad(fun,0,half,epsabs=1e-12)[0]
def estimate(profile):
    return radial_mass*integrate(lambda w: profile(w)*theta(w)**2)/(3*integrate(profile))
print(f'shape Ma={radial_mass:.12g} M_uniform={estimate(lambda w:1):.12g} M_central={estimate(lambda w:.02+theta(w)**4):.12g} K_one={p["sigma"]*area*integrate(lambda w:deriv(w)**2)/3:.12g}')
dom = ReservoirDomain(radius,half,radius,12,2,2,3,p)
dom.set_current(float(flow['reference_inward_current'])/p['m'])
dom.ac = dom.metric(dom.re,dom.ze)
static = dom.response(0.)[0].real
for nu in [.001,.0005,.00025]:
    conjugate = dom.response(-nu)[0]
    forward = dom.response(nu)[0]
    value = (.5*(conjugate.real+forward.real)-static)/nu**2
    print(f'carrier n=12 outer=4 nu={nu:.8g} chi_zero={static:.12g} induced_mass_secant={value:.12g} rt_error={abs(conjugate-forward.conjugate()):.3g}')
