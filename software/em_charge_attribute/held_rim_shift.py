"""Signed rim reaction for a translated order profile and an attached layer."""

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.linalg import solve_banded
from scipy.special import expit, spherical_jn


def emit(label, value):
    print(f"{label} = {value}", flush=True)


# No finite contact energy is added to the bulk gradient term. Translate each
# radial slice into its own frame; only the axial derivative changes.
nu, radial, axial, slope = sp.symbols("nu radial axial slope", real=True)
change = ((axial - slope*nu*radial)**2 - axial**2)/2
angular = sp.factor(2*sp.pi*sp.integrate(change, (nu,-1,1)))
emit("angular excess gradient energy per kappa r^2 dr dw", angular)
assert angular.subs(slope,0) == 0

eps, ell, length = 0.02, 0.005, 37/20
k, amp = 0.5, 0.1
b_peak = amp**2*k**2/4  # constant inertia, unit shear modulus
delta = 0.002
def radius(w):
    s = w/length
    return 1+(3*s**2-2*s**3)/2


def axial_stiffness(w):
    r0 = radius(w)
    # s=(r-R)/ell; the split resolves the true radial domain at very small ell.
    def integrand(s):
        return (r0+ell*s)**2*(expit(s)*expit(-s))**2
    bounds = [-r0/ell, -10., 0., 10., np.inf]
    moment = sum(quad(integrand,a,b,epsabs=1e-14)[0] for a,b in zip(bounds,bounds[1:]))
    return 4*np.pi*eps**2*moment/(3*ell)  # per unit kappa_B


def layer_holds(ncell, shifts):
    # Stationary displacement field from discrete energy. Interior residuals
    # vanish; the first row of its unreduced gradient is the external hold.
    edges = np.linspace(0,length,ncell+1)
    h = length/ncell
    coeff = np.array([axial_stiffness(w) for w in (edges[:-1]+edges[1:])/2])/h
    bands = np.zeros((3,ncell-1))
    bands[1] = coeff[:-1]+coeff[1:]
    bands[0,1:] = bands[2,:-1] = -coeff[1:-1]
    rhs = np.zeros((ncell-1,len(shifts)))
    rhs[-1] = coeff[-1]*shifts
    interior = solve_banded((1,1),bands,rhs)
    near = -coeff[0]*interior[0]
    far = coeff[-1]*(shifts-interior[-1])
    return near, np.max(abs(near+far))


shifts = np.array([delta,-delta,delta/2,-delta/2,0])
slopes = []
for ncell in (16,64,256):
    near, residual = layer_holds(ncell, shifts)
    emit(f"{ncell} cells: holds per kappa for +d,-d,+d/2,-d/2,0", near.tolist())
    slopes.append(-(near[0]-near[1])/(2*delta))
    assert residual < 1e-13 and abs(near[0]-2*near[2]) < 1e-15
emit("pinned stiffness per kappa by grid", slopes)
emit("straight cylinder stiffness per kappa", axial_stiffness(0)/length)

# A distinct same-section probe. For the shear-energy term alone compute the
# volume overlap of a smooth depletion with the unperturbed standing fringe.
r0 = 1.
def form_integrand(r):
    return 4*np.pi*r**2*eps*expit((r0-r)/ell)*spherical_jn(0,2*k*r)
cuts = [0., r0-10*ell, r0, r0+10*ell, np.inf]
form = sum(quad(form_integrand,a,b,epsabs=1e-14)[0] for a,b in zip(cuts,cuts[1:]))
emit("midplane order-depletion form factor", form)
# Evaluate the reaction by independently integrating the translated shear
# density derivative on cross-sectional disks (radial/angle quadrature).
nodes, weights = np.polynomial.legendre.leggauss(24)
def wave_hold(y, phase):
    def integrand(r):
        phase_at_node = k*(r*nodes+y)+phase
        grad_b = k*b_peak*np.sin(2*phase_at_node)
        return 2*np.pi*eps*expit((r0-r)/ell)*r**2*np.dot(weights,grad_b)
    return sum(quad(integrand,a,b,epsabs=1e-16)[0] for a,b in zip(cuts,cuts[1:]))


phases = [0,np.pi/4,-np.pi/4,np.pi/2]
for phase in phases:
    hold = np.array([wave_hold(y,phase) for y in shifts])
    expected = b_peak*k*form*np.sin(2*(k*shifts+phase))
    emit(f"rim-carried fringe phase/pi={phase/np.pi:g}: baseline and odd slopes at d,d/2",
         [hold[-1], (hold[0]-hold[1])/(2*delta), (hold[2]-hold[3])/delta])
    assert np.max(abs(hold-expected)) < 1e-14

K, contact, x, z, y = sp.symbols("K contact x z y",real=True)
energy = (K*(y-z)**2 + contact*(z-x)**2)/2
optimum = sp.solve(sp.diff(energy,z),z)[0]
holding = sp.factor(sp.diff(energy,x).subs(z,optimum).subs(x,0))
emit("relaxed near hold for shifted far section with local contact", holding)
emit("free trace and pinned trace limits", [holding.subs(contact,0), sp.limit(holding,contact,sp.oo)])
