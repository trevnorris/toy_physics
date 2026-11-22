"""One coherent crystalline-field check on a regular curved order front.

Carrier phase and crystalline density-wave phase are independent. A complex
crystal amplitude A=chi*f*exp(i*m*phi) may lose registry in a smooth core.
At matched curvature, compare free registration with one unit of externally
prepared dislocation winding. Gravity is absent from the material fit.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad, simpson, solve_bvp
from scipy.optimize import brentq

from mouth_free_front_wave_head import emit

R = 100.0
PATCH = 8.0
HALF_SPAN = 16.0
XI = 1.0


def circumference(r):
    return R*np.sin(r/R)


def volume(r):
    # Transverse radial measure per radian and per unit length on the
    # great circle of the radius-R three-sphere.
    return circumference(r)*np.cos(r/R)


def registered_profile(grid, tolerance):
    # Write f=C*g, C=R sin(r/R). This resolves the smooth axis of the vortex
    # without excising a hole. solve_bvp handles -3g'/r at the origin.
    def equation(r, y):
        g, dg = y
        ang = r/R
        regular = -r/(3*R**2)
        use = abs(ang)>1.e-4
        regular[use] = 1/(R*np.tan(ang[use]))-1/r[use]
        rhs = (3/R**2 + ((circumference(r)*g)**2-1)/XI**2)*g
        return np.array([dg, rhs-(3*regular-np.tan(ang)/R)*dg])

    def boundary(inner, outer):
        return np.array([inner[1], circumference(PATCH)*outer[0]-1])

    r = np.linspace(0., PATCH, grid)
    guess = 1/np.sqrt(2*XI**2+r*r)
    sol = solve_bvp(equation, boundary, r, [guess, -r*guess**3],
                    S=np.array([[0.,0.],[0.,-3.]]), tol=tolerance, max_nodes=12000)
    if not sol.success:
        raise ArithmeticError(sol.message)

    def eval_profile(x):
        g,dg = sol.sol(x)
        return circumference(x)*g, np.cos(x/R)*g+circumference(x)*dg
    return eval_profile, max(sol.rms_residuals)


def column_conductance(f_squared):
    # Conductance per reference mid-area under unit drive, K0=delta=1.
    # K_perp=(1-f^2)+f^2*tanh(s/2)^2. At exact registry the inverse
    # stiffness has a nonintegrable singularity at s=0.
    if f_squared == 1.:
        return 0.
    def resistance(s):
        k = (1-f_squared)+f_squared*np.tanh(s/2)**2
        return (R/(R+s))**3/k
    return 1./quad(resistance,-HALF_SPAN,HALF_SPAN,points=[0.],
                   epsabs=1e-9,epsrel=1e-10)[0]


def measure(grid, tolerance):
    profile, residual = registered_profile(grid, tolerance)
    r = np.linspace(0., PATCH, 3001)
    f,df = profile(r)
    if abs(f[0])>1e-12 or abs(f[-1]-1.)>1e-11 or min(df)<-2e-6:
        raise ArithmeticError('Registry amplitude does not resolve the prepared core and far lattice')
    phase_gradient = np.divide(f,circumference(r),out=np.full_like(r,df[0]),where=r!=0)
    terms = [.5*df**2, .5*phase_gradient**2, .25*(1-f*f)**2/XI**2]
    # GL energy on the mid-layer, per unit extent along the dislocation, in
    # units of its stiffness. Same nonnegative action on each preparation.
    energies = [2*np.pi*simpson(volume(r)*term,x=r) for term in terms]
    half_radius=brentq(lambda x: profile(x)[0]-.5,0.,PATCH)
    x,wt = leggauss(96)
    # r=PATCH*(1-y^2) resolves square-root conductance near the fixed rim.
    y=(x+1.)/2
    radii=PATCH*(1.-y*y)
    weights=wt*PATCH*y*volume(radii)
    amplitude=profile(radii)[0]
    mean_f2=weights@(amplitude**2)/sum(weights)
    lower=weights@np.array([column_conductance(v**2) for v in amplitude])/sum(weights)
    upper=column_conductance(mean_f2)
    reference=column_conductance(0.)
    angle=np.linspace(0.,2*np.pi,257)
    emit(test='preparation',xi=XI,delta=1.,curvature_radius=R,patch_radius=PATCH,
         normal_half_span=HALF_SPAN,trial_K='1-4*chi*(1-chi)*f^2',
         envelope_bulk_constraint='f(patch_radius)=1',
         smooth_core_radius_at_half_amplitude=half_radius,
         amplitude_at_axis=f[0],amplitude_at_edge=f[-1],rms_residual=residual,
         reference_transparent_conductance=reference)
    for ensemble, winding in [('intact_free',0),('smooth_exit_free',0),
                              ('exit_imposed_positive',1),('exit_imposed_negative',-1),
                              ('intact_imposed_positive',1)]:
        phase_loop=np.exp(1j*winding*angle)
        degree=np.sum(np.angle(phase_loop[1:]*phase_loop[:-1].conjugate()))/(2*np.pi)
        defect=winding!=0
        emit(test='registry',ensemble=ensemble,winding_measured=degree,
             core_energy_components=energies if defect else [0.,0.,0.],
             energy_per_unit_defect_length=sum(energies) if defect else 0.,
             mean_amplitude_squared=mean_f2 if defect else 1.,
             relative_conductance_lower=lower/reference if defect else 0.,
             relative_conductance_upper=upper/reference if defect else 0.)
    if not (0<lower<=upper<=reference and min(energies)>0):
        raise ArithmeticError('Variational resistance bounds or positive defect energy failed')


if __name__=='__main__':
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--grid',type=int,required=True)
    arg.add_argument('--tol',type=float,default=1e-7)
    options=arg.parse_args()
    measure(options.grid,options.tol)
