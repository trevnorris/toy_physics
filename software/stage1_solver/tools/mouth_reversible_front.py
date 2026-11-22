"""Reversible coherent carrier and crystalline order on an unpinned planar patch.

The conserved order momentum pi enters both the phase gradient and the physical
mass current. Equal-density port phases specify the current through the BVP.
This is a longitudinal material test with no geometric support or shear wave.
"""
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid, simpson, solve_bvp
from scipy.linalg import eigh
from scipy.special import expit

BIAS = .4510084282857874
KAPPA = 2.


def potential(n):
    return (n - 1)**2 * (n**3 + 2*n**2 + 3*n + 4) / 20


def derivative(n):
    return (n**4 - 1)/4


def frequencies(half_solution, span, inertia, j, elements):
    """Fixed-phase second variation in all density/order modes with fixed end fields.

    Minimize the positive quadratic pi and phase terms first. Phase stiffness
    gives a nonlocal rank-one contribution to the resulting energy Hessian.
    Report static eigenvalues in the joint L2 norm, not temporal frequencies.
    """
    nodes = np.linspace(-span/2, span/2, elements+1)
    h = span/elements
    points, weights = leggauss(4)
    x = nodes[:-1,None] + h*(1+points)/2
    f, _, chi, dchi, _ = half_solution.sol(abs(x).ravel()).reshape(5,*x.shape)
    chi = np.where(x>=0, chi, 1-chi)
    n = f*f
    chi_pot = 2*(1-6*chi*(1-chi))
    # C is the phase resistance. E_phase = (Delta theta)^2 / (2*C).
    resistance = np.sum((1/n+inertia*dchi**2/n**2)*weights)*h/2
    res_f = -2/f**3 - 4*inertia*dchi**2/f**5
    res_k = 2*inertia*dchi/f**4
    density_pot = 2*derivative(n)+4*n**4 - .5*j*j*(6/f**4+20*inertia*dchi**2/f**6)
    cross = 4*j*j*inertia*dchi/f**5
    order_stiff = KAPPA - j*j*inertia/f**4
    dimension = elements+1
    H = np.zeros((2*dimension,2*dimension))
    mass = np.zeros_like(H)
    g = np.zeros(2*dimension)
    for a in range(2):
        phi_a = (1 + (2*a-1)*points)/2
        da = (2*a-1)/h
        ia = np.arange(elements)+a
        g[ia] += (res_f*phi_a)@weights*h/2
        g[dimension+ia] += (res_k*da)@weights*h/2
        for b in range(2):
            phi_b = (1 + (2*b-1)*points)/2
            db = (2*b-1)/h
            ib = np.arange(elements)+b
            v = (phi_a*phi_b)@weights*h/2
            mass[ia,ib] += v
            mass[dimension+ia,dimension+ib] += v
            H[ia,ib] += ((density_pot*phi_a*phi_b+da*db)@weights)*h/2
            H[dimension+ia,dimension+ib] += ((chi_pot*phi_a*phi_b+order_stiff*da*db)@weights)*h/2
            mixed = ((cross*phi_a*db)@weights)*h/2
            H[ia,dimension+ib] += mixed
            H[dimension+ib,ia] += mixed
    H += j*j/resistance*np.outer(g,g)
    free = np.r_[1:elements,dimension+1:2*dimension-1]
    take = np.ix_(free,free)
    vals, vecs = eigh(H[take], mass[take], subset_by_index=(0,2))
    # The order translation mode in the interior. Dirichlet ends and a P1
    # trial space lift the exact full-line zero at finite resolution.
    center_gradient = half_solution.sol(abs(nodes))[3]
    trial = np.r_[np.zeros(dimension),center_gradient][free]
    metric = mass[take]
    overlap = abs(trial @ metric @ vecs[:,0])/np.sqrt(trial @ metric @ trial)
    return dict(elements=elements, static_eigenvalues=vals.tolist(),
                lowest_mode_overlap_with_order_translation=overlap)


def calculate(span, bias, inertia, check_hessian=False):
    def equations(x, y, pars):
        j, = pars
        f, df, c, dc, _ = y
        n = f*f
        u = j/n
        stiffness = KAPPA-inertia*u*u
        if np.min(stiffness) <= 0 or np.min(f) <= .5:
            raise ArithmeticError('Trial left regular subcritical branch')
        ddf = 2*f*(derivative(n)+.5*u*u+inertia*u*u*dc*dc/n)
        ddc = (2*c*(1-c)*(1-2*c)-4*inertia*u*u*df*dc/f)/stiffness
        phase_slope = u*(1+inertia*dc*dc/n)
        return np.array([df,ddf,dc,ddc,phase_slope])

    def boundary(left, right, pars):
        return np.array([left[1],left[2]-.5,left[4],right[0]-1,right[2],right[4]-bias/2])

    x = np.linspace(0, span/2, int(span*2)+1)
    c = expit(-x)
    guess = np.array([np.ones(x.size), np.zeros(x.size), c, -c*(1-c), bias*x/span])
    sol = solve_bvp(equations,boundary,x,guess,p=[bias/span],tol=1e-10,max_nodes=14000)
    if not sol.success:
        raise ArithmeticError(sol.message)
    x = np.linspace(0,span/2,6001)
    f,df,c,dc,phase=sol.sol(x)
    j, = sol.p
    n = f*f
    u = j/n
    pi = inertia*u*dc
    # Stationary four-dimensional unidirectional Noether stress per unit area.
    # Both the absence of order dissipation and port energy j*mu follow from
    # the canonical action (mu=0). Stress is independently checked with
    # the coordinate derivative of the collocated solution.
    normal = .5*df**2 + .5*(KAPPA-inertia*u**2)*dc**2 + .5*n*u*u - potential(n) - c*c*(1-c)**2
    rhs = equations(x,sol.sol(x),sol.p)
    d2f = sol.sol(x,1)[1]
    # Momentum flux expressed as the carrier pressure, density-gradient
    # stress, convective flux, and the reciprocal order stress.
    stress = (n*u*u + n*derivative(n) - potential(n) + (df*df-f*d2f)/2
              + (KAPPA+inertia*u*u)*dc**2/2-c*c*(1-c)**2)
    entrainment = 2*simpson(pi*dc/n,x=x)
    velocity_integral = 2*simpson(u,x=x)
    row = dict(span=span, spatial_phase_bias=bias, order_inertia=inertia, current_per_area=j,
               plain_velocity_integral=velocity_integral, entrained_phase_integral=entrainment,
               phase_error=velocity_integral+entrainment-bias, minimum_density=n.min(),
               center_order_slope=dc[0], minimum_effective_order_stiffness=min(KAPPA-inertia*u*u),
               total_normal_flux_mean=normal.mean(), normal_flux_range=np.ptp(normal),
               two_stress_max_error=max(abs(stress-normal)),
               max_bvp_rms_residual=max(sol.rms_residuals),
               endpoint_order_stress=KAPPA*dc[-1]**2/2, phase_chemical_frequency=0,
               net_steady_energy_input=0, irreversible_power=0)
    if check_hessian:
        row['static_spectrum']=[frequencies(sol,span,inertia,j,e) for e in (160,320,640)]
    if abs(row['phase_error'])>1e-9 or row['normal_flux_range']>2e-10 or row['two_stress_max_error']>3e-9:
        raise ArithmeticError(f'Independent phase or stress closure failed: {row}')
    return row


def run():
    results=[calculate(32.,BIAS,1.,True), calculate(32.,-BIAS,1.),
             calculate(32.,0.,1.), calculate(32.,BIAS,0.), calculate(48.,BIAS,1.)]
    print(json.dumps(dict(tests=results),indent=2,allow_nan=False),flush=True)
    ref, reverse, rest, uncoupled, longer=results
    checks=(abs(ref['current_per_area']+reverse['current_per_area'])<1e-10,
            abs(rest['current_per_area'])<1e-12,
            abs(rest['center_order_slope']+.25)<1e-9,
            abs(uncoupled['entrained_phase_integral'])<1e-14,
            ref['current_per_area']<uncoupled['current_per_area'],
            longer['endpoint_order_stress']<ref['endpoint_order_stress'],
            all(line['static_eigenvalues'][1]>0 and line['lowest_mode_overlap_with_order_translation']>.995
                for line in ref['static_spectrum']))
    if not all(checks):
        raise ArithmeticError(f'Conservative front controls: {checks}')


if __name__=='__main__':
    run()
