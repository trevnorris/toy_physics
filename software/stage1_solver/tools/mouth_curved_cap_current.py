"""Check whether the Route A planar layer connects identical quiet, flowing baths.

Above unit bath sound speed the stationary density eigenvalues are imaginary.
Starting on the unique decaying order tail, integrate inward to the half-order
level and compare with the reflected trajectory from the other bath. The
independent variable is log(chi) to resolve the asymptotic tail without a
finite reflecting reservoir. This test presumes a monotone single interface.
"""
import numpy as np
from scipy.integrate import solve_ivp

from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA, WELL, stationary_front
from mouth_poly_shallow_release import emit


def tail_to_center(j, eps, tol):
    decay = np.sqrt(2*WELL/(KAPPA-ORDER_INERTIA*j*j))
    density_term = 2*ORDER_INERTIA*j*j*decay**2/(4*decay**2-4*(1-j*j))
    def derivative(log_c, state):
        f, df, g = state
        c = np.exp(log_c)
        u = j/f**2
        pi = ORDER_INERTIA*u*g
        ddf = 2*f*((f**8-1)/4+(u*u-j*j)/2+u*pi*g/f**2)
        ddc = (2*WELL*c*(1-c)*(1-2*c)-4*ORDER_INERTIA*u*u*g*df/f)/(KAPPA-ORDER_INERTIA*u*u)
        return np.array([df,ddf,ddc])*c/g
    starting = [1+density_term*eps**2,-2*decay*density_term*eps**2,-decay*eps*(1-eps)]
    sol = solve_ivp(derivative,np.log([eps,.5]),starting,method='DOP853',rtol=tol,atol=tol*eps)
    if not sol.success:
        raise ArithmeticError(sol.message)
    f, df, g = sol.y[:,-1]
    # Relative to the uniform bath. Integrated normal stationarity in the
    # fixed-current Routhian and its conserved normal stress on the centerline.
    potential = f**10/20-f**2/4+.2-.5*j*j*(f-1/f)**2
    potential += WELL/16
    stationary_stress = (.5*df**2+.5*(KAPPA-ORDER_INERTIA*j*j/f**4)*g**2-potential)
    return dict(center_density=f*f,receiver_df=df,feed_df=-df,
                density_slope_jump=2*df,order_slope=g,first_integral_error=stationary_stress)


def run():
    # A subsonic positive control uses the native free density-decay amplitude.
    for j in (.28334069168539316,.6,.8):
        try:
            fields = stationary_front(j,20.)
        except ArithmeticError as exc:
            emit(kind='subsonic_native',current=j,quiet=False,reason=str(exc))
        else:
            f, ddf, _, g = fields(0.)
            emit(kind='subsonic_native',current=j,quiet=True,center_density=float(f*f),
                 center_order_slope=float(g),center_density_amplitude_curvature=float(ddf))
    for j in (1.001,1.02,1.1,1.25,1.4,1.6,1.65):
        rows = []
        for eps,tol in [(1.e-4,1.e-9),(1.e-6,2.e-11),(1.e-8,2.e-12)]:
            result = tail_to_center(j,eps,tol)
            emit(kind='quiet_supersonic_tail',current=j,initial_order=eps,rtol=tol,
                 propagating_stationary_density_wavenumber=2*np.sqrt(j*j-1),**result)
            rows.append(result)
        if abs(rows[2]['receiver_df']-rows[1]['receiver_df'])>1.e-8 or abs(rows[2]['first_integral_error'])>1.e-8:
            raise ArithmeticError('Tail convergence or conserved stationary normal stress failed')


if __name__=='__main__':
    run()
