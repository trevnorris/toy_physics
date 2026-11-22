"""Independent depth-domain tangent calculation for a test population in a draining radial well.

An initially uniform sphere of tags enters with a radial speed distribution.
The first variation of its Cartesian position under a bend q H e is
a (n.e) n + b [e - (n.e) n]. All forces below are normalized by site mass.
The density dipole is the angular crossing count about the displaced slice
center. The axial guide fixes the depth clock independently of the bend.
"""
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import solve_ivp

from mouth_calibrated_flare_coating import read_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, TEMPERATURE, emit


def quadratic(nodes, u):
    a, b, c = nodes
    linear = -3*a + 4*b - c
    quadratic_part = 2*(a - 2*b + c)
    return a + u*(linear + u*quadratic_part), linear + 2*u*quadratic_part


def run(mesh=34):
    trace, intake = cubic_trace(mesh)
    calibration, _, _ = read_inputs(mesh)
    scale = 2*calibration['shared_scale']
    well_depth, alpha, _ = BETWEEN
    entry = 1.175
    half = trace.crossing(.5)
    gates = [entry, (entry+half)/2, half, *(trace.crossing(c) for c in (.25,.1,.01))]
    knots = np.unique(np.r_[gates, trace.edges[(trace.edges>entry)&(trace.edges<gates[-1])]])
    nodes, thermal_weights = hermgauss(12)
    thermal_weights /= np.sqrt(np.pi)
    launches = np.r_[0., np.sqrt(2*scale*TEMPERATURE)*nodes]
    weights = {'cold': np.r_[1., np.zeros(12)], 'thermal': np.r_[0., thermal_weights]}
    radius, radius_slope = trace.geometry.radius(entry)
    state = np.zeros((len(launches), 9))
    state[:,0] = radius
    state[:,1] = radius_slope * trace.speed(entry) + launches
    initial = state.copy()
    emit(test='tangent_setup', mesh=mesh, intake=intake, gates=gates, potential_scale=scale,
         well_width=PITCH, radial_offsets=launches.tolist(), gaussian_weights=thermal_weights.tolist(),
         ode='DOP853 in depth, split at archived Q2 element edges and reporting cuts',
         rtol=3e-11, atol=3e-13)

    def report(w, state):
        r, v, a, u, b, z, impulse, shape_impulse, integration_by_parts = state.T
        complement, d_complement = taper(w,entry,half)
        H = 1-complement
        H_dot = -d_complement*trace.speed(w)
        momentum = (u+2*z)/3
        dipole = 2*(b-H)/r
        angular_shift = (b-H)/r
        R, dRdw = trace.geometry.radius(w)
        chi = trace.sample(trace.chi,w)[0]
        exponent = np.exp(-alpha*(r-R)/PITCH)
        specific_energy = .5*(v-dRdw*trace.speed(w))**2 + scale*chi*well_depth*exponent*(exponent-2)
        for label, weight in weights.items():
            emit(test='tangent_crossing', supply=label, depth=w, order=chi, profile=H,
                 angular_density_dipole_slope=float(weight@dipole),
                 lab_angular_density_dipole_slope=float(weight@(2*b/r)),
                 mean_lateral_flux_slope=float(weight@momentum),
                 net_impulse_slope=float(weight@impulse),
                 bend_holder_impulse_slope=float(weight@shape_impulse),
                 bend_material_on_wall_impulse_slope=float(-weight@shape_impulse),
                 matched_shape_boundary_residual=float(weight@(H*momentum-integration_by_parts-shape_impulse)),
                 travel_angle_rate_slope=float(weight@((z-H_dot)/r-angular_shift*v/r)),
                 instantaneous_negative_radial_energy_fraction=float(weight@(specific_energy < 0.)),
                 max_particle_momentum_balance_error=float(np.max(abs(momentum-impulse))))
        np.testing.assert_allclose(impulse,momentum,atol=2e-8,rtol=1e-9)
        np.testing.assert_allclose(shape_impulse,H*momentum-integration_by_parts,atol=2e-8,rtol=1e-9)

    report(entry, initial)
    for low, high in zip(knots[:-1],knots[1:],strict=True):
        mid=(low+high)/2
        j=np.searchsorted(trace.edges,mid)-1
        left, right=trace.edges[j:j+2]
        width=right-left
        def geometry(w):
            t=(w-left)/width
            chi=quadratic(trace.chi[j],t)[0]
            psi,psi_t=quadratic(trace.psi[j],t)
            R, slope=trace.geometry.radius(w)
            axial_speed=(psi_t/(width*psi)).imag/(1+slope*slope)
            return R,chi,axial_speed
        _,order,speed=geometry(mid)
        np.testing.assert_allclose([order,speed],[trace.sample(trace.chi,mid)[0],trace.speed(mid)], atol=2e-13)
        def rhs(w, flat):
            r,v,a,u,b,z,_,_,_=flat.reshape(-1,9).T
            R,chi,speed=geometry(w)
            complement, d_complement=taper(w,entry,half)
            H=1-complement
            exp=np.exp(-alpha*(r-R)/PITCH)
            coeff=2*scale*chi*well_depth*alpha/PITCH
            f=coeff*exp*(exp-1)
            f_r=coeff*(alpha/PITCH)*exp*(1-2*exp)
            radial=f_r*(a-H)
            theta=f*(b-H)/r
            mean=(radial+2*theta)/3
            d_dt=np.column_stack((v,f,u,radial,z,theta,mean,H*mean,
                                  -speed*d_complement*(u+2*z)/3))
            return (d_dt/speed).ravel()
        sol=solve_ivp(rhs,(low,high),state.ravel(),method='DOP853',
                      rtol=3e-11,atol=3e-13,t_eval=[high])
        if not sol.success or not np.isfinite(sol.y).all():
            raise ArithmeticError(sol.message)
        state=sol.y[:,-1].reshape(len(launches),9)
        if high in gates: report(high,state)


if __name__=='__main__':
    run()
