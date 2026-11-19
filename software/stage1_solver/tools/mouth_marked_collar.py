"""Tagged material on the feeding side of a moving spherical mouth.

Forward paths use Cartesian potential flow and an actual boundary crossing.
Compare two independent normal boundary responses at the same inward volume
rate. Units use a carrier number density and a tagged-site number density of
one. The passive orientation-site association with those carriers is a trial.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

INNER, OUTER = 1., 1.2
INTAKE = .25
SINK = INTAKE/(4*np.pi)
TRANSIT = (OUTER**3-INNER**3)/(3*SINK)
AMPLITUDE = .012
DIRECTION = np.array([0., 1.])  # cylindrical transverse distance, z


def speed(t, amplitude):
    if 0. < t < TRANSIT:
        return amplitude*np.pi/(2*TRANSIT)*np.sin(np.pi*t/TRANSIT)
    return 0.


def velocity(t, pos, amplitude, beta):
    r = np.linalg.norm(pos, axis=-1, keepdims=True)
    n = pos/r
    translation = speed(t, amplitude)
    # Lab-frame gradient of the sink and decaying dipole potentials, then
    # transform to the displaced center. Positive inward rate on both faces.
    image = beta*(INNER/r)**3*translation*(3*n*n[..., 1:2]-DIRECTION)/2
    return -SINK*n/r**2+image-translation*DIRECTION


def transit_path(r, mu, amplitude, beta, stop, observe=None):
    initial = r*np.array([np.sqrt(1-mu**2),mu])
    def exit_boundary(t, pos):
        return np.linalg.norm(pos)-INNER
    exit_boundary.terminal = True
    exit_boundary.direction = -1
    path = solve_ivp(lambda t,pos: velocity(t,pos,amplitude,beta), (0.,stop), initial,
                     events=exit_boundary, dense_output=observe is not None,
                     method="DOP853", rtol=2e-11, atol=2e-13)
    assert path.success
    if observe is None:
        return path.t_events[0][0] if len(path.t_events[0]) else np.inf
    assert observe<path.t[-1]
    pos = path.sol(observe)
    return np.arctanh(pos[1]/np.linalg.norm(pos)), path.t_events[0][0]


def footprint(r, amplitude, beta, t):
    nodes, weights = leggauss(16)
    normals = np.stack((np.sqrt(1-nodes**2),nodes), axis=1)
    rate = -np.sum(normals*velocity(t,r*normals,amplitude,beta),axis=1)
    # Actual per physical area with unit density, and dipole coefficient.
    return 2*np.pi*r*r*(weights@rate), 1.5*(weights@(nodes*rate)), min(rate)


def initial_survival(mu, t, amplitude, beta):
    # For each initial latitude find the innermost starting radius surviving t
    # by integrating backward from the inner cut, then invert its latitude.
    def previous(end_mu):
        start = INNER*np.array([np.sqrt(1-end_mu**2),end_mu])
        path = solve_ivp(lambda tt,p: velocity(tt,p,amplitude,beta), (t,0.), start,
                         method="DOP853", rtol=2e-11, atol=2e-13)
        assert path.success
        pos = path.y[:,-1]
        return pos[1]/np.linalg.norm(pos), np.linalg.norm(pos)
    best = brentq(lambda end: previous(end)[0]-mu, -1.,1., xtol=1e-12)
    threshold = previous(best)[1]
    return np.clip((OUTER**3-threshold**3)/(OUTER**3-INNER**3),0.,1.)


def weak_prediction(mu, radius, amplitude, beta, t):
    t0 = (radius**3-INNER**3)/(3*SINK)
    def rest(tt):
        return (radius**3-3*SINK*tt)**(1/3)
    def shift(tt):
        r = rest(tt)
        return -speed(tt,amplitude)/r*(1+beta*(INNER/r)**3/2)
    def timing(tt):
        r = rest(tt)
        return mu*speed(tt,amplitude)/SINK*(beta*INNER**3/r-r*r)
    return quad(shift,0,t)[0], quad(timing,0,t0)[0]


def main():
    radius = np.cbrt((INNER**3+3*OUTER**3)/4)
    mu, read = .5, TRANSIT/3
    t0 = .75*TRANSIT
    print(f"model intake={INTAKE:g} R={INNER:g} outer={OUTER:g} transit={TRANSIT:.10g} "
          f"uniform_projected_sites_per_area={(OUTER**3-INNER**3)/(3*INNER**2):.10g} "
          f"tag_radius={radius:.10g} initial_mu={mu:g}")
    for beta in (0.,1.):
        readings = {}
        for multiplier in (0.,1.,-1.,.5,-.5):
            amp = AMPLITUDE*multiplier
            inner, inner_dip, inside_min = footprint(INNER,amp,beta,TRANSIT/2)
            outer, outer_dip, outside_min = footprint(OUTER,amp,beta,TRANSIT/2)
            # Invert only the FULL positive measured inner normal footprint.
            if multiplier==1.:
                inferred = 1-inner_dip/speed(TRANSIT/2,amp)
                print(f"inner_full_inferred_beta={inferred:.10g} configured={beta:g}")
                assert abs(inferred-beta)<1e-11
            y, t_exit = transit_path(radius,mu,amp,beta,1.3*TRANSIT,read)
            delta_y = y-np.arctanh(mu)
            readings[multiplier] = np.array([delta_y,t_exit-t0])
            survive = [initial_survival(latitude,TRANSIT/2,amp,beta) for latitude in (mu,-mu)]
            print(f"beta={beta:g} X_over_peak={multiplier:g} intake_inner={inner:.12g} "
                  f"intake_outer={outer:.12g} inner_dip={inner_dip:.10g} outer_dip={outer_dip:.10g} "
                  f"min_inward_rate={min(inside_min,outside_min):.8g} ring_delta_atanh_mu={delta_y:.10g} "
                  f"exit_shift={t_exit-t0:.10g} half_transit_initial_front_back={survive[0]:.9f},{survive[1]:.9f}")
            assert inside_min>0 and outside_min>0 and max(abs(inner-INTAKE),abs(outer-INTAKE))<1e-12
            if multiplier==0:
                assert max(abs(delta_y),abs(t_exit-t0),max(abs(np.array(survive)-.5)))<1e-9
                assert abs(transit_path(OUTER,mu,amp,beta,1.1*TRANSIT)-TRANSIT)<1e-8
        errors=[]
        for size in (1.,.5):
            observed = (readings[size]-readings[-size])/2
            predicted = np.array(weak_prediction(mu,radius,size*AMPLITUDE,inferred,read))
            error = abs(observed/predicted-1)
            errors.append(error)
            print(f"gate beta={beta:g} amplitude={size:g} measured_odd_slip={observed[0]:.10g} "
                  f"predicted={predicted[0]:.10g} measured_odd_exit_shift={observed[1]:.10g} "
                  f"predicted={predicted[1]:.10g} relative_errors={error[0]:.3g},{error[1]:.3g}")
        assert max(errors[0])<.007 and np.all(errors[1]<.35*errors[0])


if __name__ == "__main__":
    main()
