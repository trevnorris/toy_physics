"""A cold freely cohesive sheet crossing a held cubic flare as a dilute parcel.

Sites share prescribed w(t) but all three visible coordinates are dynamical.
The only external visible force is the smooth radial well. Identical parcels
can be injected successively without interacting. A Gaussian density peak of
one parcel is reconstructed from positions at each accepted integration step.
"""
import argparse
import json

import numpy as np
from scipy.optimize import minimize

from mouth_calibrated_flare_coating import read_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_conversion import NEUTRAL
from neutral_patch_sites import BETWEEN, COUNT, WITHIN, emit, morse

START, GAUGE, END = 1.175, 1.55, 1.65
Q = np.array([0., 2e-5, -2e-5, 1e-5, -1e-5])
I, J = np.triu_indices(COUNT, 1)
# For a bond oriented i->j, gradient in the relative coordinate is force on i.
INCIDENCE = (np.eye(COUNT)[I] - np.eye(COUNT)[J]).T
FD = 1e-7


def planar_minimum():
    frame = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']
    seed = np.array(frame['rested_xy'][:COUNT])

    def energy(flat):
        sites = flat.reshape(COUNT, 2)
        d = sites[J]-sites[I]
        r = np.linalg.norm(d, axis=1)
        e, derivative = morse(r, *WITHIN)
        return float(e.sum()), (-INCIDENCE @ (d*(derivative/r)[:, None])).ravel()

    sol = minimize(energy, seed.ravel(), jac=True, method='L-BFGS-B',
                   options={'ftol': 5e-16, 'gtol': 1e-10, 'maxiter': 3000, 'maxcor': 24})
    points = sol.x.reshape(COUNT, 2)
    points = points-points.mean(axis=0)
    # Finite-difference every unrestrained gradient. Two translations and one
    # rotation are the null modes of a noncollinear planar minimum.
    eps = 2e-6
    hessian = np.column_stack([(energy(sol.x+d)[1]-energy(sol.x-d)[1])/(2*eps)
                               for d in eps*np.eye(sol.x.size)])
    levels = np.linalg.eigvalsh((hessian+hessian.T)/2)
    e, grad = energy(points.ravel())
    residual = np.max(np.abs(grad))
    emit(test='sheet_preparation', source=str(NEUTRAL), frame='frames1321.json/central/rested_xy',
         sheet='first 12 sites', optimizer=str(sol.message), original_energy=energy(seed.ravel())[0],
         minimum_energy=e, max_abs_free_force_native=residual, lowest_four_hessian=levels[:4].tolist(),
         centered_minimum_native=points.tolist())
    if residual > 8e-6 or np.max(abs(levels[:3])) > 2e-5 or levels[3] <= 1e-5:
        raise ArithmeticError('Unconstrained two-dimensional sheet failed the force or minimum test')
    return points


def inlet(planar, trace, args):
    z = np.array([START])
    radius, slope = trace.geometry.radius(z)
    x, y = (PITCH*planar).T
    if args.reflect:
        x = -x
    # Orthographic lift of the two tangent coordinates onto a sphere centered
    # on the x axis. The displayed launch pole at azimuth zero is +lab y.
    radial = np.sqrt(radius[0]**2-x*x-y*y)
    c, s = np.cos(args.azimuth), np.sin(args.azimuth)
    pos = np.column_stack((x, radial*c-y*s, radial*s+y*c))
    vel = pos*(slope[0]*trace.speed(z)[0]/radius[0])
    pos[:, 0] += args.origin
    return tuple(np.repeat(a[None], len(Q), axis=0) for a in (pos, vel))


def density_geometry(points, peak):
    """Gradient and minus Hessian of log density in Gaussian bandwidth units."""
    offset = points-peak[:, None, :]
    logweight = -.5*np.sum(offset**2, axis=-1)
    w = np.exp(logweight-logweight.max(axis=1, keepdims=True))
    p = w/w.sum(axis=1, keepdims=True)
    mean = np.sum(p[:, :, None]*offset, axis=1)
    centered = offset-mean[:, None, :]
    curvature = np.eye(3)-(centered*p[:, :, None]).transpose(0, 2, 1) @ centered
    density = np.exp(logweight).sum(axis=1)
    return mean, curvature, density


def find_peaks(points, starts):
    """Follow a nearby stationary density maximum, with mean shift away from its basin center."""
    peak = starts.copy()
    for _ in range(160):
        gradient, curvature, density = density_geometry(points, peak)
        error = np.linalg.norm(gradient, axis=1)
        if np.max(error) < 2e-12:
            break
        positive = np.linalg.eigvalsh(curvature)[:, 0] > 1e-4
        step = gradient.copy()
        if np.any(positive):
            step[positive] = np.linalg.solve(curvature[positive], gradient[positive, :, None])[:, :, 0]
        # A bounded Newton trial must ascend. Mean shift is the monotone fallback.
        step /= np.maximum(1., np.linalg.norm(step, axis=1)/.45)[:, None]
        trial = peak+step
        delta = points-trial[:, None, :]
        improved = np.exp(-.5*np.sum(delta**2, axis=-1)).sum(axis=1) >= density*(1.-5e-15)
        peak += np.where(improved[:, None], step, gradient)
    gradient, curvature, density = density_geometry(points, peak)
    return peak, np.linalg.norm(gradient, axis=1), curvature, density


class Marker:
    def __init__(self, pos, bandwidth):
        self.width = bandwidth*PITCH
        self.peak = np.zeros((len(Q), 3))
        self.min_curve = np.full(len(Q), np.inf)
        self.max_step = np.zeros(len(Q))
        self.max_gradient = np.zeros(len(Q))
        self.failed_nodes = np.zeros(len(Q), dtype=int)
        self.nearest_changes = np.zeros(len(Q), dtype=int)
        self.nearest = None
        self.track(pos, initial=True)
        self.tag = self.nearest.copy()

    def track(self, pos, initial=False):
        self.centroid = pos.mean(axis=1)
        self.points = (pos-self.centroid[:, None, :])/self.width
        peak, g, curvature, density = find_peaks(self.points, self.peak)
        lowest = np.linalg.eigvalsh(curvature)[:, 0]
        jump = np.linalg.norm(peak-self.peak, axis=1) if not initial else np.zeros(len(Q))
        self.failed_nodes += (g > 1e-8) | (lowest <= 1e-5) | (jump > .3)
        self.min_curve = np.minimum(self.min_curve, lowest)
        self.max_gradient = np.maximum(self.max_gradient, g)
        self.max_step = np.maximum(self.max_step, jump)
        self.peak, self.density, self.curvature = peak, density, curvature
        distances = np.linalg.norm(self.points-peak[:, None, :], axis=2)
        nearest = np.argmin(distances, axis=1)
        if not initial:
            self.nearest_changes += nearest != self.nearest
        self.nearest = nearest

    def read(self, pos):
        # Independent starts at every site plus centroid test for competing
        # basins. This is a finite search, not a proof of global uniqueness.
        n = COUNT+1
        starts = np.concatenate((self.points, np.zeros((len(Q), 1, 3))), axis=1)
        modes, residual, curv, _ = find_peaks(np.repeat(self.points, n, axis=0), starts.reshape(-1, 3))
        converged = (residual < 1e-8) & (np.linalg.eigvalsh(curv)[:, 0] > 1e-5)
        d = np.linalg.norm(modes.reshape(len(Q), n, 3)-self.peak[:, None, :], axis=-1)
        levels, axes = np.linalg.eigh(self.curvature)
        contrasts = {}
        for radius in (.25, 1.):
            directions = np.concatenate((axes.transpose(0, 2, 1), -axes.transpose(0, 2, 1)), axis=1)
            probes = self.peak[:, None, :]+radius*directions
            offset = self.points[:, None, :, :]-probes[:, :, None, :]
            boundary = np.exp(-.5*np.sum(offset**2, axis=-1)).sum(axis=-1)
            contrasts[str(radius)] = (1-boundary.max(axis=1)/self.density).tolist()
        absolute = self.centroid+self.width*self.peak
        original_site = pos[np.arange(len(Q)), self.tag]
        return absolute, original_site, dict(
            density_in_unit_site_masses=(self.density/(np.sqrt(2*np.pi)*self.width)**3).tolist(),
            log_density_minus_hessian_eigenvalues_sigma_units=levels.tolist(),
            fractional_density_drop_six_principal_directions=contrasts,
            starts_returning_to_peak=(converged.reshape(len(Q), n) & (d < 1e-5)).sum(axis=1).tolist(),
            other_peak_starts=(converged.reshape(len(Q), n) & (d >= 1e-5)).sum(axis=1).tolist(),
            unresolved_starts=(~converged).reshape(len(Q), n).sum(axis=1).tolist(),
            min_curvature_any_node=self.min_curve.tolist(), failed_tracking_nodes=self.failed_nodes.tolist(),
            max_gradient_any_node=self.max_gradient.tolist(), max_centroid_frame_jump_sigma=self.max_step.tolist(),
            nearest_inlet_id=self.tag.tolist(), nearest_current_id=self.nearest.tolist(),
            nearest_id_switches=self.nearest_changes.tolist(), peak_minus_centroid_in_pitches=(
                self.width*self.peak/PITCH).tolist())


def matched(x):
    change = x-x[0]
    slopes = [(x[i]-x[i+1])/(2*Q[i]) for i in (1, 3)]
    return dict(raw=x.tolist(), minus_zero=change.tolist(), central_slopes=np.asarray(slopes).tolist())


class Passage:
    def __init__(self, trace, planar, scale, args):
        self.pos, self.vel = inlet(planar, trace, args)
        self.v0 = self.vel[:, :, 0].mean(axis=1)
        self.marker = Marker(self.pos, args.bandwidth)
        self.scale, self.law, self.origin = scale, args.law, args.origin
        self.chi0 = trace.sample(trace.chi, np.array([START]))[0][0]
        self.ledger = np.zeros((3, len(Q)))
        self.elapsed = self.max_dt = 0.
        self.count = 0

    def rhs(self, pos, vel, at):
        inv_speed, radius, chi, shape, shape_prime = at
        d = pos[:, J, :]-pos[:, I, :]
        r = np.linalg.norm(d, axis=2)
        strength = self.scale*chi/self.chi0 if self.law == 'order' else self.scale
        derivative = morse(r/PITCH, *WITHIN)[1]
        acceleration = INCIDENCE @ (d*(strength*derivative/(PITCH*r))[:, :, None])
        direct = np.zeros(len(Q))
        if self.law != 'wall-off':
            relative = pos.copy()
            relative[:, :, 0] -= (self.origin+Q*shape)[:, None]
            distance = np.linalg.norm(relative, axis=2)
            gradient = morse(BETWEEN[2]+(distance-radius)/PITCH, *BETWEEN)[1]
            acceleration -= relative*(self.scale*chi*gradient/(PITCH*distance))[:, :, None]
            # Virtual wall translation at fixed absolute positions. Use the
            # factored difference of two Morse energies to retain the small
            # x signal of an almost tangent launch. Both energies use the same
            # depth, order and R. Pair energies cancel under wall translation.
            dx = relative[:, :, 0]
            transverse_sq = np.sum(relative[:, :, 1:]**2, axis=2)
            rp = np.sqrt((dx-FD)**2+transverse_sq)
            rm = np.sqrt((dx+FD)**2+transverse_sq)
            ep = np.exp(-BETWEEN[1]*(rp-radius)/PITCH)
            em = np.exp(-BETWEEN[1]*(rm-radius)/PITCH)
            delta_e = em*np.expm1(BETWEEN[1]*4*dx*FD/((rp+rm)*PITCH))
            direct = self.scale*chi*BETWEEN[0]*(delta_e*(ep+em-2)).mean(axis=1)/(2*FD)
        vx = vel[:, :, 0].mean(axis=1)
        ledger = np.array((direct*inv_speed, shape*direct*inv_speed, shape_prime*vx))
        return vel*inv_speed, acceleration*inv_speed, ledger

    def advance(self, nodes, coeff):
        for j, h in enumerate(np.diff(nodes)):
            k = 2*j
            p, v = self.pos, self.vel
            a1, b1, l1 = self.rhs(p, v, coeff[k])
            a2, b2, l2 = self.rhs(p+.5*h*a1, v+.5*h*b1, coeff[k+1])
            a3, b3, l3 = self.rhs(p+.5*h*a2, v+.5*h*b2, coeff[k+1])
            a4, b4, l4 = self.rhs(p+h*a3, v+h*b3, coeff[k+2])
            self.pos = p+h/6*(a1+2*a2+2*a3+a4)
            self.vel = v+h/6*(b1+2*b2+2*b3+b4)
            self.ledger += h/6*(l1+2*l2+2*l3+l4)
            dt = h/6*(coeff[k, 0]+4*coeff[k+1, 0]+coeff[k+2, 0])
            self.elapsed += dt
            self.max_dt = max(self.max_dt, dt)
            self.count += 1
            if not (np.isfinite(self.pos).all() and np.isfinite(self.vel).all() and np.isfinite(self.ledger).all()):
                raise FloatingPointError(f'Nonfinite parcel at depth {nodes[j+1]}')
            self.marker.track(self.pos)

    def read(self, w, radius, chi, shape):
        peak, tag, diagnostics = self.marker.read(self.pos)
        vx = self.vel[:, :, 0].mean(axis=1)
        net = vx-self.v0
        direct, weighted, by_parts_term = self.ledger
        by_parts = shape*vx-by_parts_term  # H(START)=0
        centroid = self.pos.mean(axis=1)
        center = self.origin+Q*shape
        relative = self.pos.copy()
        relative[:, :, 0] -= center[:, None]
        separation = np.linalg.norm(self.pos[:, I]-self.pos[:, J], axis=2)
        emit(test='cut', depth=w, elapsed=self.elapsed, n_steps=self.count, max_actual_dt=self.max_dt,
             order=chi, radius=radius, shape=shape, wall_center_x=center.tolist(),
             peak_xyz=peak.tolist(), peak_x=matched(peak[:, 0]), centroid_xyz=centroid.tolist(),
             centroid_x=matched(centroid[:, 0]), original_nearest_site_x=matched(tag[:, 0]),
             incoming_mean_vx=self.v0.tolist(), outgoing_mean_vx=matched(vx),
             peak_shift_minus_wall_shift=((peak[:, 0]-peak[0, 0])-Q*shape).tolist(),
             holding_per_mass_flux=matched(net), direct_holding_per_mass_flux=matched(direct),
             shape_holding_per_mass_flux=matched(weighted), shape_by_parts=matched(by_parts),
             net_minus_direct=(net-direct).tolist(), shape_minus_by_parts=(weighted-by_parts).tolist(),
             min_pair_gap_pitches=(separation.min(axis=1)/PITCH).tolist(),
             rms_patch_radius_pitches=(np.sqrt(np.mean(np.sum((self.pos-centroid[:, None, :])**2,
                                                             axis=2), axis=1))/PITCH).tolist(),
             rms_site_wall_gap_pitches=(np.sqrt(np.mean((np.linalg.norm(relative, axis=2)-radius)**2,
                                                        axis=1))/PITCH).tolist(), **diagnostics)


def run(args):
    if not all(map(np.isfinite, (args.dt, args.origin, args.azimuth, args.bandwidth))) or min(args.dt, args.bandwidth) <= 0:
        raise ValueError('Time resolution and Gaussian bandwidth must be positive and finite')
    trace, rate = cubic_trace(args.mesh)
    patch, _, _ = read_inputs(args.mesh)
    planar = planar_minimum()
    obj = Passage(trace, planar, 2*patch['shared_scale'], args)
    cuts = np.unique([START, 1.2, trace.crossing(.5), trace.crossing(.1), trace.crossing(.01), GAUGE, END])
    emit(test='setup', mesh=args.mesh, intake=rate, mouth_radius=1., q=Q.tolist(), bend_interval=[1.2, 2.],
         cuts=cuts.tolist(), law=args.law, origin=args.origin, azimuth=args.azimuth, reflected=args.reflect,
         max_dt_requested=args.dt, site_count=COUNT, within=list(WITHIN), normal=list(BETWEEN),
         microscopic_length=PITCH, specific_energy_scale=obj.scale, gaussian_bandwidth_pitches=args.bandwidth,
         finite_difference_translation=FD, equal_mass_weights=1/COUNT,
         marker='3D Gaussian KDE of unit site masses; centroid seed, then centroid-advected previous maximum',
         maximum_track_jump_sigma=.3, minimum_certified_curvature=1e-5, maximum_certified_gradient=1e-8,
         holding='external force on held fixture per dilute parcel mass flow; material recoil has opposite sign',
         launch='spherical lift of free planar minimum; local radial velocities Rprime times carrier speed',
         supply='independent identical packets with common prescribed w(t); no interpacket contacts')
    if abs(rate-.25)>1e-9 or cuts[0] != START or cuts[-1] != END:
        raise ValueError('Unmatched carrier or sections')
    def read(w):
        z = np.array([w])
        obj.read(w, trace.geometry.radius(z)[0][0], trace.sample(trace.chi, z)[0][0], 1-taper(z, 1.2, 2.)[0][0])
    read(START)
    edges = np.unique(np.r_[cuts, trace.edges[(trace.edges > START) & (trace.edges < END)]])
    for a, b in zip(edges[:-1], edges[1:], strict=True):
        n = max(1, int(np.ceil(trace.transit(a, b, 8)/args.dt)))
        while True:
            nodes = np.linspace(a, b, n+1)
            z = np.linspace(a, b, 2*n+1)
            z[0], z[-1] = np.nextafter(a, b), np.nextafter(b, a)
            inv = 1/trace.speed(z)
            max_dt = ((b-a)/(6*n)*(inv[:-1:2]+4*inv[1::2]+inv[2::2])).max()
            if max_dt <= args.dt*(1+1e-12):
                break
            n = max(n+1, int(np.ceil(n*max_dt/args.dt)))
        radius = trace.geometry.radius(z)[0]
        chi = trace.sample(trace.chi, z)[0]
        if np.min(chi)<0 or np.max(chi)>1:
            raise ValueError('Order outside the physical range')
        taper_value, taper_derivative = taper(z, 1.2, 2.)
        obj.advance(nodes, np.column_stack((inv, radius, chi, 1-taper_value, -taper_derivative)))
        if b in cuts:
            read(b)
    transit = trace.transit(START, END, 12)
    emit(test='endpoint', transit_from_quadrature=transit, transit_from_steps=obj.elapsed,
         absolute_time_error=abs(transit-obj.elapsed), finite=True)
    np.testing.assert_allclose(obj.elapsed, transit, atol=1e-9, rtol=1e-9)


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--mesh', type=int, required=True, choices=(26, 34))
    cli.add_argument('--dt', type=float, required=True, help='maximum integration step in carrier time')
    cli.add_argument('--law', choices=('full', 'order', 'wall-off'), default='full')
    cli.add_argument('--origin', type=float, default=0., help='common absolute x displacement')
    cli.add_argument('--azimuth', type=float, default=0., help='inlet rotation about the lab x axis, radians')
    cli.add_argument('--reflect', action='store_true', help='reflect the inlet sheet and velocity in x')
    cli.add_argument('--bandwidth', type=float, default=.85, help='density instrument sigma in native Morse lengths')
    run(cli.parse_args())
