"""Read one junction's identity, zero-load slip and ordered-end work through breakup."""
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve_banded

from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import digest, emit


PRESTRAIN = 1e-3
SPAN = .04
CUT = .25
LINE_LENGTH = 8.


def live_attachment_stiffness(chi, law):
    if law == 'order':
        return chi/2
    if law == 'constant':
        return np.full_like(chi, .5)
    raise ValueError('Unknown live contact stiffness law')


def force_on_ordered_end(k, offsets, slip):
    return np.sum(k*(offsets-np.asarray(slip)[..., None]), axis=-1)


def force_probe(k, offsets):
    """Displace a frozen copy. Probe motion never enters the work history."""
    slips = np.array([-2., 2.])*PRESTRAIN
    forces = force_on_ordered_end(k, offsets, slips)
    if not np.any(k):
        return forces, None
    if not forces[0] > 0 > forces[1]:
        raise ValueError('Frozen force probe does not bracket zero')
    root = slips[0] - forces[0]*np.diff(slips)[0]/np.diff(forces)[0]
    np.testing.assert_allclose(k @ (offsets-root), 0., atol=1e-16)
    return forces, float(root)


def depth_steps(trace, gates, shifts, dx):
    # Split at every force discontinuity and FE change; integrate the carrier clock.
    breaks = np.unique(np.r_[gates, *(trace.edges-s for s in [0., *shifts])])
    breaks = breaks[(breaks >= min(gates)) & (breaks <= max(gates))]
    roots, weights = leggauss(5)
    depth = [breaks[0]]
    for left, right in zip(breaks[:-1], breaks[1:], strict=True):
        count = max(1, int(np.ceil(trace.transit(left, right, 5)/(.22*dx))))
        depth.extend(np.linspace(left, right, count+1)[1:])
    depth = np.asarray(depth)
    width = np.diff(depth)[:, None]
    samples = depth[:-1, None]+width*(roots+1)/2
    dt = np.sum(width*weights/(2*trace.speed(samples)), axis=1)
    return depth, dt


def run_case(trace, mesh, dx, registry, offsets, name, folder, frozen=False, span=SPAN,
             stiffness_law='order'):
    shifts = np.array([-span/2, span/2])
    start, end = (trace.crossing(c) for c in (.95, .01))
    # Fixed material features occupy an angular row of pitch one. The transverse
    # shear and the registry are perpendicular. Each branch keeps its partner ID
    # until the shrinking window no longer covers its marked feature.
    threshold = 2*abs(registry)
    if not np.isfinite(threshold) or np.any(trace.sample(trace.chi, start+shifts)[0] <= threshold):
        raise ValueError('Choose a registry offset with both attachments engaged at entry')
    breaks = np.full(len(shifts), np.inf)
    reaches_edge = trace.sample(trace.chi, end+shifts)[0] <= threshold
    if np.any(reaches_edge) and not frozen:
        breaks[reaches_edge] = trace.crossing(threshold)-shifts[reaches_edge]
    gates = [start, end, trace.crossing(.5), *breaks[np.isfinite(breaks)]]
    depth, dt = depth_steps(trace, gates, shifts, dx)
    chi = trace.sample(trace.chi, depth[:, None]+shifts)[0]
    center_order = trace.sample(trace.chi, depth)[0]
    if frozen:
        chi[:] = chi[0]
        center_order[:] = center_order[0]
    gap = chi/2-abs(registry)
    # Compare two distinct hypotheses for the stiffness of an individual live
    # attachment at identical geometric registry and connectivity events.
    live_k = live_attachment_stiffness(chi, stiffness_law)
    intact = np.ones_like(shifts, dtype=bool)
    entry_k = live_k[0]
    _, baseline = force_probe(entry_k, offsets)
    initial_energy = .5*entry_k @ (offsets-baseline)**2
    times = np.r_[0., np.cumsum(dt)]

    # Unit-area shear strip, rho=mu=1. The endpoint owns half a cell of ordered
    # material. Its velocity is a dynamic state; it is never assigned from force.
    n = round(LINE_LENGTH/dx)
    cut = round(CUT/dx)
    if abs(cut*dx-CUT)>1e-12 or times[-1] >= LINE_LENGTH:
        raise ValueError('Receiver must resolve the cut and outlast the record before reflection')
    mass = np.full(n, dx)
    mass[0] /= 2
    u, v = np.zeros(n), np.zeros(n)  # displacements relative to uniform entry offset
    diag = np.full(n, 2/dx)
    diag[0] = 1/dx
    band = np.zeros((3, n))
    events, probes, midpoints = [], [], []
    work = cut_work = broken_energy = softening = 0.
    rupture = None

    def read(tag, j):
        k = np.where(intact, live_k[j], 0.)
        forces, slip = force_probe(k, offsets)
        x = baseline+u[0]
        return dict(stage=tag, depth=depth[j], time=times[j], center_order=center_order[j],
                    lobe_order=chi[j].tolist(), overlap_margin=gap[j].tolist(),
                    live_lobes=intact.tolist(), wall_end_motion=0., ordered_end_slip=x,
                    ordered_end_velocity=v[0], force_on_ordered=force_on_ordered_end(k, offsets, x),
                    probe_forces=forces.tolist(), zero_force_slip=slip,
                    endpoint_work=work, cut_work=cut_work,
                    receiver_energy=.5*(mass @ v**2 + np.diff(np.r_[u, 0.]) @ np.diff(np.r_[u, 0.])/dx),
                    softening=softening, rupture_loss=broken_energy)

    probes.append(read('entry', 0))
    for j, step in enumerate(dt, start=1):
        k_left = np.where(intact, live_k[j-1], 0.)
        k_right = np.where(intact, live_k[j], 0.)
        k = (k_left+k_right)/2
        # Implicit midpoint. Use the force on the receiver from its external
        # contact and read its velocity from the integrated material field.
        f = np.diff(np.r_[u, 0.])/dx
        f[1:] -= f[:-1].copy()
        x_left = baseline+u[0]
        f[0] += k @ (offsets-x_left)
        scale = step**2/4
        band[0, 1:] = band[2, :-1] = -scale/dx
        band[1] = mass+scale*diag
        band[1, 0] += scale*k.sum()
        velocity = solve_banded((1, 1), band, mass*v+step*f/2, check_finite=False)
        next_u = u+step*velocity
        x_right = baseline+next_u[0]
        force = k @ (offsets-(x_left+x_right)/2)
        cut_force = ((u[cut-1]-u[cut])+(next_u[cut-1]-next_u[cut]))/(2*dx)
        work += step*force*velocity[0]
        cut_work += step*cut_force*velocity[cut]
        softening += .25*(k_left-k_right) @ ((offsets-x_left)**2+(offsets-x_right)**2)
        kprobes, intercept = force_probe(k, offsets)
        midpoints.append([times[j-1]+step/2, step, (depth[j]+depth[j-1])/2,
                          *intact, *k, *kprobes, np.nan if intercept is None else intercept,
                          (x_left+x_right)/2, force, velocity[0], cut_force, velocity[cut]])
        u, v = next_u, 2*velocity-v
        lost = (depth[j] == breaks) & intact & (not frozen)
        if np.any(lost):
            np.testing.assert_allclose(gap[j, lost], 0., atol=1e-12)
            events.append(read('before_loss', j))
            broken_energy += .5*k_right[lost] @ (offsets[lost]-x_right)**2
            intact[lost] = False
            events.append(read('after_loss', j))
            if not np.any(intact):
                rupture = events[-1]
        if np.any(gap[j, intact] < -1e-12) or np.any(gap[j, ~intact] > 1e-12):
            raise ValueError('Contact history inconsistent with observed overlap')
        if abs(center_order[j]-.5)<1e-12:
            probes.append(read('half_order', j))
    finish = read('end_of_observation', len(depth)-1)
    probes.append(finish)
    cut_energy = .5*(mass[cut:] @ v[cut:]**2+np.sum(np.diff(np.r_[u[cut:], 0.])**2)/dx)
    stored = .5*(live_k[-1]*intact) @ (offsets-baseline-u[0])**2
    eps_energy = PRESTRAIN**2
    for sample in [*events, *probes]:
        np.testing.assert_allclose(sample['receiver_energy'], sample['endpoint_work'], atol=eps_energy*2e-10, rtol=0)
    np.testing.assert_allclose(cut_work, cut_energy, atol=eps_energy*2e-10, rtol=0)
    np.testing.assert_allclose(work+softening+broken_energy+stored, initial_energy, atol=eps_energy*2e-10, rtol=0)
    if (rupture is None) != bool(np.any(np.isinf(breaks))):
        raise ValueError('Unexpected break or censoring')
    if rupture is not None and work != rupture['endpoint_work']:
        raise ValueError('Extinct tag performed work')
    if frozen or span == 0 or np.all(offsets == 0.):
        np.testing.assert_allclose(work, 0., atol=eps_energy*1e-15)

    columns = ('time dt depth_midpoint live_upstream live_downstream k_upstream k_downstream '
               'probe_force_minus probe_force_plus zero_force_slip ordered_end_slip '
               'force_on_ordered ordered_end_velocity cut_force cut_velocity').split()
    output = folder / f'{name}.npz'
    np.savez_compressed(output, columns=columns, samples=midpoints)
    row = dict(test='tagged_contact', case=name, mesh=mesh, dx=dx, registry=registry, span=span,
               head=trace.head, intake=trace.intake, frozen=frozen, offsets=offsets.tolist(),
               live_stiffness=stiffness_law,
               wall_site='wall-0', ordered_site='ordered-0', lobes=['upstream','downstream'],
               probe_slips=[-2*PRESTRAIN, 2*PRESTRAIN], cut=CUT, initial_energy=initial_energy,
               status='still_connected' if rupture is None else 'broken', events=events, probes=probes,
               raw_path=str(output), raw_sha256=digest(output), fields=trace.hashes)
    emit(**row)
    return row


def main(folder, registry_offset=None, stiffness_law='order'):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    offsets = PRESTRAIN*np.array([-1., 1.])
    if registry_offset is not None:
        trace = WallTrace(48, .5)
        run_case(trace, 48, 1/128, registry_offset, offsets, 'selected_tag', folder,
                 stiffness_law=stiffness_law)
        return
    for mesh in (34, 48):
        trace = WallTrace(mesh, .5)
        for dx in (1/32, 1/64, 1/128):
            run_case(trace, mesh, dx, .125, offsets, f'n{mesh}_line{round(1/dx)}', folder)
    def control(name, registry=.125, offset_scale=1., **kwargs):
        return run_case(trace, 48, 1/128, registry, offsets*offset_scale, name, folder, **kwargs)
    control('alternate_tag', registry=.25)
    control('nonconverting', frozen=True)
    control('coincident', span=0.)
    control('no_prestrain', offset_scale=0.)
    reverse, half = control('reversed', offset_scale=-1.), control('half', offset_scale=.5)
    np.testing.assert_allclose(reverse['events'][-1]['endpoint_work'],
                               4*half['events'][-1]['endpoint_work'], rtol=2e-11)
