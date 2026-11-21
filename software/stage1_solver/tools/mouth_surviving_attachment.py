"""Probe the upstream attachment on matched frozen copies of a tagged 0PN mouth."""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_contact_crossing_force import REGISTRY_K, registry_force
from mouth_order_wave_pulse import digest, emit
from mouth_tagged_contact import force_on_ordered_end, live_attachment_stiffness


ARCHIVE = Path('research/pde_ledger_v3/_measurements/mouth_tagged_contact.json')
STEPS = (1e-5, 5e-6)


def source_case(archive, run, case):
    source = archive['verified_runs'][run]
    rows = source['literal_stdout']
    if source['worker_exit'] or source['guarded_exit'] or not source['zero_swap']:
        raise ValueError('Unguarded or incomplete input')
    row = next(r for r in map(json.loads, rows) if r['case'] == case)
    if (row['lobes'] != ['upstream', 'downstream'] or row['registry'] != .125 or
            row['mesh'] != 48 or row['status'] != 'broken'):
        raise ValueError('Unexpected tag preparation')
    return row


def stages(case):
    return dict(zip(('entry', 'half_center_order', 'before_neighbor_loss', 'after_neighbor_loss',
                     'before_tag_loss', 'after_tag_loss'),
                    [*case['probes'][:2], *case['events']], strict=True))


def isolate_from_archived_forces(frame, probe_positions, rest):
    # This channel uses only the two archived total force readings and geometry.
    # Continue the affine force to the sibling's stress-free slip to unload it.
    positions = np.asarray(probe_positions)
    forces = np.asarray(frame['probe_forces'])
    if not positions[0] <= rest[1] <= positions[1] or rest[0] == rest[1]:
        raise ValueError('Lobes cannot be isolated by this measurement')
    unloaded_neighbor_force = np.interp(rest[1], positions, forces)
    return unloaded_neighbor_force/(rest[0]-rest[1])


def test_attachment(case, law):
    rest = np.asarray(case['offsets'])
    history = []
    for name, f in stages(case).items():
        live = np.asarray(f['live_lobes'])
        order = np.asarray(f['lobe_order'])
        k = live_attachment_stiffness(order, law)*live
        x = f['ordered_end_slip']
        base = force_on_ordered_end(k, rest, x)
        np.testing.assert_allclose(base, f['force_on_ordered'], atol=2e-18, rtol=0.)
        transducer = isolate_from_archived_forces(f, case['probe_slips'], rest)
        probes = []
        for delta in STEPS:
            # Move the upstream wall anchor only, without changing either bond's
            # intrinsic rest offset. Hold material position, registry and order fixed.
            motion = np.array([delta, 0.])
            plus = force_on_ordered_end(k, rest+motion, x)
            minus = force_on_ordered_end(k, rest-motion, x)
            slope = (plus-minus)/(2*delta)
            probes.append(dict(step=delta, force_plus=plus, force_minus=minus,
                               change_plus=plus-base, change_minus=minus-base, slope=slope))
        emit(test='single_attachment_probe', law=law, stage=name, time=f['time'],
             depth=f['depth'], local_order=order[0], center_order=f['center_order'],
             live_lobes=live.tolist(), ordered_end_position=x, baseline_force=base,
             signed_probes=probes, slope_from_archived_forces=transducer)
        for probe in probes:
            np.testing.assert_allclose(probe['slope'], transducer, atol=5e-14, rtol=0.)
            np.testing.assert_allclose(probe['change_plus'], -probe['change_minus'], atol=2e-18)
        if not live[0] and any(probe['change_plus'] or probe['change_minus'] for probe in probes):
            raise ValueError('Disconnected tag contributed incremental force')
        history.append(transducer)
    return history


def sideways_force(registry, displacement, chi):
    """Orthogonal component of the trial central well at fixed material registry."""
    radius = np.hypot(registry, displacement)
    return registry_force(radius, chi, REGISTRY_K) * displacement/radius


def well_test(staged):
    # The old scalar force is periodic along a row. The two-coordinate continuation
    # uses only the tagged partner's well in its local cell (radius < half pitch).
    registry = .125
    for name, chi in staged:
        if chi/2 <= np.hypot(registry, max(STEPS)):
            raise ValueError('A finite probe would cross the geometric edge')
        force_ratio = -registry_force(registry, chi, REGISTRY_K)/registry
        probes = []
        for step in STEPS:
            plus, minus = [sideways_force(registry, sign*step, chi) for sign in (1., -1.)]
            probes.append(dict(step=step, force_plus=plus, force_minus=minus, slope=-(plus-minus)/(2*step)))
        extrapolated = (4*probes[1]['slope']-probes[0]['slope'])/3
        emit(test='central_well', stage=name, chi=chi, registry=registry, probes=probes,
             zero_step_slope=extrapolated, predicted_from_registry=force_ratio)
        np.testing.assert_allclose(extrapolated, force_ratio, rtol=0., atol=5e-15)
    # Uniformly spaced sites across a periodic pitch. The two probes use the same
    # fixed sites, including dead sites, and therefore sample the edge of the well.
    nodes = (np.arange(65536)+.5)/65536-.5
    roots, weights = leggauss(8)
    for chi in (.95, .75, .5, .3, .25):
        step = STEPS[-1]
        radius = np.hypot(nodes, step)
        # Explicit single-cell radial force: the two sites nearest the periodic
        # edge are already outside this well for every order sampled here.
        lateral = -REGISTRY_K*step*np.maximum(1-(2*radius/chi)**2, 0.)
        patch_slope = -lateral.mean()/step
        live_mean = weights @ (REGISTRY_K*(1-roots**2))/2
        emit(test='uniform_registry_population', chi=chi,
             active_fraction=float(np.mean(2*abs(nodes)<chi)), patch_slope=patch_slope,
             continuum_patch_slope=chi*live_mean, mean_live_slope=live_mean)
        np.testing.assert_allclose(patch_slope, chi*live_mean, rtol=0., atol=2e-9)


def main():
    archive = json.loads(ARCHIVE.read_text())
    order = source_case(archive, 'contact_tagged_work01', 'n48_line128')
    fixed = source_case(archive, 'contact_tagged_constant01', 'selected_tag')
    for key in ('fields', 'offsets', 'wall_site', 'ordered_site', 'intake', 'head', 'registry', 'span'):
        if order[key] != fixed[key]:
            raise ValueError('Unmatched preparations')
    for a, b in zip(stages(order).values(), stages(fixed).values(), strict=True):
        if any(a[key] != b[key] for key in ('time', 'depth', 'lobe_order', 'live_lobes')):
            raise ValueError('States sampled at different stages of conversion')
    emit(test='preparation', prior=str(ARCHIVE), prior_sha256=digest(ARCHIVE),
         wall_site=order['wall_site'], ordered_site=order['ordered_site'], tag='upstream',
         registry=order['registry'], intake=order['intake'], head=order['head'])
    for case, law in ((order, 'order'), (fixed, 'constant')):
        history = test_attachment(case, law)
        emit(test='tag_history', law=law, retained_at_last_live=history[-2]/history[0])
        np.testing.assert_allclose(history[2], history[3], atol=5e-14)
    sampled = stages(order)
    well_test([(key, sampled[key]['lobe_order'][0])
               for key in ('entry', 'half_center_order', 'after_neighbor_loss')] +
              [('upstream_order_0.30', .30), ('upstream_order_0.26', .26)])


if __name__ == '__main__':
    main()
