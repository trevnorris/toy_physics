"""Passive two-attachment test along a measured wall transit at 0PN.

Lengths are tangential displacements in units of a small imposed mismatch.
Energies per contact area are in units of full-order stiffness times that
mismatch squared. A no-incoming-wave termination receives the sideways work.
"""
import numpy as np
from scipy.integrate import solve_ivp

from mouth_contact_tracking import WallTrace
from mouth_order_wave_pulse import emit


def measure(trace, span, impedance, offsets):
    shifts = np.array([-span/2, span/2])
    start = trace.crossing(.95)
    cuts = {level: trace.crossing(level) for level in (.5, .01, .0001, .000001)}
    end = cuts[.000001]
    # Break at element boundaries for the carrier and for both stiffness traces.
    shifted_edges = [trace.edges-offset for offset in [0., *shifts]]
    all_edges = np.unique(np.concatenate([*shifted_edges, [start], list(cuts.values())]))
    edges = all_edges[(all_edges >= start) & (all_edges <= end)]

    def stiffness(w):
        chi, gradient = trace.sample(trace.chi, w+shifts)
        return chi/2, gradient/2

    def unloaded(k):
        # Read the zero of the holding-force curve from two independent trial slips.
        probes = np.array([-2., 2.])
        loads = (k*(probes[:, None]-offsets)).sum(axis=1)
        return probes[0]-loads[0]*np.diff(probes)[0]/np.diff(loads)[0]

    initial, _ = stiffness(start)
    slip = unloaded(initial)
    initial_slip = slip
    stored = .5*np.dot(initial, (slip-offsets)**2)
    state = np.array([slip, 0., 0.])  # slip, work into the line, energy lost to softening
    readings = []
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        width = right-left
        # The two gradients are linear within this interval. Verify both one-sided
        # endpoint values, without crossing into the neighboring finite element.
        d1 = stiffness(left+width*.25)[1]
        d3 = stiffness(left+width*.75)[1]
        if np.max([1.5*d1-.5*d3, 1.5*d3-.5*d1]) > 1e-11:
            raise ValueError('An attachment stiffens during the measured transit')
        if min(stiffness(right)[0]) <= 0:
            raise ValueError('Loss of a resolved attachment before the terminal reading')

        def rates(w, y):
            # Take the interior limit of the three FE traces at each subinterval end.
            interior = np.clip(w, left+width*1e-10, right-width*1e-10)
            k, dk = stiffness(interior)
            extension = y[0]-offsets
            force = -np.dot(k, extension)
            drift = trace.speed(interior)
            return [force/(impedance*drift), force**2/(impedance*drift),
                    -.5*np.dot(dk, extension**2)]

        solution = solve_ivp(rates, (left, right), state, method='DOP853', atol=2e-13, rtol=2e-10)
        if not solution.success:
            raise ArithmeticError(solution.message)
        state = solution.y[:, -1]
        for level, cut in cuts.items():
            if right != cut:
                continue
            k, _ = stiffness(right)
            remaining_energy = .5*np.dot(k, (state[0]-offsets)**2)
            readings.append(dict(order=level, depth=cut, unloaded_offset=unloaded(k),
                                 actual_offset=state[0], line_work=state[1],
                                 softening_loss=state[2], remaining_energy=remaining_energy))
    emit(span=span, impedance=impedance, offsets=offsets.tolist(),
         initial_unloaded=initial_slip, initial_energy=stored, readings=readings)
    for row in readings:
        np.testing.assert_allclose(row['line_work']+row['softening_loss']+row['remaining_energy'],
                                   stored, rtol=0, atol=3e-10)
        if row['line_work'] < 0 or row['softening_loss'] < -1e-13:
            raise ArithmeticError('Negative passive budget')
    if span == 0 or offsets[0] == offsets[1]:
        np.testing.assert_allclose([r['line_work'] for r in readings], 0, atol=1e-18)
    return readings


def main():
    for mesh in (34, 48):
        for head in (0., .5):
            trace = WallTrace(mesh, head)
            emit(test='junction_setup', mesh=mesh, head_fraction=head, intake=trace.intake,
                 field_sha256=trace.hashes)
            offsets = np.array([-1., 1.])
            measure(trace, 0., 1., offsets)
            main_case = measure(trace, .04, 1., offsets)
            if (mesh, head) == (48, .5):
                measure(trace, .04, 1., np.zeros(2))
                reverse = measure(trace, .04, 1., -offsets)
                for row, mirror in zip(main_case, reverse, strict=True):
                    np.testing.assert_allclose(mirror['unloaded_offset'], -row['unloaded_offset'], atol=1e-12)
                    np.testing.assert_allclose(mirror['line_work'], row['line_work'], atol=1e-12)
                for impedance in (.5, 2.):
                    measure(trace, .04, impedance, offsets)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tagged-readouts', metavar='DIRECTORY',
                        help='track one site pair through breakup and save signed force and velocity traces')
    parser.add_argument('--registry-offset', type=float, metavar='PITCHES',
                        help='select a single initially connected tag on the fine mesh')
    parser.add_argument('--live-stiffness', choices=('order', 'constant'), default='order',
                        help='stiffness of a surviving attachment for a selected tag')
    args = parser.parse_args()
    if args.registry_offset is not None and not args.tagged_readouts:
        parser.error('--registry-offset requires --tagged-readouts DIRECTORY')
    if args.live_stiffness != 'order' and args.registry_offset is None:
        parser.error('--live-stiffness constant requires --registry-offset PITCHES')
    if args.tagged_readouts:
        from mouth_tagged_contact import main as track_tag
        track_tag(args.tagged_readouts, args.registry_offset, args.live_stiffness)
    else:
        main()
