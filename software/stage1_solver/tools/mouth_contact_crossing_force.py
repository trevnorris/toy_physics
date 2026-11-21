"""Two independent instruments on a candidate material patch in the draining throat.

A line of neighboring sites registers along the radial row. A separate shear
component acts along an angular tangent. Both components use pitch units.
All new microscopic constants are trials, held fixed across all drain paths.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_carrier_contact import follow_cell
from mouth_contact_age import PRIOR, emit, trace_path
from mouth_free_contact import Carrier


# The radial periodic well and anchor have stiffnesses .25 and 1. The
# transverse link and its elastic anchor have unit stiffness; drag is unit
# in the transverse direction and .5 in registry. Units are local model units.
REGISTRY_K = 0.25
REGISTRY_DRAG = 0.5
PUSHES = np.array([0., 1.e-4, -1.e-4, 5.e-5, -5.e-5])


def registry_force(coordinate, width, strength):
    gap = coordinate - np.floor(coordinate + .5)
    return -strength * gap * np.maximum(1 - (2 * gap / width)**2, 0.)


def frames(t, width, shift, count, phase, strength):
    """Let actual features move relative to carrier-mounted elastic anchors."""
    anchor = (np.arange(count) + phase) / count
    offset = np.zeros(count)
    for _ in range(64):
        offset = registry_force(anchor + offset, width[0], strength)
    if np.max(abs(offset - registry_force(anchor + offset, width[0], strength))) > 1.e-12:
        raise ArithmeticError("Initial material row is out of equilibrium")
    entry_offset = offset.copy()
    for j, time in enumerate(t):
        if j:
            dt = time - t[j-1]
            initial_force = registry_force(anchor+shift[j-1]+offset, width[j-1], strength) - offset
            middle = offset + .5 * dt * initial_force / REGISTRY_DRAG
            mid_anchor = anchor + .5 * (shift[j-1] + shift[j])
            mid_width = .5 * (width[j-1] + width[j])
            offset += dt * (registry_force(mid_anchor+middle, mid_width, strength) - middle) / REGISTRY_DRAG
        coord = anchor + shift[j] + offset
        partner = np.floor(coord + .5).astype(int)
        active = abs(coord - partner) < width[j] / 2
        yield active, partner, offset - entry_offset


def image_crossing(t, width, shift, count, strength):
    # This instrument observes an unforced material row and never reads a force.
    slip = occupancy_error = 0.
    for j, (active, partner, delta) in enumerate(frames(t, width, shift, count, .173, strength)):
        if j == 0:
            retained = active.copy()
            entry_count = int(active.sum())
        else:
            retained &= active & (partner == previous)
        previous = partner
        slip = max(slip, np.max(abs(delta)))
        occupancy_error = max(occupancy_error, abs(active.mean() - width[j]))
    return dict(sites=count, entry=entry_count, exit=int(active.sum()),
                original=int(retained.sum()), exit_fraction=active.mean(),
                survival=retained.sum()/active.sum(), max_relative_excursion=slip,
                max_occupancy_order_difference=occupancy_error)


def read_force(t, width, shift, count, mobile, strength):
    # Different material phases from the imager. The force transducer reads
    # the fixed support. It receives no original-link IDs or survival counts.
    y = np.zeros((count, len(PUSHES)))
    birth_offset = np.zeros_like(y)
    for j, (active, partner, _) in enumerate(frames(t, width, shift, count, .739, strength)):
        if j:
            dt = t[j] - t[j-1]
            previous_mean = y.mean(axis=0)
            if mobile:
                # Exact transverse motion over an interval with fixed partners:
                # dy/dt = -y + active*(h + birth_offset - y).
                rate = 1. + old_active
                eq = old_active[:, None] * (PUSHES + birth_offset) / rate[:, None]
                y += (eq-y) * (-np.expm1(-rate*dt))[:, None]
            changed = active & (~old_active | (partner != old_partner))
            birth_offset[changed] = y[changed] - PUSHES
        old_active, old_partner = active, partner
    link_reading = ((PUSHES+birth_offset-y)*active[:, None]).mean(axis=0)
    # For a compliant support, read the elastic and viscous reaction on the
    # opposite face using a backward difference of its mean feature motion.
    if mobile and len(t) > 1:
        support = y.mean(axis=0) + (y.mean(axis=0)-previous_mean)/dt
    else:
        support = link_reading.copy()
    return support, link_reading, int(active.sum())


def normalized(force, occupied_fraction):
    if force[0] != 0.:
        raise ArithmeticError("Unforced transverse reaction")
    return force[1:] / (PUSHES[1:] * occupied_fraction)


def flat_calibration(strength):
    # Three readings in a stationary flat patch, made before the throat tests.
    times = [.5, 1.]
    flat = np.full(2048, .75)
    _, initial, occupied = read_force(np.array([0.]), flat[:1], np.zeros(1), 256, True, strength)
    values = [normalized(initial, occupied/256)[0]]
    for stop in times:
        t = np.linspace(0., stop, len(flat))
        force, _, occupied = read_force(t, flat, np.zeros_like(flat), 256, True, strength)
        values.append(normalized(force, occupied/256)[0])
    ratio = (values[2]-values[1])/(values[1]-values[0])
    tau = -times[0]/np.log(ratio)
    plateau = (values[1]-ratio*values[0])/(1-ratio)
    emit(test="stationary_patch", times=[0.] + times, normalized_forces=values,
         inferred_tau=tau, inferred_plateau=plateau)
    return tau, plateau


def renewal_control(tau, plateau):
    # Translate one carrier-mounted row out and back at constant order. Keep
    # the transverse step applied throughout. Only registration is cycled.
    for steps, count in ((512, 2048), (2048, 8192)):
        t = np.linspace(0., 1.5, steps+1)
        width = np.full_like(t, .25)
        shift = np.interp(t, [0., .5, .75, 1., 1.5], [0., 0., 1., 0., 0.])
        imaged = image_crossing(t, width, shift, count, 0.)
        for mobile in (False, True):
            force, link, occupied = read_force(t, width, shift, count, mobile, 0.)
            gain = normalized(force, occupied/count)
            predicted = imaged['survival'] * (plateau+(1-plateau)*np.exp(-t[-1]/tau)) if mobile else imaged['survival']
            emit(test="flat_return", sites=count, steps=steps, mobile=mobile,
                 imaged=imaged, predicted=predicted, normalized_gains=gain.tolist(),
                 support_minus_link_gain=max(abs(normalized(force-link, occupied/count))))


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matched-width-control', action='store_true',
                        help='carry registry with the flow, preserving the measured order fraction')
    control = parser.parse_args().matched_width_control
    strength = 0. if control else REGISTRY_K
    tau, plateau = flat_calibration(strength)
    if control:
        renewal_control(tau, plateau)
    records = json.loads(PRIOR.read_text())["runs"]["mouth_mixed_order_ripple02"]["stdout_lines"]
    for line in records:
        case = json.loads(line)
        if case["test"] != "archive_patch" or case["mesh"] != 48:
            continue
        path_time, r, w, order = trace_path(case, 2000)
        carrier = Carrier(Path(case["carrier_file"]).parent)
        shift, _ = follow_cell(carrier, path_time, np.array([r[0], w[0]]), .005)
        for count, steps in ((2048, 256), (8192, 1024)):
            t = np.linspace(0., path_time[-1], steps+1)
            width = np.interp(t, path_time, order)
            registration = np.interp(t, path_time, shift)
            imaged = image_crossing(t, width, registration, count, strength)
            prediction = imaged["survival"]
            calibrated = prediction * (plateau + (1-plateau)*np.exp(-t[-1]/tau))
            for mobile in (False, True):
                force, link, occupied = read_force(t, width, registration, count, mobile, strength)
                gain = normalized(force, occupied/count)
                estimate = calibrated if mobile else prediction
                emit(test="crossing", profile=case["profile"], state=case["state"],
                     intake=case["intake"], carrier=case["carrier_file"], order=case["order_file"],
                     registry_stiffness=strength, traversal=t[-1], steps=steps, mobile=mobile, imaged=imaged,
                     pushes=PUSHES.tolist(), force=force.tolist(), link=link.tolist(),
                     normalized_gains=gain.tolist(), readout_occupied=occupied,
                     predicted_survival_only=prediction, predicted_flat_calibrated=estimate,
                     gain_minus_calibrated=gain[0]-estimate,
                     support_minus_link_gain=max(abs(normalized(force-link, occupied/count))))
                if mobile and max(abs(normalized(force-link, occupied/count))) > .015:
                    raise ArithmeticError("Opposite-face force imbalance on weak probe")
                if np.ptp(gain) > 1.e-8:
                    raise ArithmeticError("Probe sign/strength disagreement")


if __name__ == "__main__":
    run()
