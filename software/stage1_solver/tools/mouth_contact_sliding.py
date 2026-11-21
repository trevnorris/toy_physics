"""Test a weak shear contact defined by overlap of sliding material features.

The relative tangential position is periodic and uniform in the unforced
ensemble. Its unwrapped velocity is constant over a parcel flight, with equal
populations of either sign. The occupied interval has width equal to the saved
local order. Leaving it breaks the material reference, even on later re-entry.
"""
import json

import numpy as np

from mouth_contact_age import FREQUENCIES, PRIOR, emit, readings, trace_path


def cutoffs(t, q):
    span = 2*(t[-1]-t[:-1])
    return np.min((q[:-1]-q[-1])/span), np.min((q[:-1]+q[-1])/span)


def continuum(t, q, speed):
    # Intersect all overlap intervals, on one unwrapped copy, backwards from
    # the final section. A full period traveled includes loss of the reference.
    missed_width = speed*(t[-1]-t)-(q-q[-1])/2
    largest = np.maximum.accumulate(missed_width[::-1])[::-1]
    kernel = np.maximum(0., 1-largest/q[-1])
    drive_nodes = np.unique(np.linspace(0, len(t)-1, 17).round().astype(int))
    drive = np.exp(-1j*t[drive_nodes, None]*FREQUENCIES)
    jumps = np.vstack([drive[0], np.diff(drive, axis=0)])
    return kernel[0], kernel[drive_nodes] @ jumps/drive[-1], drive_nodes


def follow_pairs(t, q, speed, drive_nodes, sites=16384):
    # Direct force reading from material pair extensions. Changes of partner
    # and losses of contact are detected on the actual periodic geometry.
    step_force = np.zeros(len(FREQUENCIES), dtype=complex)
    retained = occupied = 0
    probes = set(drive_nodes)
    offsets = (np.arange(sites)+.314159)/sites
    for direction in (-1, 1):
        old = np.zeros(sites, dtype=bool)
        last_partner = np.zeros(sites, dtype=int)
        reference = np.zeros((sites, len(FREQUENCIES)), dtype=complex)
        last_strain = np.ones(len(FREQUENCIES), dtype=complex)
        for j, (time, width) in enumerate(zip(t, q, strict=True)):
            displacement = offsets+direction*speed*time
            partner = np.floor(displacement+.5).astype(int)
            engaged = abs(displacement-partner) < width/2
            if j == 0:
                old = engaged.copy()
            else:
                continuous = last_engaged & engaged & (partner == last_partner)
                new = engaged & ~continuous
                reference[new] = last_strain
                old &= continuous
            if j in probes:
                last_strain = np.exp(-1j*time*FREQUENCIES)
            last_engaged, last_partner = engaged, partner
        occupied += engaged.sum()
        retained += old.sum()
        step_force += (last_strain-reference[engaged]).sum(axis=0)/last_strain
    return retained/occupied, step_force/occupied, occupied/(2*sites)


def run():
    data = json.loads(PRIOR.read_text())["runs"]["mouth_mixed_order_ripple02"]
    cases = [json.loads(s) for s in data["stdout_lines"]]
    for fields in cases:
        if fields["test"] != "archive_patch":
            continue
        t, r, w, q = trace_path(fields, 2000)
        v_full, v_zero = cutoffs(t, q)
        emit(test="cutoffs", profile=fields["profile"], state=fields["state"], mesh=fields["mesh"],
             intake=fields["intake"], q_entry=q[0], q_exit=q[-1], traversal=t[-1],
             maximum_speed_for_full_retention=v_full, minimum_speed_for_zero_retention=v_zero,
             carrier=fields["carrier_file"], carrier_sha256=fields["carrier_sha256"],
             order=fields["order_file"], order_sha256=fields["order_sha256"])
        if fields["mesh"] != 48:
            continue
        for speed in (0., .25, .5):
            fraction, gain, nodes = continuum(t, q, speed)
            sampled, force, active = follow_pairs(t, q, speed, nodes)
            emit(test="held_and_driven", profile=fields["profile"], state=fields["state"],
                 speed=speed, predicted_retention=fraction, counted_retention=sampled,
                 observed_order=active, omega=FREQUENCIES.tolist(), drive_times=t[nodes].tolist(),
                 predicted_gain=readings(gain), direct_force=readings(force),
                 max_force_difference=max(abs(gain-force)))
            if max(abs(gain-force)) > .004 or abs(sampled-fraction) > .001:
                raise ArithmeticError("Contact force disagrees with uninterrupted overlap")


if __name__ == "__main__":
    run()
