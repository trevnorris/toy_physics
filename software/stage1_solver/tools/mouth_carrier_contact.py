"""Try the saved four-dimensional carrier itself as the source of contact slip.

Two initially orthogonal material segments in the radial/depth plane are
tracked at common material times. The projection of one onto the other in
current repeat-length units gives a registry trajectory for the trial contact.
"""
import json
from pathlib import Path

import numpy as np

from mouth_contact_age import FREQUENCIES, PRIOR, emit, readings, trace_path
from mouth_free_contact import Carrier


def follow_cell(carrier, t, center, spacing):
    displacements = np.array([[0., 0.], [-.5, 0.], [.5, 0.], [0., -.5], [0., .5]])
    markers = center + spacing*displacements

    def velocity(points):
        r, w = points.T
        radius = carrier.radius(w)[0]
        if np.any((w <= .5) | (r <= 0) | (r >= radius)):
            raise ValueError("Marker left the conduit chart")
        _, _, _, ur, uw = carrier.eval(r/radius, w)
        return np.column_stack((ur, uw))

    def observables():
        across = markers[2]-markers[1]
        down = markers[4]-markers[3]
        return (across @ down)/(across @ across), markers[0].copy()

    initial, position = observables()
    shifts, positions = [initial], [position]
    for dt in np.diff(t):
        k1 = velocity(markers)
        k2 = velocity(markers+.5*dt*k1)
        k3 = velocity(markers+.5*dt*k2)
        k4 = velocity(markers+dt*k3)
        markers += dt*(k1+2*k2+2*k3+k4)/6
        shift, position = observables()
        shifts.append(shift)
        positions.append(position)
    return np.asarray(shifts)-initial, np.array(positions)


def overlap(q, shifts):
    lower = -q/2-shifts
    upper = q/2-shifts
    left = np.maximum.accumulate(lower[::-1])[::-1]
    right = np.minimum.accumulate(upper[::-1])[::-1]
    return np.maximum(right-left, 0.)/q[-1]


def count_forces(t, q, shifts, drive_nodes):
    # Keep individual partner identity through wraparound and re-entry. Count
    # contacts on a separate uniformly spaced set of material positions.
    phases = (np.arange(32768)+.3167)/32768
    references = np.zeros((len(phases), len(FREQUENCIES)), dtype=complex)
    strain = np.ones(len(FREQUENCIES), dtype=complex)
    for j in range(len(t)):
        coord = phases+shifts[j]
        partner = np.floor(coord+.5).astype(int)
        active = abs(coord-partner) < q[j]/2
        if j:
            replaced = active & (~old_active | (partner != old_partner))
            references[replaced] = strain
        if j in drive_nodes:
            strain = np.exp(-1j*t[j]*FREQUENCIES)
        old_active, old_partner = active, partner
    return np.mean(strain-references[active], axis=0)/strain


def run():
    data = json.loads(PRIOR.read_text())["runs"]["mouth_mixed_order_ripple02"]
    cases = [json.loads(line) for line in data["stdout_lines"]]
    for p in cases:
        if p["test"] != "archive_patch":
            continue
        t, r, w, q = trace_path(p, 2000)
        carrier = Carrier(Path(p["carrier_file"]).parent)
        drive_nodes = np.unique(np.linspace(0,len(t)-1,17).round().astype(int))
        drive = np.exp(-1j*t[drive_nodes, None]*FREQUENCIES)
        jumps = np.vstack([drive[0], np.diff(drive, axis=0)])
        for spacing in (.02, .01, .005):
            shifts, centers = follow_cell(carrier,t,np.array([r[0],w[0]]),spacing)
            distance = abs(shifts[:-1]-shifts[-1])
            finite = distance > 1e-12
            onset = min((q[:-1]-q[-1])[finite]/(2*distance[finite]))
            result = overlap(q,shifts)
            gain = result[drive_nodes] @ jumps/drive[-1]
            measured = count_forces(t,q,shifts,set(drive_nodes))
            center_error = np.hypot(centers[:,0]-r,centers[:,1]-w).max()
            emit(profile=p["profile"], state=p["state"], mesh=p["mesh"], spacing=spacing,
                 intake=p["intake"], shift_range=[shifts.min(), shifts.max()], exit_shift=shifts[-1],
                 full_survival_aspect_threshold=onset, equal_spacing_retention=result[0],
                 double_aspect_retention=overlap(q,2*shifts)[0], center_error=center_error,
                 omega=FREQUENCIES.tolist(), normalized_force=readings(gain),
                 direct_force=readings(measured), maximum_force_difference=max(abs(measured-gain)))
            if center_error > 3e-4 or max(abs(gain-measured)) > .004:
                raise ArithmeticError("Marker path or summed contact force mismatch")


if __name__ == "__main__":
    run()
