"""Test open, reversible shear response against held-mouth gauge histories.

The input is the measured material velocity on the near gauge of a source-free
sublayer. The output is its signed traction on the receiving side. The trial
extends the order of a single transiting wall parcel uniformly behind that
plane and has no incoming wave. Both assumptions concern the receiver.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from mouth_contact_tracking import WallTrace
from neutral_patch_conversion import LEVELS, clock
from neutral_patch_pulse_inverse import read_pulse


RECORD = Path('research/pde_ledger_v3/_measurements/neutral_patch_shear_pulse.json')
COUNTS = (1201, 2401, 4801)


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def report(**data):
    print(json.dumps(data, allow_nan=False), flush=True)


def wave_stress(time, velocity, order, speed_squared):
    """Traction / Z1 from X_tt = f(t) X_yy with y = distance / c1.

    Initial displacement and velocity vanish. The remote boundary lies beyond
    the largest physical characteristic distance for the duration of the run.
    The spatial coordinate increases into the receiver. Positive traction
    transfers positive sideways momentum into it.
    """
    dt = time[1] - time[0]
    dy = dt / 0.8
    if np.max(speed_squared) > 1 + 1e-10 or np.min(speed_squared) <= 0:
        raise ValueError('Speed outside the tested characteristic cone')
    # For a surface reading only the round trip matters. Taking the entire
    # numerical stencil cone here also makes the initial far closure immaterial.
    shape = (velocity.shape[0], len(time) + 2)
    previous, current, following = (np.zeros(shape) for _ in range(3))
    drive = cumulative_trapezoid(velocity, time, axis=1, initial=0.)
    stress = np.zeros_like(velocity)
    current[:, 0] = drive[:, 1]
    for i in range(1, len(time)):
        stress[:, i] = order[i] * (3*current[:, 0] - 4*current[:, 1] + current[:, 2]) / (2*dy)
        if i == len(time) - 1:
            break
        # Only cells at most i from the port can be reached by this stencil.
        middle = slice(1, i+1)
        following[:, middle] = (2*current[:, middle] - previous[:, middle] +
                                (dt/dy)**2 * speed_squared[i] *
                                (current[:, :i] - 2*current[:, middle] + current[:, 2:i+2]))
        following[:, 0] = drive[:, i+1]
        previous, current, following = current, following, previous
    return stress


def gain(predictor, observed):
    return max(0., float(np.sum(predictor*observed)/np.sum(predictor**2)))


def rel_l2(predicted, observed):
    return float(np.linalg.norm(predicted-observed)/np.linalg.norm(observed))


def evaluate(folder_paths, time, chi, count, control):
    full, half, stresses, half_stresses = [], [], [], []
    for path in folder_paths:
        times, v, s, hv, hs, *_ = read_pulse(path, count)
        np.testing.assert_allclose(times, time, rtol=0, atol=1e-11)
        full.append(v[0])
        half.append(hv[0])
        stresses.append(s[0])
        half_stresses.append(hs[0])
    velocities = np.array(full + half)
    np.testing.assert_array_equal(velocities[:, 0], 0.)
    if control:
        # Constant order with fixed inertia has the exact outgoing stress
        # sqrt(chi) v in the same normalization.
        flat = np.full(count, .36)
        observed = np.sqrt(flat)*velocities[:1]
        numerical = wave_stress(time, velocities[:1], flat, flat)
        error = rel_l2(numerical, observed)
        report(test='flat_manufactured_trace', count=count, relative_l2=error)
        if error > .001:
            raise ArithmeticError('Constant-order propagation failed')
    return {
        'constant_participating_mass': wave_stress(time, velocities, chi, chi),
        'participating_mass_proportional_to_order': velocities*chi,
    }, np.asarray(stresses), np.asarray(half_stresses)


def main():
    hashes = json.loads(RECORD.read_text())['gauge_file_sha256']
    paths = [Path(key) for key in hashes]
    if not all(p.parent.joinpath('run_'+p.stem).is_dir() and checksum(p) == hashes[str(p)] for p in paths):
        raise ValueError('Missing or changed guarded gauge file')
    trace = WallTrace(48, .5)
    clock_time, depths, orders, _, marks = clock(trace, 2)
    if np.any(np.diff(orders) > 1e-9):
        raise ValueError('Order must decrease on this clock')
    report(test='material_clock', intake=trace.intake, clock_sources=trace.hashes,
           gauge_sources={str(p): hashes[str(p)] for p in paths},
           order=orders[marks].tolist(), time=clock_time[marks].tolist(), w=depths[marks].tolist(),
           constant_mass_speed_relative_to_entry=np.sqrt(np.asarray(LEVELS)/LEVELS[0]).tolist())
    prior = None
    for count in COUNTS:
        time = np.linspace(0, clock_time[-1], count)
        chi = np.interp(time, clock_time, orders)
        predictions, stress, half_stress = evaluate(paths, time, chi, count, True)
        if prior is not None:
            discrepancy = max(rel_l2(p[::2], q) for key in predictions
                              for p, q in zip(predictions[key], prior[key], strict=True))
            report(test='temporal_spatial_refinement', count=count, worst_relative_l2=discrepancy)
            if discrepancy > .001:
                raise ArithmeticError('Pulse refinement failed')
        prior = predictions
    train = time <= .65*time[-1]
    rows = []
    for law, predicted in predictions.items():
        full, half = np.split(predicted, 2)
        shared = gain(full[:, train], stress[:, train])
        for i, path in enumerate(paths):
            # Each individual fit gives the elastic model extra freedom.
            early_gain = gain(full[i, train], stress[i, train])
            optimal = gain(full[i], stress[i])
            rows.append(dict(law=law, source=path.stem, shared_gain=shared, own_early_gain=early_gain,
                             early_rel_l2=rel_l2(early_gain*full[i, train], stress[i, train]),
                             late_rel_l2=rel_l2(early_gain*full[i, ~train], stress[i, ~train]),
                             late_half_rel_l2=rel_l2(early_gain*half[i, ~train], half_stress[i, ~train]),
                             shared_late_rel_l2=rel_l2(shared*full[i, ~train], stress[i, ~train]),
                             minimum_whole_trace_rel_l2=rel_l2(optimal*full[i], stress[i])))
    report(test='signed_surface_stress', rows=rows)


if __name__ == '__main__':
    main()
