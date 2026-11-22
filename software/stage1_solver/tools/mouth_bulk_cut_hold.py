"""Hold a tangential population displacement for one transit in a four-periodic cell."""
import argparse
import ctypes as ct
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace

import numpy as np

from mouth_bulk_cut_grip import REFERENCE, STEP, load_state, sha
from mouth_bulk_isotropic_4d import radial
from mouth_bulk_relaxed_shape import PeriodicMaterial
from mouth_bulk_temperature_steps import Contacts, FRAMES, TRANSIT
from mouth_bulk_thermal_cell import ThermalCell
from neutral_patch_sites import emit


FROZEN = REFERENCE.with_name('mouth_bulk_cut_grip.json')
OFFSETS = dict(zero=0., plus=STEP, minus=-STEP, half_plus=STEP/2, half_minus=-STEP/2)


def library(work, archived_kernel_sha):
    """Build this integrator with the same source bytes used by the frozen-cut experiment."""
    work.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_suffix('.cpp')
    base = source.with_name('mouth_bulk_thermal_cell.cpp')
    if sha(base) != archived_kernel_sha:
        raise ValueError('Thermal force kernel differs from the frozen-cut reference')
    digest = hashlib.sha256(source.read_bytes()+base.read_bytes()).hexdigest()
    dest = (work/f'cut_hold_{digest}.so').resolve()
    if not dest.exists():
        with tempfile.TemporaryDirectory(dir=work) as tmp:
            compiled = Path(tmp)/'cut.so'
            subprocess.run(['g++', '-std=c++17', '-O3', '-march=native', '-fPIC', '-shared',
                            str(source), '-o', str(compiled)], check=True)
            os.replace(compiled, dest)
    lib = ct.CDLL(str(dest))
    arr = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    lib.thermal_make.argtypes = [ct.c_int, arr, arr, ct.c_double, arr]
    lib.thermal_make.restype = ct.c_void_p
    lib.thermal_free.argtypes, lib.thermal_free.restype = [ct.c_void_p], None
    lib.thermal_read.argtypes, lib.thermal_read.restype = [ct.c_void_p, arr, arr, arr], ct.c_double
    lib.cut_advance.argtypes = [ct.c_void_p, ct.c_int, ct.c_double, arr, arr, arr, ct.c_double,
                                ct.c_double, arr, arr]
    lib.cut_advance.restype = ct.c_int
    return lib, digest


def frozen_reading(record, source, history):
    rows = [r for line in record['literal_stdout'] if (r := json.loads(line))['state'] == str(source)]
    if len(rows) != 1 or rows[0]['sha256'] != sha(source) or rows[0]['history'] != history:
        raise ValueError('Archive has no matching preparation for the initial force comparison')
    return rows[0]


class Cut:
    def __init__(self, start, box, cutoff):
        self.start, self.box, self.cutoff = start, box, cutoff
        self.inverse = np.linalg.inv(box)
        self.normal = self.inverse[:, 0]/np.linalg.norm(self.inverse[:, 0])
        self.tangent = np.ascontiguousarray(box[1]/np.linalg.norm(box[1]))
        self.area = np.linalg.det(box)*np.linalg.norm(self.inverse[:, 0])
        np.testing.assert_allclose(self.normal@self.tangent, 0., atol=1e-13)
        np.testing.assert_allclose(self.area, np.sqrt(np.linalg.det(box[1:]@box[1:].T)), rtol=1e-13)
        self.labels = ((start@self.inverse)[:, 0] % 1. >= .5).astype(int)
        if len(np.unique(self.labels)) != 2:
            raise ValueError('Need two populated material halves')
        self.d = self.labels-self.labels.mean()
        self.w = self.d@self.d
        self.metric = PeriodicMaterial(len(start), cutoff)
        self.ever_crossed = np.zeros(len(start), bool)

    def coordinate(self, pos):
        # Use group means independently of the native weighted projection.
        s = pos@self.tangent
        return s[self.labels==1].mean()-s[self.labels==0].mean()

    def at_frame(self, pos):
        center_drift = (pos-self.start).mean(axis=0)
        origin = center_drift@self.inverse[:, 0]
        y = pos@self.inverse
        u = y[:, 0]-origin
        spatial = (u % 1. >= .5).astype(int)
        crossed = self.labels != spatial
        self.ever_crossed |= crossed
        i, j, image = self.metric.pairs(y, self.box)
        bonds = (y[j]-y[i]+image)@self.box
        radius = np.linalg.norm(bonds, axis=1)
        _, du, _ = radial(radius, self.cutoff)
        tangential = (du/radius)*(bonds@self.tangent)
        tagged = tangential@(self.labels[j]-self.labels[i])/(2*self.area)
        # Both periodically repeated planes, with opposing orientations. An
        # image bond that crosses both contributes to each with opposite sign.
        v = u[j]+image[:, 0]
        crossings = np.stack((np.floor(u[i])-np.floor(v), np.floor(v-.5)-np.floor(u[i]-.5)))
        np.testing.assert_array_equal(crossings.sum(axis=0), spatial[j]-spatial[i])
        face_loads = crossings@tangential/self.area
        return dict(pair_tagged_traction=float(tagged), plane_face_tractions=face_loads.tolist(),
                    plane_mean_traction=float(face_loads.mean()), changed_side=int(crossed.sum()),
                    ever_changed_side_at_frames=int(self.ever_crossed.sum()), upper_in_plane=int(spatial.sum()),
                    plane_origin_fraction=float(origin), rigid_drift=center_drift.tolist())


def check_projection(lib, x, box, cutoff, tangent, dt, temp):
    # An unequal partition stresses both the relative-coordinate and zero-total-force factors.
    # Use the analytic Python potential and separately take the two ordinary group means.
    ids = np.arange(len(x)) % 3 == 0
    counts = np.array([int(ids.sum()), int((~ids).sum())])
    d = ids.astype(float)-ids.mean()
    w = counts.prod()/len(d)
    noise = np.random.default_rng(7547).normal(size=(1, len(x), 4))*np.sqrt(2*temp*dt)
    oracle = PeriodicMaterial(len(x), cutoff)
    _, grad, _ = oracle.read(x@np.linalg.inv(box), box)
    free = x+noise[0]-dt*grad

    def relative(pos):
        return (pos[ids].mean(axis=0)-pos[~ids].mean(axis=0))@tangent

    target = relative(x)
    correction = target-relative(free)
    projected = free.copy()
    projected[ids] += correction*counts[1]/len(x)*tangent
    projected[~ids] -= correction*counts[0]/len(x)*tangent
    force = (grad[ids]@tangent).sum()  # Internal gradients sum to zero.
    stochastic = -w*relative(noise[0])/dt
    reference = np.array([w*correction/dt, stochastic, force, force])
    actual, stats = x.copy(), np.empty(8)
    mat = ThermalCell(lib, box, len(x), cutoff)
    status = lib.cut_advance(mat.handle, 1, dt, noise, d, tangent, w, target, actual, stats)
    if status:
        raise ArithmeticError(f'One-step validation failed, status {status}')
    np.testing.assert_allclose(actual, projected, atol=2e-12, rtol=0)
    np.testing.assert_allclose(stats[:4], reference, atol=3e-8, rtol=0)
    np.testing.assert_allclose((actual-free).sum(axis=0), 0., atol=2e-12)
    np.testing.assert_allclose(relative(actual), target, atol=2e-13, rtol=0)
    emit(test='projection_check', group_counts=counts.tolist(), max_position_error=float(abs(actual-projected).max()),
         actual_generalized_force=float(stats[0]), reference_generalized_force=float(reference[0]),
         thermal_impulse=float(dt*stats[1]), reference_thermal_impulse=float(dt*stochastic),
         conditional_force=float(stats[2]), analytic_tagged_force=float(force))


def advance_frame(mat, x, cut, target, rng, temp, dt, steps):
    # Same sequence, batching, dt and temperature for all five signed branches.
    # Native step sums include every pre-step force throughout the interval.
    sums, peak = np.zeros(4), np.zeros(2)
    noise_square = 0.
    first = None
    for offset in range(0, steps, 128):
        count = min(128, steps-offset)
        noise = rng.normal(size=(count, len(x), 4))*np.sqrt(2*temp*dt)
        stats = np.empty(8)
        failed = mat.lib.cut_advance(mat.handle, count, dt, noise, cut.d, cut.tangent,
                                     cut.w, target, x, stats)
        if failed:
            raise ArithmeticError(f'Constraint integration failed at batch offset {offset}, status {failed}')
        if first is None:
            first = stats[6]
        sums += stats[:4]
        peak = np.maximum(peak, stats[4:6])
        noise_square += stats[7]
    scale = 2*cut.area
    applied, thermal, conditional, gradient = sums/(steps*scale)
    np.testing.assert_allclose(peak[0]/scale, 0., atol=5e-8)
    np.testing.assert_allclose(peak[1], 0., atol=1e-11)
    np.testing.assert_allclose([conditional, applied-thermal], gradient, atol=5e-9, rtol=0)
    return dict(mean_actuator_traction=float(applied), mean_thermal_traction=float(thermal),
                mean_conditional_traction=float(conditional), mean_site_gradient_traction=float(gradient),
                rms_thermal_traction=float(np.sqrt(noise_square/steps)/scale),
                thermal_impulse_per_two_areas=float(thermal*steps*dt),
                max_force_balance_error_per_two_areas=float(peak[0]/scale),
                max_centroid_error=float(peak[1]), first_conditional_traction=float(first/scale))


def run(args):
    record = json.loads(FROZEN.read_text())
    prior = json.loads(REFERENCE.read_text())
    source, start, box, setup, end = load_state(prior, args.rung, args.history, 'prepare')
    frozen = frozen_reading(record, source, args.history)
    core_path = 'software/stage1_solver/tools/mouth_bulk_thermal_cell.cpp'
    lib, kernel_sha = library(args.work, record['program_and_source_pins'][core_path])
    temp, dt, duration, cutoff = (setup[k] for k in ('temperature', 'dt', 'transit', 'cutoff'))
    steps = round(duration/(FRAMES*dt))
    np.testing.assert_allclose([FRAMES*steps*dt, duration], TRANSIT, rtol=1e-13)
    cut = Cut(start, box, cutoff)
    q = OFFSETS[args.branch]
    target = cut.coordinate(start)+q
    x = np.ascontiguousarray(start+q*cut.d[:, None]*cut.tangent)
    np.testing.assert_allclose((x-start).mean(axis=0), 0., atol=1e-14)
    if args.check_only:
        check_projection(lib, x, box, cutoff, cut.tangent, dt, temp)
        return
    mat = ThermalCell(lib, box, len(x), cutoff)
    contacts = Contacts(SimpleNamespace(pos=start, box=box))
    sample = SimpleNamespace(pos=x, box=box)
    cross_labels = cut.labels[:, None] != cut.labels
    initial_cross = contacts.initial & cross_labels
    density = len(x)/mat.volume
    seed = [31049, args.history, args.rung]
    rng = np.random.default_rng(seed)
    pin = sha(source)
    emit(test='setup', branch=args.branch, rung=args.rung, history=args.history, source=str(source),
         source_sha256=pin, frozen_sha256=sha(FROZEN), kernel_sha256=kernel_sha, temperature=temp,
         transit=duration, dt=dt, mobility=1., steps_per_frame=steps, frames=FRAMES, cutoff=cutoff,
         q=q, baseline_coordinate=float(target-q), target=float(target), fixed_box=box.tolist(),
         fixed_density=float(density), sites=len(x), upper_count=int(cut.labels.sum()),
         normal=cut.normal.tolist(), tangent=cut.tangent.tolist(), area_per_face=float(cut.area),
         noise_seed=seed, total_initial_cross_contacts=int(initial_cross.sum()),
         preparation_end_pressure=end['readings'][0]['pressure'])
    expected = frozen['spatial']['baseline']
    if q:
        level = int(abs(q)<STEP)
        expected += frozen['spatial']['signed_force_increments'][level][int(q<0)]
    snapshots = []
    previous = None
    for frame in range(FRAMES+1):
        interval = None
        if frame:
            interval = advance_frame(mat, x, cut, target, rng, temp, dt, steps)
            # Measured constraint reaction at the first step of this interval
            # versus an analytic all-image force sum at the preceding snapshot.
            np.testing.assert_allclose(interval['first_conditional_traction'], previous, atol=5e-8, rtol=0)
        snapshots.append(x.copy())
        energy, grad, virial = mat.read(x)
        direct = cut.at_frame(x)
        previous = direct['pair_tagged_traction']
        native = np.sum(cut.d*(grad@cut.tangent))/(2*cut.area)
        np.testing.assert_allclose(native, previous, atol=5e-9, rtol=1e-11)
        centroid = cut.coordinate(x)
        np.testing.assert_allclose(centroid, target, rtol=0, atol=1e-11)
        trace_pressure = -np.trace(virial)/(4*mat.volume)
        pressure = density*temp+trace_pressure  # Identical material diagnostic to the preparation.
        if frame == 0:
            np.testing.assert_allclose(native, expected, atol=5e-9, rtol=0)
            if q == 0:
                np.testing.assert_allclose([energy/len(x), pressure],
                                           [end['readings'][0][k] for k in ('energy_per_site','pressure')], atol=5e-9)
        edges = contacts.read(sample)
        emit(test='frame', frame=frame, elapsed=frame*steps*dt, interval=interval, centroid=float(centroid),
             energy_per_site=float(energy/len(x)), material_pressure=float(pressure),
             pair_virial_pressure=float(trace_pressure), ideal_pressure=float(density*temp),
             projected_thermal_pressure=float((len(x)-.25)*temp/mat.volume),
             site_gradient_traction=float(native), initial_cross_contacts_missing=int((initial_cross&~contacts.current).sum()),
             new_cross_contacts=int((contacts.current&cross_labels&~initial_cross).sum()), **direct, **edges)
    np.savez(args.work/'frames.npz', pos=np.asarray(snapshots), box=box, labels=cut.labels, start=start,
             tangent=cut.tangent, target=target, q=q, dt=dt, frame_steps=steps, temperature=temp,
             cutoff=cutoff, source_sha256=pin, branch=args.branch, rung=args.rung, history=args.history)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('branch', choices=OFFSETS)
    parser.add_argument('--rung', required=True, choices=(0, 4), type=int)
    parser.add_argument('--history', required=True, choices=(0, 1, 2), type=int)
    parser.add_argument('--work', required=True, type=Path)
    parser.add_argument('--check-only', action='store_true', help='Check one unequal-population constraint step and exit')
    run(parser.parse_args())
