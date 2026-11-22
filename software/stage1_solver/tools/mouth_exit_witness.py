"""Signed virtual work in the bonded positive-radius chart at a saved exit-release endpoint."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_live_waist_reopen import reopen
from mouth_poly_shallow_release import CUT, ROOT, Trial, emit, source
from mouth_sealed_shape_relax import Meridian
from mouth_waist_outer_order import exterior


def sha(data):
    return hashlib.sha256(data).hexdigest()


def inputs(path, radius, mesh, trial):
    guard = json.loads((path.parent / 'outcome.json').read_text())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError('Endpoint guard or child failed')
    digest = sha(path.read_bytes())
    rows = map(json.loads, (path.parent / 'stdout').read_text().splitlines())
    parent, = [r for r in rows if r.get('kind') == 'endpoint'
               and digest in (r.get('snapshot_sha256'), r.get('sha256'))]
    with np.load(path, allow_pickle=False) as f:
        if (f['radius'].item(), f['mesh'].item()) != (radius, mesh):
            raise ValueError('Snapshot radius or mesh mismatch')
        if f['inputs'].tolist() != trial.fingerprints:
            raise ValueError('Snapshot uses a different material preparation')
        shifts, chi, seed = f['shifts'], f['chi'], (complex(f['pole'].item()), f['vector'])
    if (shifts.shape != (6,) or chi.shape != trial.chi.shape or seed[1].shape != trial.seed[1].shape
            or not all(np.isfinite(x).all() for x in (shifts, chi, *seed))
            or chi.min() < 0 or chi.max() > 1):
        raise ValueError('Malformed endpoint')
    previous = ROOT / 'research/pde_ledger_v3/_measurements/mouth_poly_shallow_release.json'
    archived = json.loads(previous.read_text())['runs'][f'balance_r{radius}n{mesh}']
    raw = archived['literal_stdout']
    if archived['guard_and_child_exit'] or sha(raw.encode()) != archived['source_stdout_sha256']:
        raise ValueError('Unverified previous equilibrium')
    rows = list(map(json.loads, raw.splitlines()))
    if rows[0]['inputs'] != trial.fingerprints or rows[-1]['kind'] != 'balance':
        raise ValueError('Preparation changed since the three-coordinate balance')
    return shifts, chi, seed, digest, parent['status'], rows[-1], archived['source_stdout_sha256']


def validate_fields(info):
    # Both CoupledLoads and Trial supply the current and pole tests under these names.
    outflow = info['outflow'] if 'outflow' in info else info['terminal_current']
    if (info['mode_overlap'] < .94 or info['pole_residual'] > 1e-7
            or abs(info['current']-outflow) > 1e-6*max(abs(info['current']), 1e-12)):
        raise ArithmeticError('Current or continued occupied mode unresolved')


def shape_work(trial, shifts, chi, seed, center):
    # The receiver is minimized at each shape. Its separate order and elastic energy
    # changes include redistribution. For elastic virtual work hold its mapped nodal
    # order at the central value and differentiate the integrated reference density.
    shape = Meridian(trial.radius, trial.x + shifts[:4])
    spacing = trial.cr.copy()
    spacing[0] += shifts[5]

    def receiver_load(w, x):
        return source(shape, spacing, w/shape.length, x, True)

    ext, (w, x, side_chi) = exterior(shape, trial.n, CUT, return_side_field=True,
                                      side_order_work=receiver_load)
    if abs(ext['static_work'] - center['energy']['static_outside']) > 2e-7:
        raise ArithmeticError('Failed to recover receiving order for elastic virtual work')
    t = w/shape.length
    for axis, label in ((1, 'log_exit_radius'), (3, 'upper_minus_lower'),
                        (4, 'inner_log_spacing'), (5, 'outer_log_spacing')):
        for h in (.001, .0005):
            probes = []
            for sign in (-1, 1):
                shifted = shifts.copy()
                shifted[axis] += sign*h
                other = Meridian(trial.radius, trial.x + shifted[:4])
                try:
                    if other.wall()['minimum_radius'] <= 0:
                        raise ValueError('Positive-radius chart ended')
                    _, _, _, _, info, _ = trial.evaluate(shifted, chi, seed)
                    validate_fields(info)
                    if axis == 1:
                        coeff = source(other, spacing, t, x, True)
                        info['fixed_order_elastic'] = (info['energy']['static_inside']
                                                      + np.sum(coeff * side_chi))
                    probes.append(info)
                except (ArithmeticError, RuntimeError, ValueError) as err:
                    emit(kind='shape_work_unresolved', axis=label, offset=sign*h, reason=str(err))
                    probes.append(None)
            if any(info is None for info in probes):
                continue
            minus, plus = probes
            slopes = {k: (plus['energy'][k]-minus['energy'][k])/(2*h) for k in center['energy']}
            changes = [sum(v-center['energy'][k] for k, v in p['energy'].items()) for p in probes]
            extra = {}
            if axis == 1:
                elastic = (plus['fixed_order_elastic']-minus['fixed_order_elastic'])/(2*h)
                extra = dict(fixed_order_elastic=elastic,
                             order_by_envelope=sum(slopes[k] for k in (
                                 'inner_order', 'outer_order', 'static_inside', 'static_outside'))-elastic)
            emit(kind='shape_work', axis=label, h=h, sector_slopes=slopes, total_slope=sum(slopes.values()),
                 changes_minus_plus=changes, curvature=sum(changes)/h**2,
                 minimum_overlap=min(p['mode_overlap'] for p in probes), **extra)


def order_work(trial, chi, grad, owner, loads, info, elastic):
    wave = loads.wave
    omega, v = wave.seed
    density = loads.carrier.f.copy()
    base = info['energy']

    def change(direction, h, vstart):
        moved = chi + h*direction
        if np.min(moved) < 0 or np.max(moved) > 1 or np.any(moved[~owner.free] != chi[~owner.free]):
            raise ValueError('Order variation crosses an imposed value')
        wave.seed, wave.target = (omega, vstart.copy()), omega.real
        loads.carrier.f = density.copy()
        _, _, readout = loads.evaluate(moved, 1., BIAS)
        validate_fields(readout)
        return (owner.energy_gradient(moved)[0]-base['inner_order']
                + np.sum(elastic*(moved-chi))
                + (readout['carrier_energy']-base['two_sided_carrier'])
                + (readout['wave_energy']-base['strained_wave']))

    def sample(direction, steps, vstart):
        samples = []
        for h in steps:
            try:
                samples.append(dict(offset=h, change=change(direction, h, vstart)))
            except (ArithmeticError, RuntimeError, ValueError) as err:
                samples.append(dict(offset=h, reason=str(err)))
        return samples

    # This two-sided, fixed-width perturbation has the same support profile on both
    # grids. The chi(1-chi) factor preserves the active set in either direction.
    t = np.linspace(0, 1, chi.shape[0])[:, None]
    x = np.linspace(0, 1, chi.shape[1])[None, :]
    patch = chi*(1-chi)*np.exp(-((t-.7)**2+(x-.55)**2)/(2*.12**2))
    patch[~owner.free] = 0.
    emit(kind='local_order_patch', center_t_x=[.7, .55], sigma=.12, slope=np.sum(grad*patch),
         max_direction=float(np.max(patch)), samples=sample(patch, (-.004, -.002, .002, .004), v))

    # An independent nodal direction measures a true finite change of the same
    # functional. Restrict the signed test to nodes with room on both sides.
    h = .001
    interior = owner.free & (chi > h) & (chi < 1-h)
    if interior.any():
        index = np.unravel_index(np.argmax(np.where(interior, abs(grad), -1.)), chi.shape)
        direction = np.zeros_like(chi)
        direction[index] = 1.
        emit(kind='local_order', tested_for_selection=int(np.sum(interior)), index=list(map(int, index)),
             fractions=[index[0]/(2*trial.n), index[1]/trial.n], order=chi[index], slope=grad[index],
             samples=sample(direction, (-h, -h/2, h/2, h), v))
    else:
        emit(kind='local_order', tested_for_selection=0)

    upper = owner.free & (chi == 1)
    if upper.any():
        i, j = np.unravel_index(np.argmax(np.where(upper, grad, -np.inf)), chi.shape)
        direction = np.zeros_like(chi)
        direction[i, j] = -1.
        emit(kind='upper_bound_order', nodes=int(upper.sum()), index=[int(i), int(j)],
             melt_slope=-grad[i, j], samples=sample(direction, (.001, .0005), v))
    else:
        emit(kind='upper_bound_order', nodes=0)

    zeros = np.argwhere(owner.free & (chi == 0))
    if not len(zeros):
        emit(kind='reopened_order', zeros=0, resolved=0, unresolved=0)
        return
    project = wave.project.tocsc()
    weights = project @ chi.T.ravel()
    mass = wave.cut_m + wave.matrix(weights[:, None, None]*wave.mc)
    active = mass.diagonal() > 0
    if np.any(v[~active] != 0):
        raise ArithmeticError('Saved central mode extends outside its active domain')
    edge = v[wave.port.ids]
    port_derivative = (wave.port.block(omega+2e-6)-wave.port.block(omega-2e-6))/4e-6
    denominator = 2*omega*(v@(mass@v)) + edge@(port_derivative@edge)
    cell_form = wave.kc - omega**2*wave.mc
    negative = resolved = activated = 0
    best = None
    unresolved = []
    for i, j in zeros:
        col = j*chi.shape[0]+i
        cells = slice(project.indptr[col], project.indptr[col+1])
        try:
            _, correction, ids, extension, cond = reopen(
                project.indices[cells], project.data[cells], wave, cell_form, v, active)
            # grad already contains the direct bilinear wave term at the zero node.
            slope = grad[i, j] + (correction/denominator).real
            if not np.isfinite(slope) or not np.isfinite(cond):
                raise ArithmeticError('Nonfinite one-sided order slope')
        except (ArithmeticError, np.linalg.LinAlgError) as err:
            unresolved.append(dict(index=[int(i), int(j)], reason=str(err)))
            continue
        resolved += 1
        activated += int(bool(len(ids)))
        negative += int(slope < -1e-6)
        if best is None or slope < best[0]:
            best = (slope, i, j, ids, extension, cond)
    details = {}
    if best is not None:
        slope, i, j, ids, extension, cond = best
        extended = v.copy()
        extended[ids] = extension
        direction = np.zeros_like(chi)
        direction[i, j] = 1.
        details = dict(min_slope=slope, index=[int(i), int(j)], fractions=[i/(2*trial.n), j/trial.n],
                       extension_nodes=len(ids), extension_condition=cond,
                       samples=sample(direction, (.001, .0005, .00025), extended))
    emit(kind='reopened_order', zeros=len(zeros), resolved=resolved, unresolved=len(unresolved),
         unresolved_sites=unresolved, reopening_wave_nodes=activated, below_minus_1e_minus_6=negative,
         **details)


def run(radius, mesh, path):
    trial = Trial(radius, mesh, bonded_core=True)
    shifts, chi, seed, digest, status, start, start_sha = inputs(path, radius, mesh, trial)
    _, grad, owner, loads, info, elastic = trial.evaluate(shifts, chi, seed)
    validate_fields(info)
    if (abs(loads.wave.seed[0]-seed[0]) > 2e-7 or
            not np.array_equal(chi[~owner.free], trial.chi[~owner.free])):
        raise ArithmeticError('Could not reproduce the saved pole or order boundary')
    entrance = Meridian(radius, trial.x+shifts[:4])
    emit(kind='endpoint', radius=radius, mesh=mesh, parent_status=status, snapshot_sha256=digest,
         previous_balance_stdout_sha256=start_sha, exit_radius=entrance.exit,
         previous_cavity_fraction=start['cavity_fraction'],
         cavity_fraction_change=info['cavity_fraction']-start['cavity_fraction'],
         cavity_fraction_ratio=info['cavity_fraction']/start['cavity_fraction'],
         previous_wave_energy=start['energy']['strained_wave'],
         fixed_action_energy_ratio=info['energy']['strained_wave']/start['energy']['strained_wave'],
         **info)
    shape_work(trial, shifts, chi, loads.wave.seed, info)
    order_work(trial, chi, grad, owner, loads, info, elastic)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius', type=int, choices=(1, 2), required=True)
    parser.add_argument('--mesh', type=int, choices=(48, 72), required=True)
    parser.add_argument('--snapshot', type=Path, required=True)
    args = parser.parse_args()
    run(args.radius, args.mesh, args.snapshot)
