"""Frozen comparison of the two determinant stencils and one-sided order work at R4 step 108."""
import argparse
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np

from mouth_r4_aligned_stencil import barrier_means
from mouth_r4_newton_stop import PIN_SHA, read_state, sha


def force_history(study):
    journal = study/'primary_guard/stdout'
    work = study.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_newton_stop_work.json'
    expected = json.loads(work.read_bytes())['source_sha256']['primary_guard/stdout']
    if sha(journal) != expected:
        raise ValueError('Changed accepted force history')
    history = []
    for line in journal.read_text().splitlines():
        event = json.loads(line)
        if event.get('outer', 0) >= 90:
            history.append({k:event[k] for k in ('kind','outer','max_force','active_force','alpha','field_residual')
                            if k in event})
    return history


def stencil_change(gauge, labels, dq):
    """Return new minus old unit gauge energy at fixed label displacement on the moving grid."""
    from mouth_cp10_label_gauge import depth_grid, profile
    from mouth_r4_smooth_static import BITS, RULE
    dtype = np.result_type(dq, np.longdouble)
    L = np.longdouble(gauge.shape.length) * np.exp(dq[0])
    axes = [np.asarray(a, dtype=dtype) for a in gauge.axes]
    axes[1] = depth_grid(L, gauge.reference)
    lo = np.column_stack([axis[gauge.cells[:,k]] for k, axis in enumerate(axes)])
    widths = np.column_stack([np.diff(axis)[gauge.cells[:,k]] for k, axis in enumerate(axes)])
    off = labels.reshape(-1,3)[gauge.ids].astype(dtype) - gauge.base.reshape(-1,3)[gauge.ids]
    points = off + lo[:,None,:] + BITS*widths[:,None,:]
    old, aligned = barrier_means(points, widths)
    measure = np.zeros(len(lo), dtype=dtype)
    for q in RULE:
        x, t, angle = (lo + widths*q).T
        R, _ = profile(gauge.shape, t, dq)
        y = (t-.06)/.56
        bump = np.where((y.real>0)&(y.real<1), 16*y*y*(1-y)**2, 0)
        perturbation = bump * (gauge.polar + dq[5]) * (3*np.cos(angle)**2-1)/2
        radius = x*R+x**3*perturbation
        radial_derivative = R+3*x*x*perturbation
        phi = np.prod(np.where(BITS, q, 1-np.array(q)), axis=1)
        chi = gauge.order.astype(dtype) @ phi.astype(dtype)
        weight = (1-chi)**2 * (1+2*chi)
        measure += .5*np.pi * widths.prod(axis=1)*radial_derivative*radius**2*L*np.sin(angle)*weight
    delta = aligned-old
    return np.sum(delta*measure), delta


def run(study, oracle_path):
    from mouth_r4_newton_stop import STATES
    pin = study/'pin/software/stage1_solver/tools'
    manifest = study/'pinned.json'
    if sha(manifest) != PIN_SHA or any(sha(pin/k) != v for k,v in json.loads(manifest.read_bytes()).items()):
        raise ValueError('Changed physical source')
    sys.path.insert(0, str(pin))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_cp10_spatial_crystal import MODES
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_slaved_melt import FLOW
    from mouth_r4_smooth_fields import GaugeGeometry
    from mouth_r4_smooth_held import SmoothHeld
    from mouth_r4_smooth_static import Gauge

    reference_path = study.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_newton_held_stop.json'
    reference = json.loads(reference_path.read_bytes())
    with np.load(oracle_path, allow_pickle=False) as archive:
        if str(archive['sha']) != STATES['last_active'][1]:
            raise ValueError('Wrong order energy sample')
        ids, physical_work = archive['ids'], archive['physical']
    with TemporaryDirectory(dir=study.parent) as tmp:
        driver = SmoothHeld(study.parent, Path(tmp), 1e-8)
        cp, fields, labels = read_state(study/'primary', 'last_active', driver, CHECKPOINT_SHA)
        shape = driver.shape(driver.x)
        trial = ShapeTrial(driver.base, driver.prep, driver.x, 8,24,6,3,2, precondition=False,
                           axes=axes_for_depth(shape.length, driver.reference_depth))
        trial.model.extended_energy_sum = True
        order = order_grid(trial, fields)
        solid = SpatialCrystal(trial.b)
        elastic, forces, source, _ = solid.evaluate(labels, order)
        load = trial.ext + np.bincount(trial.b.ids.ravel(), source.transpose(0,2,1).ravel(),
                                       minlength=trial.model.size)
        trial.model.material_load = load
        physical, _ = trial.model.evaluate(fields)
        gauge = Gauge(trial, driver.reference_depth, order)
        zero = np.zeros(6, dtype=np.longdouble)
        difference, per_cell = stencil_change(gauge, labels, zero)
        # Use the physical quadrature and joined order registry independently of the shape calculation.
        geometry = GaugeGeometry(trial)
        old_density = geometry.work_density(labels).astype(np.longdouble)
        density = old_density + np.concatenate([volume * per_cell for volume in geometry.volume])
        if np.any(density < 0):
            raise ArithmeticError('Negative gauge density in physical-volume work')
        chi = fields[2, trial.first]
        if not np.array_equal(ids, np.flatnonzero((chi == 0) & trial.free)):
            raise ValueError('Order sweep omitted or repeated free zero nodes')
        N = geometry.V.astype(np.longdouble)
        c = np.clip(N @ chi.astype(np.longdouble), 0, 1)
        weights = (1-c)**2*(1+2*c)
        old = np.sum(old_density*weights)
        new = np.sum(density*weights)
        if abs(old-reference['material']['gauge_energy']/1e-8) > 1e-11 or abs((new-old)-difference)>1e-12:
            raise ArithmeticError('Order work and moving-grid gauge means disagree')
        delta = np.longdouble(1e-9)
        # Exact cubic difference for each joined Q1 coefficient; summed coalesced entries on symmetry axes.
        a = N.T @ (density*(-6*c*(1-c)))
        b = N.power(2).T @ (density*(-3+6*c))
        d = N.power(3).T @ (2*density)
        gauge_work = delta*(a[ids]+delta*(b[ids]+delta*d[ids]))
        if np.any(gauge_work > 0):
            raise ArithmeticError('Positive gauge work on an increasing zero-order coefficient')
        shape_delta = []
        shape_fd = []
        for vector in np.eye(6,dtype=np.longdouble):
            tangent = vector.astype(np.clongdouble)*1e-24j
            shape_delta.append(float(-stencil_change(gauge,labels,tangent)[0].imag / 1e-24))
            width = np.longdouble(2)**-26
            plus = stencil_change(gauge,labels,width*vector)[0]
            minus = stencil_change(gauge,labels,-width*vector)[0]
            shape_fd.append(float((minus-plus)/(2*width)))
        second_elastic, second_forces, _, _ = solid.evaluate(labels, order)
        currents = trial.measure(fields, load)
        invariance = dict(elastic_vs_archived=float(elastic-reference['material']['elastic_energy']),
                          physical_total_vs_archived=float(physical-reference['physical_energy']),
                          elastic_before_minus_after=float(elastic-second_elastic),
                          all_physical_crystal_forces_max_abs=float(np.max(abs(second_forces-forces))),
                          currents_vs_archived={k:float(currents[k]-reference['current'][k]) for k in FLOW})
        if max(map(abs, [v for k,v in invariance.items() if k!='currents_vs_archived']
                   + list(invariance['currents_vs_archived'].values()))) > 1e-10:
            raise ArithmeticError(f'Physical invariance readback failed: {invariance}')
        results = {}
        for eps in (1e-8, 1e-10):
            g = np.longdouble(eps)*gauge_work
            total = physical_work+g
            pos = physical_work>0
            failures = [dict(node=int(ids[j]),physical=float(physical_work[j]),gauge=float(g[j]),
                             total=float(total[j])) for j in np.flatnonzero(total<=0)]
            results[str(eps)] = dict(old_gauge=float(eps*old),new_gauge=float(eps*new),
                    change_gauge=float(eps*difference), signed_shape_force_changes=(eps*np.array(shape_delta)).tolist(),
                    central_shape_change_errors=(eps*(np.array(shape_fd)-shape_delta)).tolist(),
                    all_node_positive=bool(np.all(total>0)), positive=int(np.sum(total>0)),
                    zero=int(np.sum(total==0)), negative=int(np.sum(total<0)), failures=failures,
                    max_abs_gauge_over_positive_physical=float(np.max(abs(g[pos])/physical_work[pos])))
        return dict(step=108,state_sha256=cp['sha256'], physical_work_sha256=sha(oracle_path), shape_modes=MODES,
                    free_zero=len(ids), physical_total_native_extended=float(physical), elastic_energy=float(elastic),
                    unit_gauge_energy_change=float(difference),unit_shape_force_changes=shape_delta,
                    central_unit_shape_force_changes=shape_fd,invariance=invariance, comparisons=results,
                    crystal_and_field_rematches=force_history(study))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('study', type=Path)
    p.add_argument('physical_order_work', type=Path)
    args = p.parse_args()
    print(json.dumps(run(args.study.resolve(), args.physical_order_work.resolve()), indent=2, allow_nan=False))
