"""Held-state finite differences after solving the melt gauge separately on every trial."""
import argparse
from pathlib import Path

import numpy as np

from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_shape_relax import ShapeTrial
from mouth_cp10_spatial_crystal import MODES
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, LabelRelaxation
from mouth_r4_eulerian_release import order_grid
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_harmonic_shape import HarmonicShapeDriver
from mouth_r4_stage2_outcome import verify_guard
from mouth_whole_r4_relax import atomic_json


def differences(evaluate, zero, analytic, steps, sign):
    """Keep signed partial-energy differences above the large common bath energy."""
    rows = []
    for h in steps:
        banks = {}
        for side in (-1,1):
            try:
                values, data = evaluate(side*h)
            except (InvalidSpatialCrystal,InfeasibleLabels) as error:
                witness = error.diagnostics if isinstance(error,InvalidSpatialCrystal) else error.witness
                banks[str(side)] = dict(valid=False,reason=str(error),witness=witness)
            else:
                rate = sign*np.sum(values-zero,dtype=np.longdouble)/(side*h)
                banks[str(side)] = dict(valid=True,one_sided=float(rate),**data)
        row = dict(h=h,banks=banks)
        if all(bank['valid'] for bank in banks.values()):
            row['central'] = .5*(banks['-1']['one_sided']+banks['1']['one_sided'])
            row['central_minus_analytic'] = row['central']-analytic
        rows.append(row)
    return dict(analytic=analytic,steps=rows)


def active_probes(solid, harmonic, labels, chi, limiting_vertex):
    trace = harmonic._pack(labels)
    projected = harmonic._lift(trace)
    energy, gradient, _ = harmonic._point(solid,projected,chi)
    start = np.array([energy],dtype=np.longdouble)
    components = harmonic._components()
    starts = np.cumsum([0]+[cm.free.sum() for cm in components])
    probes = {}
    for c,cm in enumerate(components):
        # An interior stiffness variation can be two-sided even when the remote
        # harmonic corner has vanishingly little clearance. Also vary a label
        # which directly controls the harmonic extension at that corner.
        which = np.flatnonzero(cm.free)
        uncoupled = np.asarray(cm.cross.getnnz(axis=0)).ravel()==0
        active_gradient = gradient[starts[c]:starts[c+1]]
        index = np.argmax(np.where(uncoupled,np.abs(active_gradient),-1.))
        if not uncoupled[index]:
            raise ArithmeticError('No interior ordered coordinate')
        selected = [('interior',index)]
        limiting_group = cm.groups[np.ravel_multi_index(tuple(limiting_vertex),chi.shape)]
        if cm.free[limiting_group]:
            selected.append(('limiting_vertex',int(np.searchsorted(which,limiting_group))))
        for name,local in selected:
            absolute = int(starts[c]+local)
            def bank(inc):
                perturbed = trace.copy()
                perturbed[absolute] += inc
                extended = harmonic._lift(perturbed)
                total, _, info = harmonic._point(solid,extended,chi)
                return np.array([total],dtype=np.longdouble),dict(
                        min_chart=info['min_chart_Jinv_bound'],harmonic_residual=info['max_unit_harmonic_residual'])
            result = differences(bank,start,float(gradient[absolute]),(2.e-6,1.e-6),1.)
            group = which[local]
            result.update(component=c,free_coordinate=absolute,
                          node=list(map(int,np.unravel_index(cm.source[group],chi.shape))),
                          incident_melt_cross_entries=int(cm.cross.getnnz(axis=0)[local]))
            probes[f'{name}_{c}'] = result
    return probes


def run(scratch, case, output):
    area = scratch/'r4_eulerian_stage2_01'
    verify_guard(area/(case+'_guard'),1)
    endpoint = read(area/case/'label_stall.json')
    state = area/case/endpoint['file']
    signed_name = 'secondary_fixed_field' if case=='secondary' else case
    signed = scratch/'r4_spatial_stall_force01'/signed_name
    verify_guard(signed.with_name(signed_name+'_guard'),0)
    original = read(signed/'readout.json')
    if not original['complete'] or digest(state)!=original['state_sha256'] or original['state_sha256']!=endpoint['sha256']:
        raise ValueError('Analytic endpoint source mismatch')
    with np.load(state,allow_pickle=False) as f:
        if str(f['parent_sha256'])!=CHECKPOINT_SHA or int(f['radius'])!=4 or int(f['iteration'])!=10:
            raise ValueError('Wrong endpoint')
        driver = HarmonicShapeDriver(scratch,output.parent,float(f['epsilon']))
        driver.fields,driver.labels,driver.x = [f[k].copy() for k in ('field','labels','x')]
    shape = driver.shape(driver.x)
    trial = ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
                       axes=axes_for_depth(shape.length,driver.reference_depth))
    solid = SpatialCrystal(trial.b)
    chi = order_grid(trial,driver.fields)
    harmonic = LabelRelaxation(trial.b,driver.labels,chi,driver.epsilon)
    reduced = driver.parts_at(trial,solid)
    zero = reduced(np.zeros(6))
    original_labels = driver.labels
    stats = dict(case=case,epsilon=driver.epsilon,state_sha256=digest(state),
                 analytic_archive_sha256=digest(signed/'readout.json'),modes={},
                 base_reextension_max=float(np.max(np.abs(reduced.last_labels-original_labels))),
                 base_active_reextension_max=float(np.max(np.abs((reduced.last_labels-original_labels)[harmonic.stencil]))),
                 base_energies=zero.astype(float).tolist(),physical_mode_order=MODES)
    atomic_json(output,stats)
    corner = read(scratch/'r4_eulerian_corner01'/(case+'_support.json'))
    stats['active_probes'] = active_probes(solid,harmonic,original_labels,chi,corner['vertex'])
    atomic_json(output,stats)
    for i,(name,unit) in enumerate(zip(MODES,np.eye(6))):
        # At a stationary slaved coordinate the envelope partial is precisely
        # the previously checked fixed-offset derivative at the SAME state.
        # Its bath, crystal and K(q) derivatives were evaluated independently.
        analytic = original['modes'][name]['force']
        def bank(inc):
            values = reduced(inc*unit)
            return values,dict(min_chart=reduced.last_info['min_chart_Jinv_bound'],
                               min_chart_active=reduced.last_info['min_chart_Jinv_active_bound'])
        h = 1.e-5 if i==0 else 1.e-4
        stats['modes'][name] = differences(bank,zero,analytic,(h,h/2),-1.)
        atomic_json(output,stats)
    stats['complete'] = True
    atomic_json(output,stats)
    print('Active central checks:', {k:r['steps'][-1].get('central_minus_analytic') for k,r in stats['active_probes'].items()})
    print('Shape central checks:', {k:r['steps'][-1].get('central_minus_analytic') for k,r in stats['modes'].items()})


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('case',choices=('primary','secondary'))
    p.add_argument('output',type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    run(a.scratch.resolve(),a.case,a.output.resolve())
