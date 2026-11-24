"""One bounded order relaxation at fixed R4 labels, phase, amplitude and geometry."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np


def read(file):
    return json.loads(file.read_bytes())


def digest(file):
    with file.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def run(area, out, result):
    pin = area / 'pin/software/stage1_solver/tools'
    if any(digest(pin/name) != h for name,h in read(area/'sources.json').items()):
        raise ValueError('Pinned physical source changed')
    sys.path.insert(0, str(pin))
    from mouth_continuous_order import KAPPA, WELL, mat
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_r4_eulerian_release import Driver, order_grid, melted_cells
    from mouth_r4_eulerian_relax_labels import LabelRelaxation, absmax, chart_orientation
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_eulerian_corner_audit import point
    from mouth_r4_refinement_review import _collar_reference
    from mouth_r4_melt_barrier import MeltBarrier
    from mouth_r4_stage2_outcome import verify_guard
    from mouth_r4_fold_order import all_low_after

    out.mkdir(parents=True, exist_ok=False)
    verify_guard(area/'primary_guard', 1)
    stop = read(area/'primary/held_stop.json')
    original = area/'primary'/stop['file']
    if digest(original) != '418376c3eb1922f6335543f43dfee11e531cebbb6958d37afd73e844571b653a' or digest(original) != stop['sha256']:
        raise ValueError('Incorrect 4c3c0a59 held state')
    with np.load(original) as saved:
        fields, labels, q = (saved[k] for k in ('field','labels','x'))
        epsilon = float(saved['epsilon'])
        if (int(saved['radius']),int(saved['winding']),epsilon) != (4,1,1.e-8):
            raise ValueError('Incorrect physical parameters')
        parent = saved['parent_sha256']
    driver = Driver(area.parent, out, epsilon)
    if not np.array_equal(q, driver.x):
        raise ValueError('Changed fixed shape')
    length = driver.shape(q).length
    trial = ShapeTrial(driver.base,driver.prep,q,8,24,6,3,2,precondition=False,
                       axes=axes_for_depth(length,driver.reference_depth))
    model = trial.model
    model.extended_energy_sum = True
    if not np.all(fields[1,model.feed]==0.) or not np.all(fields[1,model.bulk]==.4510084282857874):
        raise ValueError('Incorrect far reservoir')
    trial.K = sum(mat(g,KAPPA*model.hq) for g in trial.G).tocsr()
    trial.diag = trial.K.diagonal() + trial.V.power(2).T @ (2*WELL*model.hq)
    crystal = SpatialCrystal(trial.b)
    order = order_grid(trial,fields)
    _,_,src,_ = crystal.evaluate(labels,order,False)
    inner = np.bincount(trial.b.ids.ravel(),weights=src.transpose(0,2,1).ravel(), minlength=model.size)
    model.material_load = trial.ext + inner
    joined = trial.pull.T @ model.material_load
    cell = [15,9,1]
    def measure(u,name):
        chi = order_grid(trial,u)
        mask = melted_cells(chi)
        opt = LabelRelaxation(trial.b,labels,chi,epsilon)
        folder = out/name
        folder.mkdir()
        gauge_model = MeltBarrier(trial.b,opt,folder)
        unit,gg,parts,_ = gauge_model.sample(labels)
        elastic,_,local,_ = crystal.evaluate(labels,chi,False)
        if not np.array_equal(local,src):
            raise ArithmeticError('Physical crystal order deposition depends on order at fixed labels')
        energy,g = model.evaluate(u)
        c = u[2,trial.first]
        grd = trial.pull.T@g[2]
        kkt = np.where(c==0,np.minimum(grd,0),np.where(c==1,np.maximum(grd,0),grd))
        with mp.workdps(65):
            sites = {key:point(trial.b,opt.axes,labels,chi,cell,v,_collar_reference) for key,v in
                     [('limiting_corner',[1,1,0]),('opposite_vertex',[0,0,1])]}
        return dict(physical_energy=float(energy),crystal_energy=float(elastic),gauge_energy=float(epsilon*unit),
                    total_energy=float(energy+epsilon*unit),gauge_parts=parts,
                    fixed_mask_active_gauge_grad=epsilon*absmax(gauge_model.extract(gg,'free')),
                    fixed_mask_melt_gauge_grad=epsilon*absmax(gauge_model.extract(gg,'dead')),
                    chart=chart_orientation(opt.axes,labels,opt.active),melt_cells=int(mask.sum()),
                    mixed_corner_order_vertices=chi[15:17,9:11,1:3].tolist(),sites=sites,
                    sampled_front=all_low_after(opt.axes[1]*length,chi.transpose(0,2,1),0.),
                    order_max_by_depth=chi.max(axis=(0,2)).tolist(),
                    depth_knots=(opt.axes[1]*length).tolist(),order_raw=absmax(grd[trial.free]),
                    physical_order_kkt=absmax(kkt[trial.free]),projected=absmax(trial.projected(c,grd)),
                    amplitude_residual=absmax(g[0,model.masks[0]]),phase_residual=absmax(g[1,model.masks[1]]))
    before = measure(fields,'before')
    old_mask = melted_cells(order)
    old_field = fields.copy()
    iterations = []
    def record(i,energy,resid,c,scale,damping):
        iterations.append(dict(index=i,order_part_energy=float(energy),previous_projected_residual=float(resid),
                               fraction=float(scale),damping=damping))
    count,resid = trial.solve_order(fields,joined,False,2.e-9,log=record)
    after = measure(fields,'after')
    if not np.array_equal(old_field[:2],fields[:2]) or not np.array_equal(old_field[2,~trial.free[trial.number]], fields[2,~trial.free[trial.number]]):
        raise ArithmeticError('Fixed fields or distant order trace changed')
    new_mask = melted_cells(order_grid(trial,fields))
    if any(after['sites'][k]['Jinv'] != before['sites'][k]['Jinv'] for k in before['sites']):
        raise ArithmeticError('Fixed label or geometry changed')
    delta = (after['physical_energy']-before['physical_energy']) + (after['gauge_energy']-before['gauge_energy'])
    test = count>0 and resid<2.e-9 and after['physical_order_kkt']<2.e-7 and delta<0
    snapshot = out/'fixed_label_order.npz'
    np.savez_compressed(snapshot, field=fields,labels=labels,x=q,epsilon=epsilon,parent_sha256=parent,
                        source_sha256=digest(original),radius=4,winding=1,iteration=10)
    data = dict(source_state_sha256=digest(original),frozen_source_manifest_sha256=digest(area/'sources.json'),
                measurement_source_sha256=digest(Path(__file__)),saved_sha256=digest(snapshot),
                epsilon=epsilon,phase_drop=.4510084282857874,radius=4,winding=1,order_steps=count,
                passed=test,before=before,after=after,iterations=iterations,total_energy_change=delta,
                physical_energy_change=after['physical_energy']-before['physical_energy'],
                gauge_energy_change=after['gauge_energy']-before['gauge_energy'],
                cells_entered_melt=int((new_mask&~old_mask).sum()),cells_left_melt=int((old_mask&~new_mask).sum()),
                old_new_different_order_nodes=int(np.count_nonzero(old_field[2]!=fields[2])))
    with result.open('x') as f:
        json.dump(data,f,indent=2,allow_nan=False)
        f.write('\n')
    print(count, before['physical_order_kkt'],after['physical_order_kkt'],delta,after['sampled_front'],
          {k:(s['physical']['w'],s['chi'],s['Jinv']) for k,s in after['sites'].items()},flush=True)
    if not test:
        raise ArithmeticError('Physical obstacle stationarity or same-functional energy gate failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('area','out','measurement'):
        p.add_argument(key,type=Path)
    args=p.parse_args()
    run(args.area.resolve(), args.out.resolve(), args.measurement.resolve())
