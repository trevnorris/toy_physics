"""Read the accepted held R4 step through the archived smooth solver before Newton descent."""
import argparse
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np


def load_state(driver, folder):
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
    saved = folder/'last_active.npz'
    record = read(saved.with_suffix('.json'))
    if digest(saved) != record['sha256'] or record['outer'] < 29:
        raise ValueError('Missing accepted step or checkpoint hash mismatch')
    with np.load(saved, allow_pickle=False) as f:
        if (int(f['radius']), int(f['winding']), str(f['parent_sha256'])) != (4,1,CHECKPOINT_SHA):
            raise ValueError('Different radius, winding or reference state')
        if not np.array_equal(f['x'],driver.x) or float(f['epsilon']) != 1.e-8:
            raise ValueError('Changed held shape or gauge')
        driver.fields,driver.labels = f['field'].copy(),f['labels'].copy()
    return record


def inspect(study):
    from mouth_r4_smooth_held import SmoothHeld, SpatialCrystal, Gauge, order_grid, melted_cells
    from mouth_r4_smooth_fields import read_fields
    from mouth_r4_slaved_melt import FLOW
    from mouth_r4_projected_fields import location
    from mouth_r4_fold_order import all_low_after
    from mouth_whole_r4_relax import atomic_json
    result = study/'baseline'
    result.mkdir()
    with TemporaryDirectory(dir=result) as directory:
        driver = SmoothHeld(study.parent,Path(directory),1.e-8)
        cp = load_state(driver,study/'paused')
        trial = driver.trial()
        driver.current_trial = trial
        solid = SpatialCrystal(trial.b)
        order = order_grid(trial,driver.fields)
        barrier = driver.context(trial,driver.labels,order,melted_cells(order))
        driver.full_gauge = Gauge(trial,driver.reference_depth,order)
        _,gradient,state = driver.point(driver.labels,order,solid,barrier)
        physical,load,fields = read_fields(driver,trial,solid,driver.labels,driver.fields)
        flows = trial.measure(driver.fields,load)
        grad_by_node = np.zeros_like(driver.labels)
        channels=[]
        j=0
        for c,cm in enumerate(barrier.components):
            groups = np.flatnonzero(cm.free)
            signed = -gradient[j:j+len(groups)]
            j += len(groups)
            nodes = cm.source[groups]
            grad_by_node[...,c].flat[nodes] = signed
            for k in np.argsort(np.abs(signed))[-5:]:
                node=np.unravel_index(nodes[k],order.shape)
                channels.append(dict(node=list(map(int,node)),group=int(groups[k]),component=c,
                                     signed_force=float(signed[k]),chi=float(order[node]),
                                     **{k:v for k,v in location(trial,barrier.opt.axes,node).items() if k!='node'}))
        assert j==len(gradient)
        norms = np.linalg.norm(grad_by_node,axis=-1)
        top=[]
        for k in np.argsort(norms.ravel())[-5:][::-1]:
            node=np.unravel_index(k,order.shape)
            top.append(dict(**location(trial,barrier.opt.axes,node),chi=float(order[node]),
                            signed_components=grad_by_node[node].tolist(),l2=float(norms[node])))
        mat={k:state[k] for k in ('gauge_energy','elastic_energy','max_free_label_force','min_Jinv_active','min_Jinv_global')}
        diffs={k:mat[k]-cp['material'][k] for k in mat}
        if any(x!=0 for x in diffs.values()):
            raise ArithmeticError(f'Archived crystal force/energy changed on readback: {diffs}')
        result_data=dict(state_sha256=cp['sha256'],old_step=cp['outer'],
                    old_logged_differences=diffs,material=mat,
                    total_energy=physical+state['gauge_energy'],
                    physical_energy=physical,field=fields,current={k:flows[k] for k in FLOW},
                    full_melt_front=all_low_after(trial.model.z[1],order.transpose(0,2,1),0.),
                    active_gradient_l2=float(np.linalg.norm(gradient)),active_gradient_count=len(gradient),
                    top5_nodes_by_vector_norm=top,top5_scalar_forces=sorted(channels,
                                  key=lambda rec:abs(rec['signed_force']),reverse=True)[:5])
        atomic_json(result/'record.json',result_data)
        np.save(result/'active_gradient.npy',gradient,allow_pickle=False)
        print(result_data,flush=True)


if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument('study',type=Path)
    study=args.parse_args().study.resolve()
    # Re-run the immutable *readout* only; the relaxation is never invoked.
    sys.path.insert(0,str(study.parent/'r4_smooth_paired01/pin/software/stage1_solver/tools'))
    inspect(study)
