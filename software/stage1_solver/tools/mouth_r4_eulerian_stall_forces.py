"""Six frozen-field shape forces at an R4 label stall, with feasible signed probes."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_continuous_depth_mesh import axes_for_depth
from mouth_continuous_shape_relax import ShapeTrial
from mouth_cp10_field_shape import Bath, COMPONENTS
from mouth_cp10_label_gauge import work as gauge_work
from mouth_cp10_nested_forces import integral, kernel
from mouth_cp10_spatial_crystal import MODES
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_eulerian_relax_labels import chart_orientation, InfeasibleLabels
from mouth_r4_eulerian_release import Driver, block_axes, melted_cells, order_grid
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_whole_r4_relax import atomic_json


def audit(scratch, source, output):
    output.mkdir(parents=True,exist_ok=False)
    record = read(source/'label_stall.json')
    state = source/'label_stall.npz'
    if digest(state)!=record['sha256']:
        raise ValueError('Saved stall state and its guard metadata disagree')
    guard=read(source.parent/(source.name+'_guard')/'outcome.json')
    if not guard['limitsVerified'] or read(source/'issue.json')['error']!='LabelStall':
        raise ValueError('Missing bounded original LabelStall endpoint')
    pins=read(source.parent/'sources.json')
    here=Path(__file__).resolve().parent
    for name,want in pins.items():
        if digest(here/name)!=want:
            raise ValueError(f'Original solver source changed: {name}')
    with np.load(state) as f:
        if str(f['parent_sha256'])!=CHECKPOINT_SHA or int(f['radius'])!=4 or int(f['winding'])!=1 or int(f['iteration'])!=10:
            raise ValueError('Wrong held state')
        eps=float(f['epsilon'])
        driver=Driver(scratch.resolve(),output,eps)
        driver.labels,driver.x,driver.fields=(f[k].copy() for k in ('labels','x','field'))
    # Fixed-field work requires geometry and stored fields only. In particular,
    # constructing a sparse descent preconditioner would consume several GiB.
    base_shape=driver.shape(driver.x)
    trial=ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
                     axes=axes_for_depth(base_shape.length,driver.reference_depth))
    axes=block_axes(trial.b)
    c=order_grid(trial,driver.fields)
    melt=melted_cells(c)
    labels=driver.labels
    independent_orientation=chart_orientation(axes,labels,~melt)
    spatial=SpatialCrystal(trial.b)
    energy,gradient,_,diagnostics=spatial.evaluate(labels,c,True)
    parts=driver.parts_at(trial,spatial)
    stored=parts(np.zeros(6))
    if abs(energy-stored[-2])>1e-12 or abs(energy-record['label']['elastic_energy'])>1e-12:
        raise ArithmeticError('Crystal state was not reproduced')
    result=dict(state_sha256=digest(state),original_guard=guard,source=source.name,scope='fixed saved fields and unrelaxed labels at entrance R4; not a stationary-state force',
                epsilon=eps,components=list(COMPONENTS)+['welded_crystal','harmonic_gauge'],energy=stored.astype(float).tolist(),
                accepted_orientation=independent_orientation,orientation_min_sampled=diagnostics['min_Jinv_global'],modes={})
    def save():
        atomic_json(output/'readout.json',result)
    save()
    # Fixed nodal-label analytic shape derivatives. The production signed probes
    # also move the normalized identity depth; include that work explicitly.
    depth=trial.shape.length
    t=axes[1]
    tidot=np.where(t*depth<depth-.12,(t*depth/(depth-.12))**25,1.)-t
    baseline=integral(kernel(output),trial.b,axes,tidot,labels,c,(spatial.points,spatial.weights),True)
    if abs(energy-baseline['energy'])>1e-12:
        raise ValueError('Fixed-label and production crystal evaluators disagree')
    label_work=-np.sum(gradient[...,1]*tidot[None,:,None],dtype=np.longdouble)
    result['solid_label_depth_transport_force']=float(label_work)
    bath=Bath(trial,driver.base,driver.fields,driver.reference_depth)
    for i,unit in enumerate(np.eye(6)):
        inc=1.e-24j*unit
        cs=bath.energy(inc)
        moved=labels.astype(complex)
        if i==0:
            moved[...,1]+=1.e-24j*tidot[None,:,None]
        cg=gauge_work(trial.shape,trial.polar_shift,driver.reference_depth,axes,moved,melt,inc)
        analytic=np.r_[-cs.imag/1.e-24, baseline['analytic'][i]+(float(label_work) if i==0 else 0.), -eps*cg.imag/1.e-24]
        mode=result['modes'][MODES[i]]=dict(analytic=analytic.astype(float).tolist(),force=float(analytic.sum()),signed=[])
        save()
        for scale in (1.,.5):
            h=scale*(1e-5 if i==0 else 1e-4)
            item=dict(h=h,banks={})
            for s in (-1,1):
                try:
                    values=parts(s*h*unit)
                except (InvalidSpatialCrystal,InfeasibleLabels) as e:
                    item['banks'][str(s)]=dict(feasible=False,issue=str(e),diagnostics=getattr(e,'diagnostics',getattr(e,'witness',None)))
                else:
                    item['banks'][str(s)]=dict(feasible=True,energies=values.astype(float).tolist(),one_sided=float((stored-values).sum()/(s*h)))
            if all(x['feasible'] for x in item['banks'].values()):
                minus,plus=[np.array(item['banks'][str(s)]['energies'],dtype=np.longdouble) for s in (-1,1)]
                dforce=(minus-plus)/(2*h)
                item['force']=float(dforce.sum())
                item['analytic_minus_fd']=float(analytic.sum()-dforce.sum())
                item['fd_components']=dforce.astype(float).tolist()
                item['valid_central']=True
            else:
                item['valid_central']=False
            mode['signed'].append(item)
            save()
    result['complete']=True
    save()


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('scratch',type=Path)
    ap.add_argument('source',type=Path)
    ap.add_argument('output',type=Path)
    args=ap.parse_args()
    audit(args.scratch,args.source,args.output)
