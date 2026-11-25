"""Reconstruct the two rejected step-167 search directions without advancing any field."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def load_state(folder, name):
    record = read(folder/(name+'.json'))
    path = folder/record['file']
    if digest(path) != record['sha256']:
        raise ValueError(f'Invalid saved-state digest at {path}')
    with np.load(path,allow_pickle=False) as data:
        if (int(data['radius']),int(data['winding']),float(data['epsilon'])) != (4,1,1e-8):
            raise ValueError('Wrong preparation')
        return record, data['labels'], data['field'], data['x']


def scan_history(folder):
    tight = []
    anchor = None
    for n in range(149,167):
        record = read(folder/f'step_{n:05}.json')
        if record['outer'] != n:
            raise ValueError(f'Missing step {n}')
        search = record['search']
        if search['jacobi']['rebuilt']:
            anchor = n-1
        for trial in search['rejections']:
            site = trial['site']
            location = site.get('chart_fractions',site.get('local'))
            if location is not None and 'cell' in site:
                point = (tuple(site['cell']),tuple(location))
                if point not in tight:
                    tight.append(point)
    if anchor is None:
        raise ValueError('Missing Hessian diagonal anchor')
    return anchor,tight


def reconstruct(study, out):
    pin = study/'held_pin'
    source = pin/'software/stage1_solver/tools'
    if not all(digest(source/k)==v for k,v in read(pin/'manifest.json').items()):
        raise ValueError('Historical solver pin changed')
    sys.path.insert(0,str(source))
    from mouth_continuous_depth_mesh import axes_for_depth
    from mouth_continuous_shape_relax import ShapeTrial
    from mouth_r4_aligned_crystal import AlignedCrystal
    from mouth_r4_corner_projection import project_halfspaces
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_eulerian_solid import SpatialCrystal
    from mouth_r4_smooth_fields import GaugeGeometry
    from mouth_r4_smooth_held import SmoothHeld

    out.mkdir(parents=True,exist_ok=False)
    folder = study/'primary'
    anchor,tight = scan_history(folder)
    anchor_record, previous, old_fields, old_x = load_state(folder,f'step_{anchor:05}')
    last, labels, fields, shape = load_state(folder,'step_00166')
    if last['sha256']!='4b84a0d2ad53889a1d619e68960df70233ceb589fc4fe6c336c6e59e832a2694':
        raise ValueError('Wrong step-166 state')
    if not (np.array_equal(old_fields,fields) and np.array_equal(old_x,shape)):
        raise ValueError('The reconstructed diagonal belongs to different fields')
    issue = read(folder/'ISSUE.json')['witness']
    if issue['material']!=last['material'] or (issue['cg_reason'],issue['cg_iterations'])!=('nonpositive curvature',2):
        raise ValueError('Inconsistent step-167 terminal witness')

    driver = SmoothHeld(study.parent,out,1e-8)
    if not np.array_equal(driver.x,shape):
        raise ValueError('Changed held shape')
    trial = ShapeTrial(driver.base,driver.prep,shape,8,24,6,3,2,precondition=False,
                       axes=axes_for_depth(driver.shape(shape).length,driver.reference_depth))
    trial.model.extended_energy_sum=True
    trial.smooth_order_geometry=GaugeGeometry(trial,aligned=True)
    solid=SpatialCrystal(trial.b)
    order=order_grid(trial,fields)
    joint=AlignedCrystal(driver,trial,solid,labels,order)
    _, grad, material = joint.evaluate(labels)
    if material!=last['material']:
        raise ArithmeticError('Native frozen-state force differs from the archived step')
    _, history_diag = joint.inverse_diagonal(previous,anchor_record['material']['chart']['chart_det_global_lower'])
    archived = last['search']['jacobi']
    if any(v!=archived[k] for k,v in history_diag.items() if k!='rebuilt'):
        raise ArithmeticError('Reconstructed Jacobi assembly differs from recorded curvature')
    diag, cached = joint.inverse_diagonal(labels,material['chart']['chart_det_global_lower'])
    if cached['rebuilt']:
        raise ArithmeticError('Expected to reuse the historical inverse diagonal')

    # The native solve truncated on its second CG curvature test. Reproduce
    # that recurrence with the pinned HVP and check both curvature signs.
    rhs=-grad
    pj_raw=diag*rhs
    rz=float(rhs@pj_raw)
    h1=joint.hvp(labels,pj_raw)
    curvature1=float(pj_raw@h1)
    if curvature1<=0:
        raise ArithmeticError('First CG curvature differs from the trajectory')
    beta=rz/curvature1
    newton_raw=beta*pj_raw
    residual=rhs-beta*h1
    z=diag*residual
    next_direction=z+(float(residual@z)/rz)*pj_raw
    curvature2=float(next_direction@joint.hvp(labels,next_direction))
    if curvature2>=0 or float(grad@newton_raw)>=0:
        raise ArithmeticError('Changed CG stop or descent orientation')
    # The direction switches when a new certificate witness enters the
    # projection. Reconstruct the final tested projection for each bank.
    series=issue['tests']
    saved = {}
    metadata = {}
    for index,(name,raw) in enumerate((('newton',newton_raw),('pj',pj_raw))):
        family=series[14*index:14*(index+1)]
        if [v['alpha'] for v in family] != [2.**(-i) for i in range(14)]:
            raise ValueError('Unexpected line-search path')
        for v in family:
            site=v['site']
            location=site.get('chart_fractions',site.get('local'))
            if location is not None and 'cell' in site:
                entry=(tuple(site['cell']),tuple(location))
                if entry not in tight:
                    tight.append(entry)
        normals=[]
        applied=[]
        for cell,local in tight:
            normal,det=joint.corner_normal(labels,cell,local)
            rate=float(normal@raw)
            if det<.1 and rate<-.9*det and np.any(normal):
                normals.append(normal)
                applied.append(dict(cell=cell,local=local,det=det,rate_before_projection=rate))
        plan=project_halfspaces(raw,normals,lambda v:diag*v) if normals else raw
        saved['packed_'+name]=plan
        d=saved[name]=joint.lift(plan)
        slope=float(grad@plan)
        width=family[-1]['alpha']
        shifted=labels+width*d
        _,_,probe=joint.evaluate(shifted)
        crystal=probe['elastic_energy']-material['elastic_energy']
        gauge=float(joint.epsilon*joint.weighted.work(labels,shifted))
        delta=crystal+gauge
        archived_delta=family[-1]['site']['delta']
        metadata[name]=dict(slope=slope,alpha=width,crystal=crystal,gauge=gauge,
                        delta=delta,archived_delta=archived_delta,work_error=delta-archived_delta,
                        projected_sites=applied)
        if not(slope<0 and delta>0 and abs(delta-archived_delta)<2e-14):
            raise ArithmeticError(f'Failed terminal reproduction: {name}: {metadata[name]}')
    array_path=out/'directions.npz'
    with array_path.open('xb') as dest:
        np.savez_compressed(dest,**saved,packed_gradient=grad,labels=labels,order=order)
    result=dict(command=sys.orig_argv, script_sha256=digest(Path(__file__)),state_sha256=last['sha256'],
                pin_sha256=digest(pin/'manifest.json'),directions_sha256=digest(array_path),
                issue_sha256=digest(folder/'ISSUE.json'),jacobi_anchor=anchor,
                anchor_sha256=anchor_record['sha256'],inverse_diagonal=cached,
                two_cg_curvatures=[curvature1,curvature2],cg_first_step=beta,retained_corner_count=len(tight),
                slopes={n:v['slope'] for n,v in metadata.items()},terminal_banks=metadata,
                radius=4,epsilon=joint.epsilon,optimization_steps=0,passed=True)
    temp=out/'result.tmp'
    temp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    temp.replace(out/'result.json')
    print(json.dumps(result,allow_nan=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    reconstruct(args.study.resolve(),args.output.resolve())
