"""Compare both orientation entry points and the frozen held functional at step 108."""
import argparse
import ctypes as ct
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from mouth_r4_aligned_crystal import AlignedCrystal, cell_call
from mouth_r4_aligned_derivatives import EPS, load
from mouth_r4_aligned_rematch import FIELDS
from mouth_r4_aligned_stencil import EDGES, samples
from mouth_r4_depth_cert import digest
from mouth_r4_edge_descent import CELL, NODES
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, chart_orientation, cell_order
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_slaved_melt import FLOW
from mouth_r4_smooth_fields import GaugeGeometry, read_fields


def pinned_kernel(pin, directory, expected):
    names=('mouth_r4_eulerian_solid.cpp','mouth_cp10_spatial_crystal.cpp','mouth_free_shape_solid.cpp')
    if any(digest(pin/n)!=expected[n] for n in names):
        raise ValueError('Historical physical kernel differs from the reviewed binary inputs')
    so=directory/'baseline_physics.so'
    subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(pin/names[0]),'-o',str(so)],check=True)
    return ct.CDLL(str(so)).r4_eulerian_material


def check_state(study, case, directory, old_native, exact):
    driver,bare,derivative_solid,order=load(study,case,directory)
    del bare,derivative_solid
    # Use the same projected field metric as the original Gate-1 rematch.
    trial=driver.trial()
    trial.smooth_order_geometry=GaugeGeometry(trial,aligned=True)
    solid=SpatialCrystal(trial.b)
    state=json.loads((study/case/'rematched.json').read_bytes())
    frozen=json.loads((study/case/'labels.json').read_bytes())['material']
    labels=driver.labels
    old=SpatialCrystal(trial.b)
    old_native.argtypes=solid.function.argtypes
    old_native.restype=solid.function.restype
    old.function=old_native
    e0,g0,c0,_=old.evaluate(labels,order)
    e1,g1,c1,_=solid.evaluate(labels,order)
    joint=AlignedCrystal(driver,trial,solid,labels,order)
    _,_,info=joint.evaluate(labels)
    physical,load_on_order,field=read_fields(driver,trial,solid,labels,driver.fields)
    currents=trial.measure(driver.fields,load_on_order)
    control=dict(elastic_energy=e1-e0,all_elastic_label_gradient_max=float(np.max(abs(g1-g0))),
                 all_elastic_order_load_max=float(np.max(abs(c1-c0))),
                 total_energy=(physical+field['gauge_energy'])-state['total'],
                 gauge_energy=info['gauge_energy']-frozen['gauge_energy'],
                 physical_energy=physical-state['final']['physical_energy'],
                 peak_active_label_force=info['max_free_label_force']-frozen['max_free_label_force'])
    field_delta={k:field[k]-state['final'][k] for k in FIELDS}
    flow={k:float(currents[k]) for k in FLOW}
    flow_delta={k:v-state['current'][k] for k,v in flow.items()}
    with (directory/'static_readback.json').open('x') as stream:
        json.dump(dict(state=control,fields=field_delta,currents=flow_delta),stream,allow_nan=False)
    arm=json.loads((study/case/'hessian.json').read_bytes())['orientation']
    direction=arm['perturbed_radial_node']
    if (digest(study/case/'rematched.npz')!=exact['inputs'][case] or
            tuple(direction['node'])!=NODES[0] or direction['direction']!=exact['reconstruction']['radial_direction']):
        raise ValueError('The degree of freedom differs from the diagnostic checkpoint')
    axes=joint.opt.axes
    width=np.array([[np.diff(ax)[i] for ax,i in zip(axes,CELL)]],np.longdouble)
    ids=tuple(np.array(NODES).T)
    active=cell_order(order)
    positive={};negative={}
    for param in (-.95,-4.):
        y=labels.copy()
        y[tuple(direction['node'])+(0,)]+=param*direction['direction']
        p=y[ids][None].astype(np.longdouble)
        fixed,t,moving=samples(p,width)
        edge=EDGES.index((0,1,2))
        if param==-.95:
            control_sha=hashlib.sha256(y.tobytes()).hexdigest()
            if control_sha!=exact['reconstruction']['labels_sha256_after']:
                raise ValueError('Positive witness replay differs from the signed trial')
            standalone=chart_orientation(axes,y,active)
            _,g=cell_call('values',p,width,gradient=True)
            _,_,_,m=solid.evaluate(y,order,gradient=False)
            _,_,coupled=joint.evaluate(y)
            bound=m['min_chart_Jinv_bound']
            historical=exact['chart']['chart_det_global_lower']
            positive=dict(label_sha256=control_sha,physical_lower=bound,standalone=standalone,
                          coupled_lower=coupled['min_chart_Jinv_bound'],prior_python_lower=historical,
                          analytic_exact_rational_lower=exact['target_cell_exact_proof']['certified_lower'],
                          difference_from_prior_bound=bound-historical,edge_min=float(moving[0,edge]),
                          local_unit_opening_slope=float(g[0,0,0]*direction['direction']))
            positive['passed']=(standalone['chart_lowest_cell']==list(CELL)
                and coupled['min_chart_Jinv_bound']==bound==standalone['chart_det_global_lower']
                and 0<bound<=positive['analytic_exact_rational_lower'] and abs(bound-historical)<1e-14
                and positive['local_unit_opening_slope']<0)
        else:
            negative.update(vertices_and_midpoints_min=float(fixed.min()),quadratic_edge_sample=float(moving[0,edge]),
                            sample_fraction=float(t[0,edge]))
            try:
                chart_orientation(axes,y,active)
            except InfeasibleLabels as exc:
                negative['chart']=exc.witness
            try:
                solid.evaluate(y,order,gradient=False)
            except InvalidSpatialCrystal as exc:
                negative['physical']=dict(reason=str(exc),**exc.diagnostics)
            try:
                cell_call('values',p,width)
            except InfeasibleLabels as exc:
                negative['barrier']=exc.witness
            cc=negative.get('chart',{}); pc=negative.get('physical',{})
            negative['passed']=(negative['quadratic_edge_sample']<0<negative['vertices_and_midpoints_min']
                and cc.get('cell')==pc.get('cell')==list(CELL) and cc.get('code')==pc.get('code')==2
                and cc.get('chart_det',0)<0 and cc.get('chart_fractions')==pc.get('local')
                and cc.get('chart_lower')==pc.get('chart_bound')
                and negative.get('barrier',{}).get('stencil_index')==27+edge)
    return dict(passed=all(v==0 for table in (control,field_delta,flow_delta) for v in table.values())
                and positive['passed'] and negative['passed'],
                state_sha256=state['sha256'],epsilon=driver.epsilon,
                energy_and_field=dict(total=physical+field['gauge_energy'],physical=physical,elastic=e1,
                    gauge=info['gauge_energy'],peak_active_label_force=info['max_free_label_force'],
                    field_residuals={k:field[k] for k in FIELDS},currents=flow),
                invariance_delta=dict(state=control,field=field_delta,currents=flow_delta),
                positive=positive,negative=negative)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',type=Path)
    parser.add_argument('output_directory',type=Path)
    args=parser.parse_args()
    study=args.study.resolve(); out=args.output_directory.resolve()
    out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).resolve().parent
    anchor=source.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_depth_cert.json'
    exact=json.loads(anchor.read_bytes())
    old=pinned_kernel(study.parent/'r4_newton_held01/pin/software/stage1_solver/tools',out,exact['source_digests'])
    result=dict(command=sys.orig_argv,prior_exact_certificate_sha256=digest(anchor),
        source_sha256={n:digest(source/n) for n in ('mouth_r4_eulerian_solid.cpp','mouth_r4_eulerian_solid.py',
            'mouth_r4_full_label_work.cpp','mouth_r4_eulerian_relax_labels.py',Path(__file__).name)},
        gauges={})
    for case in EPS:
        directory=out/case; directory.mkdir()
        result['gauges'][case]=check_state(study,case,directory,old,exact)
        with (directory/'controls.json').open('x') as stream:
            json.dump(result['gauges'][case],stream,allow_nan=False,indent=2)
    result['passed']=all(v['passed'] for v in result['gauges'].values())
    with (out/'result.json').open('x') as f:
        json.dump(result,f,allow_nan=False,indent=2); f.write('\n')
    print(json.dumps(result,allow_nan=False),flush=True)
    if not result['passed']:
        raise ArithmeticError('Whole-cell controls or held-state invariance failed')


if __name__=='__main__':
    main()
