"""Static branchwise Eulerian energy on the imaged R4 quadrature and unchanged port fields."""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_continuous_shape_relax import ContinuousRelease
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, PORTS, digest, read
from mouth_exact_melt_trial import ExactMeltTrial
from mouth_free_shape_solid import Data, Results, Pointer, pointer
from mouth_r4_refinement_audit import inputs
from mouth_r4_refined_solid import RefinedSolid

HERE=Path(__file__).resolve().parent


def run(root,output,which):
    owner=inputs(root)
    owner.trial_type=ExactMeltTrial
    if which=='cp10':
        src=root/'r4_melted_survival01/coarse/state_00010.npz'
        expected=CHECKPOINT_SHA
    else:
        src=root/'r4_local_refinement01/first_step/material.npz'
        expected=read(root/'r4_local_refinement01/first_step/result.json')['state_sha256']
    if digest(src)!=expected:
        raise ValueError('Wrong material state')
    with np.load(src,allow_pickle=False) as old:
        if int(old['radius'])!=4 or int(old['winding'])!=1:
            raise ValueError('Requires the saved wound R4 mouth')
        x,raw,u=(old[k].copy() for k in ('x','q','field'))
    trial=ContinuousRelease.build(owner,x,False)
    body=trial.b
    is_child=which!='cp10'
    extra=raw[-3:].copy() if is_child else np.zeros(3)
    q=raw[:-3] if is_child else raw
    work=RefinedSolid(body) if is_child else body
    _,_,dep,_=work.lower_bound_load(raw,u[2])
    trial.model.material_load=dep+trial.ext
    total=float(trial.model.evaluate(u)[0])
    flow=trial.measure(u,trial.model.material_load)
    archive=owner.cp['current'] if not is_child else read(root/'r4_local_refinement01/first_step/result.json')['flow']
    output.mkdir(exist_ok=False,parents=True)
    lib=output/'image.so'
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(HERE/'mouth_inverse_material_energy.cpp'),'-o',str(lib)],check=True)
    func=ct.CDLL(str(lib)).image_energy
    func.restype=ct.c_int
    func.argtypes=[ct.POINTER(Data),ct.POINTER(Results),Pointer,ct.c_double,ct.c_int,Pointer,Pointer]
    cells=np.zeros((24,24,6,2,2,2))
    diag=np.zeros(8)
    body.work=lambda data,results:func(data,results,pointer(extra),float(work.seam) if is_child else 0.,
                          int(is_child),pointer(cells),pointer(diag))
    native=body.evaluate(q,u[2],load=False,gradient=False)[0]
    np.save(output/'cell_energy.npy',cells,allow_pickle=False)
    res=cells[...,1]-cells[...,0]
    delta=float(np.sum(res,dtype=np.longdouble))
    if not (diag[0]>0 and diag[3]>0 and diag[4]<0 and np.isfinite(cells).all()):
        raise ArithmeticError('Physical integrand or negative control failed')
    info=dict(state_sha256=expected,state=which,quadrature='Imaged material quadrature with each native reference tensor; epsilon=0.',
              lagrangian_total=total,lagrangian_crystal=native,branchwise_E_Euler=float(np.sum(cells[...,1],dtype=np.longdouble)),
              branchwise_E_Lagr=float(np.sum(cells[...,0],dtype=np.longdouble)),
              summed_cell_delta=delta,total_after_replacing_integral=total+delta,
              max_abs_element_energy_delta=float(np.max(abs(res))),element_energy_hash=digest(output/'cell_energy.npy'),
              active_quadrature_points=int(diag[0]),max_abs_J_roundtrip=diag[1],
              principal_stretch_roundtrip_error_upper_frobenius=diag[2],positive_order_rejection=dict(chi=diag[3],
              J_inverse_perturbed=diag[4],landed_chart=diag[5:].tolist()),
              field_held=True,port_values={k:flow[k] for k in PORTS},
              maximum_archived_port_delta=max(abs(flow[k]-archive[k]) for k in PORTS))
    (output/'result.json').write_text(json.dumps(info,indent=2,allow_nan=False)+'\n')
    print(json.dumps(info,allow_nan=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('output',type=Path);p.add_argument('which',choices=('cp10','post'))
    a=p.parse_args()
    run(a.root.resolve(),a.output.resolve(),a.which)
