"""Signed crystal work on the two registered R4 charts under the saved six shape variations.

The analytic observable is spatial Cauchy stress work, integrated over each
ordered branch. The signed energy probes use physical volume and inverse
Jacobian. Neither operation supplies a single global Eulerian reference field.
"""
import argparse
import ctypes as ct
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_free_shape_solid import Data, Solid, pointer
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def measure(area, folder):
    code=Path(__file__).resolve().with_suffix('.cpp')
    lib=folder/'work.so'
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(code),'-o',str(lib)],check=True)
    kernel=ct.CDLL(str(lib)).inverse_chart_work
    ptr=ct.POINTER(ct.c_double)
    kernel.argtypes=[ct.POINTER(Data),ptr,ct.c_double,ct.c_double,ptr]
    kernel.restype=ct.c_int
    parent=read(area/'r4_melted_survival01/coarse/latest.json')
    step=read(area/'r4_local_refinement01/first_step/result.json')
    inputs=[('checkpoint10',area/'r4_melted_survival01/coarse'/parent['file'],parent['sha256']),
            ('postshape',area/'r4_local_refinement01/first_step/material.npz',step['state_sha256'])]
    prep=area/'supercore_joined_front01/adopt309/state.npz'
    depth_ref=read(area/'r4_whole_mouth01/coarse/endpoint.json')['depth']
    out={}
    mode_names=('log_depth','scaled_log_exit','mean_width','depth_asymmetry','waist_profile','polar_equatorial')
    for name,path,expected in inputs:
        if digest(path)!=expected:raise ValueError(f'Archive changed in {name}')
        with np.load(path,allow_pickle=False) as state:
            if str(state['state_sha'])!=digest(prep) or int(state['winding'])!=1 or int(state['radius'])!=4:
                raise ValueError('Incorrect reference preparation')
            if 'shape_params' in state:
                params, extra=state['shape_params'],float(state['shape_extra'])
            else:
                # Shape is constant over the isolated refined crystal step.
                old=area/'r4_feasible_descent01/coarse/stalled_00011.npz'
                if digest(old)!=step['source_sha256']:raise ValueError('Prior shape changed')
                with np.load(old,allow_pickle=False) as first:
                    if not np.array_equal(first['x'],state['x']):raise ValueError('Different restart shape')
                    params, extra=first['shape_params'],float(first['shape_extra'])
            epsilon=float(state['x'][5])
            def body_at(delta):
                shape=Meridian(4.,params+np.array([delta,0,0,0]),extra=extra)
                model=axes_model(shape,depth_ref)
                body=Solid(prep,shape,epsilon,8,24,6,model,degree=3,subdivisions=2)
                body.melted_labels=True
                body.data.exact_melt_support=1
                body.buffers['base'].fill(0.)
                return body
            body=body_at(0.)
            offs=np.ascontiguousarray(state['solid'])
            field=state['field'][2].copy()
            correction=np.ascontiguousarray(field[body.ids])
            surplus=np.ascontiguousarray(state['surplus'] if 'surplus' in state else np.zeros(3))
            seam=read(area/'r4_local_refinement01/audit/preservation.json')['quadrature_seam_angle']
            body.data.offset,body.data.correction=pointer(offs),pointer(correction)
            values=np.zeros(31)
            h=1.e-4
            rc=kernel(ct.byref(body.data),pointer(surplus),seam,h,pointer(values))
            if rc:raise ArithmeticError(f'Invalid order/volume in chart work: {rc}')
            e=values[7:].reshape(6,4)
            analytic=-values[1:7]
            full=(e[:,0]-e[:,1])/(2*h)
            half=(e[:,2]-e[:,3])/h
            # At fixed fractional crystal coordinates only depth also moves
            # spatial order knots and integration cuts in the production probes.
            # Evaluate those changes separately using the original native solid.
            from mouth_r4_refined_solid import RefinedSolid
            depth=[]
            for dx in (-h,h,-h/2,h/2):
                bt=body_at(dx)
                if np.any(surplus):
                    obj=RefinedSolid(bt,folder/'compiled')
                    q=np.r_[bt.pack(offs),surplus]
                else:
                    obj=bt
                    q=bt.pack(offs)
                depth.append(obj.evaluate(q,field,load=False,gradient=False)[0])
            dt_1,dt_2=(depth[0]-depth[1])/(2*h),(depth[2]-depth[3])/h
            r=dict(state_sha256=expected,modes=mode_names,branchwise_E_crystal=values[0],
                analytic_fixed_spatial_knots=analytic.tolist(),signed_fixed_knots_minus_plus=e.tolist(),
                signed_fixed_knots_h=full.tolist(),signed_fixed_knots_half=half.tolist(),
                residual_analytic_minus_h=(analytic-full).tolist(),residual_analytic_minus_half=(analytic-half).tolist(),
                production_signed_depth_energies=depth,production_depth_h=dt_1,production_depth_half=dt_2,
                depth_order_knot_and_quadrature_work_h=dt_1-full[0],
                depth_order_knot_and_quadrature_work_half=dt_2-half[0])
            if name=='checkpoint10':
                total=np.asarray(parent['shape_force'])
                material=full.copy();material[0]=dt_1
                r.update(saved_total=total.tolist(),crystal_at_matched_kinematics=material.tolist(),
                         saved_total_minus_crystal=(total-material).tolist())
            out[name]=r
            print(json.dumps({'state':name,'analysis':r}),flush=True)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('folder',type=Path)
    args=p.parse_args()
    folder=args.folder.resolve()
    folder.mkdir(exist_ok=True,parents=True)
    if (folder/'work.json').exists():raise FileExistsError('Retain existing work')
    d=measure(args.area.resolve(),folder)
    (folder/'work.json').write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
