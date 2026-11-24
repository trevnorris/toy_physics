"""Static two-metric crystal integral and shape work on the reconstructed spatial-label mesh."""
import argparse
import ctypes as ct
import json
import subprocess
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_free_shape_solid import Solid, Data, pointer, Pointer
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian
from mouth_continuous_depth_mesh import axes_for_depth

MODES=['log_depth','scaled_log_exit','mean_width','depth_asymmetry','waist_profile','polar_equatorial']


def load_context(sc,labels):
    prep=sc/'supercore_joined_front01/adopt309/state.npz'
    depth_ref=read(sc/'r4_whole_mouth01/coarse/endpoint.json')['depth']
    cp=sc/'r4_melted_survival01/coarse/state_00010.npz'
    assert digest(cp)==CHECKPOINT_SHA
    with np.load(cp) as data:
        assert int(data['radius'])==4 and int(data['winding'])==1 and str(data['state_sha'])==digest(prep)
        shape_params=data['shape_params']; extra=float(data['shape_extra']); polar=float(data['x'][5])
        field_order=data['field'][2]
    sh=Meridian(4.,shape_params,extra=extra)
    b=Solid(prep,sh,polar,8,24,6,axes_model(sh,depth_ref),degree=3,subdivisions=2)
    with np.load(labels) as f:
        xi=np.ascontiguousarray(f['labels'])
        chi=np.ascontiguousarray(f['order'])
        for key in ('x','t','a'):
            if not np.array_equal(f[key],b.buffers['field_'+key]): raise ValueError('Changed physical axis '+key)
        if not np.array_equal(chi,field_order[b.ids].transpose(0,2,1)): raise ValueError('Spatial order differs from checkpoint')
    if xi.shape!=chi.shape+(3,):raise ValueError('Shape of common label field')
    return prep,depth_ref,b,shape_params,extra,polar,xi,chi


def signed_geometry(context,mode=0,delta=0.):
    prep,dref,b,parameters,extra,polar,*_=context
    if delta==0:return b
    p=parameters.copy()
    if mode<4:p[mode]+=delta*(1. if mode==0 else .25)
    elif mode==4:extra+=.25*delta
    else:polar+=delta
    s=Meridian(4.,p,extra=extra)
    return Solid(prep,s,polar,8,24,6,axes_model(s,dref),degree=3,subdivisions=2)


def integrate(kernel,bodies,xi,chi,degree):
    b=bodies[0]
    data=(Data*len(bodies))(*(v.data for v in bodies))
    foo,w=np.polynomial.legendre.leggauss(degree)
    pts=np.ascontiguousarray((1+foo)*.5);weights=np.ascontiguousarray(w*.5)
    t=b.buffers['field_t']; L=b.data.length
    tdot=np.where(t*L<L-.12,(t*L/(L-.12))**25,1.)-t
    totals=np.zeros(31);orientation=np.zeros(17);neg=np.zeros(8);cells=np.zeros((24,24,6,2))
    rc=kernel(data,len(bodies),pointer(np.ascontiguousarray(tdot)),pointer(xi),pointer(chi),
              pointer(pts),pointer(weights),degree,pointer(totals),pointer(orientation),pointer(neg),pointer(cells))
    if rc:raise ArithmeticError(f'Spatial crystal integration failed ({rc})')
    return totals,orientation,neg,cells


def main(sc,labels,folder):
    folder.mkdir(exist_ok=False,parents=True)
    lib=folder/'crystal.so'
    subprocess.run(['g++','-O3','-fPIC','-std=c++17','-shared',str(Path(__file__).with_suffix('.cpp')),'-o',str(lib)],check=True)
    kernel=ct.CDLL(str(lib)).cp10_spatial_integral
    kernel.argtypes=[ct.POINTER(Data),ct.c_int]+[Pointer]*5+[ct.c_int]+[Pointer]*4
    kernel.restype=ct.c_int
    ctx=load_context(sc,labels)
    _,_,b,_,_,_,xi,chi=ctx
    bodies=[b]
    for mode in range(6):
        for h in (-1.e-4,1.e-4,-5.e-5,5.e-5):
            bodies.append(signed_geometry(ctx,mode,h))
    archive=sc/'r4_inverse_material01/cp10'
    proof=read(archive/'result.json')
    assert digest(archive/'cell_energy.npy')==proof['element_energy_hash'] and proof['state_sha256']==CHECKPOINT_SHA
    cellref=np.load(archive/'cell_energy.npy')[...,0].sum(axis=-1)
    assert abs(cellref.sum()-proof['lagrangian_crystal'])<1e-12
    report=dict(checkpoint=CHECKPOINT_SHA,labels=digest(labels),semantics='Transport spatial nodal order and labels together; reference phase fixed by label, no gauge contribution here.',levels=[])
    for n in (3,6,9):
        totals,orientation,neg,cells=integrate(kernel,bodies if n==6 else bodies[:1],xi,chi,n)
        defect=cells-cellref
        index=np.unravel_index(np.argmax(np.abs(defect)),defect.shape)
        r=dict(degree=n,order_positive_nonpositive_Jinv_samples=int(orientation[13]),melt_nonpositive_Jinv_samples=int(orientation[14]),
                 min_active_Jinverse_sampled=orientation[15],min_all_inverse_chart_det_sampled=orientation[16],
                 negative_control_core_collar=neg.reshape(2,4).tolist())
        if orientation[0]:
            r['worst_inverted']=dict(label=orientation[1:4].tolist(),chi=orientation[4],Jinverse=orientation[5],
                              physical_chart=orientation[6:9].tolist(),spatial_cell=(orientation[9:12].astype(int)).tolist(),core=bool(orientation[12]))
        if orientation[13]==0:
            r.update(crystal_energy=totals[0],delta_C=totals[0]-proof['lagrangian_crystal'],
                     max_abs_element_energy_delta=float(np.max(abs(defect))),worst_material_element=list(map(int,index)),
                     worst_element_energies=dict(native=cellref[index],spatial=cells[index]))
            if n==6:
                es=totals[7:].reshape(6,4)
                h=(es[:,0]-es[:,1])/2.e-4
                h2=(es[:,2]-es[:,3])/1.e-4
                r.update(modes=MODES,analytic=totals[1:7].tolist(),signed_energies=es.tolist(),
                         forces_h=h.tolist(),forces_half=h2.tolist(),residual_analytic_minus_half=(totals[1:7]-h2).tolist())
            np.save(folder/f'cell_energy_{n}.npy',cells)
        report['levels'].append(r)
        (folder/'result.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(json.dumps(r),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('sc',type=Path);p.add_argument('labels',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args()
    main(a.sc.resolve(),a.labels.resolve(),a.out.resolve())
