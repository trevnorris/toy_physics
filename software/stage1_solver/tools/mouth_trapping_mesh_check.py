"""Follow the selected degree-six pole under wave-only refinement of a fixed n=72 material."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import Trial, emit
from mouth_shallow_puncture_fields import PORT
from mouth_trapping_modes import Sector, closed_seeds, pole


def project(previous, new, vector):
    full = np.zeros(len(previous.allpoints), complex)
    full[previous.live] = vector
    sampled = np.empty(len(new.points), complex)
    slab = previous.allpoints[:, 1] <= 0.
    r, w = previous.allpoints[slab].T
    rs, ws = np.unique(r), np.unique(w)
    if len(rs)*len(ws)!=len(r):
        raise ValueError('Unstructured feed grid')
    arr = np.empty((len(rs),len(ws)),complex)
    arr[np.searchsorted(rs,r),np.searchsorted(ws,w)] = full[slab]
    out = new.points[:, 1] <= 0
    sampled[out] = RegularGridInterpolator((rs,ws),arr)(new.points[out])
    nt = len(previous.ids)//previous.mesh
    cells = previous.ids.reshape(previous.mesh,nt,4)
    cavity = np.empty((previous.mesh+1,nt+1),complex)
    for a,b,i in ((0,0,0),(1,0,1),(0,1,2),(1,1,3)):
        cavity[a:a+previous.mesh,b:b+nt] = full[cells[:,:,i]]
    interp = RegularGridInterpolator((np.linspace(0,1,previous.mesh+1),np.linspace(0,1,nt+1)),cavity)
    r, w = new.points[~out].T
    t = w/new.shape.length
    x = r/new.shape.profile(t)[0]
    if max(np.max(t),np.max(x))>1+1e-8:
        raise ValueError('Interpolated point outside material chart')
    sampled[~out]=interp(np.clip(np.c_[x,t],0.,1.))
    return sampled


def row(radius,m,h,stats):
    omega=stats['frequency']
    return dict(radius=radius,wave_mesh=m,log_exit_step=h,frequency_minus_2pi=omega-2*np.pi,
                frequency_minus_FE_n2=omega-stats['depth_thresholds'][2],
                unit_slab_n2_continuum=2*np.pi,**stats)


def run(radius, meshes):
    n=72
    tr=Trial(radius,n,bonded_core=True)
    path=Path(f'_scratch/s11c/mouth_exit_release01/r{radius}n72/endpoint.npz')
    shifts,chi,_,digest,status,*_=inputs(path,radius,n,tr)
    if status!='balanced' or tr.mesh!=20:
        raise ValueError('Different native control')
    folder=Path(f'_scratch/s11c/mouth_trapping_modes01/r{radius}n72l6')
    guard=json.loads((folder/'outcome.json').read_bytes())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
        raise ValueError('Native degree-six spectrum guard failed')
    rows=[json.loads(z) for z in (folder/'stdout').read_text().splitlines()]
    if rows[0]['sha256']!=digest:
        raise ValueError('Different native geometry')
    known,=[s for s in rows if s['kind']=='mode' and s['label']=='clamped_1']
    previous=Sector(tr,shifts,chi,6)
    seed, observed=pole(previous,closed_seeds(previous,3)[0])
    if abs(observed['frequency']-known['frequency'])>1e-8 or abs(observed['loss']-known['loss'])>1e-9:
        raise ArithmeticError('Failed to recover saved degree-six resonance')
    emit(kind='reconstructed',mesh_order=n,endpoint_sha256=digest,**row(radius,20,0.,observed))
    for m in meshes:
        higher=copy.copy(tr)
        higher.mesh=m
        op=Sector(higher,shifts,chi,6)
        seeded=(seed[0],project(previous,op,seed[1]))
        center,stats=pole(op,seeded)
        if stats['seed_overlap']<.90 or stats['detached']:
            emit(kind='unresolved_center',**row(radius,m,0.,stats))
            return
        emit(kind='center',degrees_of_freedom=len(op.points),shear_components=op.ncomp,
             **row(radius,m,0.,stats))
        signed={}
        for h in (-.01,-.005,.005,.01):
            shift=shifts.copy()
            shift[1]+=h
            test=Sector(higher,shift,chi,6)
            if not np.array_equal(test.live,op.live):
                raise ArithmeticError('Frozen order chart changed supported wave nodes')
            try:
                other,s=pole(test,center)
                if s['seed_overlap']<.96 or s['detached']:
                    raise ArithmeticError('Unresolved signed branch')
                signed[h]=(other[0],s['loss'])
                emit(kind='signed',exit=test.shape.exit,**row(radius,m,h,s))
            except (ValueError,RuntimeError,ArithmeticError) as exc:
                emit(kind='unresolved_signed',radius=radius,wave_mesh=m,log_exit_step=h,reason=str(exc))
        for h in (.005,.01):
            if -h in signed and h in signed:
                difference=2*op.shape.exit*np.sinh(h)
                emit(kind='work',radius=radius,wave_mesh=m,h=h,
                     domega_dexit=float((signed[h][0]-signed[-h][0]).real/difference),
                     dloss_dexit=float((signed[h][1]-signed[-h][1])/difference))
        previous,seed=op,center


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,choices=(1,2),required=True)
    parser.add_argument('--wave-meshes',type=int,nargs='+',default=[28,36,48])
    args=parser.parse_args()
    if any(a>=b for a,b in zip([20,*args.wave_meshes[:-1]],args.wave_meshes)):
        raise ValueError('Wave meshes must refine')
    run(args.radius,args.wave_meshes)
