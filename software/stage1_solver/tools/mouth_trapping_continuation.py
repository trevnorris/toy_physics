"""Look for torsional poles missed by clamped seeds, and resolve the physical flux partition."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse.linalg import eigsh

from mouth_trapping_modes import Sector, form, pole
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_exit_witness import inputs
from mouth_strained_wave import cell_tangents


def free_seeds(op):
    local=SimpleNamespace(mesh=op.mesh,ids=op.ids,kc=op.mc,mc=op.mc)
    k,_=cell_tangents(local,op.shape,op.coeff,degree=op.degree,core_reference=op.reference)
    block=form(op.ids,k*op.weights[:,None,None],len(op.live))[op.live][:,op.live]
    in_cavity=(op.inside.diagonal()>0) & (op.labels==op.connected[0])
    K,M=block[in_cavity][:,in_cavity],op.inside[in_cavity][:,in_cavity]
    vals,vecs=eigsh(K,k=4,M=M,sigma=-.01,which='LM',v0=np.ones(K.shape[0]))
    starts=[]
    for index in np.argsort(vals):
        # Rigid rotation of the isolated degree-one pocket is a zero-frequency
        # field. It is outside the set of oscillations in this test.
        if vals[index]<1e-8:
            continue
        v=np.zeros(op.K.shape[0],complex)
        v[in_cavity]=vecs[:,index]
        starts.append((np.sqrt(vals[index])-1e-4j,v))
    return starts[:3]


def radiated(op,om,v):
    # Curl and symmetric strain differ by the same outer boundary term in
    # K and the DtN. The pole is unchanged; stored physical strain energy
    # must use Ksym to compare the power to the energy inside the cut.
    edge=op.port.project@v[op.port.ids]
    radius=op.port.radius
    kin=abs(om)**2*np.vdot(v,op.M@v).real
    strain=np.vdot(v,op.K@v).real-2*radius**3*np.vdot(edge,edge).real
    if strain<=0:
        raise ArithmeticError('Nonpositive physical strain work')
    power=.5*om.real*radius**4*op.port.traction(om.real).imag*abs(edge)**2
    open_ids=op.port.cutoff<om.real
    return dict(physical_power_per_energy=4*power[open_ids].sum()/(kin+strain),
                physical_channels=(4*power[open_ids]/(kin+strain)).tolist())


def family(radius,n,degree):
    trial=Trial(radius,n,bonded_core=True)
    base=ROOT/'_scratch/s11c/mouth_exit_release01'/f'r{radius}n{n}'/'endpoint.npz'
    sh,c,(freq,archived),digest,status,*_=inputs(base,radius,n,trial)
    if status!='balanced':
        raise ValueError('Unexpected frozen background')
    ops={}
    for h in (-.01,-.005,0.,.005,.01):
        shifts=sh.copy(); shifts[1]+=h
        ops[h]=Sector(trial,shifts,c,degree)
    op=ops[0.]
    if not all(np.array_equal(val.live,op.live) for val in ops.values()):
        raise ArithmeticError('Active nodes changed under signed exit variation')
    emit(kind='setup',radius=radius,mesh=n,degree=degree,sha256=digest,exit=op.shape.exit)
    starts=free_seeds(op)
    starts=[(f'free_{i}',seed) for i,seed in enumerate(starts,1)]
    if degree!=4:
        full=np.zeros(len(op.live),complex)
        full[op.allpoints[:,0]>0]=archived
        incoming=full[op.live]
        free_edge=op.port.project@incoming[op.port.ids]
        k=np.vdot(incoming,op.K@incoming).real-2*op.port.radius**3*np.vdot(free_edge,free_edge).real
        guess=np.sqrt(k/np.vdot(incoming,op.M@incoming).real)
        starts.append(('degree_four_profile',(guess-1e-4j,incoming)))
    for label,seed in starts:
        try:
            middle,stats=pole(op,seed)
        except (ArithmeticError,ValueError,RuntimeError) as exc:
            emit(kind='unresolved_start',label=label,reason=str(exc)); continue
        emit(kind='mode',label=label,start=seed[0].real,**stats,**radiated(op,*middle))
        ends={}
        for h in (-.01,-.005,.005,.01):
            try:
                continued,data=pole(ops[h],middle)
                if data['seed_overlap']<.96:
                    raise ArithmeticError('Insufficient shape-continuation overlap')
                ends[h]=continued[0]
                emit(kind='signed',label=label,log_exit_step=h,exit=ops[h].shape.exit,
                     **data,**radiated(ops[h],*continued))
            except (ValueError,ArithmeticError,RuntimeError) as exc:
                emit(kind='unresolved_signed',label=label,log_exit_step=h,reason=str(exc))
        for h in (.005,.01):
            if h in ends and -h in ends:
                den=ops[h].shape.exit-ops[-h].shape.exit
                emit(kind='exit_work',label=label,step=h,domega_dexit=(ends[h]-ends[-h]).real/den,
                     dloss_dexit=-2*(ends[h]-ends[-h]).imag/den)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name,choices in [('radius',(1,2)),('mesh',(48,72)),('degree',range(1,7))]:
        p.add_argument('--'+name,type=int,choices=choices,required=True)
    args=p.parse_args()
    family(args.radius,args.mesh,args.degree)
