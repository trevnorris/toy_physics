"""Measure matched phases on the same relaxed bonded lattice and resolved exit."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import Coupled, get_material_load
from mouth_bonded_material_motion import Bonded
from mouth_joined_angular_work import deformed_work
from mouth_joined_front_audit import Front
from mouth_joined_front_stationarity import solve
from mouth_joined_front_trace import face
from mouth_poly_shallow_release import CUT
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import emit
from mouth_supercore_terminal_order import terminal


def run(a):
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for path in (a.state,a.joint,a.flow,a.wound):check_guard(path)
    def record(path):
        return list(map(json.loads,(path.parent/'stdout').read_text().splitlines()))
    joint=record(a.joint)
    if joint[-1]['kind']!='held_equilibrium' or joint[-1]['sha256']!=sha(a.joint):
        raise ValueError('Unconfirmed coupled source')
    j=next(r for r in joint if r['kind']=='preparation')
    p=next(r for r in record(a.flow) if r['kind']=='preparation')
    if p['source_sha256']!=sha(a.state) or p['source_sha256']!=j['source_sha256'] or j['flow_sha256']!=sha(a.wound):
        raise ValueError('Mismatched geometry or mechanical reference')
    base=Front(a.state,tuple(p['grid']),p['winding'])
    with np.load(a.joint) as jf:
        field=jf['field'];d=jf['displacement']
    b=Bonded(a.state,(d.shape[1]-1)//3,d.shape[0]-1,a.wound)
    q=d[b.mask]*b.scale
    n=base.order_mesh;t=np.linspace(0,1,2*n+1);x=np.linspace(0,1,n+1)
    osrc=b.chart.evaluate(base.labels,(t,x,base.c),((1-np.cos(np.pi*t))/2,.1*x+.9*x*x,base.outer))['outer_source']
    profile=terminal(base.sh,n,CUT,lambda wa,xo:osrc)
    m=Coupled(base.sh,base.c,base.N,base.grid,terminal_profile=profile)
    m.material_load=get_material_load(b,q,m)
    m.warm_path=a.flow
    if base.N==1:
        u=field
        err=float(np.max(abs(m.norm(m.vars(m.evaluate(u)[1])))))
        if err>=2e-7:raise ArithmeticError('Reloaded joint field is not stationary')
    else:
        u,info=solve(m,False)
        err=info['stationarity']
    ur,info=solve(m,True)
    emit(kind='same_lattice',source_sha256=sha(a.state),joint_sha256=sha(a.joint),winding=base.N,
         drive_residual=err,rest_residual=info['stationarity'],**m.measure(u))
    for side in (-1,1):
        wt,n,ma,j,T=face(m,base.c,u,side,16)
        T-=face(m,base.c,ur,side,16)[-1]
        use=(n>.5)&(ma<1.);pos=np.maximum(j,0)
        emit(kind='disk',side=side,dense_share=wt[use]@pos[use]/(wt@pos),signed=wt@j,
             drive_stress=wt@T,dense_drive_stress=wt[use]@T[use])
    for key,sol in [('drive',u),('rest',ur)]:
        for h in (.0004,.0002):
            g=(deformed_work(m,sol,-h)-deformed_work(m,sol,h))/(2*h)
            emit(kind='polar',state=key,h=h,depth_span=[.06,.62],force=g.sum())
    np.savez_compressed(a.save,fields=u,rest=ur)
    emit(kind='fields',sha256=sha(a.save))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','joint','flow','wound','save'):p.add_argument(key,type=Path)
    run(p.parse_args())
