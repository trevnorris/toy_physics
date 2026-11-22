"""One-sided conserved-wave energy slopes at individually reopened zero-order nodes."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_fields import BIAS, CoupledLoads
from mouth_live_waist_order import Order, state
from mouth_sealed_shape_relax import Meridian


def reopen(cells, weight, wave, q, v, active):
    nodes=wave.ids[cells]
    vec=np.r_[v,0j][nodes]
    local_form=weight[:,None,None]*q[cells]
    direct=np.einsum('ci,cij,cj->',vec,local_form,vec)
    new=np.unique(nodes[(nodes>=0)&(~np.r_[active,True][nodes])])
    if not len(new):
        return direct,0j,new,np.empty(0,complex),1.
    form=np.zeros((len(new),len(new)),complex)
    cross=np.zeros(len(new),complex)
    for ids,mat,values in zip(nodes,local_form,vec):
        indices=np.flatnonzero(np.isin(ids,new))
        if len(indices):
            ix=np.searchsorted(new,ids[indices])
            form[np.ix_(ix,ix)]+=mat[np.ix_(indices,indices)]
            cross[ix]+=(mat@values)[indices]
    condition=np.linalg.cond(form)
    if condition>1e10:
        raise ArithmeticError('Near-degenerate opening of previously silent shear nodes')
    extension=-np.linalg.solve(form,cross)
    return direct,cross@extension,new,extension,condition


def test(point,n,m,file):
    saved=np.load(file)
    chi=saved['order']
    omega=complex(saved['pole'])
    v=saved['vector']
    shape=Meridian(1,state(1,point)['coordinates'])
    order=Order(shape,n)
    loads=ActiveWave.install(CoupledLoads(shape,n,2*n,field_mesh=m,wave_seed=(omega,v),
                                          target_frequency=omega.real))
    wave=loads.wave
    base_order,go=order.energy_gradient(chi)
    ce,gc,_=loads.carrier.response(chi,BIAS)
    base=base_order+ce+omega.real
    weights=wave.project@chi.T.ravel()
    mass=wave.cut_m+wave.matrix(weights[:,None,None]*wave.mc)
    active=mass.diagonal()>0
    if np.any(v[~active] != 0):
        raise ArithmeticError('Saved vector is nonzero on a structurally absent degree of freedom')
    boundary=v[wave.port.ids]
    dw=2e-6
    slope=(wave.port.block(omega+dw)-wave.port.block(omega-dw))/(2*dw)
    den=2*omega*(v@(mass@v))+boundary@(slope@boundary)
    proj=wave.project.tocsc()
    form=wave.kc-omega**2*wave.mc
    rows=[]
    for i,j in np.argwhere((chi==0)&order.free):
        column=j*(2*n+1)+i
        sl=slice(proj.indptr[column],proj.indptr[column+1])
        direct,correction,ids,extension,condition=reopen(proj.indices[sl],proj.data[sl],wave,form,v,active)
        material=go[i,j]+gc[i,j]
        change=(correction/den).real
        rows.append(dict(i=int(i),j=int(j),slope=material+((direct+correction)/den).real,
                         wave=((direct+correction)/den).real,material=material,correction=change,
                         ids=ids,values=extension,cond=condition))
    if not rows:
        raise ArithmeticError('Saved center has no exact lower-bound nodes')
    lowest=min(rows,key=lambda row:row['slope'])
    most=max(rows,key=lambda row:abs(row['correction']))
    def concise(r):
        return {key:r[key] for key in ('i','j','slope','material','wave','correction')}
    def difference(r):
        direction=np.zeros_like(chi)
        direction[r['i'],r['j']]=1.
        vstart=v.copy()
        vstart[r['ids']]=r['values']
        slopes=[]
        for h in (.002,.001):
            wave.seed=(omega,vstart.copy())
            wave.target=omega.real
            trial=chi+h*direction
            val,_,info=loads.evaluate(trial)
            total=order.energy_gradient(trial)[0]+val
            slopes.append((total-base)/h)
        extrapolated=2*slopes[1]-slopes[0]
        return dict(**concise(r),slope_at_002=slopes[0],slope_at_001=slopes[1],
                    extrapolated_slope=extrapolated,absolute_difference=abs(extrapolated-r['slope']))
    print(json.dumps(dict(snapshot=str(file),sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
          point=point,mesh=n,field_mesh=m,zero_nodes=len(rows),inactive_wave_nodes=int(sum(~active)),
          newly_activated_cases=sum(len(r['ids'])>0 for r in rows),
          worst_extension_condition=max(r['cond'] for r in rows),
          minimum_slope=lowest['slope'],maximum_slope=max(r['slope'] for r in rows),
          negative_nodes=int(sum(r['slope']<0 for r in rows)),below_minus_1e_minus_6=int(sum(r['slope']< -1e-6 for r in rows)),
          near_zero_positive=int(np.sum((chi>0)&(chi<1e-11)&order.free)),
          finite_difference=[difference(r) for r in ([lowest] if most is lowest else [lowest,most])])),
          flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--point',choices=['mid','tight'],required=True)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--field-mesh',type=int,required=True)
    p.add_argument('--snapshot',type=Path,required=True)
    a=p.parse_args()
    test(a.point,a.mesh,a.field_mesh,a.snapshot)
