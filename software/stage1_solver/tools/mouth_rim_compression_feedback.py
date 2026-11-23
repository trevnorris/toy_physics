"""Compare the same supported last-ring patch on legacy and consistent R4 fields.

Only the interior of the previous enlarged patch can deform. Return its changed
registered crystal load to the driven coherent fields; test the updated solid
force and field stationarity on the same saved iterate.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scipy.optimize._numdiff import group_columns

from mouth_bonded_front_adjoint import angular_weights
from mouth_last_rings_enriched import Patch
from mouth_rim_consistent_fields import from_archive, checked_pair
from mouth_rim_shared_crystal import TerminalQuadrature, assemble, deposit
from mouth_supercore_joined_front import emit,file_sha,check_guard

OFFSETS=[(0,0),(-1,0),(0,-1),(-1,-1)]


class LocalMaterial(TerminalQuadrature):
    def _radial_lines(self,t):
        old=super()._radial_lines(t)
        # New solid cells/diagonals are exactly nested in the old ones.
        step=self.tt[1]/self.mesh
        fraction=(t/step)%1
        xs=np.linspace(self.xt[30],self.xt[33],3*self.mesh+1)
        return np.unique(np.r_[old,xs,xs[:-1]+np.diff(xs)*(1-fraction)])


def collect(quad,build=False):
    dx=quad.xt[-1]-quad.xt[-2]
    dt=quad.tt[1]
    xt,tt=np.meshgrid(quad.xt,quad.tt)
    quad.first=29
    quad.tstart=quad.tt[29]
    kept={}
    original=np.zeros((len(quad.ax),len(quad.z)))
    for v in quad.blocks():
        ids=v['ids']
        x=(xt.ravel()[ids]*v['bary']).sum(axis=1)
        t=(tt.ravel()[ids]*v['bary']).sum(axis=1)
        A=(x-quad.xt[32])/dx
        T=(t-quad.tt[31])/dt
        active=np.zeros(len(x),bool)
        for i,j in OFFSETS:
            active |= np.maximum.reduce([abs(A-i),abs(T-j),abs(A+T-i-j)])<1.
        if not active.any():continue
        def get(key):return v[key][active]
        deposit(original,quad.ax,quad.z,get('x'),get('z'),get('measure')*get('W'))
        if not build:continue
        r0=x[active]*quad.sh.radius(t[active]*quad.L)[0]
        hrat=get('r')/r0
        B=get('B')
        inv=np.linalg.inv(B)
        vals=dict(rbase=r0,zbase=t[active]*quad.L,r=get('r'),z=get('z'),B=B,
                  j0=get('J')/(np.linalg.det(B)*hrat**2),refhoop=get('hoop')/hrat,
                  K=inv@get('K')@inv.transpose(0,2,1),w=get('measure'),
                  phi=np.maximum(0,1-np.maximum.reduce([abs(A[active]),abs(T[active]),abs(A[active]+T[active])])))
        for k,arr in vals.items():kept.setdefault(k,[]).append(arr)
    if not build:return original
    point=SimpleNamespace(**{k:np.concatenate(v) for k,v in kept.items()},
                          sh=quad.sh,vcell=(dx,dt),order=None)
    return original,point


def source_at(patch,U,ax,z):
    ring=patch.ring
    patch.eval(U)
    K,hoop,J=patch.last
    off=patch.V@U.reshape(-1,2)
    rr,zz=ring.r+off[:,0],ring.z+off[:,1]
    W=(np.trace(K,axis1=1,axis2=2)+2*hoop**2-4)/2-np.log(J)+.25*(J-1)**2
    plane=np.zeros((len(ax),len(z)))
    deposit(plane,ax,z,rr/ring.sh.radius(zz)[0],zz,ring.w*W)
    return plane


def stationary_patch(patch,U,metric=None):
    # Recenter the existing positive-metric, admissible material descent.
    from mouth_last_rings_enriched import descend
    if metric is None:
        basis=patch.V.astype(bool)
        pattern=sparse.kron(basis.T@basis,np.ones((2,2),dtype=bool),format='csc')
        colors=group_columns(pattern)
        H=np.zeros((len(U),len(U)))
        h=1e-7
        for color in np.unique(colors):
            step=(colors==color).astype(float)*h
            slopes=(patch.eval(U+step)[1]-patch.eval(U-step)[1])/(2*h)
            for k in np.flatnonzero(step):
                ii=pattern.indices[pattern.indptr[k]:pattern.indptr[k+1]]
                H[ii,k]=slopes[ii]
        eig,V=np.linalg.eigh((H+H.T)/2)
        metric=V/np.sqrt(np.maximum(abs(eig),1.))
    class Offset:
        def eval(self,v):
            return patch.eval(U+v)
    du,n,status=descend(Offset(),metric)
    return U+du,metric,n,status


def run(a):
    checked_pair(a.old_fields)
    base,tr,m,pair,level,sha,d,_=from_archive(a.state,a.wall,a.fields,True)
    check_guard(a.load)
    with np.load(a.fields) as f:
        if f['load_sha'].item()!=file_sha(a.load) or f['prior_fields_sha'].item()!=file_sha(a.old_fields):
            raise ValueError('Mismatched original field or full crystal integration')
    ax,z=m.rr[0],m.z[1]
    with np.load(a.load) as f:
        if not np.array_equal(f['ax'],ax) or not np.array_equal(f['z'],z):raise ValueError('Different load grid')
        parent_plane=f['plane'].copy()
        sub=int(f['subdivisions'])
    # Subtract original contributions to this domain from the original rule,
    # then replace them with their new displacement-mesh integral.
    args=(a.state,a.wall,ax,z,sub)
    old_source=collect(TerminalQuadrature(*args))
    quad=LocalMaterial(*args)
    quad.mesh=a.mesh
    corrected,point=collect(quad,True)
    external=parent_plane-old_source
    w=angular_weights(m.a)
    def set_order(field):
        mean=np.einsum('raj,a->jr',field[2,m.ids[1,0]],w)
        point.order=[quad.inner,(z/quad.L,ax,mean)]
    with np.load(a.old_fields) as f:old=f['wound']
    set_order(old)
    patch=Patch(point,a.mesh,OFFSETS)
    U0=patch.zero.ravel()
    if np.max(abs(source_at(patch,U0,ax,z)-corrected))>1e-10:raise ArithmeticError('Material energy and deposited source differ')
    emit(kind='inputs',checkpoint=file_sha(a.wall),drive=sha,old=file_sha(a.old_fields),quadrature=sub,
         repeated_patch_rule_l1=float(abs(old_source-corrected).sum()),level=level,dofs=len(U0))
    def sample(U,field,name,metric=None):
        set_order(field)
        U,metric,n,reason=stationary_patch(patch,U,metric)
        score=patch.sample(U)
        emit(kind=name,iterations=n,stopping=reason,**score)
        return U,metric,score
    oldU,metric,old_score=sample(U0,old,'legacy_field_same_quadrature')
    del old
    field=pair[0]
    U,metric,new_score=sample(U0,field,'corrected_field_same_quadrature',metric)
    for i in range(20):
        plane=external+source_at(patch,U,ax,z)
        m.material_load=assemble(m,plane)
        grad=np.max(abs(m.norm(m.vars(m.evaluate(field)[1]))))
        if grad>=2e-7:
            field,info=m.solve(rest=False,initial=field)
        final_field=np.max(abs(m.norm(m.vars(m.evaluate(field)[1]))))
        set_order(field)
        score=patch.sample(U)
        np.savez_compressed(a.output,field=field,patch=U,plane=plane,
                            checkpoint_sha=file_sha(a.wall),field_sha=sha,grid=level)
        emit(kind='feedback',sweep=i,field_residual=final_field,prior_field_residual=grad,
             old_min=old_score['min_stretch'],frozen_corrected_min=new_score['min_stretch'],
             plane_l1_change=float(abs(plane-parent_plane).sum()),**score)
        if final_field<2e-7 and score['max_force']<3e-5:
            return
        if final_field>=2e-7:raise ArithmeticError('Unresolved responsive conversion field')
        U,metric,score=sample(U,field,'solid_sweep',metric)
        if score['max_force']>=3e-5:raise ArithmeticError('Unresolved local material force')
    raise ArithmeticError('Local alternating solve did not meet simultaneous force criteria')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','wall','old_fields','fields','load','output'):p.add_argument(key,type=Path)
    p.add_argument('--mesh',type=int,default=8)
    run(p.parse_args())
