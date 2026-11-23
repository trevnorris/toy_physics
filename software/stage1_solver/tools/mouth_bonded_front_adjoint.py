"""Use the same material quadrature to send crystalline work to the front and back.

Radius, depth, meridian, background order and external skeletal registration
are held; all nodes inside the ordered core and inner collar can move together.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from mouth_bonded_material_motion import Bonded
from mouth_joined_front_audit import Front
from mouth_joined_front_stationarity import solve
from mouth_poly_shallow_release import CUT
from mouth_polycrystal_material import stable_energy
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import Joined,emit
from mouth_supercore_terminal_order import terminal


def angular_weights(angles):
    a,b=angles[:-1],angles[1:]
    right=(np.sin(b)-np.sin(a))/(b-a)-np.cos(b)
    result=np.zeros(len(angles))
    result[:-1]+=np.cos(a)-np.cos(b)-right
    result[1:]+=right
    return result


def get_material_load(body, q, model):
    u=body.free.copy();u[body.mask]=q/body.scale
    nodes=u.reshape(-1,2)[body.id]
    val=np.einsum('ij,ijk->ik',body.basis,nodes)
    B=np.stack([np.einsum('ij,ijk->ik',d,nodes) for d in (body.dr,body.dw)],axis=2)+np.eye(2)
    r=body.r+val[:,0];w=body.w+val[:,1]
    hrat=r/body.r
    jac=body.j*np.linalg.det(B)*hrat**2
    if np.min(jac)<=0:raise ValueError('Inverted material quadrature state')
    tr=np.trace(B@body.K@B.transpose((0,2,1)),axis1=1,axis2=2)+2*(body.h*hrat)**2
    weight=body.measure*stable_energy(tr,jac)
    rx=model.rr[0];zz=model.z[1];x=r/body.chart.shape.radius(w)[0]
    i=np.clip(np.searchsorted(rx,x)-1,0,len(rx)-2)
    k=np.clip(np.searchsorted(zz,w)-1,0,len(zz)-2)
    a=(x-rx[i])/(rx[i+1]-rx[i]);b=(w-zz[k])/(zz[k+1]-zz[k])
    if np.min([a.min(),b.min()]) < -1e-9 or max(a.max(),b.max())>1+1e-9:
        raise ValueError('Material quadrature outside ordered chart')
    plane=np.zeros((len(rx),len(zz)))
    for di,si in ((0,1-a),(1,a)):
        for dk,sk in ((0,1-b),(1,b)):
            np.add.at(plane,(i+di,k+dk),weight*si*sk)
    load=plane[:,None,:]*angular_weights(model.a)[None,:,None]
    return np.bincount(model.ids[1,0].ravel(),weights=load.ravel(),minlength=model.size)


class Coupled(Joined):
    def evaluate(self,u,read=False):
        if read:
            en,chi,n,j,v,T=super().evaluate(u,True)
            return en+self.material_load@u[2],chi,n,j,v,T
        en,grad=super().evaluate(u)
        grad[2]+=self.material_load
        return en+self.material_load@u[2],grad


def set_front(body, model, field):
    correction=field[2,model.ids[1,0]]
    mean=np.einsum('raj,a->jr',correction,angular_weights(model.a))
    body.front=(model.z[1]/model.shape.length,model.rr[0],mean)


def run(a):
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    for path in (a.state,a.motion,a.flow):check_guard(path)
    inp=next(row for row in map(json.loads,(a.flow.parent/'stdout').read_text().splitlines())
             if row.get('kind')=='preparation')
    if sha(a.state)!=inp['source_sha256']:raise ValueError('Unmatched reservoir fields')
    base=Front(a.state,tuple(inp['grid']),inp['winding'])
    with np.load(a.motion) as f:
        disp=f['displacement'];n=(disp.shape[1]-1)//3;k=disp.shape[0]-1
        if f['source_sha256'].item()!=sha(a.state):raise ValueError('Unmatched material state')
    b=Bonded(a.state,n,k,a.flow)
    q=disp[b.mask]*b.scale
    p=b.chart
    nr=base.order_mesh;t=np.linspace(0,1,2*nr+1);x=np.linspace(0,1,nr+1)
    out=p.evaluate(base.labels,(t,x,base.c),((1-np.cos(np.pi*t))/2,.1*x+.9*x*x,base.outer))['outer_source']
    profile=terminal(base.sh,nr,CUT,lambda wa,xo:out)
    m=Coupled(base.sh,base.c,base.N,base.grid,terminal_profile=profile)
    m.warm_path=None
    with np.load(a.flow) as f:
        u=f['fields' if not a.rest else 'rest']
    emit(kind='preparation',source_sha256=sha(a.state),displacement_sha256=sha(a.motion),flow_sha256=sha(a.flow),
         winding=base.N,grid=base.grid,rest=a.rest,material_mesh=[3*n,k],geometry=[base.sh.length,base.sh.exit])
    def objective(v):
        try:return b.energy(v)
        except ValueError:return 1e10+v@v,2*v
    for it in range(1,61):
        set_front(b,m,u)
        # Analytic adjoint test at actual overlap, signed in the full stored motion.
        def transfer(v):
            set_front(b,m,u);a=b.energy(v)[0];b.front=None
            off=b.energy(v)[0];set_front(b,m,u)
            return a-off,get_material_load(b,v,m)@u[2]
        both=np.array(transfer(q))
        h=5e-5
        ds=(np.array(transfer(q*(1+h)))-transfer(q*(1-h)))/(2*h)
        if abs(both[0]-both[1])>5e-10 or abs(ds[0]-ds[1])>4e-8:
            raise ArithmeticError(f'Work transfer failed {both} {ds}')
        previous=q.copy()
        opt=minimize(objective,q,jac=True,method='L-BFGS-B',
                     options=dict(gtol=3e-9,ftol=1e-15,maxiter=2400,maxls=60))
        q=opt.x.copy()
        m.material_load=get_material_load(b,q,m)
        m.seed=u.copy()
        # Warm prescribed values are supplied by the recorded preparation.
        f,info=solve(m,a.rest)
        set_front(b,m,f)
        E,gg=b.energy(q)
        res=np.max(abs(b.native[b.mask]))
        upd=np.max(abs(f[2]-u[2]))
        maxstep=np.max(abs(q-previous)/b.scale)
        mE=m.evaluate(f)[0]
        base_front=b.front;b.front=None
        Eref=b.energy(q)[0];b.front=base_front
        emit(kind='joint_sweep',iteration=it,material_nodal_residual=res,field_stationarity=info['stationarity'],
             peak_step=maxstep,order_change=upd,energy_varying_terms=Eref+mE,
             transfer_work=both.tolist(),signed_transfer_derivative=ds.tolist())
        u=f
        displacement=np.zeros(b.mask.shape);displacement[b.mask]=q/b.scale
        pth=a.folder/'coupled.npz';tmp=pth.with_suffix('.partial.npz')
        np.savez_compressed(tmp,field=u,displacement=displacement,radial_fraction=b.x,depth_fraction=b.t)
        tmp.replace(pth)
        # Independent absolute tolerances on the unscaled weak material force
        # and on the native front/carrier degrees of freedom.
        if res<3e-5 and info['stationarity']<2e-7 and upd<2e-5:
            result=m.measure(u)
            if result['port_error']>3e-5:raise ArithmeticError('Unclosed bath currents')
            emit(kind='held_equilibrium',sha256=sha(pth),material_nodal_residual=res,
                 field_stationarity=info['stationarity'],**b.energy(q,True),**result)
            return
    raise ArithmeticError('Coupled material and front did not jointly settle on the held shape')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','motion','flow','folder'):p.add_argument(key,type=Path)
    p.add_argument('--rest',action='store_true')
    run(p.parse_args())
