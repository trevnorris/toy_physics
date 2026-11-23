"""Check wall-load changes using the angular solid integral and a local nodal variation."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_angular_energy_review import angular_material
from mouth_bonded_front_adjoint import angular_weights,get_material_load
from mouth_bonded_material_motion import Bonded
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_exit_readout import accepted
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import emit,file_sha


def wall_load(state,point,endpoint,component):
    check_guard(endpoint)
    prep,row,x,d0,u0=accepted(state,point)
    with np.load(endpoint) as f:
        if str(f['input_sha'])!=row['sha256'] or not np.array_equal(f['x'],x):
            raise ValueError('Input shape mismatch')
        s=float(f['internal_amplitude']); d1=f['displacement'];u1=f['field']
    n,k=(d0.shape[1]-1)//3,d0.shape[0]-1
    base=Front(state,tuple(prep['field_grid']),1)
    def change(y):
        t=Trial(base,state,y,n,k,False)
        b=Bonded(state,n,k,shape=t.shape,full_depth=True)
        parts=[]
        for d,u,a in ((d0,u0,0.),(d1,u1,s)):
            E=angular_material(b,d,y[5],a)[0]
            t.model.material_load=get_material_load(b,d[b.mask]*b.scale,t.model)
            parts.append([E,t.model.evaluate(u)[0]])
        return np.sum(np.subtract(*reversed(parts)))
    for axis in (component,5):
        for step in (1e-4,5e-5):
            dx=step*np.eye(6)[axis]
            lo,hi=change(x-dx),change(x+dx)
            emit(kind='signed_load_change',checkpoint=row['sha256'],equilibrium=file_sha(endpoint),axis=axis,
                 step=step,delta=(lo-hi)/(2*step),energy_diff_minus=lo,energy_diff_plus=hi)


def node_work(state,failed):
    outcome=json.loads((failed.parent/'outcome.json').read_text())
    limits=json.loads((failed.parent/'limit-validation.json').read_text())
    if not outcome['limitsVerified'] or not limits['verified'] or int(limits['actual']['memory.swap.max']) or outcome['childOutcome']['exitCode']!=1:
        raise ValueError('Expected saved failed guarded trial')
    with np.load(failed) as f:
        x,d,u=f['x'],f['displacement'],f['field'];s=float(f['internal_amplitude'])
        meta=json.loads(str(f['readout']))
        if str(f['state_sha'])!=file_sha(state):raise ValueError('Source mismatch')
    base=Front(state,(13,10,14),1)
    t=Trial(base,state,x,(d.shape[1]-1)//3,d.shape[0]-1,False)
    b,m=t.body,t.model
    g=RegularGridInterpolator((np.linspace(0,1,b.grid.shape[1]),np.linspace(0,1,b.grid.shape[0])),b.grid.T)
    corr=(u[2,m.ids[1,0]].transpose(0,2,1))@angular_weights(m.a)
    live=RegularGridInterpolator((m.rr[0],m.z[1]),corr)
    iz=np.flatnonzero(b.t==.96875).item()
    ix=np.argmin(abs(b.x-.9893179002349632))
    mask=np.any(b.id==iz*d.shape[1]+ix,axis=1)
    def term(a):
        dd=d.copy(); dd[iz,ix,0]+=a
        val=dd.reshape(-1,2)[b.id[mask]]
        delta=np.einsum('qi,qik->qk',b.basis[mask],val)
        r,z=b.r[mask]+delta[:,0],b.w[mask]+delta[:,1]
        B=np.eye(2)+np.stack([np.einsum('qi,qik->qk',v[mask],val) for v in (b.dr,b.dw)],axis=-1)
        hoop=b.h[mask]*r/b.r[mask]
        K=B@b.K[mask]@B.transpose(0,2,1)
        J=b.j[mask]*np.linalg.det(B)*(r/b.r[mask])**2
        W=(np.trace(K,axis1=1,axis2=2)+2*hoop**2-4)/2-np.log(J)+.25*(J-1)**2
        y=r/t.shape.radius(z)[0]
        if min(z/t.shape.length)<=.62:
            raise ArithmeticError('Nodal basis overlaps the angular deformation')
        return (g((y,z/t.shape.length))+live((y,z)))*W
    z0=b.t[iz]*t.shape.length
    deformed=[b.x[ix]*t.shape.radius(z0)[0]+d[iz,ix,0],z0+d[iz,ix,1]]
    reference=term(0)
    weights=b.measure[mask]
    emit(kind='failed_step_site',checkpoint=file_sha(failed),native_sweep=meta['sweep'],
         position=deformed,fraction=deformed[1]/t.shape.length,node=[int(iz),int(ix)])
    for h in (1e-5,1e-6,1e-7,1e-8,1e-9):
        up=weights@(term(h)-reference); dn=weights@(term(-h)-reference)
        emit(kind='local_signed_force',step=h,saved_material_amplitude=s,
             radial_force=(dn-up)/(2*h),left_radial_force=dn/h,right_radial_force=-up/h)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path);p.add_argument('point',type=Path);p.add_argument('--endpoint',type=Path)
    p.add_argument('--peak-coordinate',type=int)
    a=p.parse_args()
    if a.endpoint: wall_load(a.state,a.point,a.endpoint,a.peak_coordinate)
    else: node_work(a.state,a.point)
