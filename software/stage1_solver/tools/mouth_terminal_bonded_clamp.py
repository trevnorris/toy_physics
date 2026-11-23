"""Complement the clamped-disk probe with an admissible virtual motion of both banks.

Both the field and solid complementary motions are their common enclosing
translation minus their own disk displacement. They coincide on each sealed
wall tracing the terminal region and on the disk. The primary wall stays fixed.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_ring_solid import interpolation
from mouth_last_rings_fluid import GRIDS,read_source
from mouth_rim_consistent_fields import checked_pair
from mouth_rim_shared_crystal import TerminalQuadrature,moving_energy
from mouth_terminal_disk_work import DiskMotion,STEPS
from mouth_supercore_joined_front import file_sha


def ramp(bounds,z):
    lo,lmid,hmid,hi=bounds
    return (np.clip((z-lo)/(lmid-lo),0,1)*np.clip((hi-z)/(hi-hmid),0,1),
            ((z>lo)&(z<lmid))/(lmid-lo)-((z>hmid)&(z<hi))/(hi-hmid))


def run(a):
    checked_pair(a.fields)
    c=json.loads(a.clamp.read_text())
    disk=json.loads(a.disk.read_text())
    x,_,_,_=read_source(a.state,a.wall)
    base=Front(a.state,GRIDS['native'],1)
    trial=Trial(base,a.state,x,11,32,False)
    motion=DiskMotion(trial.model)
    with np.load(a.fields) as f:
        if not (c['fields']==disk['fields']==file_sha(a.fields)):raise ValueError('Fields differ')
        # Two controls are stored in the same topology. Only the order field
        # on the duct is needed to evaluate the crystal part of the work.
        from mouth_last_ring_surface_review import FieldChart
        chart=FieldChart(trial.shape,base.c,1,f['grid'],precondition=False,terminal_profile=trial.model.terminal_profile)
        idx=chart.ids[1,0]
        delta=np.einsum('raj,a->jr',(f['wound'][2]-f['rest'][2])[idx],angular_weights(f['angles']))
        quad=TerminalQuadrature(a.state,a.wall,f['inside_x'],f['duct_z'],64)
    readings={n:np.zeros(6) for n in ('near','far')}
    max_error=0.
    for v in quad.blocks():
        r,z=v['r'],v['z']
        chi,cx,cz=interpolation(quad.z,quad.ax,delta,z,v['x'])
        R,slope=quad.sh.radius(z)
        cz-=cx*v['x']*slope/R
        hat=v['ids']//len(quad.xt)==len(quad.tt)-1
        h=(hat*v['bary']).sum(axis=1)
        dh=np.stack([(hat*v[k]).sum(axis=1) for k in ('gx','gt')],axis=1)
        dh=np.einsum('qi,qij->qj',dh,np.linalg.inv(v['B']))
        g,gr,gz=motion.profile(r,z)
        K,J,W,wt=v['K'],v['J'],v['W'],v['measure']
        for name,values in readings.items():
            A,Az=ramp(motion.ramps[name],z)
            us=A-h*g
            sr=-dh[:,0]*g-h*gr
            sz=Az-dh[:,1]*g-h*gz
            uf=A*(1-g)
            stress=wt@(chi*(K[:,0,1]*sr+(K[:,1,1]-1+.5*J*(J-1))*sz))
            relative=wt@(W*cz*(us-uf))
            values[0]+=stress
            values[1]+=relative
            for i,step in enumerate(STEPS):
                energies=[]
                for s in (-step,step):
                    target=z+s*us
                    pre=z.copy()
                    for _ in range(3):
                        Av,Ap=ramp(motion.ramps[name],pre)
                        gv,_,gp=motion.profile(r,pre)
                        pre-=(pre+s*Av*(1-gv)-target)/(1+s*(Ap*(1-gv)-Av*gp))
                    xp=r/quad.sh.radius(pre)[0]
                    max_error=max(max_error,float(np.max(xp-1)),float(np.max(pre-quad.L)))
                    q=interpolation(quad.z,quad.ax,delta,pre,xp)[0]
                    en=moving_energy(K,J,v['hoop'],sr,sz,s)
                    energies.append(wt@(q*en-chi*W))
                values[2+i]+=(energies[1]-energies[0])/(2*step)
            values[4]+=wt@(chi*(K[:,1,1]-1+.5*J*(J-1))*Az)
    if max_error>1e-12:raise ValueError(f'Bonded bank variation outside field chart {max_error}')
    result={}
    for name,v in readings.items():
        prev=c['measured'][name]
        comp=v[0]+v[1]+prev['complement_fluid']
        total=disk['reactions'][name]['physical_actuator']+comp
        external=prev['external_virtual']
        if max(abs(v[2:4]-sum(v[:2])))>5e-5 or abs(total-external-v[4])>3e-5:
            raise ArithmeticError('Unresolved complementary signed-work test')
        result[name]=dict(solid_complement=sum(v[:2]),solid_signed=v[2:4].tolist(),
                          complement_total=comp,disk=disk['reactions'][name]['physical_actuator'],
                          sum_of_virtual_works=total,common_translation=external,common_solid=v[4],
                          strong_wall_banks=prev['sloping_banks'],closed=prev['measured_closed'],
                          weighted_disk_plus_banks_and_complement=total+prev['sloping_banks'],
                          boundary_identity_error=total+prev['sloping_banks']-prev['measured_closed'])
    return dict(fields=c['fields'],common_motion='The solid and the field follow the same enclosed translation when the disk and its welded complementary virtual motion are added.',max_chart_extrapolation=max_error,measurements=result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('state','wall','fields','disk','clamp','output'):p.add_argument(n,type=Path)
    a=p.parse_args()
    a.output.write_text(json.dumps(run(a),indent=2,allow_nan=False)+'\n')
