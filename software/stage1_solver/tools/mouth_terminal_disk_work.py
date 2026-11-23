"""Direct axial work of bowing the terminal disk on the rejected one-quantum R4 wall.

Keep the common rim and both sealed wall banks fixed; carry the local reference
order with the chosen field-coordinate motion. Only the last row of solid
elements moves. Compare short, near and far extensions of the same disk trace.
"""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from mouth_bonded_front_adjoint import angular_weights
from mouth_closed_rim_surface_review import Reader
from mouth_disk_pillbox import dot, physical, stresses, transformed
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_ring_solid import interpolation
from mouth_last_ring_surface_review import FieldChart
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_rim_consistent_fields import checked_pair
from mouth_rim_shared_crystal import TerminalQuadrature, moving_energy
from mouth_supercore_joined_front import file_sha

STEPS=(1e-5,5e-6)


class DiskMotion:
    def __init__(self,ref):
        self.shape=ref.shape
        self.rx=ref.rr[0]
        t=np.clip((self.rx-.90)/.10,0.,1.)
        self.weight=1.-t*t*(3.-2.*t)
        L,h=ref.shape.length,ref.h
        lo,hi=ref.z[1],ref.z[2]
        i=np.argmin(abs(lo-(L-h)))
        j=np.argmin(abs(hi-(L+h)))
        self.ramps={'short':[L-h,L,L,L+h]}
        for name,k in [('near',2),('far',3)]:
            self.ramps[name]=[lo[i-k],lo[i-k+1],hi[j+k-1],hi[j+k]]

    def profile(self,r,z):
        R,slope=self.shape.radius(z)
        x=r/R
        i=np.clip(np.searchsorted(self.rx,x)-1,0,len(self.rx)-2)
        deriv=np.diff(self.weight)[i]/np.diff(self.rx)[i]
        g=np.interp(x,self.rx,self.weight)
        return g,deriv/R,-deriv*x*slope/R

    def lift(self,name,r,z):
        left,mid1,mid2,right=self.ramps[name]
        A=np.clip((z-left)/(mid1-left),0,1)*np.clip((right-z)/(right-mid2),0,1)
        Az=((z>left)&(z<mid1))/(mid1-left)-((z>mid2)&(z<right))/(right-mid2)
        g,gr,gz=self.profile(r,z)
        return A*g,A*gr,A*gz+Az*g,Az*g


def carrier_work(reader,motion):
    m=reader.model
    names=tuple(motion.ramps)
    # Columns: analytic physical work; part from axial shell caps; native
    # analytic work; signed physical and native derivatives at both steps.
    out={name:np.zeros(7) for name in names}
    gauss,weights=roots_legendre(2)
    gauss=(gauss+1)/2
    bits=np.array(list(product((0,1),repeat=3)))
    for slab in (1,2):
        xs,ang,zaxis=m.rr[0],m.a,m.z[slab]
        ok=(zaxis[1:]>motion.ramps['far'][0]) & (zaxis[:-1]<motion.ramps['far'][-1])
        ix,ia,iz=np.meshgrid(np.arange(len(xs)-1),np.arange(len(ang)-1),
                              np.flatnonzero(ok),indexing='ij')
        ix,ia,iz=ix.ravel(),ia.ravel(),iz.ravel()
        dx,da,dz=np.diff(xs)[ix],np.diff(ang)[ia],np.diff(zaxis)[iz]
        nodes=[v.values[ix[:,None]+bits[:,0],ia[:,None]+bits[:,1],iz[:,None]+bits[:,2]]
               for v in reader.fields[slab,0]]
        for gx,ga,gt in product(gauss,repeat=3):
            q=np.array([gx,ga,gt])
            f=np.where(bits,q,1-q)
            N=f.prod(axis=1)
            deriv=np.stack([(2*bits[:,i]-1)*f[:,(i+1)%3]*f[:,(i+2)%3] for i in range(3)])
            x=xs[ix]+gx*dx
            a=ang[ia]+ga*da
            z=zaxis[iz]+gt*dz
            R,slope=m.shape.radius(z)
            r=R*x
            w=4*np.pi*r*r*np.sin(a)*R*dx*da*dz/8.
            F,Fp=reader.vortex(r*np.sin(a))
            states=[]
            native=[]
            for n in nodes:
                val=np.einsum('i,qij->qj',N,n)
                grad=np.einsum('ki,qij->qjk',deriv,n)
                grad[:,:,0]/=(R*dx)[:,None]
                grad[:,:,1]/=(r*da)[:,None]
                grad[:,:,2]/=dz[:,None]
                grad[:,:,2]-=(slope*x)[:,None]*grad[:,:,0]
                states.append(physical(1,r,a,val.T,grad.transpose(1,0,2),F,Fp))
                y=val[:,0]**2
                un=F**10*(y-1)**2*(y**3+2*y**2+3*y+4)/20
                native.append((grad[:,0],un))
            t0,t1=(stresses(s) for s in states)
            Tr=(t0[0]-t1[0]).sum(axis=0)
            Tz=(t0[1]-t1[1]).sum(axis=0)
            for name in names:
                _,gr,gz,cap=motion.lift(name,r,z)
                values=out[name]
                ana=-w@(Tr*gr+Tz*gz)
                values[0]+=ana
                values[1]-=w@(Tz*cap)
                def native_shift(s,grad,potential,step):
                    gf,_,_,n,swirl,_=s
                    def squared(v):
                        vz=v[:,2]/(1.+step*gz)
                        return (v[:,0]-step*gr*vz)**2+v[:,1]**2+vz**2
                    return (1.+step*gz)*(.5*(F**2*squared(grad)-squared(gf)-n*swirl)
                                             +potential-(n**5-5*n+4)/20)
                # Complex step through the native change of variables is a
                # third, sign-sensitive check independent of a real secant.
                complex_step=1j*1e-25
                dc=native_shift(states[0],*native[0],complex_step)-native_shift(states[1],*native[1],complex_step)
                values[2]+=ana+w@(dc.imag/complex_step.imag)
                for j,h in enumerate(STEPS):
                    signed=[]
                    for e in (-h,h):
                        de=transformed(states[0],gr,gz,e)-transformed(states[1],gr,gz,e)
                        dn=native_shift(states[0],*native[0],e)-native_shift(states[1],*native[1],e)
                        signed.append([w@de,w@(de+dn)])
                    values[3+2*j:5+2*j]+=(np.array(signed[1])-signed[0])/(2*h)
    return out


def crystal_work(quad,delta,motion):
    """Move disk labels, keeping every other original crystal node fixed.

    The axial crystalline lift is the terminal-row material hat times the exact
    shared disk weight. For each virtual state, sample the scalar order on
    the inverse field-coordinate map at the new solid position.
    """
    out={name:np.zeros(6) for name in motion.ramps}
    for v in quad.blocks():
        r,z,K,J,W,vol=[v[key] for key in ('r','z','K','J','W','measure')]
        chi,cx,cz=interpolation(quad.z,quad.ax,delta,z,v['x'])
        R,slope=quad.sh.radius(z)
        cz-=cx*v['x']*slope/R
        hat=(v['ids']//len(quad.xt)==len(quad.tt)-1)
        A=(hat*v['bary']).sum(axis=1)
        grads=np.stack([(hat*v[g]).sum(axis=1) for g in ('gx','gt')],axis=1)
        grads=np.einsum('qi,qij->qj',grads,np.linalg.inv(v['B']))
        g,gr,gz=motion.profile(r,z)
        us=A*g
        sr=gr*A+g*grads[:,0]
        sz=gz*A+g*grads[:,1]
        elastic=vol@(chi*(K[:,0,1]*sr+(K[:,1,1]-1+.5*J*(J-1))*sz))
        for name,values in out.items():
            uf=motion.lift(name,r,z)[0]
            relative=vol@(W*cz*(us-uf))
            values[0]+=elastic
            values[1]+=relative
            # Split advection into the solid's contribution at fixed spatial
            # order and the field-coordinate contribution; retain both signs.
            values[2]+=vol@(W*cz*us)
            for i,h in enumerate(STEPS):
                signed=[]
                for e in (-h,h):
                    target=z+e*us
                    pre=z.copy()
                    for _ in range(3):
                        eta,_,etaz,_=motion.lift(name,r,pre)
                        pre-=(pre+e*eta-target)/(1+e*etaz)
                    if np.max(abs(pre+e*motion.lift(name,r,pre)[0]-target))>2e-12 or np.max(pre)>quad.L+2e-12:
                        raise ArithmeticError('Crystal did not follow the specified disk trace')
                    c=interpolation(quad.z,quad.ax,delta,pre,r/quad.sh.radius(pre)[0])[0]
                    perturbed=moving_energy(K,J,v['hoop'],sr,sz,e)
                    signed.append(vol@(c*perturbed-chi*W))
                values[3+i]+=(signed[1]-signed[0])/(2*h)
    return out


def run(state,wall,fields):
    x,d,_,_=read_source(state,wall)
    checked_pair(fields)
    base=Front(state,GRIDS['native'],1)
    tr=Trial(base,state,x,11,32,False)
    motion=DiskMotion(tr.model)
    with np.load(fields) as archive:
        level=archive['level'].item()
        m=FieldChart(tr.shape,base.c,1,GRIDS[level],precondition=False,terminal_profile=tr.model.terminal_profile)
        m.rr=tuple(archive[n] for n in ('inside_x','outside_x'))
        m.z=tuple(archive[n] for n in ('feed_z','duct_z','receiver_z'))
        pair=[archive[n] for n in ('wound','rest')]
        assert np.array_equal(m.a,archive['angles'])
        if archive['checkpoint_sha'].item()!=file_sha(wall):raise ValueError('Wrong wall')
        reader=Reader(m,base.c,*pair)
        dc=pair[0][2,m.ids[1,0]]-pair[1][2,m.ids[1,0]]
        order=np.einsum('raj,a->jr',dc,angular_weights(m.a))
        load_sha=archive['load_sha'].item() if 'load_sha' in archive else None
    fluid=carrier_work(reader,motion)
    solid=crystal_work(TerminalQuadrature(state,wall,m.rr[0],m.z[1],64),order,motion)
    rows={}
    for name,F in fluid.items():
        S=solid[name]
        mat=float(S[0]+S[1])
        err_field=max(abs(F[[3,5]]-F[0]).max(),abs(F[[4,6]]-F[2]).max())
        err_solid=max(abs(S[3:5]-mat))
        if err_field>3e-5 or err_solid>5e-5:
            raise ArithmeticError(f'Unresolved signed-work check {name} {err_field} {err_solid}')
        rows[name]=dict(field_physical=F[0],axial_cap_part=F[1],field_native=F[2],
                        field_signed_physical=F[[3,5]].tolist(),field_signed_native=F[[4,6]].tolist(),
                        solid_strain_work=S[0],solid_relative_order_work=S[1],solid_spatial_advection=S[2],
                        solid_signed_work=S[3:5].tolist(),native_actuator=F[2]+mat,physical_actuator=F[0]+mat)
    return dict(wall=file_sha(wall),fields=file_sha(fields),load=load_sha,grid=level,
                disk_rule='R4 bonded outer wall clamped, unit motion inboard of 0.90 normalized radius, cubic taper sampled on native radial nodes. All original nonterminal solid sites fixed.',signed_steps=STEPS,ramps=motion.ramps,reactions=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('state','wall','fields','output'):p.add_argument(arg,type=Path)
    a=p.parse_args()
    data=run(a.state,a.wall,a.fields)
    a.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
