"""Signed work and one-sided traction on two independently moved slit banks of the fixed R4 wall."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_closed_rim_surface_review import quadrature, bracket
from mouth_disk_pillbox import dot, physical, stresses, transformed
from mouth_last_ring_solid import interpolation
from mouth_rim_disk_work_review import Disk
from mouth_rim_shared_crystal import TerminalQuadrature, moving_energy
from mouth_supercore_joined_front import emit,file_sha


class Bank:
    def __init__(self,disk,bank,strip,width):
        self.disk,self.side,self.strip,self.edge=disk,bank,strip,width

    def motion(self,r,z):
        b,s=self.disk.sh.radius(z)
        if self.side:
            width=28-b; x=(r-b)/width; factor=1-x
            q=np.clip(1-x/self.edge,0.,1.)
            dq=np.where((x>=-1e-12)&(x<self.edge),-1/self.edge,0.)
        else:
            width=b; x=r/b; factor=x
            q=np.clip((x-self.edge)/(1-self.edge),0.,1.)
            dq=np.where((x>self.edge)&(x<=1+1e-12),1/(1-self.edge),0.)
        a,b,c=self.strip
        az=np.maximum(0,np.minimum((z-a)/(b-a),(c-z)/(c-b)))
        dz=np.where((z>a)&(z<b),1/(b-a),np.where((z>b)&(z<c),-1/(c-b),0))
        gr=az*dq/width
        return q*az,gr,q*dz-s*factor*gr

    def physical(self,x,a,z):
        model=self.disk.model
        xi,theta,ws=(p.ravel() for p in np.broadcast_arrays(x,a,z))
        b,slope=self.disk.sh.radius(ws)
        width=28-b if self.side else b
        radial=b+xi*width if self.side else xi*width
        pts=np.column_stack((xi,theta,ws))
        F,Fp=self.disk.reader.vortex(radial*np.sin(theta))
        for interp in self.disk.reader.fields[1,self.side]:
            vals=interp(pts)
            raw=[]
            for i in range(3):
                t=pts.copy()
                lo,hi=bracket(interp.grid[i],pts[:,i]);t[:,i]=lo
                low=interp(t);t[:,i]=hi
                raw.append((interp(t)-low)/(hi-lo)[:,None])
            dr=raw[0]/width[:,None]
            grads=np.stack([dr,raw[1]/radial[:,None],raw[2]-(slope*(1-xi if self.side else xi))[:,None]*dr],axis=-1)
            gs=tuple(grads[:,k] for k in range(3))
            state=physical(model.winding,radial,theta,vals.T,gs,F,Fp)
            yield state,vals[:,0],gs[0],F

    def work(self,n):
        m=self.disk.model;sh=self.disk.sh
        lo,peak,hi=self.strip
        radial=m.rr[self.side]
        xs=radial[radial<self.edge] if self.side else radial[radial>self.edge]
        x,wx=quadrature(np.unique(np.r_[xs,self.edge]),n)
        z,zw=quadrature(np.unique(np.r_[lo,peak,hi,m.z[1][(m.z[1]>lo)&(m.z[1]<hi)]]),n)
        angle,aw=quadrature(m.a,n)
        x,a=np.meshgrid(x,angle,indexing='ij')
        wx=wx[:,None]*aw
        signed=np.zeros((2,2))
        native=np.zeros_like(signed)
        sectors=np.zeros(3)
        sharp=np.zeros(2)
        boundary=0.
        for zz,wt in zip(z,zw):
            b,slope=sh.radius(zz)
            width=28-b if self.side else b
            r=b+x*width if self.side else x*width
            _,gr,gz=self.motion(r,zz)
            gr,gz=gr.ravel(),gz.ravel()
            area=4*np.pi*wt*width*(wx*r*r*np.sin(a)).ravel()
            for sign,(state,t,gv,F) in zip((1.,-1.),self.physical(x,a,zz)):
                tr,tz=stresses(state)
                sectors+=sign*((tr*gr+tz*gz)@area)
                gf=state[0];n0=state[3];swirl=state[4];y=t*t
                interaction=F**10*(y-1)**2*(y**3+2*y*y+3*y+4)/20
                for k,h in enumerate((2e-6,1e-6)):
                    for i,e in enumerate((-h,h)):
                        density=transformed(state,gr,gz,e)
                        signed[k,i]+=sign*(area@density)
                        det=1+e*gz
                        ngu=gv.copy();ngu[:,2]/=det;ngu[:,0]-=e*gr*ngu[:,2]
                        ngf=gf.copy();ngf[:,2]/=det;ngf[:,0]-=e*gr*ngf[:,2]
                        delta=(.5*F**2*dot(ngu,ngu)+interaction
                               -.5*(dot(ngf,ngf)+n0*swirl)-(n0**5-5*n0+4)/20)
                        native[k,i]+=sign*(area@(density+det*delta))
            eta=self.motion(b,zz)[0]
            tr=self.disk.reader.sample(1,self.side,b,angle,zz,np.array([1.,0.,0.]))
            tz=self.disk.reader.sample(1,self.side,b,angle,zz,np.array([0.,0.,1.]))
            wall_weight=(1-2*self.side)*4*np.pi*wt*b*b*np.sin(angle)*aw
            sharp+=(wall_weight*eta)@(tr-slope*tz)
            # d(eta_trace)/dz along the moving sealed bank; factorized vortex
            # term uses oriented radial area, not just meridional arc length.
            der=1/(peak-lo) if zz<peak else -1/(hi-peak)
            F,Fp=self.disk.reader.vortex(b*np.sin(angle))
            idx=0. if self.side else 1.
            point=np.column_stack([np.full(len(angle),idx),angle,np.full(len(angle),zz)])
            amp=[interp(point)[:,0] for interp in self.disk.reader.fields[1,self.side]]
            boundary+=.5*wall_weight@(np.sin(angle)*der*F*Fp*(amp[0]**2-amp[1]**2))
        steps=np.array([2e-6,1e-6])
        sign_work=-(signed[:,1]-signed[:,0])/(2*steps)
        native_work=-(native[:,1]-native[:,0])/(2*steps)-boundary
        if np.max(abs(sign_work-sectors.sum()))>2e-5:
            raise ArithmeticError('Signed physical field work inconsistent with analytic pullback')
        return dict(quadrature=n,canonical_sector_reactions=sectors.tolist(),direct_reaction=sign_work.tolist(),
                    native_including_vortex_bank=native_work.tolist(),vortex_bank_energy_derivative=float(boundary),
                    outward_wall_axial=float(sharp[1]),outward_wall_mass=float(sharp[0]))


def run(a):
    disk=Disk(a.state,a.wall,a.fields)
    sh=disk.sh;m=disk.model
    with np.load(a.reference) as g:
        axes=[g['inside_x'],g['outside_x']]; zd=g['duct_z']
    for coarse,fine in zip([*axes,zd],[*m.rr,m.z[1]]):
        if any(np.min(abs(v-fine))>2e-12 for v in coarse):
            raise ValueError('Probes must have identical physical boundaries')
    join=np.argmin(abs(zd-(sh.length-m.h)))
    strips={'near':zd[[-9,-5,-1]],'middle':zd[[-25,-17,-9]],
            'upstream':zd[[join-4,join-3,join-1]]}
    if strips['upstream'][0]<.62*sh.length:
        raise ValueError('Probe overlaps angular deformation')
    probes={}
    for name,strip in strips.items():
        for s in (0,1):
            for cells in (4,8):
                edge=axes[s][cells if s else -1-cells]
                key=(name,s,cells);probes[key]=Bank(disk,s,strip,edge)
    with np.load(a.fields) as f:
        mean=np.einsum('raj,a->jr',(f['wound'][2]-f['rest'][2])[m.ids[1,0]],angular_weights(m.a))
    quad=TerminalQuadrature(a.state,a.wall,m.rr[0],m.z[1],64)
    work={key:0. for key in probes}
    finite={key:np.zeros(2) for key in probes}
    for v in quad.blocks():
        chi=interpolation(quad.z,quad.ax,mean,v['z'],v['x'])[0]
        if not chi.any():continue
        K,J=v['K'],v['J']; wtchi=v['measure']*chi
        for key,b in probes.items():
            if b.side or not (v['z'].min()<b.strip[-1] and v['z'].max()>b.strip[0]):continue
            _,ur,uz=b.motion(v['r'],v['z'])
            work[key]-=wtchi@(K[:,0,1]*ur+(K[:,1,1]-1+.5*J*(J-1))*uz)
            for n,h in enumerate((-1e-6,1e-6)):
                finite[key][n]-=wtchi@(moving_energy(K,J,v['hoop'],ur,uz,h)-v['W'])
    emit(kind='source',**disk.source_hash,reference=file_sha(a.reference),load_plane=a.level)
    for key,b in probes.items():
        sign=(finite[key][1]-finite[key][0])/(2e-6)
        if abs(sign-work[key])>3e-5:
            raise ArithmeticError('Signed solid work check failed')
        for n in (2,4):
            row=b.work(n)
            full=row['native_including_vortex_bank'][-1]+work[key]
            emit(kind='bank_test',strip=key[0],bank=('duct','outside')[b.side],radial_cells=key[2],
                 depth_knots=b.strip.tolist(),normalized_radial_bound=float(b.edge),
                 solid_reaction=work[key],signed_solid=sign,total_native=full,
                 native_minus_wall=full-row['outward_wall_axial'],**row)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','wall','fields','reference'):p.add_argument(name,type=Path)
    p.add_argument('level',choices=('double','local'))
    run(p.parse_args())
