"""Whole-disk axial virtual translations enclosing all banks of the held R4 exit.

A unit FE displacement covers the entire mobile layer and the common rim.
The three boundaries of the support taper only in exterior cells. The weak
surface flux is an average over closed surfaces, with outward orientation.
A separately signed native-energy pullback includes the two sealed wall terms.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_disk_pillbox import dot, physical, stresses, lag, transformed
from mouth_joint_polar_shape import Trial
from mouth_joined_front_audit import Front
from mouth_polar_exit_readout import accepted
from mouth_bonded_front_adjoint import set_front

from mouth_supercore_joined_front import check_guard, emit, file_sha, ORDER_INERTIA


class Closed:
    def __init__(self,m,b,d,pair,window_stride=1):
        self.window_stride=window_stride
        self.m,self.b,self.d,self.pair=m,b,d,pair
        wt,r,a,z,F,Fp,*_=m.D.T
        self.states=[]
        self.order_grad=[]
        self.solid_grad=[]
        q=d[b.mask]*b.scale
        self.natural=[]
        for U in pair:
            t,p,c=(m.V@v for v in U)
            gv,gp,gc=(np.column_stack([g@u for g in m.G]) for u in U)
            gc+=m.ref[:,1:]
            state=physical(m.winding,r,a,(t,p,m.ref[:,0]+c),(gv,gp,gc),F,Fp)
            self.states.append(state)
            self.natural.append((t,gv))
            _,g=m.evaluate(U)
            self.order_grad.append(g[2]);set_front(b,m,U);b.energy(q)
            self.solid_grad.append(b.native.copy())
        self.Tr,self.Tz=[np.sum(v,axis=0) for v in (
            a-bb for a,bb in zip(stresses(self.states[0]),stresses(self.states[1])))]
        self.j=lag(self.states[0])[1]*self.states[0][3][:,None]-lag(self.states[1])[1]*self.states[1][3][:,None]
        self.mass_lump=m.V.T@wt
        self.znodes=np.empty(m.size)
        self.bank=[]
        for (slab,side),ids in m.ids.items():
            self.znodes[ids]=np.broadcast_to(m.z[slab],ids.shape)
        self.clip=(abs(z-m.shape.length)<m.h+1.e-12)
        self.shape=m.shape
        self.wt=wt

    def lifts(self,k,j):
        m=self.m
        L=m.shape.length
        iz=np.flatnonzero(np.isclose(m.z[1],L-m.h,atol=1e-11,rtol=0.))[0]
        jo=np.flatnonzero(np.isclose(m.z[2],L+m.h,atol=1e-11,rtol=0.))[0]
        stride=self.window_stride
        lo1,lo2=m.z[1][[iz-k*stride,iz-(k-1)*stride]]
        hi1,hi2=m.z[2][[jo+(k-1)*stride,jo+k*stride]]
        x1,x2=m.rr[1][[j*stride,(j+1)*stride]]
        def ax(z):
            return np.clip((z-lo1)/(lo2-lo1),0,1)*np.clip((hi2-z)/(hi2-hi1),0,1)
        az=ax(self.znodes)
        rad=np.ones_like(az)
        for (slab,side),ids in m.ids.items():
            if side:
                rad[ids]=np.broadcast_to(np.clip((x2-m.rr[1])/(x2-x1),0,1)[:,None,None],ids.shape)
        return ax,az,rad,dict(lower_half_thickness=[L-lo2,L-lo1],upper_half_thickness=[hi1-L,hi2-L],
                              lateral_radius_at_disk=(L*0+m.shape.exit+(28.-m.shape.exit)*np.array([x1,x2])).tolist())

    def native_wall(self,az):
        m=self.m
        slope=np.diff(az[m.ids[1,0][-1,0]])/np.diff(m.z[1])
        use=abs(slope)>0
        if np.any(m.z[1][:-1][use]<.62*m.shape.length):
            raise ValueError('Polar wall distortion overlaps boundary-work shell')
        work=0.
        for c in ((1.-1/np.sqrt(3))/2,(1.+1/np.sqrt(3))/2):
            a=(m.a[:-1]+np.diff(m.a)*c)[:,None]
            for h in ((1.-1/np.sqrt(3))/2,(1.+1/np.sqrt(3))/2):
                z=m.z[1][:-1]+np.diff(m.z[1])*h
                r,_=m.shape.radius(z)
                F,Fp=m._vortex(r[None,:]*np.sin(a))
                weights=4*np.pi*r**2*np.sin(a)**2*F*Fp*np.diff(m.a)[:,None]*np.diff(m.z[1])*slope/4.
                for side in (0,1):
                    ids=m.ids[1,side][0 if side else -1]
                    def val(u):
                        v=u[0,ids]
                        return (1-c)*((1-h)*v[:-1,:-1]+h*v[:-1,1:])+c*((1-h)*v[1:,:-1]+h*v[1:,1:])
                    work+=(1-2*side)*.5*np.sum(weights*(val(self.pair[0])**2-val(self.pair[1])**2))
        return work

    def audit(self,k,j):
        m=self.m
        wt,r,ang,z,F,*_=m.D.T
        ax,aznode,rnode,geometry=self.lifts(k,j)
        A=m.V@aznode;rad=m.V@rnode
        eta=aznode*rnode
        grad=np.stack([G@eta for G in m.G],axis=1)
        za=(m.G[2]@aznode)*rad
        gr=(m.G[0]@rnode)*A
        gz=(m.G[2]@rnode)*A
        if np.max(abs(grad-np.column_stack((gr,0*gr,gz+za))))>1.e-10:
            raise ArithmeticError('Closed test gradient does not factor')
        if min(m.V[self.clip & (r<m.shape.radius(z)[0])]@eta)<1-1e-12:
            raise ArithmeticError('Not all of the mobile rim lies inside closed support')
        if np.max(abs(grad[(z<.62*m.shape.length)&(z>.06*m.shape.length)]))>2e-11:
            raise ValueError('Weight shell intersects polar displacement')
        # The receiving solid labels and order are identical in drive and rest.
        # The bonded solid has drive-dependent order only in the mobile band.
        # There eta=1 up to and beyond the side; its strain variation is zero.
        if np.max(abs(eta[np.r_[m.feed,m.bulk]]))!=0:
            raise ValueError('Local test touches far-bath constraint')
        below=z<m.shape.length
        def faces(vr,vz):
            return np.array([-wt[below]@(vz*za)[below],-wt[~below]@(vz*za)[~below],
                             -wt@(vr*gr+vz*gz)])
        momentum=faces(self.Tr,self.Tz)
        current=faces(self.j[:,0],self.j[:,2])
        # Hold material reference coordinates and convect all primary fields
        # under the same axial map. Normal material elasticity cancels here.
        exact=np.sum(momentum)
        signed={}
        native={}
        wall=self.native_wall(aznode)
        for step in (2e-5,1e-5):
            pert=[]
            factored=[]
            for s in (-step,step):
                dE=transformed(self.states[0],gr,za+gz,s)-transformed(self.states[1],gr,za+gz,s)
                pert.append(wt@dE)
                det=1+s*(za+gz)
                terms=[]
                for state,(t,gv) in zip(self.states,self.natural):
                    # Native factorization differs only in the amplitude channel.
                    gf,gp,gc,n,swirl,chi=state
                    old=gf.copy();old[:,2]=gf[:,2]/det;old[:,0]-=s*gr*old[:,2]
                    new=gv.copy();new[:,2]=gv[:,2]/det;new[:,0]-=s*gr*new[:,2]
                    y=t*t
                    interaction=F**10*(y-1)**2*(y**3+2*y**2+3*y+4)/20
                    raw=(.5*(dot(old,old)+n*swirl)+(n**5-5*n+4)/20)
                    terms.append(.5*F**2*dot(new,new)+interaction-raw)
                factored.append(pert[-1]+wt@(det*(terms[0]-terms[1]))+s*wall)
            signed[str(step)]=(pert[1]-pert[0])/(2*step)
            native[str(step)]=(factored[1]-factored[0])/(2*step)
        if max(abs(v-exact) for v in signed.values())>1e-4:
            raise ArithmeticError(f'Signed material geometry work disagrees {exact}, {signed}')
        if abs(sum(current))>3e-5:raise ArithmeticError(f'All-bank weak mass balance failed {current}')
        work=self.constraint_diagnostics(ax,eta)
        return dict(kind='closed_work',axial_band=k,lateral_band=j,**geometry,
                    flux_order=['bottom_outward','top_outward','side_outward'],mass=current.tolist(),
                    outward_mass=sum(current),momentum=momentum.tolist(),outward_momentum=exact,
                    signed_common_translation_work=signed,native_factored_translation=native,
                    sealed_vortex_boundary_term=wall,crystal_flux_drive_minus_rest=0.,
                    **work)

    def constraint_diagnostics(self,ax,eta):
        m,b=self.m,self.b
        xx,tt=np.meshgrid(b.x,b.t)
        z=tt*m.shape.length+self.d[...,1]
        test=ax(z)  # all solid nodes lie inside the radial plateau
        delta=self.solid_grad[0]-self.solid_grad[1]
        held=~b.mask[...,1]
        # Fixed-order reaction quoted for one explicit Eulerian interpolation
        # of the convective derivative. At the rim each bank has its own slope.
        # Evaluate two alternative extensions only at shared rim nodes to expose
        # the uncertainty; this partition must not be added to common work.
        fixed=~m.masks[2]
        cap=m.ids[1,0][-1,:,-1]
        zz=self.znodes
        fields=[s[2][:,2] for s in self.states]
        def proj(v,mask=None):
            if mask is None:return (m.V.T@(self.wt*v))/self.mass_lump
            w=self.wt*mask
            mass=m.V.T@w
            return np.divide(m.V.T@(w*v),mass,out=np.zeros_like(mass),where=mass>0)
        vals=[]
        _,r,ang,z,*_=m.D.T
        below=z<m.shape.length
        inside=r<m.shape.radius(z)[0]
        for case in ('all','incident','receiving','outside'):
            work=0.
            for i,v in enumerate(fields):
                den=proj(v)
                if case!='all':
                    take={'incident':below&inside,'receiving':~below&inside,'outside':~inside}[case]
                    den[cap]=proj(v,take)[cap]
                work+=(1-2*i)*np.sum(self.order_grad[i][fixed]*(-eta*den)[fixed])
            vals.append(work)
        return dict(bonded_axial_pin_work_fixed_eulerian_order=float((delta[...,1]*test)[held].sum()),
                    free_solid_axial_work_fixed_eulerian_order=float((delta[...,1]*test)[~held].sum()),
                    outer_crystal_pin_difference=0.,bath_phase_and_amplitude_pin_work=0.,
                    advected_fixed_order_work_extensions=dict(zip(('all_volume','duct_rim','upper_disk_rim','outer_bank_rim'),vals)),
                    convected_order_value_change_on_fixed_labels=0.)


def run(a):
    prep,rec,x,d,baseline=accepted(a.state,a.checkpoint)
    check_guard(a.fields)
    from json import loads
    previous=[loads(s) for s in (a.fields.parent/'stdout').read_text().splitlines()]
    if not any(row.get('sha256')==file_sha(a.fields) for row in previous):
        raise ValueError('Field digest missing')
    with np.load(a.fields) as z:
        if z['checkpoint_sha'].item()!=rec['sha256']:raise ValueError('Unmatched states')
        pair=[z['wound'],z['rest']]
    if not np.array_equal(pair[0],baseline):raise ValueError('Input drive fields disagree')
    base=Front(a.state,tuple(prep['field_grid']),1)
    t=Trial(base,a.state,x,(d.shape[1]-1)//3,d.shape[0]-1,True)
    t.model.material_load=t.body.get_material_load(d[t.body.mask]*t.body.scale,t.model)
    audits(t.model,t.body,d,pair,rec['sha256'],file_sha(a.fields))


def audits(m,b,d,pair,checkpoint,fields_hash):
    study=Closed(m,b,d,pair)
    emit(kind='preparation',checkpoint_sha256=checkpoint,controls_sha256=fields_hash,
         active_halfwidth=m.h,solid_cells=[d.shape[1]-1,d.shape[0]-1],joint_geometry=[m.shape.length,m.shape.exit],
         drive_current=m.measure(pair[0])['feed_current'],rest_current=m.measure(pair[1])['feed_current'])
    for axial in (2,3):
        for side in (4,5):emit(**study.audit(axial,side))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('state','checkpoint','fields'):p.add_argument(key,type=Path)
    run(p.parse_args())
