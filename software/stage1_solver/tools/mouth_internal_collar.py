"""Two-component meridional material maps for the finite supersolid collars.

The piecewise-linear triangular unknowns are inverse material labels on a chart fitted to the wall. This
keeps the order field Eulerian as material sites move. The reference measure
is varied along with the inverse label map; the weak equation in a
direct map also has the occupancy force W grad(chi).
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss, legder, legval

from mouth_poly_shallow_release import CUT, exterior, source
from mouth_sealed_shape_relax import Meridian


def interpolation(t,x,grid):
    tz,rx,chi=grid
    j=np.clip(np.searchsorted(tz,t)-1,0,len(tz)-2)
    i=np.clip(np.searchsorted(rx,x)-1,0,len(rx)-2)
    h=(t-tz[j])/(tz[j+1]-tz[j]); k=(x-rx[i])/(rx[i+1]-rx[i])
    if np.min(h)<-1e-9 or np.max(h)>1+1e-9 or np.min(k)<-1e-9 or np.max(k)>1+1e-9:
        raise ValueError('Collar extends outside its Eulerian order field')
    n=len(rx)
    ids=(j*n+i)[:,None]+np.array([0,1,n,n+1])
    weights=np.stack(((1-h)*(1-k),(1-h)*k,h*(1-k),h*k),-1)
    return (chi.ravel()[ids]*weights).sum(axis=1),ids,weights


class CollarPair:
    def __init__(self,shape,co,na,nt,quad=3):
        self.shape,self.na,self.nt=shape,na,nt
        self.a=np.linspace(-.5,.5,2*na+1); self.t=np.linspace(0.,1.,nt+1)
        et,tt=np.meshgrid(self.a,self.t)
        self.identity=np.stack((et,tt),axis=-1)
        self.mask=np.ones_like(self.identity,dtype=bool)
        self.mask[:,[0,na,-1],0]=False
        self.mask[:,-1,1]=False
        self.mask[[0,-1],:,1]=False
        self.size=et.size
        nx=2*na+1
        nid=(np.arange(nt)[:,None]*nx+np.arange(2*na)).ravel()
        corners=np.array([0,1,nx,nx+1])
        self.cell_ids=nid[:,None]+corners
        q,qw=leggauss(quad);q=(1+q)/2;qw=qw/2
        # Two positive-area triangles per chart rectangle; the label
        # determinant is constant over an entire triangle.
        v,u=np.meshgrid(q,q,indexing='ij');v,u=v.ravel(),u.ravel()
        h=u; k=(1-u)*v
        qweight=np.tile(np.outer(qw,qw).ravel()*(1-u),2)
        corner0=np.array([0,1,nx]); corner1=np.array([1+nx,nx,1])
        triangles=np.stack((nid[:,None]+corner0,nid[:,None]+corner1),axis=1)
        self.triangles=triangles.reshape(-1,3)
        self.id=np.repeat(self.triangles,len(u),axis=0)
        bfs=np.stack((1-h-k,h,k),axis=1)
        self.B=np.tile(bfs,(len(self.triangles),1))
        sgn=np.repeat(np.tile([1.,-1.],len(nid)),len(h))
        self.E=sgn[:,None]*np.array([-1.,1.,0.])*(2*na)
        self.T=sgn[:,None]*np.array([-1.,0.,1.])*nt
        self.eta=(self.B*et.ravel()[self.id]).sum(axis=1)
        self.z=(self.B*tt.ravel()[self.id]).sum(axis=1)
        self.qw=4*np.pi*np.tile(qweight,len(nid))/(2*na*nt)
        self.co=co
        b,db,_=shape.profile(self.z)
        registration=np.exp(legval(2*self.z-1,co))
        pt=2*legval(2*self.z-1,legder(co))
        self.r=b*np.exp(-registration*self.eta)
        self.re=-registration*self.r
        self.rt=self.r*(shape.length*db/b-self.eta*registration*pt)
        self.L=shape.length
        self.inside=self.eta>0
        self.coordinate=np.where(self.inside,self.r/b,(self.r-b)/(CUT-b))
        self.inverse_nodes=self.id.ravel()
        self.chart_shape=et.shape

    def _measure(self,inner,outer):
        chi=np.empty_like(self.r); probes=[]
        for mask,grid in ((self.inside,inner),(~self.inside,outer)):
            c,ids,weights=interpolation(self.z[mask],self.coordinate[mask],grid)
            chi[mask]=c
            probes.append((mask,ids,weights,grid[-1].shape))
        if chi.min()<-1.e-11 or chi.max()>1.+1.e-11:
            raise ValueError('Invalid occupancy')
        return chi,probes

    def _corners(self,l):
        verts=l.reshape(-1,2)[self.triangles]
        h=(verts[:,1]-verts[:,0])*2*self.na
        v=(verts[:,2]-verts[:,0])*self.nt
        return np.min(h[:,0]*v[:,1]-h[:,1]*v[:,0])

    def energy_density(self,labels,r,w):
        """Eulerian elastic energy W/J independent of occupancy, zero outside either collar."""
        r,w=np.broadcast_arrays(np.asarray(r,dtype=float),np.asarray(w,dtype=float))
        out=np.zeros(r.shape)
        t=np.clip(w/self.L,0.,1.)
        b,slope,_=self.shape.profile(t)
        lp=legval(2*t-1,self.co)
        p=np.exp(lp)
        e=-np.log(np.maximum(r/b,1.e-100))/p
        live=(e>=-.5)&(e<=.5)&(w>=0)&(w<=self.L)
        if not np.any(live):
            return out
        e,t,p,b,r,slope=e[live],t[live],p[live],b[live],r[live],slope[live]
        xx=(e+.5)*2*self.na; yy=t*self.nt
        i=np.minimum(xx.astype(int),2*self.na-1);j=np.minimum(yy.astype(int),self.nt-1)
        u,v=xx-i,yy-j
        upper=u+v>1
        ids=self.triangles[2*(j*2*self.na+i)+upper]
        h,k=np.where(upper,1-u,u),np.where(upper,1-v,v)
        nodes=labels.reshape(-1,2)[ids]
        A=np.sum(np.column_stack((1-h-k,h,k))*nodes[:,:,0],axis=1)
        sign=np.where(upper,-1.,1.)
        Ae,Te=((nodes[:,1]-nodes[:,0])*(sign*2*self.na)[:,None]).T
        At,Tt=((nodes[:,2]-nodes[:,0])*(sign*self.nt)[:,None]).T
        det=Ae*Tt-At*Te
        re=-p*r;rt=r*(self.L*slope/b-e*p*2*legval(2*t-1,legder(self.co)))
        hoop=r/(self.shape.neck+A)
        J=-re*self.L*hoop**2/det
        trace=((re*Tt-rt*Te)**2+(rt*Ae-re*At)**2+self.L**2*(Te**2+Ae**2))/det**2+2*hoop**2
        out[live]=(.5*(trace-4)-np.log(J)+.25*(J-1)**2)/J
        return out

    def _value(self,l,chi,read=False):
        nodes=l.reshape(-1,2)[self.id]
        pos=np.einsum('qi,qi->q',self.B,nodes[:,:,0])+self.shape.neck
        Ae,Te=np.einsum('qi,qik->kq',self.E,nodes)
        At,Tt=np.einsum('qi,qik->kq',self.T,nodes)
        d=Ae*Tt-At*Te
        min_det=self._corners(l)
        if min_det<=0 or np.min(pos)<=0:
            raise ValueError('Folded material map')
        # Fbar = [re,rt; 0,L] D^{-1}, changing the sign of its last column
        # recovers F for reference Z=1-T without changing either singular value.
        f11=(self.re*Tt-self.rt*Te)/d; f12=(-self.re*At+self.rt*Ae)/d
        f21=-self.L*Te/d; f22=self.L*Ae/d
        h=self.r/pos
        J=-self.re*self.L*h*h/d
        xx=f11*f11+f21*f21
        yy=f12*f12+f22*f22
        xy=f11*f12+f21*f22
        W=.5*(xx+yy+2*h*h-4)-np.log(J)+.25*(J-1)**2
        Q=-1+.5*J*(J-1)
        weight=self.qw*pos**2*d
        measure=weight*chi
        energy=measure@W
        gS=2*measure/pos*(W-h*h-Q)
        h11=W-Q-xx;h22=W-Q-yy;h12=-xy
        c=measure/d
        gae=c*(h11*Tt-h12*At); gat=c*(-h11*Te+h12*Ae)
        gte=c*(h12*Tt-h22*At); gtt=c*(-h12*Te+h22*Ae)
        def scatter(v):
            return np.bincount(self.inverse_nodes,weights=v.ravel(),minlength=self.size).reshape(self.chart_shape)
        g=np.stack((scatter(gS[:,None]*self.B+gae[:,None]*self.E+gat[:,None]*self.T),
                    scatter(gte[:,None]*self.E+gtt[:,None]*self.T)),-1)
        if not read:
            return energy,g
        tr=xx+yy
        top=np.sqrt(.5*(tr+np.sqrt(np.maximum(0.,(xx-yy)**2+4*xy*xy))))
        bottom=(-self.re*self.L/d)/top
        extrema=dict(min_J=float(J.min()),max_J=float(J.max()),min_stretch=float(np.minimum(h,bottom).min()),
                     max_stretch=float(np.maximum(h,top).max()),min_inverse_corner_det=float(min_det),
                     zero_order_quadrature_sites=int(np.sum(chi==0)),min_order=float(chi.min()))
        return energy,g,weight*W,extrema

    def evaluate(self,labels,inner,outer,diagnostics=False):
        chi,probes=self._measure(inner,outer)
        energy,grad,coeff,extrema=self._value(labels,chi,True)
        result=dict(energy=energy,gradient=grad,residual=float(np.max(abs(grad[self.mask]))))
        for key,(mask,ids,w,shape) in zip(('inner_source','outer_source'),probes):
            result[key]=np.bincount(ids.ravel(),weights=(w*coeff[mask,None]).ravel(),
                                     minlength=np.prod(shape)).reshape(shape)
        if diagnostics:
            # These conjugates act on boundary reference labels. A net
            # displacement of an edge in the spatial meridian needs separate
            # global shape work, including the exterior order/carrier work.
            outer_edge=grad[:,0,0]
            result.update(extrema,receiver_edge_label_work_sum=float(outer_edge.sum()),
                          seam_label_work_sum=float(grad[:,self.na,0].sum()),
                          core_join_label_work_sum=grad[:,-1].sum(axis=0).tolist(),
                          max_free_radial=float(np.max(abs(grad[:,:,0][self.mask[:,:,0]]))),
                          max_free_depth=float(np.max(abs(grad[:,:,1][self.mask[:,:,1]]))))
        return result

    def solve(self,inner,outer,initial=None,gauge=1.e-6):
        native,_=self._measure(inner,outer)
        # Null elements have no lattice. This extension only selects an
        # orientation-preserving chart there. Test both the native residual
        # and the dependence on its vanishing coefficient.
        nel=len(self.triangles)
        empty=np.repeat(native.reshape(nel,-1).max(axis=1)==0,len(native)//nel)
        chi=native+gauge*empty
        weight=self.qw*self.shape.neck**2*chi
        diag=np.bincount(self.inverse_nodes,weights=(weight[:,None]*(self.E**2+self.T**2+self.B**2)).ravel(),minlength=self.size)
        # The floor conditions the search metric only. Energy and weak force
        # use native chi on every quadrature site, including exact zeros.
        scale=np.repeat(np.sqrt(np.maximum(diag,1.e-6*diag.max())).reshape((*self.chart_shape,1)),2,axis=2)
        mask=self.mask & (np.broadcast_to(diag.reshape((*self.chart_shape,1)),self.mask.shape)>0)
        labels=self.identity.copy() if initial is None else initial.copy()
        labels[~self.mask]=self.identity[~self.mask]
        energy,g=self._value(labels,chi)
        g=g[mask]/scale[mask]
        history=[]
        status="iteration allowance reached"
        for iteration in range(3500):
            if np.linalg.norm(g,np.inf)<2.e-8:
                status="native weak residual converged"
                break
            vec=g.copy(); alpha=[]
            for step,change,rho in reversed(history):
                a=rho*np.dot(step,vec);alpha.append(a);vec-=a*change
            if history:
                step,change,rho=history[-1]
                vec*=np.dot(step,change)/np.dot(change,change)
            for (step,change,rho),a in zip(history,reversed(alpha)):
                vec+=step*(a-rho*np.dot(change,vec))
            if g@vec<=0:
                vec=g.copy();history=[]
            direction=-vec
            slope=g@direction
            step_size=1.
            accepted=False
            for _ in range(55):
                trial=labels.copy()
                trial[mask]+=step_size*direction/scale[mask]
                if self._corners(trial)>0:
                    next_energy,next_g=self._value(trial,chi)
                    if next_energy<=energy+1.e-4*step_size*slope+2.e-15*max(1.,abs(energy)):
                        accepted=True;break
                step_size*=.5
            if not accepted:
                status="no representable orientation-preserving descent"
                break
            delta=step_size*direction
            next_g=next_g[mask]/scale[mask]
            dg=next_g-g
            prod=dg@delta
            if prod>1.e-12*np.linalg.norm(dg)*np.linalg.norm(delta):
                history.append((delta,dg,1/prod))
                history=history[-16:]
            energy,g,labels=next_energy,next_g,trial
            if np.max(abs(delta))<3.e-14:
                status="orientation or precision limited before stationarity"
                break
        result=self.evaluate(labels,inner,outer,True)
        result.update(iterations=iteration,optimizer_status=status,
                      scaled_residual=float(np.max(abs(g))),active_dofs=int(mask.sum()),
                      gauge=float(gauge),empty_triangles=int(np.sum(empty)//(len(native)//nel)))
        return labels,result


def transfer(labels,na,nt):
    t=np.linspace(0,1,labels.shape[0]);a=np.linspace(-.5,.5,labels.shape[1])
    grids=np.meshgrid(np.linspace(0,1,nt+1),np.linspace(-.5,.5,2*na+1),indexing='ij')
    y=grids[0]*(len(t)-1);x=(grids[1]+.5)*(len(a)-1)
    j=np.minimum(y.astype(int),len(t)-2);i=np.minimum(x.astype(int),len(a)-2)
    u=x-i;v=y-j;top=u+v>1
    lo=(1-u-v)[...,None]*labels[j,i]+u[...,None]*labels[j,i+1]+v[...,None]*labels[j+1,i]
    hi=(u+v-1)[...,None]*labels[j+1,i+1]+(1-u)[...,None]*labels[j+1,i]+(1-v)[...,None]*labels[j,i+1]
    return np.where(top[...,None],hi,lo)


def fields(shape,order,co):
    t=np.linspace(0,1,order.shape[0]);x=np.linspace(0,1,order.shape[1])
    kw=dict(side_order_work=lambda w,y:source(shape,co,w/shape.length,y,True),return_side_field=True)
    out,(w,y,chi)=exterior(shape,len(x)-1,CUT,**kw)
    return (t,x,order),(w/shape.length,y,chi)


def main(args):
    guard=json.loads((args.state.parent/'outcome.json').read_bytes())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError('Unvalidated material state')
    with np.load(args.state,allow_pickle=False) as z:
        shape=Meridian(4.,z['shape_params'],float(z.get('shape_extra',0.)))
        if int(z['radius'])!=4:
            raise ValueError('Expected entrance four')
        co=z['registration']
        inner,outer=fields(shape,z['order'],co)
    labels=None
    for na in args.mesh:
        nt=na*2
        chart=CollarPair(shape,co,na,nt,args.quad)
        before=chart.evaluate(chart.identity,inner,outer,True)
        initial=transfer(labels,na,nt) if labels is not None else None
        labels,after=chart.solve(inner,outer,initial,gauge=args.gauge)
        np.savez_compressed(args.out/f'labels_{na}.npz',labels=labels,shape_params=shape.params,
                            shape_extra=shape.extra,registration=co,order=inner[-1],outer_order=outer[-1],
                            outer_t=outer[0],outer_x=outer[1],radius=4)
        chosen=lambda x:{k:v for k,v in x.items() if k not in ('gradient','inner_source','outer_source')}
        print(json.dumps(dict(kind='collar',source_sha256=hashlib.sha256(args.state.read_bytes()).hexdigest(),
                              na=na,nt=nt,quad=args.quad,before=chosen(before),after=chosen(after)),
                         allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('out',type=Path)
    p.add_argument('--mesh',type=int,nargs='+',default=[4,8,12])
    p.add_argument('--quad',type=int,default=3)
    p.add_argument('--gauge',type=float,default=1.e-6)
    main(p.parse_args())
