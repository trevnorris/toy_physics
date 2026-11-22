"""Stationary reversible carrier and outgoing transverse mode loading a responsive ordered core."""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from numpy.polynomial.legendre import leggauss

from mouth_exchange_carrier import assemble, mapped_nodes, BIAS
from mouth_live_waist_wave import OccupiedWave, transfer

WIDTH = .04
ORDER_INERTIA = WIDTH**2


def core_forms(shape, n):
    """Physical Q1 radial, cross, axial integrals on the original carrier mesh."""
    radial, axial=mapped_nodes(n,shape.length,28.,22.)
    x,t=radial[0],axial[1]/shape.length
    dx, dw=np.diff(x)[None,:],np.diff(axial[1])[:,None]
    rr=np.zeros((n,n,4,4)); rw=np.zeros_like(rr); ww=np.zeros_like(rr)
    points,wt=leggauss(3)
    for b,wb in zip((points+1)/2,wt/2):
        radius,slope=shape.radius(axial[1][:-1,None]+b*dw)
        for a,wa in zip((points+1)/2,wt/2):
            s=x[None,:-1]+a*dx
            r=radius*s
            dr=np.array([b-1,1-b,-b,b])/(dx*radius)[...,None]
            dz=np.array([a-1,-a,1-a,a])/dw[...,None]-(s*slope)[...,None]*dr
            vol=(4*np.pi*wa*wb*r*r*radius*dx*dw)[...,None,None]
            rr+=vol*dr[...,None]*dr[...,None,:]
            ww+=vol*dz[...,None]*dz[...,None,:]
            rw+=vol*(dr[...,None]*dz[...,None,:]+dz[...,None]*dr[...,None,:])
    return x,t,tuple(g.reshape(-1,4,4) for g in (rr,rw,ww))


class Carrier:
    def __init__(self,shape,nr,nz,n):
        self.shape,self.nr,self.nz=shape,nr,nz
        # Only the unrepresented receiving-side half of the planar cap phase delay.
        (pieces,self.volumes,_,self.stiff,self.feed,self.bulk,self.free)=assemble(
            0.,n,WIDTH/12,28.,22.,lambda step:shape,return_system=True,sealed_side=True)
        index=np.arange(self.stiff.shape[0])
        for name,ids,_,_ in pieces:
            if name=='cap_layer':
                index[ids[:,2:].ravel()]=ids[:,:2].ravel()
        _,did=np.unique(index,return_inverse=True)
        self.density_free=np.ones(max(did)+1,bool)
        self.density_free[did[np.r_[self.feed,self.bulk]]]=False
        self.parts=[]
        rows=[];cols=[]
        for name,ids,K,_ in pieces:
            densityids=did[ids]
            self.parts.append((name,ids,K,densityids))
            rows.append(np.broadcast_to(ids[:,:,None],K.shape).ravel())
            cols.append(np.broadcast_to(ids[:,None,:],K.shape).ravel())
        self.rows,self.cols=np.concatenate(rows),np.concatenate(cols)
        self.drows,self.dcols=did[self.rows],did[self.cols]
        x,t,self.tensor=core_forms(shape,n)
        bx,bt=transfer(nr,x),transfer(nz,t)
        dx,dt=transfer(nr,x,True),transfer(nz,t,True)
        self.means=sparse.kron(bt,bx,format='csr')
        X=sparse.kron(bt,dx,format='csr')
        T=sparse.kron(dt,bx,format='csr')
        tm=(t[:-1]+t[1:])/2
        xm=(x[:-1]+x[1:])/2
        r,slope=shape.radius(tm[:,None]*shape.length)
        self.GradR=sparse.diags(np.repeat(1/r.ravel(),n))@X
        self.GradW=T/shape.length-sparse.diags((xm[None,:]*slope/r).ravel())@X
        self.core_x,self.core_t=x,t
        self.f=np.ones(max(did)+1)
        self.waist_t=np.linspace(0,1,4097)[np.argmin(shape.profile(np.linspace(0,1,4097))[0])]

    def response(self, order, bias):
        gr, gw=self.GradR@order.ravel(),self.GradW@order.ravel()
        r2=gr*gr+gw*gw
        Qrr,Qrw,Qww=self.tensor
        Q=gr[:,None,None]**2*Qrr+(gr*gw)[:,None,None]*Qrw+gw[:,None,None]**2*Qww
        phase=np.zeros(self.stiff.shape[0]);phase[self.bulk]=bias
        f=self.f.copy()
        I=ORDER_INERTIA
        for it in range(40):
            matrices=[]
            for name,ids,K,did in self.parts:
                avg=f[did].mean(axis=1)
                rho=avg**2
                if name=='core':
                    P=rho[:,None,None]*(K-(I/(rho+I*r2))[:,None,None]*Q)
                else:
                    P=avg[:,None,None]**(4 if name=='cap_layer' else 2)*K
                matrices.append(P.ravel())
            A=sparse.coo_matrix((np.concatenate(matrices),(self.rows,self.cols)),
                                shape=self.stiff.shape).tocsr()
            phase[self.free]=0.
            phase[self.free]=spsolve(A[self.free][:,self.free],-(A@phase)[self.free])
            grad=np.zeros_like(f); hvalues=[]; energy=0.
            for name,ids,K,did in self.parts:
                theta=phase[ids]-phase[ids][:,:1]
                avg=f[did].mean(axis=1)
                rho=avg**2
                def half(mat):
                    return .5*np.einsum('ni,nij,nj->n',theta,mat,theta)
                e0=half(K)
                if name=='cap_layer':
                    e=avg**4*e0
                    g=4*avg**3*e0
                    H=12*rho*e0
                    cap_current=-np.sum(avg[:,None]**4*np.einsum('nij,nj->ni',K[:,:2,:],theta))
                    other=0.
                else:
                    if name=='core':
                        denom=rho+I*r2
                        eq=half(Q)
                        dN=e0-I**2*r2*eq/denom**2
                        e=rho*(e0-I*eq/denom)
                        g=2*avg*dN
                        H=2*dN+8*rho*I**2*r2*eq/denom**3
                        ar=half(2*gr[:,None,None]*Qrr+gw[:,None,None]*Qrw)
                        aw=half(2*gw[:,None,None]*Qww+gr[:,None,None]*Qrw)
                        dgr=-rho*I/denom*(ar-2*I*gr*eq/denom)
                        dgw=-rho*I/denom*(aw-2*I*gw*eq/denom)
                        self.kinetic_density=rho
                        self.frozen_mass=I*rho/denom
                        axis=f[did].reshape(len(self.core_t)-1,len(self.core_x)-1,4)[:,0][:,[0,2]]
                        cs=np.r_[axis[:,0],axis[-1,1]]
                        axis_n=np.interp(self.waist_t,self.core_t,cs)**2
                    else:
                        e=rho*e0
                        g=2*avg*e0
                        H=2*e0
                    v0,_=self.volumes[name]
                    ff=f[did]-f[did][:,:1]
                    e+=.5*np.einsum('ni,nij,nj->n',ff,K,ff)
                    e+=v0*(rho-1)**2*(rho**3+2*rho**2+3*rho+4)/20
                    g+=v0*avg*(rho**4-1)/2
                    H+=v0*(9*rho**4-1)/2
                    np.add.at(grad,did,np.einsum('nij,nj->ni',K,ff))
                    other=K
                energy+=e.sum()
                np.add.at(grad,did,np.broadcast_to(g[:,None]/4,did.shape))
                hvalues.append((other+H[:,None,None]*np.ones((4,4))/16).ravel())
            H=sparse.coo_matrix((np.concatenate(hvalues),(self.drows,self.dcols)),
                                shape=(len(f),len(f))).tocsr()
            free=self.density_free
            ch=spsolve(H[free][:,free],-grad[free])
            if max(abs(ch))<1e-9:
                break
            f[free]+=ch
            if min(f)<.55:
                raise ArithmeticError('Finite-density carrier left near-unit branch')
        else:
            raise ArithmeticError('Density relaxation exceeded iteration count')
        self.f=f
        currents=A@phase
        intake=currents[self.bulk].sum()
        if max(abs(intake+currents[self.feed].sum()),abs(intake-cap_current))>8e-7*max(abs(intake),1e-6):
            raise ArithmeticError('Sealed two-sided carrier currents disagree')
        order_gradient=(self.GradR.T@dgr+self.GradW.T@dgw).reshape(order.shape)
        axis_d=np.interp(self.waist_t,np.linspace(0,1,self.nz+1),
                         np.gradient(order[:,0],self.shape.length/self.nz))
        return energy,order_gradient,dict(carrier_energy=energy,current=intake,terminal_current=cap_current,
                            axis_density=axis_n,axis_order_gradient_squared=axis_d**2,
                            H_pipi_axis=1/I+axis_d**2/axis_n,density_iterations=it+1,
                            density_min=min(f)**2,carrier_gradient_norm=np.linalg.norm(order_gradient))


class CoupledLoads:
    def __init__(self,shape,nr,nz,field_mesh=14,port=4.2,wave_seed=None,target_frequency=None):
        self.carrier=Carrier(shape,nr,nz,field_mesh)
        self.wave=OccupiedWave(shape,nr,nz,field_mesh,port,wave_seed,target_frequency)

    def evaluate(self,chi,action=1.,bias=BIAS):
        order=np.asarray(chi)
        if not np.isfinite(order).all() or order.shape != (self.carrier.nz+1,self.carrier.nr+1):
            raise ValueError('Invalid material field on the normalized meridian')
        if action<0:
            raise ValueError('Action must be nonnegative')
        ce,cg,info=self.carrier.response(order,bias)
        if not action:
            info['wave_energy']=0.
            return ce,cg,info
        frequency,wg,readout=self.wave.response(order)
        info.update(readout,wave_energy=action*frequency,wave_gradient_norm=np.linalg.norm(wg))
        return ce+action*frequency,cg+action*wg,info
