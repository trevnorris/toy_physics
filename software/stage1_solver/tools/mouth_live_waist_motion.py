"""Projected Route A order momentum and inertia at fixed carrier density on a prepared shape."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_live_waist_fields import ORDER_INERTIA, BIAS


def integrated_basis(shape, x, t):
    w=t*shape.length
    dx,dw=np.diff(x)[None,:],np.diff(w)[:,None]
    drint=np.zeros((len(t)-1,len(x)-1,4))
    dwint=np.zeros_like(drint)
    roots,weights=leggauss(3)
    for b,wb in zip((1+roots)/2,weights/2):
        radius,slope=shape.radius(w[:-1,None]+b*dw)
        for a,wa in zip((1+roots)/2,weights/2):
            u=x[None,:-1]+a*dx
            dr=np.array([b-1,1-b,-b,b])/(dx*radius)[...,None]
            dz=np.array([a-1,-a,1-a,a])/dw[...,None]-(u*slope)[...,None]*dr
            vol=(4*np.pi*wa*wb*(u*radius)**2*radius*dx*dw)[...,None]
            drint+=vol*dr
            dwint+=vol*dz
    return drint.reshape(-1,4),dwint.reshape(-1,4)


class CarrierMotion:
    """The prepared density must first be set by Carrier.response(order, bias)."""
    def __init__(self, carrier, order, eta, bias=BIAS):
        self.c=carrier
        self.gr=carrier.GradR@order.ravel()
        self.gw=carrier.GradW@order.ravel()
        self.er=carrier.GradR@eta.ravel()
        self.ew=carrier.GradW@eta.ravel()
        self.eta=carrier.means@eta.ravel()
        self.vol=carrier.volumes['core'][0]
        self.grad_basis=integrated_basis(carrier.shape,carrier.core_x,carrier.core_t)
        matrices=[]
        for name,ids,K,dids in carrier.parts:
            rho=carrier.f[dids].mean(axis=1)**2
            if name=='core':
                self.ids=ids
                self.rho=rho
                rr,rw,ww=carrier.tensor
                Q=self.gr[:,None,None]**2*rr+(self.gr*self.gw)[:,None,None]*rw+self.gw[:,None,None]**2*ww
                M=rho[:,None,None]*(K-self.coeff(self.gr,self.gw)[:,None,None]*Q/rho[:,None,None])
            else:
                M=rho[:,None,None]**(2 if name=='cap_layer' else 1)*K
            matrices.append(M.ravel())
        self.A=sparse.coo_matrix((np.concatenate(matrices),(carrier.rows,carrier.cols)),shape=carrier.stiff.shape).tocsr()
        self.phase=np.zeros(carrier.stiff.shape[0])
        self.phase[carrier.bulk]=bias
        self.phase[carrier.free]=spsolve(self.A[carrier.free][:,carrier.free],-(self.A@self.phase)[carrier.free])
        theta=self.phase[self.ids]-self.phase[self.ids][:,:1]
        self.theta_moments=[np.einsum('ci,cij,cj->c',theta,g,theta) for g in carrier.tensor]
        self.theta_grad_integral=[np.einsum('ci,ci->c',theta,b) for b in self.grad_basis]
        self.base=self.kinetic(self.gr,self.gw)

    def coeff(self,gr,gw):
        return ORDER_INERTIA*self.rho/(self.rho+ORDER_INERTIA*(gr*gr+gw*gw))

    def kinetic(self,gr,gw):
        rr,rw,ww=self.theta_moments
        return -.5*np.sum(self.coeff(gr,gw)*(gr*gr*rr+gr*gw*rw+gw*gw*ww))

    def held(self, a):
        gr=self.gr+a*self.er
        gw=self.gw+a*self.ew
        coef=self.coeff(gr,gw)
        vr,vw=self.theta_grad_integral
        return dict(energy_change=self.kinetic(gr,gw)-self.base,
                    frozen_mass=np.sum(coef*self.vol*self.eta**2),
                    momentum_offset=np.sum(coef*self.eta*(gr*vr+gw*vw)))

    def matched_mass(self):
        out=self.held(0.)
        dr,dw=self.grad_basis
        local=(self.coeff(self.gr,self.gw)*self.eta)[:,None]*(self.gr[:,None]*dr+self.gw[:,None]*dw)
        h=np.zeros_like(self.phase)
        np.add.at(h,self.ids,local)
        free=self.c.free
        dp=spsolve(self.A[free][:,free],h[free])
        induced=h[free]@dp
        if induced < -1.e-15:
            raise ArithmeticError('Nonpositive passive phase-elimination increment')
        out.update(induced_phase_mass=induced,matched_mass=out['frozen_mass']+induced)
        return out
