"""Simultaneous shape, lattice and order descent with both native carrier baths.

A line search may move the wall before the internal label problem has settled.
A failed line search or a small nodal determinant never qualifies as equilibrium.
"""
import argparse
import contextlib
import hashlib
import io

import numpy as np
from scipy.optimize import minimize
from pathlib import Path

from mouth_joint_collar_release import Release
from mouth_exchange_carrier import BIAS
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import Order
from mouth_supercore_shape import emit, order_matrices


class Simultaneous(Release):
    def initialize(self):
        sh=self.model(self.start)
        pair=self.pair(sh)
        c=self.initial_chi.copy()
        owner=Order(sh,self.a.order)
        c[~owner.free]=owner.seed[~owner.free]
        y=pair.identity.copy() if self.initial_labels is None else self.initial_labels
        _,(_,_,out)=self.outside(sh,np.zeros_like(c),True)
        self.base_labels=y.copy()
        self.base_chi=c.copy()
        self.label_mask=pair.mask.copy()
        self.cm=owner.free
        # Positive diagonal scaling sets numerical units only; retain each
        # unconstrained label, including those instantaneously in liquid.
        self.yscale=np.full(y.shape,np.sqrt(4*np.pi*sh.neck**2))
        self.yscale[:,:,0]*=np.sqrt(self.a.na/self.a.nt)
        self.yscale[:,:,1]*=np.sqrt(self.a.nt/self.a.na)
        self.cscale=np.sqrt(2*owner.stiff.diagonal()+2*owner.lump).reshape(c.shape)
        self.sscale=np.array([80.,80.,90.,90.,90.])
        self.slice=5+self.label_mask.sum()
        self.dcache=None
        self.amp=None
        zz=np.r_[np.zeros(5),np.zeros(self.label_mask.sum()),c[self.cm]*self.cscale[self.cm]]
        return zz

    def unpack(self,z):
        v=self.start+z[:5]/self.sscale
        labels=self.base_labels.copy()
        labels[self.label_mask]+=z[5:self.slice]/self.yscale[self.label_mask]
        c=self.base_chi.copy()
        c[self.cm]=z[self.slice:]/self.cscale[self.cm]
        return v,labels,c

    def objective(self,z):
        if self.dcache is not None and np.array_equal(z,self.dcache[0]):
            return self.dcache[1:3]
        v,labels,c=self.unpack(z)
        sh=self.model(v)
        pair=self.pair(sh)
        if pair._corners(labels)<=0 or sh.wall()['minimum_radius']<.4 or sh.exit>8 or sh.length>6:
            # Reject infeasible line-search candidates. No elastic modulus is
            # added to liquid cells; no such step is accepted as equilibrium.
            step=z-self.last_z
            return 1.e10+step@step,2*step
        owner=Order(sh,self.a.order)
        core=core_source(sh,self.co,self.ref,self.co,self.t,self.x)
        coeff=pair.evaluate(labels,self.inner(c),self.outer(np.zeros_like(c)))
        out,(_,_,outer)=self.outside(sh,coeff['outer_source'],True)
        work=pair.evaluate(labels,self.inner(c),self.outer(outer),diagnostics=True)
        flow=self.wound(sh,c)
        with contextlib.redirect_stdout(io.StringIO()):
            amp,phase,carrier=flow.solve(BIAS,self.amp)
        order_E,og=owner.energy_gradient(c)
        _,fg=order_matrices(flow,c,amp,phase)(c)
        cg=og+fg+core+work['inner_source']
        energy=order_E+(core*c).sum()+out['energy']+work['energy']+carrier['excess_energy']
        h=.0003
        sg=np.array([(self.fixed_fields(v+h*d,c,labels,amp,phase)-
                      self.fixed_fields(v-h*d,c,labels,amp,phase))/(2*h) for d in np.eye(5)])
        jac=np.r_[sg/self.sscale,work['gradient'][self.label_mask]/self.yscale[self.label_mask],
                  cg[self.cm]/self.cscale[self.cm]]
        stats=dict(shape_force=(-sg).tolist(),lattice_nodal_residual=work['residual'],
                   min_corner_inverse_J=work['min_inverse_corner_det'],min_J=work['min_J'],
                   max_J=work['max_J'],min_stretch=work['min_stretch'],max_stretch=work['max_stretch'],
                   order_projected_residual=float(max(abs((c-np.clip(c-cg,0,1))[self.cm]))),
                   bath_order_residual=out['residual'],carrier_residual=carrier['residual'],
                   sharp_current=carrier['total_current'],depth=sh.length,exit=sh.exit,
                   waist=sh.wall()['minimum_radius'],seam_work=work['seam_label_work_sum'],
                   receiver_edge_work=work['receiver_edge_label_work_sum'])
        self.dcache=z.copy(),energy,jac,stats,outer,amp
        return energy,jac

    def accept(self,z):
        energy,_=self.objective(z)
        if self.dcache is None or not np.array_equal(z,self.dcache[0]):
            raise ArithmeticError('Optimizer accepted an infeasible candidate')
        _,_,_,stats,outer,self.amp=self.dcache
        self.last_z=z.copy()
        v,y,c=self.unpack(z)
        sh=self.model(v)
        path=self.a.folder/'state.npz'
        tmp=path.with_suffix('.partial.npz')
        np.savez_compressed(tmp,radius=4,winding=1,shape_params=sh.params,shape_extra=sh.extra,
                            labels=y,order=c,outer_order=outer,outer_t=self.ot,outer_x=self.ox,
                            registration=self.co)
        tmp.replace(path)
        self.count+=1
        emit(kind='accepted',iteration=self.count,energy=energy,
             state_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),**stats)

    def descend(self):
        z=self.initialize()
        self.last_z=z.copy()
        self.accept(z)
        limits=np.array([[-3.8,1.6],[-7.,2.7],[-2.,2.],[-2.,2.],[-2.,2.]])
        shape_bounds=list(zip(*((limits-self.start[:,None])*self.sscale[:,None]).T))
        bounds=shape_bounds+[(None,None)]*self.label_mask.sum()+list(zip(
            np.zeros(self.cm.sum()),self.cscale[self.cm]))
        opt=minimize(self.objective,z,jac=True,bounds=bounds,method='L-BFGS-B',callback=self.accept,
                     options=dict(gtol=1.e-7,ftol=3.e-16,maxls=55,maxiter=2500,maxcor=20))
        self.accept(opt.x)
        vals=self.dcache[3]
        v=self.unpack(opt.x)[0]
        active=bool(np.any((v-limits[:,0]<1e-5)|(limits[:,1]-v<1e-5)))
        passed=bool(max(abs(np.array(vals['shape_force'])))<.01 and vals['lattice_nodal_residual']<3.e-4
                    and vals['order_projected_residual']<5.e-5 and not active)
        emit(kind='endpoint',stationary_in_released_space=passed,active_shape_bound=active,
             optimizer=str(opt.message),accepted_count=self.count,**vals)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('folder',type=Path)
    for name in ('order','na','nt','flow','polar'):
        p.add_argument('--'+name,type=int,required=True)
    Simultaneous(p.parse_args()).descend()
