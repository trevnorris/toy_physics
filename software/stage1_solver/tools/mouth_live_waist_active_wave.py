"""Exact zero-shear cells of the outgoing waist operator, preserving its global material chart."""
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from mouth_live_waist_wave import OccupiedWave
from mouth_exchange_shear import solve


class ActiveWave(OccupiedWave):
    @classmethod
    def install(cls, loads):
        """Adopt the existing pole and all pinned operator data; no changes to old launches."""
        replacement=cls.__new__(cls)
        replacement.__dict__.update(loads.wave.__dict__)
        feed=sparse.diags((replacement.uniform[0][:,1]<=0).astype(float))
        replacement.cut_k=(feed@replacement.cut_k@feed).tocsc()
        replacement.cut_m=(feed@replacement.cut_m@feed).tocsc()
        loads.wave=replacement
        return loads

    def pole(self, weights):
        if weights.min()>0:
            self.active_info=dict(zero_cells=0,inactive_nodes=0,wave_components=1,detached_mode=False)
            return super().pole(weights)
        k=self.matrix(weights[:,None,None]*self.kc)
        m=self.matrix(weights[:,None,None]*self.mc)
        mass=self.cut_m+m; stiffness=self.cut_k+k
        mass.eliminate_zeros()
        active=mass.diagonal()>0
        ids=np.flatnonzero(active)
        ncomp,labels=connected_components(mass[ids][:,ids],directed=False)
        prior,vector=self.seed
        past=np.array([np.vdot(vector[ids[labels==i]],mass[ids[labels==i]][:,ids[labels==i]]
                       @vector[ids[labels==i]]).real for i in range(ncomp)])
        selected=int(np.argmax(past))
        if past[selected] <= 0 or past[selected] < .88*sum(past):
            raise ArithmeticError('Ambiguous occupied wave between disconnected shear regions')
        ids=ids[labels==selected]
        K=stiffness[ids][:,ids]; M=mass[ids][:,ids]; inside=m[ids][:,ids]
        pts=self.uniform[0][ids]
        detached=not np.any(abs(pts[:,0]-self.radius)<1e-9)
        self.active_info=dict(zero_cells=int(sum(weights==0)),inactive_nodes=int(sum(~active)),
                              wave_components=ncomp,detached_mode=detached)
        if detached:
            if len(ids)<2:
                raise ArithmeticError('No finite-dimensional oscillating pocket to continue')
            vals,vec=eigsh(K,k=1,M=M,sigma=prior.real**2,which='LM',v0=vector[ids].real)
            omega=np.sqrt(vals[0])+0j
            v=vec[:,0].astype(complex)
            residual=np.linalg.norm(K@v-omega**2*(M@v))/np.linalg.norm(K@v)
            if residual>1e-8:
                raise ArithmeticError('Detached pocket eigensystem not resolved')
            data=dict(frequency=omega.real,amplitude_loss=0.,confined_kinetic_fraction=1.,
                      outward_power_over_interior_energy=0.,pole_residual=residual)
        else:
            data,v,_,=solve(self.shape,self.mesh,self.radius,0.,prior=prior,incoming_vector=vector[ids],
                            shape_work=False,supplied_system=(pts,K,M,inside),target_frequency=self.target)
            omega=complex(data['frequency'],-data['amplitude_loss'])
        full=np.zeros(len(active),complex)
        full[ids]=v
        overlap=abs(np.vdot(vector,mass@full))**2/(np.vdot(vector,mass@vector).real)
        if overlap<.88:
            raise ArithmeticError(f'Occupied zero-shear continuation overlap is {overlap}')
        self.seed=omega,full
        self.target=omega.real
        return data,full,mass,omega,overlap

    def response(self, order):
        weight=np.asarray(self.project@order.T.ravel()).ravel()
        if min(weight)>0:
            value,grad,data=super().response(order)
            data.update(self.active_info)
            return value,grad,data
        if min(weight)<0:
            raise ValueError('Negative shear stiffness on a physical cell')
        overlap_floor=1.
        if self.seed is None:
            self.pole(np.ones_like(weight))
            for s in (.25,.5,.75):
                *_,ov=self.pole(1+s*(weight-1))
                overlap_floor=min(ov,overlap_floor)
        out,v,mass,omega,ov=self.pole(weight)
        vec=np.r_[v,0j][self.ids]
        kcell=np.einsum('ci,cij,cj->c',vec.conj(),self.kc,vec).real
        mcell=np.einsum('ci,cij,cj->c',vec.conj(),self.mc,vec).real
        ktotal=np.vdot(v,self.cut_k@v).real+weight@kcell
        self.strain_weights=np.asarray(self.project.T@kcell).reshape(self.nr+1,self.nz+1).T/ktotal
        self.inertia_weights=np.asarray(self.project.T@mcell).reshape(self.nr+1,self.nz+1).T/np.vdot(v,mass@v).real
        diff=np.einsum('ci,cij,cj->c',vec,self.kc-omega**2*self.mc,vec)
        boundary_v=v[self.port.ids]
        derivative=(self.port.block(omega+2e-6)-self.port.block(omega-2e-6))/(4e-6)
        den=2*omega*(v@(mass@v))+boundary_v@(derivative@boundary_v)
        grad=np.asarray(self.project.T@(diff/den).real).reshape(self.nr+1,self.nz+1).T.copy()
        return out['frequency'],grad,dict(wave_frequency=out['frequency'],
                   wave_pocket_fraction=out['confined_kinetic_fraction'],
                   wave_loss_rate=out['outward_power_over_interior_energy'],
                   pole_amplitude_loss=out['amplitude_loss'],mode_overlap=min(overlap_floor,ov),
                   pole_residual=out['pole_residual'],wave_minimum_cell_order=float(min(weight)),**self.active_info)
