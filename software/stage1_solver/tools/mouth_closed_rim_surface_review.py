"""Independent sharp enclosing-surface flux of a radius-four rim.

Differentiate the physical carrier and order energy under a local affine
axial displacement by complex step.  All integrations use piecewise Q1
fields reconstructed from nodal values.  The closed surface is outside the
free order band and beyond the polar deformation at both caps and the wall
within its axial interval.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq

from mouth_joint_polar_shape import Trial
from mouth_joined_front_audit import Front
from mouth_polar_exit_readout import accepted
from mouth_supercore_joined_front import KAPPA, ORDER_INERTIA, WELL, check_guard, emit, file_sha
from mouth_wound_carrier import straight_vortex


def quadrature(knots, n):
    knots = np.unique(np.asarray(knots))
    x, w = leggauss(n)
    widths = np.diff(knots)
    return ((knots[1:]+knots[:-1])[:,None]/2+widths[:,None]*x/2).ravel(), (widths[:,None]*w/2).ravel()


def bracket(axis, pt):
    i = np.clip(np.searchsorted(axis, pt, side='right')-1, 0, len(axis)-2)
    return axis[i], axis[i+1]


def stress_energy(amp_grad, phi_grad, order_grad, n, chi, swirl2):
    """Local density eliminated with respect to the order momentum."""
    sq = lambda a: np.sum(a*a, axis=1)
    projection = np.sum(phi_grad*order_grad, axis=1)
    return (.5*(sq(amp_grad)+n*swirl2) + (n**5-5*n+4)/20
            + .5*n*(sq(phi_grad)-ORDER_INERTIA*projection**2/(n+ORDER_INERTIA*sq(order_grad)))
            + .5*KAPPA*sq(order_grad)+WELL*(chi*(1-chi))**2)


class Reader:
    def __init__(self, model, base_order, driven, resting):
        self.model = model
        self.vortex = straight_vortex(model.winding,28.)
        meridian = RegularGridInterpolator(tuple(np.linspace(0,1,k) for k in base_order.shape),base_order)
        self.fields = {}
        for slab in (1,2):
            for outer in (0,1):
                axes = (model.rr[outer],model.a,model.z[slab])
                x,_,z = np.meshgrid(*axes,indexing='ij')
                b = model.shape.radius(z)[0]
                if outer:
                    ref = model.terminal_profile.outside(b+x*(28-b),z,slab==2)
                elif slab == 1:
                    ref = meridian((z/model.shape.length,x))
                else:
                    ref = model.terminal_profile(x,z-model.shape.length)
                self.fields[slab,outer] = []
                ids = model.ids[slab,outer]
                for u in (driven,resting):
                    nodal = u[:,ids].copy()
                    nodal[2] += ref
                    self.fields[slab,outer].append(RegularGridInterpolator(axes,np.moveaxis(nodal,0,-1)))

    def sample(self, slab, outer, radius, angle, depth, normal):
        r,a,z = [q.ravel() for q in np.broadcast_arrays(radius,angle,depth)]
        sh = self.model.shape
        b,slope = sh.radius(z)
        width = 28-b if outer else b
        x = (r-b)/width if outer else r/b
        coords = np.stack((x,a,z),axis=1)
        F,Fp = self.vortex(r*np.sin(a))
        rows = []
        for interp in self.fields[slab,outer]:
            vals = interp(coords)
            delta = []
            for i,axis in enumerate(interp.grid):
                left,right=bracket(axis,coords[:,i])
                end=coords.copy()
                end[:,i]=left
                begin=interp(end)
                end[:,i]=right
                delta.append((interp(end)-begin)/(right-left)[:,None])
            dr = delta[0]/width[:,None]
            gradients = np.stack((dr,delta[1]/r[:,None],
                                  delta[2]-dr*(slope*(1-x if outer else x))[:,None]),axis=-1)
            amp,phase,chi = vals.T
            amp_grad = F[:,None]*gradients[:,0] + (Fp*amp)[:,None]*np.stack((np.sin(a),np.cos(a),0*a),axis=1)
            phi_grad,chi_grad = gradients[:,1],gradients[:,2]
            n = (F*amp)**2
            c = np.einsum('ij,ij->i',phi_grad,chi_grad)*ORDER_INERTIA/(n+ORDER_INERTIA*np.sum(chi_grad**2,axis=1))
            j = n[:,None]*(phi_grad-c[:,None]*chi_grad)
            # F = I + i*eps e_z tensor n_hat, Cauchy flux = -d(comoving energy)/d eps.
            s = 1j*1e-25
            jac = 1+s*normal[2]
            pull = lambda grad: grad - (s*grad[:,2]/jac)[:,None]*normal
            fluid = stress_energy(pull(amp_grad),pull(phi_grad),pull(chi_grad),n,chi,(1/(r*np.sin(a)))**2)
            output = np.stack((j@normal,-np.imag(jac*fluid)/s.imag),axis=1)
            rows.append(output)
        return rows[0]-rows[1]

    def sealed_increment(self,h0,h1,n,ramp_start=None):
        # Outward normals of the two material banks at the same sloped wall.
        # The orientation length cancels the area normalization, so evaluate
        # radial and axial fluxes separately and contract with (1,-db/dz).
        m=self.model
        L=m.shape.length
        z0,z1=L-h1,L-h0
        axis=m.z[1]
        z,wz=quadrature([z0,*axis[(axis>z0)&(axis<z1)],z1],n)
        a,wa=quadrature(m.a,n)
        zz,aa=np.meshgrid(z,a,indexing='ij')
        b,s=m.shape.radius(zz)
        weight=4*np.pi*b**2*np.sin(aa)*wz[:,None]*wa[None,:]
        if ramp_start is not None:
            weight *= np.clip((zz-z0)/(ramp_start-z0),0,1)
        sides=[]
        for outer in (0,1):
            tr=self.sample(1,outer,b,aa,zz,np.array([1.,0.,0.]))
            tz=self.sample(1,outer,b,aa,zz,np.array([0.,0.,1.]))
            sides.append((1-2*outer)*weight.ravel()@(tr-s.ravel()[:,None]*tz))
        return np.sum(sides,axis=0)

    def integrate(self, R, h, n):
        m = self.model
        sh = m.shape
        L = sh.length
        if L-h <= .62*L or h<=m.h:
            raise ValueError('Enclosure crosses polar chart or terminal-order support')
        a,wa=quadrature(m.a,n)
        results={}
        for slab,z,sign in ((1,L-h,-1),(2,L+h,1)):
            b=float(sh.radius(z)[0])
            for outer in (0,1):
                grid=b+(28-b)*m.rr[1] if outer else b*m.rr[0]
                if outer:
                    grid=np.r_[grid[grid<R],R]
                r,wr=quadrature(grid,n)
                rr,aa=np.meshgrid(r,a,indexing='ij')
                values=self.sample(slab,outer,rr,aa,z,np.array([0.,0.,sign]))
                weight=4*np.pi*rr**2*np.sin(aa)*wr[:,None]*wa[None,:]
                results[f'cap_{slab}_{outer}']=weight.ravel()@values
        for slab,z0,z1 in ((1,L-h,L),(2,L,L+h)):
            if np.min(R-sh.radius(np.linspace(z0,z1,501))[0])<=0:
                raise ValueError('Side intersects the bent wall')
            axis=m.z[slab]
            cuts=np.r_[z0,axis[(axis>z0)&(axis<z1)],z1]
            roots=[]
            for frac in m.rr[1][1:-1]:
                def diff(z):
                    b=sh.radius(z)[0]
                    return b+(28-b)*frac-R
                for lo,hi in zip(cuts[:-1],cuts[1:]):
                    if diff(lo)*diff(hi)<0:
                        roots.append(brentq(diff,lo,hi))
            z,wz=quadrature([*cuts,*roots],n)
            zz,aa=np.meshgrid(z,a,indexing='ij')
            vals=self.sample(slab,1,R,aa,zz,np.array([1.,0.,0.]))
            weights=4*np.pi*R**2*np.sin(aa)*wz[:,None]*wa[None,:]
            results[f'side_{slab}']=weights.ravel()@vals
        return {k:v.tolist() for k,v in results.items()}, np.sum(list(results.values()),axis=0)


def run(args):
    check_guard(args.controls)
    prep,row,x,d,u = accepted(args.state,args.snapshot)
    base=Front(args.state,tuple(prep['field_grid']),1)
    trial=Trial(base,args.state,x,(d.shape[1]-1)//3,d.shape[0]-1,False)
    expected=row['sha256']
    if args.relaxed:
        check_guard(args.relaxed)
        with np.load(args.relaxed) as f:
            if str(f['input_sha'])!=expected or str(f['state_sha'])!=file_sha(args.state) or not np.array_equal(f['x'],x):
                raise ValueError('Angular equilibrium is from a different fixed wall or preparation')
            u=f['field'].copy()
            expected=file_sha(args.relaxed)
    with np.load(args.controls) as f:
        if str(f['checkpoint_sha'])!=expected:
            raise ValueError('Mismatched drive−rest controls')
        mode=trial.model
        if 'grid' in f:
            from mouth_polar_joined import PolarCoupled
            rows=list(map(json.loads,(args.controls.parent/'stdout').read_text().splitlines()))
            if (not args.relaxed or rows[0]['endpoint']!=expected or
                    not any(r.get('kind')=='matched_states' and r['sha256']==file_sha(args.controls) for r in rows)):
                raise ValueError('Unverified nested controls')
            mode=PolarCoupled(trial.shape,base.c,1,tuple(f['grid']),polar_shift=x[5],
                              precondition=False,terminal_profile=mode.terminal_profile)
            u=f['wound'].copy()
        elif not np.array_equal(u,f['wound']):
            raise ValueError('Wound field changed on the original carrier grid')
        rest=f['rest'].copy()
    reader=Reader(mode,base.c,u,rest)
    ids=mode.ids[1,0]
    delta=u[2]-rest[2]
    if (np.max(abs(delta[ids[:,:,mode.z[1] < mode.shape.length-mode.h-1e-13]])) or
            np.max(abs(delta[mode.ids[2,0][:,:,mode.z[2] > mode.shape.length+mode.h+1e-13]])) or
            any(np.any(delta[mode.ids[j,1]]) for j in (1,2))):
        raise ArithmeticError('Direct drive-induced crystal surface stress has support beyond the front')
    emit(kind='sources',state=file_sha(args.state),snapshot=row['sha256'],controls=file_sha(args.controls),
         relaxed=None if not args.relaxed else file_sha(args.relaxed),face_model='physical Q1 strong form and signed affine energy')
    if args.slit_only:
        L=mode.shape.length
        for n in (7,10):
            for h in (.15,.22,.27):
                emit(kind='sharp_slit',h=h,quadrature=n,
                     outward_into_slit=reader.sealed_increment(0.,h,n).tolist())
            iz=np.searchsorted(mode.z[1],L-mode.h-1e-12)
            for k in (2,3):
                lo1,lo2=mode.z[1][iz-k:iz-k+2]
                emit(kind='weak_slit',k=k,quadrature=n,bottom_ramp=[lo1,lo2],
                     outward_into_slit=reader.sealed_increment(0.,L-lo1,n,lo2).tolist())
        return
    for h0,h1 in ((0.,.15),(0.,.22),(0.,.27),(.15,.22),(.22,.27),(.15,.27)):
        emit(kind='incremental_sealed_wall',h0=h0,h1=h1,inner_out_minus_outer_in=reader.sealed_increment(h0,h1,8).tolist())
    for radius in (4.5,5.5):
        for h in (.15,.22,.27):
            vlow,low=reader.integrate(radius,h,4)
            pieces,v=reader.integrate(radius,h,8)
            emit(kind='closed_cylinder',radius=radius,h=h,columns=['outward_mass','outward_axial_momentum'],
                 all_faces=pieces,net=v.tolist(),degree4_to_8=(v-low).tolist())


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','snapshot','controls'):
        p.add_argument(name,type=Path)
    p.add_argument('--relaxed',type=Path)
    p.add_argument('--slit-only',action='store_true')
    run(p.parse_args())
