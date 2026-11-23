"""Release a fixed-wall interior angular deformation of the bonded R4 solid.

The two angular factors are G = h x^3[e+s(1-x)] and P2(cos(theta)).
e deforms the wall and its Eulerian order chart. s only moves the material.
The correction has compact depth support away from the mobile terminal front.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize

from mouth_bonded_front_adjoint import set_front
from mouth_bonded_material_motion import Bonded
from mouth_bonded_material_residual import polish
from mouth_joined_angular_work import bump
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_polar_bonded import PolarBonded, order_at, require_map
from mouth_polar_exit_readout import accepted
from mouth_supercore_joined_front import emit, file_sha


class AngularBonded(PolarBonded):
    internal_amplitude = 0.

    def _strains(self, q):
        u = self.free.copy()
        u[self.mask] = q/self.scale
        nd = u.reshape(-1,2)[self.id]
        pos = np.einsum('qi,qik->qk', self.basis, nd)
        B = np.eye(2) + np.stack([np.einsum('qi,qik->qk', v,nd)
                                   for v in (self.dr,self.dw)],axis=-1)
        r,z = self.r+pos[:,0],self.w+pos[:,1]
        hoop = self.h*r/self.r
        J = self.j*np.linalg.det(B)*(r/self.r)**2
        BK = B@self.K
        K = BK@B.transpose(0,2,1)
        return r,z,B,hoop,J,BK,K

    def _chart(self,r,z):
        sh = self.chart.shape
        b,bw,bww = sh.profile(z/self.L)
        x = r/b
        h,hw = bump(z,self.L)
        y = (z/self.L-.06)/.56
        hww = np.where((y>0)&(y<1),32*(1-6*y+6*y*y)/(.56*self.L)**2,0.)
        e,s = self.epsilon,self.internal_amplitude
        f = x**3*(e+s*(1-x))
        fx = 3*x**2*(e+s)-4*s*x**3
        fxx = 6*x*(e+s)-12*s*x*x
        l = bw/b
        G,Gr,Gz = h*f,h*fx/b,hw*f-h*fx*x*l
        Grr = h*fxx/b**2
        Grz = (hw*fx-h*l*(fx+x*fxx))/b
        Gzz = hww*f-2*hw*fx*x*l+h*(fxx*x*x*l*l+fx*x*(2*l*l-bww/b))
        return b,l,x,h,hw,G,Gr,Gz,Grr,Grz,Gzz

    def energy(self,q,diagnose=False):
        # In the compact material window the front correction vanishes; the
        # ordinary weak work elsewhere is included exactly once by Bonded.
        spherical = Bonded.energy(self,q,diagnose)
        if self.front is not None:
            tt,_,inc = self.front
            window = (tt[:-1]<.62)&(tt[1:]>.06)
            if np.any(inc[:-1][window]) or np.any(inc[1:][window]):
                raise ValueError('Interior angular deformation overlaps mobile order')
        r,z,B,hoop,J,BK,K = self._strains(q)
        if diagnose:
            return self._diagnose(q,spherical,r,z,hoop,J,K)
        active = (z>.06*self.L)&(z<.62*self.L)
        r,z,B,hoop,J,BK,K = (f[active] for f in (r,z,B,hoop,J,BK,K))
        rr,rz,zz = K[:,0,0],K[:,0,1],K[:,1,1]
        b,l,x,h,hw,G,Gr,Gz,Grr,Grz,Gzz = self._chart(r,z)
        angular = G/r
        require_map(Gr,angular)
        W0 = (rr+zz+2*hoop**2-4)/2 - np.log(J) + .25*(J-1)**2
        p0 = -1+.5*J*(J-1)
        chi,cr,cw = order_at(self,r,z)
        H2 = hoop**2
        invt = np.linalg.inv(B).transpose(0,2,1)
        A0 = BK + p0[:,None,None]*invt
        dA = np.zeros_like(B)
        dpos = np.zeros((len(r),2))
        wt = self.measure[active]
        energy_delta = deriv = 0.
        # Scalar derivatives of G, G_r, G_z for material amplitude.
        k = x**3*(1-x)
        kx = x**2*(3-4*x)
        mode,mode_r,mode_z = k*h,kx*h/b,k*hw-h*kx*x*l
        for p,dp,w in zip(*self.quad):
            a,g,v = 1+Gr*p,1+angular*p,Gz*p
            u = angular*dp
            jac = J*a*g*g
            press = -1+.5*jac*(jac-1)
            W = (a*a*rr+2*a*v*rz+(1+v*v)*zz+H2*(u*u+2*g*g)-4)/2 - np.log(jac)+.25*(jac-1)**2
            # Invert the fixed spatial order chart of the angularly deformed
            # wall at the new landing position of this material element.
            dest = r+G*p
            pre = r.copy()
            for _ in range(4):
                xp = pre/b
                wall_a = 1+3*self.epsilon*h*xp*xp*p/b
                pre -= (pre+self.epsilon*h*xp**3*p-dest)/wall_a
            xp = pre/b
            wall_a = 1+3*self.epsilon*h*xp*xp*p/b
            wall_v = self.epsilon*xp**3*(hw-3*h*l)*p
            if np.max(abs(pre+self.epsilon*h*xp**3*p-dest))>3e-13:
                raise ValueError('Fixed order chart inversion failed')
            cp,cpr,cpw = order_at(self,pre,z)
            chi_r = cpr*a/wall_a
            chi_z = cpw+cpr*(v-wall_v)/wall_a
            # W derivatives for the four independent entries of the angular map.
            wa = a*rr+v*rz+press/a
            wv = a*rz+v*zz
            wg = 2*H2*g+2*press/g
            wu = H2*u
            wg_total = wg*p+wu*dp
            wr = wa*p*Grr+wv*p*Grz+wg_total*(Gr-angular)/r + (H2*(u*u+2*g*g)+2*press)/r
            wz = wa*p*Grz+wv*p*Gzz+wg_total*Gz/r
            dpos[:,0] += w*(W*chi_r+cp*wr-W0*cr-chi*2*(H2+p0)/r)
            dpos[:,1] += w*(W*chi_z+cp*wz-W0*cw)
            D = np.empty_like(B)
            D[:,0,:] = (a*a)[:,None]*BK[:,0,:] + (a*v)[:,None]*BK[:,1,:]
            D[:,1,:] = (a*v)[:,None]*BK[:,0,:]+(1+v*v)[:,None]*BK[:,1,:]
            dA += w*(cp[:,None,None]*(D+press[:,None,None]*invt)-chi[:,None,None]*A0)
            energy_delta += w*(wt @ (cp*(W-W0)+(cp-chi)*W0))
            deriv += w*(wt @ (cp*(wa*p*mode_r+wv*p*mode_z+wg_total*mode/r)
                                  + W*cpr*mode*p/wall_a))
        ids = self.id[active].ravel()
        for axis in (0,1):
            row = wt[:,None]*(dpos[:,axis,None]*self.basis[active]+dA[:,axis,0,None]*self.dr[active]
                               +dA[:,axis,1,None]*self.dw[active])
            self.native[...,axis] += np.bincount(ids,weights=row.ravel(),minlength=self.npts).reshape(self.mask.shape[:2])
        self.angular_gradient = float(deriv)
        self.angular_energy_delta = float(energy_delta)
        return spherical[0]+energy_delta, self.native[self.mask]/self.scale

    def _diagnose(self,q,base,r,z,hoop,J,K):
        b,_,x,_,_,G,Gr,Gz,*_ = self._chart(r,z)
        require_map(Gr,G/r)
        cut = self.x[self.cut_node]
        regions = dict(all=np.ones(len(r),bool),terminal_collar=(z/self.L>.80)&(x>cut),
                       near_tip=(z/self.L>.95)&(x>.90))
        extrema = {k:dict(min_stretch=np.inf,max_stretch=0.,min_J=np.inf,max_J=0.) for k in regions}
        corners = [(p,dp) for p,dp in zip(*self.quad[:2])]+[(1.,0.),(-.5,0.)]
        for p,dp in corners:
            a,g,v = 1+Gr*p,1+G/r*p,Gz*p
            C = np.zeros((len(r),4,4))
            C[:,:2,:2] = K
            C[:,2,2] = C[:,3,3] = hoop**2
            A = np.zeros_like(C)
            A[:,0,0],A[:,0,1],A[:,0,2] = a,v,G/r*dp
            A[:,1,1],A[:,2,2],A[:,3,3] = 1.,g,g
            ev = np.linalg.eigvalsh(A@C@A.transpose(0,2,1))
            jac = J*a*g*g
            for name,mask in regions.items():
                ids = np.flatnonzero(mask)
                i = ids[np.argmin(ev[mask,0])]
                out = extrema[name]
                new_min = np.sqrt(ev[i,0])
                if new_min<out['min_stretch']:
                    out.update(min_stretch=float(new_min),node_site=[float(x[i]),float(z[i]/self.L),float(p)])
                out.update(max_stretch=max(out['max_stretch'],float(np.sqrt(ev[mask,-1].max()))),
                           min_J=min(out['min_J'],float(jac[mask].min())),
                           max_J=max(out['max_J'],float(jac[mask].max())))
        energy,_ = self.energy(q)
        return {**base,**extrema['all'],'energy':energy,'regions':extrema}


def equilibrate_scalar(body,q):
    def deriv(s):
        body.internal_amplitude = s
        body.energy(q)
        return body.angular_gradient
    center = body.internal_amplitude
    if abs(deriv(center)) < 1.e-10:
        return center
    width = 1.e-3
    for _ in range(12):
        lo,hi = center-width,center+width
        if deriv(lo) < 0 < deriv(hi):
            result = brentq(deriv,lo,hi,xtol=8.e-13)
            body.internal_amplitude = result
            return result
        width *= 2
    raise ArithmeticError('No bracketing minimum in the fixed-wall angular coordinate')


def solve_fixed(state_path,x,d,field,folder=None,*,angular_free=True):
    state_path = Path(state_path)
    base = Front(state_path,(13,10,14),1)
    n,k = (d.shape[1]-1)//3,d.shape[0]-1
    trial = Trial(base,state_path,x,n,k,True)
    b = AngularBonded(state_path,n,k,shape=trial.shape,full_depth=True,epsilon=x[5])
    m = trial.model
    q = d[b.mask]*b.scale
    u = field.copy()
    m.material_load = b.get_material_load(q,m)
    initial = trial.energy(d,u)
    set_front(b,m,u)
    before = b.energy(q,True)
    emit(kind='initial_angular',energy=initial,material_residual=float(np.max(abs(b.native[b.mask]))),
         angular_slope=b.angular_gradient,crystal=before,grid=[3*n,k],receiver=trial.receiver)

    def save(name,d,u,info):
        if folder is None:
            return
        path = Path(folder)/name
        tmp = path.with_suffix('.partial.npz')
        np.savez_compressed(tmp,x=x,field=u,displacement=d,internal_amplitude=b.internal_amplitude,
                            shape_params=trial.shape.params,shape_extra=trial.shape.extra,
                            state_sha=file_sha(state_path),input_sha=getattr(solve_fixed,'input_sha',''),
                            readout=json.dumps(info))
        tmp.replace(path)

    # The bond can rearrange. A tangent to the small angular correction is
    # used only as a preconditioner; acceptance uses the full nonlinear law.
    select = equilibrate_scalar if angular_free else lambda body,vector: 0.
    for sweep in range(30):
        set_front(b,m,u)
        s = select(b,q)
        E,qg = b.energy(q)
        e0,g0 = Bonded.energy(b,q)
        tangent = qg-g0
        old = q.copy()
        def local(y):
            try:
                en,grad = Bonded.energy(b,y)
            except ValueError as exc:
                if str(exc) != 'Nonpositive trial map or outside crystalline order grid':
                    raise
                return 1.e16+y@y,2*y
            return en+tangent@(y-old),grad+tangent
        opt = minimize(local,old,method='L-BFGS-B',jac=True,
                       options=dict(ftol=3e-16,gtol=2e-10,maxiter=3500,maxls=60))
        search = opt.x-old
        for half in range(26):
            q = old+2.**(-half)*search
            try:
                en,grad = b.energy(q)
            except ValueError as exc:
                if str(exc) != 'Nonpositive trial map or outside crystalline order grid':
                    raise
                continue
            if en <= E+8*np.spacing(abs(E)):
                break
        else:
            q = old
        s = select(b,q)
        b.energy(q)
        native = float(abs(b.native[b.mask]).max())
        if sweep >= 2 and native >= 3e-5:
            q,info = polish(b,q,target=3e-5)
            s = select(b,q)
            emit(kind='force_polish',**info)

        m.material_load = b.get_material_load(q,m)
        f,fi = m.solve(initial=u)
        set_front(b,m,f)
        e,g = b.energy(q)
        native = float(abs(b.native[b.mask]).max())
        gs = b.angular_gradient
        change = float(np.max(abs(f[2]-u[2])))
        b.front = None
        er,_ = b.energy(q)
        exchange = float((e-er)-m.material_load@f[2])
        energy = trial.static+er+m.evaluate(f)[0]
        info = dict(sweep=sweep+1,internal_amplitude=s,shape_fixed=True,angular_free=angular_free,
                    angular_slope=gs,native_meridional_residual=native,
                    field_residual=fi['stationarity'],order_step=change,exchange_error=exchange,
                    weak_receiver_residual=trial.receiver['residual'],energy=energy,
                    energy_change=energy-initial,meridional_step=float(np.max(abs(q-old)/b.scale)))
        emit(kind='fixed_wall_sweep',**info)
        u = f
        if abs(exchange)>1.e-9:
            raise ArithmeticError('Material and front work differ')
        dnew = np.zeros_like(d);dnew[b.mask] = q/b.scale
        save('latest.npz',dnew,u,info)
        if native < 3e-5 and (not angular_free or abs(gs)<1e-8) and change<2e-5 and fi['stationarity']<2e-7:
            set_front(b,m,u)
            solid = b.energy(q,True)
            flow = m.measure(u)
            if solid['min_J']<=0 or flow['port_error']>=3e-5:
                raise ArithmeticError('Orientation or reservoir-current closure failed')
            # Signed bracketing controls at fixed meridional solution.
            es = b.energy(q)[0]
            arcs = []
            for delta in ((-2e-4,-1e-4,1e-4,2e-4) if angular_free else ()):
                b.internal_amplitude = s+delta
                arcs.append((delta,b.energy(q)[0]-es,b.angular_gradient))
            b.internal_amplitude = s
            if not all(dE>0 and off*slope>0 for off,dE,slope in arcs):
                raise ArithmeticError('Reported internal mode lacks a restoring slope')
            info.update(crystal=solid,flow=flow,signed_angular_steps=arcs,
                        constrained_angular_force=None if angular_free else -gs,
                        physical_material_change=float(np.max(abs(dnew-d))))
            save('equilibrium.npz',dnew,u,info)
            emit(kind='fixed_wall_equilibrium',**info)
            return dnew,u,s,info
    raise ArithmeticError(f'Unsettled fixed wall at material residual {native}, angular derivative {gs}')


def run(args):
    prep,row,x,d,u = accepted(args.state,args.snapshot)
    emit(kind='preparation',state=file_sha(args.state),snapshot=file_sha(args.snapshot),
         source=prep['joint'],coordinates=x.tolist(),material_cells=prep['material_cells'],
         initial_iteration=row['iteration'],field_grid=prep['field_grid'])
    solve_fixed.input_sha = file_sha(args.snapshot)
    solve_fixed(args.state,x,d,u,args.folder,angular_free=not args.clamp_angular)
    emit(kind='saved',sha256=file_sha(args.folder/'equilibrium.npz'))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('state','snapshot','folder'):
        parser.add_argument(name,type=Path)
    parser.add_argument('--clamp-angular',action='store_true')
    run(parser.parse_args())
