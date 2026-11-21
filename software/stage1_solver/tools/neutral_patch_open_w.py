"""Remove plane constraints and measure relaxed in-plane shear at zero external w load."""
import argparse

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

from neutral_patch_sites import BETWEEN, GAP, WITHIN, emit, morse
from neutral_patch_structural_endpoints import Cell, preparation


def mix(z, power):
    t = np.minimum(abs(z)/GAP, 1.)
    f = t**3*(10-15*t+6*t*t)
    df = 30*t*t*(1-t)**2*np.sign(z)/GAP
    ddf = 60*t*(1-t)*(1-2*t)/GAP**2
    if power == 2:
        return f*f, 2*f*df, 2*(df*df+f*ddf)
    return f, df, ddf


class OpenCell(Cell):
    def __init__(self, side, cutoff, power):
        super().__init__(side, cutoff)
        self.power = power

    def evaluate(self, pos, shear=0., tangent=False):
        d = pos[self.j]-pos[self.i]
        d[:,:2] += self.shifts-self.box*np.rint(d[:,:2]/self.box)
        d[:,0] += shear*d[:,1]
        rr = np.einsum('ij,ij->i',d,d)
        keep = rr < self.cutoff**2
        i,j,d,r = self.i[keep],self.j[keep],d[keep],np.sqrt(rr[keep])
        f,df,ddf = mix(d[:,2],self.power)
        ea,da = morse(r,*WITHIN)
        eb,db = morse(r,*BETWEEN)
        diff = eb-ea
        g = (da+f*(db-da))/r
        pair_grad = g[:,None]*d
        pair_grad[:,2] += df*diff
        def incidence(values):
            return np.bincount(j,values,minlength=self.n)-np.bincount(i,values,minlength=self.n)
        grad = np.stack([incidence(v) for v in pair_grad.T],axis=1)
        grad[:,1] += shear*grad[:,0]
        stress = pair_grad[:,0]@d[:,1]/self.area
        if not tangent:
            return np.sum(ea+f*diff),grad,stress
        if shear != 0:
            raise ValueError('Tangent uses the unsheared reference')
        def second(params):
            depth,alpha,length = params
            z = np.exp(alpha*(length-r))
            return 2*depth*alpha**2*z*(2*z-1)
        s = ((1-f)*second(WITHIN)+f*second(BETWEEN)-g)/r**2
        normal = df*(db-da)/r
        hess = np.zeros((3*self.n,3*self.n))
        b = np.zeros_like(pos)
        for a in range(3):
            for k in range(3):
                h = s*d[:,a]*d[:,k] + (a==k)*g
                if a==2:
                    h += normal*d[:,k]
                if k==2:
                    h += normal*d[:,a]
                if a==k==2:
                    h += ddf*diff
                for u,v,sign in ((i,i,1),(j,j,1),(i,j,-1),(j,i,-1)):
                    np.add.at(hess,(3*u+a,3*v+k),sign*h)
                if k==0:
                    b[:,a] = incidence(h*d[:,1]+(a==1)*pair_grad[:,0])
                    if a==0:
                        born = h@(d[:,1]**2)
        return hess,b.ravel(),born

    def relax(self, initial, shear=0.):
        def unpack(flat):
            return np.vstack((initial[0],flat.reshape(-1,3)))
        def value(flat):
            energy,grad,_ = self.evaluate(unpack(flat),shear)
            return energy,grad[1:].ravel()
        sol = minimize(value,initial[1:].ravel(),jac=True,method='L-BFGS-B',
                       options=dict(ftol=2e-16,gtol=1e-10,maxiter=5000,maxls=60,maxcor=15))
        pos = unpack(sol.x)
        if shear==0:
            for _ in range(4):
                e,grad,_ = self.evaluate(pos)
                if np.max(abs(grad)) < 1e-9:
                    break
                h = self.evaluate(pos,tangent=True)[0][3:,3:]
                delta = np.linalg.solve(h,grad[1:].ravel()).reshape(-1,3)
                if np.max(abs(delta)) > .005:
                    break
                pos[1:] -= delta
        e,grad,stress = self.evaluate(pos,shear)
        err = np.max(abs(grad))
        if err > (1e-7 if shear==0 else 3e-5):
            raise ArithmeticError(f'Unresolved forces {err} at strain {shear}, {sol.message}')
        return pos,e,stress,err


def modulus(h,b,born,area,translations):
    k,modes = eigh(h)
    if np.max(abs(k[:translations]))>1e-7 or k[translations]<1e-7:
        raise ArithmeticError(f'Unresolved internal mode: {k[:translations+2]}')
    np.testing.assert_allclose(modes[:,:translations].T@b,0,atol=2e-7)
    values = k[translations:]
    projection = modes[:,translations:].T@b
    correction = -modes[:,translations:]@(projection/values)
    return (born-np.sum(projection**2/values))/area,correction,values[0]


def run(side,seed,power,kick_seed):
    flat = Cell(side,6.)
    xy = flat.settle(preparation(flat,seed))[0]
    pos0 = np.column_stack((xy,np.repeat([-.5,.5],flat.n//2)*GAP))
    cell = OpenCell(side,6.,power)
    e0,g0,s0 = cell.evaluate(pos0)
    old_e,old_g,old_s = flat.read(xy)
    np.testing.assert_allclose([e0,s0],[old_e,old_s],rtol=0,atol=1e-9)
    np.testing.assert_allclose(g0[:,:2],old_g,rtol=0,atol=1e-9)
    planar = np.tile([True,True,False],cell.n)
    h0,b0,born0 = cell.evaluate(pos0,tangent=True)
    hc,bc,ac = flat.read(xy,tangent=True)
    np.testing.assert_allclose(h0[planar][:,planar],hc,atol=1e-8,rtol=1e-10)
    np.testing.assert_allclose(b0[planar],bc,atol=1e-8,rtol=1e-10)
    np.testing.assert_allclose(born0,ac,atol=1e-8)
    confined = modulus(hc,bc,ac,cell.area,2)[0]
    emit(test='confined',side=side,seed=seed,power=power,kick_seed=kick_seed,
         modulus=confined,upper_holder_normal_force_per_area=g0[cell.n//2:,2].sum()/cell.area,
         peak_site_normal_force=np.max(abs(g0[:,2])))
    rng = np.random.default_rng(kick_seed)
    noise = rng.normal(size=pos0.shape)
    pos = pos0+.012*(noise-noise.mean(axis=0))
    for escapes in range(8):
        pos,energy,stress,residual = cell.relax(pos)
        h,b,born = cell.evaluate(pos,tangent=True)
        k,modes = eigh(h,subset_by_index=(0,4))
        if k[0] >= -1e-7:
            break
        mode = modes[:,0].reshape(-1,3)
        mode /= abs(mode).max()
        candidates = [pos+sign*.003*mode for sign in (-1,1)]
        pos = min(candidates,key=lambda p:cell.evaluate(p)[0])
        if cell.evaluate(pos)[0] >= energy:
            raise ArithmeticError('No downhill step along the negative mode')
    relaxed,dx,least = modulus(h,b,born,cell.area,3)
    fixed_w = modulus(h[planar][:,planar],b[planar],born,cell.area,2)[0]
    dx = dx.reshape(-1,3)
    dx -= dx[0].copy()
    checks=[]
    for strain in (5e-5,2.5e-5):
        results = [cell.relax(pos+sign*strain*dx,sign*strain) for sign in (-1,1)]
        _,em,sm,_ = results[0]
        _,ep,sp,_ = results[1]
        checks.append(dict(strain=strain, stress_slope=(sp-sm)/(2*strain),
                           energy_curvature=(ep+em-2*energy)/(cell.area*strain**2)))
    emit(test='released',side=side,seed=seed,power=power,kick_seed=kick_seed,escapes=escapes,
         energy_change_per_site=(energy-e0)/cell.n,affine=born/cell.area,relaxed=relaxed,
         frozen_w_at_new_state=fixed_w,residual=residual,least_mode=least,baseline_shear=stress,
         w_span=np.ptp(pos[:,2]),w_rms=np.std(pos[:,2]),
         sorted_half_mean_separation=np.diff([np.mean(t) for t in np.split(np.sort(pos[:,2]),2)])[0],
         sorted_half_w_std=[np.std(t) for t in np.split(np.sort(pos[:,2]),2)],checks=checks)
    if energy > e0+1e-8:
        raise ArithmeticError('Released state has higher energy than constrained baseline')
    for c in checks:
        np.testing.assert_allclose([c['stress_slope'],c['energy_curvature']],relaxed,
                                   rtol=.001,atol=.004)
    extended = OpenCell(side,7.,power)
    e,g,s = extended.evaluate(pos)
    np.testing.assert_allclose([e,s],[energy,stress],rtol=0,atol=2e-8)
    np.testing.assert_allclose(g,0,rtol=0,atol=2e-7)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('side',type=int,choices=[4,6])
    cli.add_argument('seed',type=int)
    cli.add_argument('power',type=int,choices=[1,2])
    cli.add_argument('--kick-seed',type=int,default=319)
    args=cli.parse_args()
    run(args.side,args.seed,args.power,args.kick_seed)
