"""Condense the mobile bonded material and joined field at each radius-four shape.

The six coordinates move depth, outlet, mean width, width asymmetry, waist and
the polar-equatorial wall. Physical material displacement and field coefficients
are held in the signed envelope derivative on a moving geometry chart.
"""
import argparse
from hashlib import sha256
from io import BytesIO
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from mouth_bonded_front_adjoint import set_front
from mouth_bonded_material_residual import polish
from mouth_joined_angular_work import bump
from mouth_joined_front_audit import Front
from mouth_polar_bonded import PolarBonded, StaticPolar
from mouth_polar_joined import PolarCoupled
from mouth_poly_shallow_release import CUT, Order
from mouth_receiver_order_control import side_source
from mouth_receiver_order_equilibrium import receiver
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_front import check_guard
from mouth_supercore_joined_front import emit, file_sha


class Trial:
    def __init__(self, base, path, x, n, k, precondition):
        self.x = x.copy()
        sh = Meridian(4., base.sh.params + np.array([1, .25, .25, .25])*x[:4],
                      extra=base.sh.extra + .25*x[4])
        self.shape = sh
        self.epsilon = x[5]
        order = base.order_mesh
        source = side_source(base, sh)
        profile, self.receiver = receiver(sh, order, CUT, lambda w, r: source,
                                          polar_shift=x[5], material_path=path)
        polar = StaticPolar(path, shape=sh, outer_order=profile.side)
        self.static = (Order(sh, order).energy_gradient(base.c)[0] +
                       self.receiver['energy']+self.receiver['static_work'] +
                       polar.inner_energy(x[5])+polar.tail_energy(x[5]))
        self.model = PolarCoupled(sh, base.c, 1, base.grid, polar_shift=x[5],
                                  terminal_profile=profile, precondition=precondition)
        self.body = PolarBonded(path, n, k, shape=sh, full_depth=True, epsilon=x[5])

    def energy(self, d, u):
        b, m = self.body, self.model
        b.front = None
        q = d[b.mask]*b.scale
        m.material_load = b.get_material_load(q, m)
        return self.static + b.energy(q)[0]+m.evaluate(u)[0]

    def settle(self, d, u, polish_target=3e-5):
        b, m = self.body, self.model
        q = d[b.mask]*b.scale
        for sweep in range(70):
            set_front(b, m, u)
            def objective(qtrial):
                try:
                    return b.energy(qtrial)
                except ValueError as exc:
                    # Inversion is an infeasible line-search step. Other errors
                    # (especially a mismatch of coupled bands) are genuine.
                    if str(exc) != 'Nonpositive trial map or outside crystalline order grid':
                        raise
                    return 1e16+qtrial@qtrial, 2*qtrial
            sol = minimize(objective, q, jac=True, method='L-BFGS-B',
                           options=dict(gtol=1e-9, ftol=1e-15, maxiter=3500, maxls=60))
            q = sol.x
            b.energy(q)
            if polish_target is not None and sweep and np.max(abs(b.native[b.mask])) >= 3e-5:
                q, polish_info = polish(b,q,target=polish_target)
                emit(kind='force_polish', coordinates=self.x.tolist(), sweep=sweep,
                     **polish_info)
            m.material_load = b.get_material_load(q, m)
            f, field_info = m.solve(initial=u)
            set_front(b, m, f)
            e, grad = b.energy(q)
            residual = abs(b.native[b.mask]).max()
            b.front = None
            e_ref, _ = b.energy(q)
            transfer = (e-e_ref)-m.material_load@f[2]
            if abs(transfer)>1e-9:
                raise ArithmeticError(f"Solid and interface work differ by {transfer}")
            update = np.max(abs(f[2]-u[2]))
            u = f
            info = dict(material_residual=residual, field_residual=field_info['stationarity'],
                        exchange_error=transfer, order_step=update, sweeps=sweep+1,
                        receiver_residual=self.receiver['residual'],
                        receiver_newton_steps=self.receiver['newton_steps'])
            if residual<3e-5 and field_info['stationarity']<2e-7 and update<2e-5:
                displacement = np.zeros_like(d)
                displacement[b.mask] = q/b.scale
                return displacement, u, self.static+e_ref+m.evaluate(u)[0], info
        raise ArithmeticError(f"Bonded solid and front did not settle: {info}")


class Release:
    def __init__(self, args):
        self.args = args
        for p in (args.state, args.joint):
            check_guard(p)
        self.base = Front(args.state, (13, 10, 14), 1)
        from json import loads
        history = (args.joint.parent/'stdout').read_text().splitlines()
        first, last = loads(history[0]), loads(history[-1])
        if (last['kind'] != 'held_equilibrium' or last['sha256'] != file_sha(args.joint)
                or first['rest'] or tuple(first['grid']) != self.base.grid):
            raise ValueError('Guarded material/field pair and held equilibrium must agree')
        if first['source_sha256'] != file_sha(args.state) or first['winding'] != 1:
            raise ValueError('Unmatched R4 state and historical joint field')
        with np.load(args.joint) as f:
            self.d = f['displacement']
            self.u = f['field']
        self.n, self.k = (self.d.shape[1]-1)//3, self.d.shape[0]-1
        self.x = np.zeros(6)
        transfer = None
        if args.from_accepted:
            from scipy.interpolate import RegularGridInterpolator
            path = args.from_accepted
            contents = path.read_bytes()
            digest = sha256(contents).hexdigest()
            limit = loads((path.parent/'limit-validation.json').read_text())
            source_rows = [loads(l) for l in (path.parent/'stdout').read_text().splitlines() if l.startswith('{')]
            if (not limit['verified'] or int(limit['actual']['memory.swap.max'])
                    or source_rows[0]['state'] != file_sha(args.state)
                    or source_rows[0]['field_grid'] != list(self.base.grid)
                    or not any(r.get('kind')=='accepted_shape' and r.get('sha256')==digest for r in source_rows)):
                raise ValueError('Unverified mesh-transfer seed')
            with np.load(BytesIO(contents)) as seed:
                self.x = seed['x'].copy()
                self.u = seed['field'].copy()
                old = seed['displacement']
            old_n, old_k = (old.shape[1]-1)//3, old.shape[0]-1
            coarse = PolarBonded(args.state,old_n,old_k,shape=self.base.sh,full_depth=True)
            fine = PolarBonded(args.state,self.n,self.k,shape=self.base.sh,full_depth=True)
            tt,xx = np.meshgrid(fine.t,fine.x,indexing='ij')
            self.d = RegularGridInterpolator((coarse.t,coarse.x),old)((tt,xx))
            transfer = dict(kind='mesh_transfer',state=digest,coordinates=self.x.tolist(),
                 source_material_mesh=[3*old_n,old_k],target_material_mesh=[3*self.n,self.k])
        self.folder = args.folder
        emit(kind='preparation', state=file_sha(args.state), joint=file_sha(args.joint),
             material_cells=[3*self.n,self.k], field_grid=self.base.grid, winding=1)
        if transfer is not None:
            emit(**transfer)

    def trial(self, x, precondition):
        return Trial(self.base, self.args.state, x, self.n, self.k, precondition)

    def work(self, x, d, u, e, h):
        values = []
        for axis in np.eye(6):
            values.append([self.trial(x+sign*h*axis, False).energy(d,u)
                           for sign in (-1.,1.)])
        pairs = np.asarray(values)
        grad = (pairs[:,1]-pairs[:,0])/(2*h)
        diag = (pairs.sum(axis=1)-2*e)/h**2
        return grad, diag

    def save(self, iteration, x, d, u, trial, energy, gradient, info):
        p = self.folder/'shape.npz'
        temp = self.folder/'shape.partial.npz'
        np.savez_compressed(temp, x=x, field=u, displacement=d, iteration=iteration,
                            shape_params=trial.shape.params, shape_extra=trial.shape.extra,
                            state_sha=file_sha(self.args.state), energy=energy)
        temp.replace(p)
        (self.folder / f'accepted_{iteration:03d}.npz').hardlink_to(p)
        t = np.linspace(0,1,1001)
        h, _ = bump(t*trial.shape.length,trial.shape.length)
        radius = trial.shape.radius(t*trial.shape.length)[0]
        pole = radius+x[5]*h
        equator = radius-.5*x[5]*h
        body = trial.body
        set_front(body, trial.model, u)
        read = body.energy(d[body.mask]*body.scale, True)
        emit(kind='accepted_shape', iteration=iteration, coordinates=x.tolist(),
             depth=trial.shape.length, exit=trial.shape.exit, pole_minimum=pole.min(),
             equator_minimum=equator.min(), widest_depth_pole=t[np.argmax(pole)]*trial.shape.length,
             widest_depth_equator=t[np.argmax(equator)]*trial.shape.length, energy=energy,
             global_forces=(-gradient).tolist(), internal=info, crystal=read,
             flow=trial.model.measure(u), sha256=file_sha(p))

    def run(self):
        x = self.x.copy()
        d, u = self.d, self.u
        h = .00010
        trial = self.trial(x,True)
        d,u,e,info = trial.settle(d,u,self.args.polish_target)
        g,diagonal = self.work(x,d,u,e,h)
        self.save(0,x,d,u,trial,e,g,info)
        metric = np.diag(1/np.maximum(abs(diagonal),1.))
        for it in range(1,251):
            if np.max(abs(g)) < 4e-4:
                finer, _ = self.work(x,d,u,e,h/2)
                emit(kind='signed_half_step', forces=(-finer).tolist(),
                     force_step_disagreement=float(abs(finer-g).max()))
                if max(abs(finer)) > 5e-4 or max(abs(finer-g)) > 1e-4:
                    raise ArithmeticError('Signed shape derivatives have not converged')
                held = trial.model.measure(u)
                if held['port_error'] >= 3e-5 or held['crystalline_range_in_band'][0] < 0 or held['crystalline_range_in_band'][1] > 1:
                    raise ArithmeticError('Endpoint fails current or mobile terminal order gate')
                emit(kind='projected_equilibrium', iteration=it-1, state=file_sha(self.folder/'shape.npz'))
                return
            step = -metric@g
            step *= min(1., .05/np.max(abs(step)))
            accepted = False
            for half in range(22):
                dx = step*2.**-half
                new_x = x + dx
                # Fixed-field screening avoids re-solving at an uphill trial.
                probe = self.trial(new_x,False)
                fixed = probe.energy(d,u)
                if fixed > e and half<20:
                    continue
                new = self.trial(new_x,True)
                nd,nu,ne,ni = new.settle(d,u,self.args.polish_target)
                if ne < e + 1e-4*(g@dx):
                    accepted = True
                    break
            if not accepted:
                raise ArithmeticError(f'No accepted shape step at forces {(-g).tolist()}')
            ng, _ = self.work(new_x,nd,nu,ne,h)
            y = ng-g
            curv = y@dx
            if curv>1e-12*np.linalg.norm(y)*np.linalg.norm(dx):
                B = np.eye(6)-np.outer(dx,y)/curv
                metric = B@metric@B.T + np.outer(dx,dx)/curv
            x,d,u,e,g,trial = new_x,nd,nu,ne,ng,new
            self.save(it,x,d,u,new,e,g,ni)
        raise ArithmeticError("Shape-force convergence unresolved")


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('state','joint','folder'):
        p.add_argument(key,type=Path)
    p.add_argument('--from-accepted',type=Path)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--lbfgs-only',action='store_const',const=None,dest='polish_target',
                      help='Reproduce the original coarse and fine descent.')
    mode.add_argument('--polish-target',type=float,default=3e-5,
                      help='Native force tolerance for directional-work polishing; the transfer used 1e-5.')
    p.set_defaults(polish_target=3e-5)
    Release(p.parse_args()).run()
