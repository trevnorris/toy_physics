"""Stationary order under a prepared wave and a through-flow on held meridians."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import factorized

from mouth_exchange_carrier import BIAS
from mouth_live_waist_fields import CoupledLoads
from mouth_neck_order_check import order_grid
from mouth_sealed_shape_relax import DELTA, Meridian


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json'
WELL = 3/DELTA


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def state(radius, point):
    rows = json.loads(DATA.read_text())['training'][str(radius)]
    return {'broad': rows['states'][-2], 'mid': rows['extension_states'][1],
            'tight': rows['extension_states'][-1]}[point]


class Order:
    def __init__(self, shape, mesh):
        self.shape = shape
        self.w, self.r, self.ids, lump, stiff = order_grid(mesh, shape)
        self.lump, self.stiff = WELL*lump, WELL*stiff
        self.free = np.ones_like(self.r, dtype=bool)
        self.free[[0, -1], :] = False
        self.free[:, -1] = False
        self.flatfree = self.free.ravel()
        gap = np.minimum(self.r[:, -1:] - self.r, shape.length-self.w[:,None])
        self.seed = .5 + .5*np.tanh(gap/(2*DELTA))
        self.seed[0] = .5 + .5*np.tanh((shape.neck-self.r[0])/(2*DELTA))
        self.seed[:, -1] = self.seed[-1] = .5
        self.norm_solver = factorized((2*self.stiff[self.flatfree][:, self.flatfree]
                                      + diags(1.5*self.lump[self.flatfree])).tocsc())
        ts = np.linspace(0., 1., 4097)
        self.t_waist = ts[np.argmin(shape.profile(ts)[0])]

    def energy_gradient(self, chi):
        y = chi.ravel() - .5
        strain = self.stiff@y
        return y@strain + self.lump@((y*y-.25)**2), (
            2*strain + self.lump*(4*y**3-y)).reshape(chi.shape)

    def axis(self, chi):
        return np.interp(self.t_waist*self.shape.length, self.w, chi[:, 0])

    def condition(self, chi, load, action, bias, tolerance=1e-5):
        """Preconditioned L-BFGS. Count a stationary field only after checking its residual."""
        chi = chi.copy()
        chi[~self.free] = self.seed[~self.free]
        def fields(c):
            eo, go = self.energy_gradient(c)
            ef, gf, info = load.evaluate(c, action=action, bias=bias)
            return eo+ef, (go+gf)[self.free], eo, info
        value, g, eo, info = fields(chi)
        memory = []
        for iteration in range(130):
            native = self.norm_solver(g)
            err = np.max(abs(native))
            if err < tolerance:
                break
            q = g.copy()
            coefficients = []
            for dx, dy, scale in reversed(memory):
                a = scale*(dx@q)
                q -= a*dy
                coefficients.append(a)
            direction = self.norm_solver(q)
            for (dx, dy, scale), a in zip(memory, reversed(coefficients)):
                direction += dx*(a - scale*(dy@direction))
            direction *= -1
            if g@direction >= 0:
                memory.clear()
                direction = -native
            slope = g@direction
            step = min(1., .18/max(abs(direction)))
            while step > 1.e-7:
                trial = chi.copy()
                trial[self.free] += step*direction
                if np.min(trial) <= .03:
                    step *= .5
                    continue
                newval, newg, neweo, newinfo = fields(trial)
                if newval <= value + .0001*step*slope + 5.e-13:
                    break
                step *= .5
            else:
                raise ArithmeticError(f'Smooth ordered branch stalls with residual {err}, minimum {chi.min()} and energy {value}')
            dx, dy = step*direction, newg-g
            curvature = dx@dy
            if curvature > 1.e-13*np.linalg.norm(dx)*np.linalg.norm(dy) and curvature>0:
                memory.append((dx,dy,1/curvature))
                if len(memory)>10:
                    memory.pop(0)
            value, g, eo, info, chi = newval, newg, neweo, newinfo, trial
        else:
            raise ArithmeticError(f'Coupled order minimization failed with residual {err}')
        m = {k:float(v) for k,v in info.items() if np.isscalar(v) and isinstance(v,(int,float,np.number))}
        m.update(total_inside_and_drives=value, order_energy=eo, axis_order=self.axis(chi),
                 min_order=float(chi.min()), max_order=float(chi.max()),iterations=iteration,
                 order_residual=err, waist_radius=float(self.shape.radius(self.shape.length*self.t_waist)[0]),
                 axial_depth=self.shape.length)
        return chi, m


def load_for(shape, n, field_mesh, target, wave_seed=None):
    return CoupledLoads(shape,n,2*n,field_mesh=field_mesh,port=4.2,
                        target_frequency=target,wave_seed=wave_seed)


def run(a):
    s = state(a.radius,a.point)
    p = np.array(s['coordinates'])
    shape = Meridian(a.radius,p)
    order = Order(shape,a.mesh)
    loads = load_for(shape,a.mesh,a.field_mesh,s['wave_frequency'])
    emit(kind='preparation', radius=a.radius, point=a.point, coordinates=p.tolist(),
         mesh=a.mesh, field_mesh=a.field_mesh, J=1.,phase_bias=BIAS,order_inertia=DELTA**2,
         well_scale=WELL, shape_width_coordinate='4t(1-t) length units', axial_coordinate='log(depth)')
    chi, blank = order.condition(order.seed,loads,0.,0.)
    emit(kind='zero',**blank)
    chi, flow = order.condition(chi,loads,0.,BIAS)
    emit(kind='flow',**flow)
    for action in (.25,.6,1.):
        chi, center = order.condition(chi,loads,action,BIAS)
        emit(kind='coupled',action=action,**center)
    if a.save is not None:
        np.savez_compressed(a.save,order=chi,frequency=center['wave_frequency'],
                            vector=loads.wave.seed[1],pole=loads.wave.seed[0])
    if a.pilot:
        return
    for coord in ('waist','depth'):
        for h in (-.002,-.001,.001,.002):
            delta = p.copy()
            if coord=='waist':
                delta[2] += h/a.radius
            else:
                delta[0] += h
            other = Meridian(a.radius,delta)
            fields = load_for(other,a.mesh,a.field_mesh,s['wave_frequency'],loads.wave.seed)
            shifted = Order(other,a.mesh)
            _, reading = shifted.condition(chi,fields,1.,BIAS)
            emit(kind='signed',coordinate=coord,offset=h,**reading)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radius',type=int,choices=(1,2),required=True)
    parser.add_argument('--point',choices=('broad','mid','tight'),required=True)
    parser.add_argument('--mesh',type=int,required=True)
    parser.add_argument('--field-mesh',type=int,required=True)
    parser.add_argument('--save',type=Path)
    parser.add_argument('--pilot',action='store_true')
    run(parser.parse_args())
