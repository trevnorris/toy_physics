"""Release terminal radius, depth, waist and the previously balanced order amplitude."""
import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import expit, logit

from mouth_poly_shallow_release import Trial, emit
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import ROOT


@dataclass
class Reading:
    energy: float
    order_slope: float
    seed: tuple
    info: dict
    chi: np.ndarray
    shifts: np.ndarray


class Release:
    def __init__(self, radius, mesh):
        self.trial = Trial(radius,mesh,bonded_core=True)
        record = json.loads((ROOT/'research/pde_ledger_v3/_measurements/mouth_poly_shallow_release.json').read_text())
        key = f'r{radius}n{mesh}'
        def accepted(family):
            row = record['runs'][f'{family}_{key}']
            if (row['guard_and_child_exit'] or
                    hashlib.sha256(row['literal_stdout'].encode()).hexdigest()!=row['source_stdout_sha256']):
                raise ValueError('Prior measurement failed validation')
            return [json.loads(line) for line in row['literal_stdout'].splitlines()]
        line = accepted('bonded')
        balance = accepted('balance')
        if (balance[-1]['kind']!='balance' or line[-1]['kind']!='directional_minimum'
                or balance[0]['inputs']!=self.trial.fingerprints):
            raise ValueError('Different parent material or missing equilibrium')
        self.held = line[-1]['step']*np.array(line[1]['displacement'])
        self.prior = balance[-1]
        self.x0 = np.array([self.prior['coordinates'][0],1.,*self.prior['coordinates'][1:]])
        e0,g,order,material,_,_ = self.trial.evaluate(np.zeros(6),self.trial.chi)
        self.reference_energy = e0
        self.initial_seed = material.wave.seed
        self.active = order.free & (self.trial.chi>0.) & (self.trial.chi<1.)
        chi0 = self.trial.chi[self.active]
        self.logits = logit(chi0)
        diag = (2*order.stiff.diagonal()+1.5*order.lump).reshape(g.shape)
        self.rate = -g[self.active]/diag[self.active]
        self.rate /= np.max(abs(self.rate*chi0*(1-chi0)))
        self.prior_hessian = np.array(self.prior['projected_hessian'])

    def measure(self, x, seed, *, preparation_fraction=1.):
        if x[1]<=0.:
            raise ValueError('Melted or closed topology outside open branch')
        shifts = self.held*preparation_fraction
        shifts[0], shifts[2] = x[0], x[2]
        shifts[1] += np.log(x[1])
        shape = Meridian(self.trial.radius,self.trial.x+shifts[:4])
        if shape.wall()['minimum_radius']<=0.:
            raise ValueError('Meridian crossed axis')
        chi = self.trial.chi.copy()
        vals = expit(self.logits+x[3]*self.rate)
        chi[self.active] = vals
        e,g,_,material,info,_ = self.trial.evaluate(shifts,chi,seed)
        if (not np.isfinite(e) or info['current']<=0. or info['mode_overlap'] < .94 or
                info['pole_residual']>1e-7 or abs(info['current']-info['outflow'])>1e-6*info['current']):
            raise ArithmeticError(f'Unresolved flow or mode at {x}: {info}')
        info['exit_radius'] = shape.exit
        return Reading(e,g[self.active]@(vals*(1-vals)*self.rate),material.wave.seed,info,chi,shifts)

    def start(self):
        seed = self.initial_seed
        for f in (.25,.5,.75,1.):
            x = self.x0*f
            x[1] = 1.
            center = self.measure(x,seed,preparation_fraction=f)
            seed = center.seed
        if (abs(center.energy-self.reference_energy-self.prior['change'])>2e-6 or
                max(abs(center.info['energy'][k]-v) for k,v in self.prior['energy'].items())>2e-6):
            raise ArithmeticError('Could not reconstruct archived balance and wave branch')
        return center

    def derivative(self,x,center,step):
        gradient = []
        for axis in (0,1,2):
            dx = np.zeros(4)
            dx[axis] = min(step,.05*x[1]) if axis==1 else step
            plus = self.measure(x+dx,center.seed)
            minus = self.measure(x-dx,center.seed)
            gradient.append((plus.energy-minus.energy)/(2*dx[axis]))
        return np.array([*gradient,center.order_slope])


def save(folder, name, problem, x, reading):
    file = folder/(name+'.npz')
    temp = folder/(name+'.tmp.npz')
    np.savez_compressed(temp,shifts=reading.shifts, chi=reading.chi,
                        pole=reading.seed[0],vector=reading.seed[1],radius=problem.trial.radius,
                        mesh=problem.trial.n,coordinates=x,inputs=problem.trial.fingerprints)
    temp.replace(file)
    return hashlib.sha256(file.read_bytes()).hexdigest()


def stability(problem,x,center):
    """Signed energy polarizations in four independent coordinates at two step sizes."""
    results = []
    for h in (.002,.001):
        basis = np.eye(4)*h
        basis[1,1] = min(h,.05*x[1])
        H = np.zeros((4,4))
        signed = []
        pair_energy = []
        for axis,dx in enumerate(basis):
            p = problem.measure(x+dx,center.seed).energy-center.energy
            m = problem.measure(x-dx,center.seed).energy-center.energy
            pair_energy.append(p+m)
            signed.append((p-m)/(2*dx[axis]))
            H[axis,axis] = (p+m)/dx[axis]**2
        for a in range(4):
            for b in range(a):
                mixed = basis[a]+basis[b]
                p = problem.measure(x+mixed,center.seed).energy-center.energy
                m = problem.measure(x-mixed,center.seed).energy-center.energy
                H[a,b] = H[b,a] = (p+m-pair_energy[a]-pair_energy[b])/(2*basis[a,a]*basis[b,b])
        vals = np.linalg.eigvalsh(H)
        emit(kind='four_coordinate_curvature',step=h,slopes=signed,hessian=H.tolist(),eigenvalues=vals.tolist())
        results.append((np.array(signed),H,vals))
    coarse, fine = results
    sensitivity = np.linalg.norm(fine[1]-coarse[1],2)
    return max(abs(fine[0]))<.002 and fine[2][0]>sensitivity


def run(radius,mesh,folder):
    folder.mkdir(exist_ok=True, parents=True)
    problem = Release(radius,mesh)
    x = problem.x0.copy()
    center = problem.start()
    energy0 = center.energy
    grad = problem.derivative(x,center,.0004)
    prior_slope = problem.prior['independent_log_exit_slopes'][-1]
    if abs(grad[1]-prior_slope)>.0015*max(1.,abs(prior_slope)):
        raise ArithmeticError('Exit slope at handoff disagrees with signed archived work')
    emit(kind='origin',radius=radius,mesh=mesh,inputs=problem.trial.fingerprints,
         axes=['log_depth','exit_over_start','middle','order_amplitude'],
         previous_log_exit_slope=prior_slope,initial_slopes=grad.tolist(),**center.info)
    indices = [0,2,3]
    metric = np.eye(4)/max(1.,2*abs(grad[1]))
    metric[np.ix_(indices,indices)] = np.linalg.inv((problem.prior_hessian+problem.prior_hessian.T)/2)
    iteration = 0
    status = 'stationarity_unresolved'
    while True:
        save(folder,'latest',problem,x,center)
        emit(kind='release',iteration=iteration,coordinates=x.tolist(),gradient=grad.tolist(),
             energy_change=center.energy-energy0,**center.info)
        if max(abs(grad))<.0008:
            status = 'balanced' if stability(problem,x,center) else 'curvature_or_force_unresolved'
            break
        step = -metric@grad
        if not np.isfinite(step).all() or step@grad>=0.:
            raise ArithmeticError('Invalid descent metric')
        step /= max(1.,np.max(abs(step))/.18, -step[1]/(.35*x[1]))
        factor = 1.
        last_error = ''
        while True:
            dx = factor*step
            if max(abs(dx))<5e-7:
                emit(kind='unresolved_step',last_probe_error=last_error)
                break
            try:
                probe = problem.measure(x+dx,center.seed)
                if probe.energy <= center.energy+1e-4*(dx@grad):
                    new_grad = problem.derivative(x+dx,probe,.0004)
                    break
            except (ArithmeticError,ValueError,RuntimeError) as err:
                last_error = str(err)
            factor /= 2.
        if max(abs(dx))<5e-7:
            break
        delta_grad = new_grad-grad
        curvature = dx@delta_grad
        if curvature>1e-10*np.linalg.norm(dx)*np.linalg.norm(delta_grad):
            left = np.eye(4)-np.outer(dx,delta_grad)/curvature
            metric = left@metric@left.T+np.outer(dx,dx)/curvature
        x,center,grad = x+dx,probe,new_grad
        iteration += 1
    digest = save(folder,'endpoint',problem,x,center)
    emit(kind='endpoint',status=status,coordinates=x.tolist(),snapshot_sha256=digest,
         energy_change=center.energy-energy0,gradient=grad.tolist(),**center.info)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--radius',type=int,choices=(1,2),required=True)
    cli.add_argument('--mesh',type=int,choices=(48,72),required=True)
    cli.add_argument('--folder',type=Path,required=True)
    a=cli.parse_args()
    run(a.radius,a.mesh,a.folder)
