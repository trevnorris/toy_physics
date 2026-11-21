"""Energy captured by an opening radius-one scalar well, measured in its final modes."""
import argparse
import json
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from scipy.sparse.linalg import eigsh

from mouth_wave_flare_work import ARCHIVE, EXIT, Wave


BENCHMARK = Path(__file__).resolve().parents[3] / 'research/pde_ledger_v3/_measurements/mouth_occupied_companion.json'


def emit(kind, **fields):
    print(json.dumps(dict(test=kind, **fields), allow_nan=False), flush=True)


def opening(t, start):
    if start is None:
        return 1.
    x = np.clip((t-start)/3, 0., 1.)
    return x*x*(3-2*x)


def incident(wave):
    z = wave.w + 2.
    window = np.where(abs(z) < 1, np.cos(np.pi*z/2)**4, 0.)
    derivative = np.where(abs(z) < 1, -2*np.pi*np.sin(np.pi*z/2)*np.cos(np.pi*z/2)**3, 0.)
    k = 2.
    omega = np.hypot(wave.p['mu_out'], wave.p['c_w']*k)
    speed = wave.p['c_w']**2*k/omega
    radial = wave.r*np.exp(-wave.r**2/2)
    q = radial*window*np.cos(k*z)
    v = radial*(-speed*derivative*np.cos(k*z) + omega*window*np.sin(k*z))
    return q.ravel(), v.ravel()


def threshold():
    run = json.loads(BENCHMARK.read_text())['guarded_runs']['occupied_companion_n34_01']
    return next(row['critical_action_at_zero_tension'] for line in run['literal_stdout'].splitlines()
                if (row := json.loads(line))['test'] == 'static')


class Experiment:
    def __init__(self, args):
        self.dt = args.dt
        p = yaml.safe_load(ARCHIVE.read_text())['parameters']
        wave = Wave(1., 1., args.step, args.pad, p)
        self.uniform = wave.free.tocsr()
        self.uniform.setdiag(self.uniform.diagonal() + p['mu_out']**2)
        self.depth = wave.contrast*wave.gate(1.)[0].ravel()
        final = self.uniform-sparse.diags(self.depth, format='csr')
        # Request modes from the bottom up until the next level clears the open exterior gap.
        count = 3
        while True:
            values, modes = eigsh(final, k=count, sigma=0., which='LM', v0=wave.seed, tol=1e-10)
            if values[-1] > p['mu_out']**2:
                break
            count *= 2
        keep = values < p['mu_out']**2
        self.omega = np.sqrt(values[keep])
        self.modes = np.ascontiguousarray(modes[:, keep].T)
        self.target = threshold()
        residuals = np.linalg.norm(final@modes-modes*values, axis=0)
        # Remote sponge starts six units beyond either axial end or the widest radius.
        # Report actual energy lost there. It is switched on at t=8, after both ramps.
        exterior = np.maximum(np.maximum(-wave.w-6, wave.w-EXIT-6), wave.r-1.5-6)
        gamma = 1.2*np.clip(exterior/3, 0, 1)**2
        self.sponge = gamma.ravel() > 0
        self.atten = np.exp(-gamma.ravel()[self.sponge]*self.dt/2)
        self.loss_coeff = -np.expm1(-gamma.ravel()[self.sponge]*self.dt)
        self.q, self.v = incident(wave)
        normalizer = np.sqrt(self.energy(self.q, self.v, 0.))
        self.q /= normalizer
        self.v /= normalizer
        emit('settings', step=args.step, padding=args.pad, dt=args.dt, grid=[wave.r.size,wave.w.size],
             frequencies=np.sqrt(values).tolist(), gap=p['mu_out'], eigen_residuals=residuals.tolist(),
             bound_modes=len(self.omega), bound_sponge_norms=np.sum(self.modes[:,self.sponge]**2, axis=1).tolist(),
             initial_unit_energy=self.energy(self.q, self.v, 0.), earlier_ground_action_threshold=self.target)
        if residuals.max() > 1e-7 or not len(self.omega):
            raise ValueError('unresolved final spectrum')

    def energy(self, q, v, s):
        return float((v@v + q@(self.uniform@q-s*self.depth*q))/2)

    def mode_energy(self, q, v):
        return ((self.modes@v)**2+(self.omega*(self.modes@q))**2)/2

    def damp(self, v):
        local = v[self.sponge]
        v[self.sponge] = self.atten*local
        return .5*np.dot(self.loss_coeff, local*local)

    def run(self, name, start, amplitude, pure=False):
        q, v = self.q*amplitude, self.v*amplitude
        if pure:
            q, v = np.sqrt(2)*self.modes[0]/self.omega[0], np.zeros_like(v)
        s = opening(0., start)
        e0 = self.energy(q,v,s)
        acceleration = -(self.uniform@q)+s*self.depth*q
        work = loss = 0.
        snapshots = []
        def record(t):
            total = self.energy(q,v,s)
            modes = self.mode_energy(q,v)
            row = dict(t=t, energy=total, work=work, absorbed=loss,
                       balance_error=total+loss-e0-work, final_mode_energies=modes.tolist())
            snapshots.append(row)
        record(0.)
        dt = self.dt
        end = 48 if amplitude==1 else 8
        times = {round(t/dt): t for t in (3,6,8,16,32,48)}
        for n in range(round(end/dt)):
            t = (n+1)*dt
            new_s = opening(t, start)
            if n >= round(8/dt):
                loss += self.damp(v)
            v += dt/2*acceleration
            q_next = q + dt*v
            acceleration = -(self.uniform@q_next)+new_s*self.depth*q_next
            v += dt/2*acceleration
            if new_s != s:
                work -= (new_s-s)/4 * np.dot(self.depth, q*q+q_next*q_next)
            q, s = q_next, new_s
            if n >= round(8/dt):
                loss += self.damp(v)
            if n+1 in times:
                record(times[n+1])
        capture = next(row for row in snapshots if row['t']==(6 if start==3 else 3))
        ground = capture['final_mode_energies'][0]
        required = self.target*self.omega[0]*e0/ground if ground>0 and e0>0 and not pure else None
        emit('trajectory', name=name, start=start, amplitude=amplitude, snapshots=snapshots,
             ground_action_at_ramp_end=ground/self.omega[0],
             incident_energy_for_conditional_ground_threshold=required)
        if max(abs(row['balance_error']) for row in snapshots)>0.008*max(e0,1.e-15):
            raise ValueError('field work and remote absorption fail to account for energy')
        return ground


def main(args):
    if args.step<=0 or args.dt<=0 or args.pad<9.5 or args.dt>.4*args.step:
        raise ValueError('use a positive, resolved step and a padded domain')
    if not np.isclose(1/args.dt, round(1/args.dt)):
        raise ValueError('requested readouts must fall on full time steps')
    test = Experiment(args)
    early = test.run('opens_early', 0, 1.)
    late = test.run('opens_during_passage', 3, 1.)
    static = test.run('already_open', None, 1.)
    half = test.run('half_packet', 3, .5)
    quiet = test.run('quiet', 3, 0.)
    test.run('static_ground', None, 1., pure=True)
    emit('controls', early=early, passage=late, static=static,
         four_times_half=4*half, quiet=quiet)
    np.testing.assert_allclose(4*half, late, rtol=1e-11)
    if quiet != 0.:
        raise ValueError('empty homogeneous state changed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('step','pad','dt'):
        parser.add_argument('--'+key, type=float, required=True)
    main(parser.parse_args())
