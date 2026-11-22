"""Resolve the fast modes of a saved fixed-cubic 4D Morse sample under a warm hold."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_bulk_isotropic_4d import Cell, library, strain_basis, structure
from neutral_patch_sites import TEMPERATURE, emit


FRAMES = 48
TRANSIT = 48.88171243122866
RESOLUTION = .6  # maximum cold decay rate times the coarse Euler step
MODES = {'undriven': None, 'bulk': 0, 'shear': 1, 'half-shear': 1}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source(path):
    parent = path.parent.parent
    outcome = json.loads((parent/'outcome.json').read_text())
    if not outcome['limitsVerified'] or outcome['exitCode'] or outcome['childOutcome']['exitCode']:
        raise ValueError(f'Unverified input at {path}')
    return np.load(path, allow_pickle=False)


def displacement(a, b):
    d = a-b
    d -= d.mean(axis=0)
    return float(np.sqrt((d*d).sum(axis=1).mean()))


class Brownian:
    def __init__(self, n, highest, refine, seed):
        self.n = n
        self.refine = refine
        self.count = int(np.ceil(TRANSIT*highest/(FRAMES*RESOLUTION)))
        self.coarse_dt = TRANSIT/(FRAMES*self.count)
        self.dt = self.coarse_dt/refine
        self.rng = np.random.default_rng(1701+seed)
        self.bridge = np.random.default_rng(9071+seed)
        if highest*self.coarse_dt > RESOLUTION:
            raise ArithmeticError('Under-resolved cold mode')

    def frame(self, states):
        scale = np.sqrt(2*TEMPERATURE*self.coarse_dt)
        for i in range(0, self.count, 128):
            k = min(128, self.count-i)
            coarse = self.rng.normal(size=(k, self.n, 4))*scale
            if self.refine == 1:
                noise = coarse
            else:
                bridge = self.bridge.normal(size=coarse.shape)*scale
                noise = np.stack(((coarse+bridge)/2, (coarse-bridge)/2), axis=1).reshape(-1,self.n,4)
            for cell, pos in states:
                if cell.lib.integrate(cell.handle, len(noise), self.dt, noise, pos):
                    raise ArithmeticError('Invalid Euler increment')


class Neighbors:
    def __init__(self, pos, box):
        self.initial = self.bonds = structure(pos,box)[0]
        self.start = pos.copy()
        self.box = box

    def read(self, pos):
        fresh, r = structure(pos, self.box)
        self.bonds = np.triu((self.bonds & (r < 1.30)) | fresh, 1)
        return dict(initial_bonds=int(self.initial.sum()),
                    lost=int((self.initial & ~self.bonds).sum()),
                    added=int((self.bonds & ~self.initial).sum()),
                    rms_since_start=displacement(pos,self.start))


def quench(cell, pos, original):
    minimum = cell.settle(pos.copy())
    old_bonds = structure(original, cell.box)[0]
    new_bonds = structure(minimum, cell.box)[0]
    return dict(energy_per_site=cell.read(minimum)[0]/cell.n,
                rms_to_original_minimum=displacement(minimum,original),
                bonds_lost_from_original=int((old_bonds & ~new_bonds).sum()),
                bonds_added_from_original=int((new_bonds & ~old_bonds).sum()),
                maximum_force=float(abs(cell.read(minimum)[1]).max()))


def run(args):
    lib = library(args.work)
    with source(args.cold) as cold:
        original, box = cold['pos'], cold['box']
        cutoff, highest, seed = float(cold['rc']), float(cold['highest']), int(cold['seed'])
    n = len(original)
    brownian = Brownian(n,highest,args.refine,seed)
    cell = Cell(lib,box,n,cutoff)
    emit(test='setup', phase=args.phase, resolution=args.refine, seed=seed, cold_sha256=digest(args.cold),
         highest=highest, dt=brownian.dt, highest_times_dt=highest*brownian.dt, transit=TRANSIT,
         steps_per_frame=brownian.count*args.refine, T=TEMPERATURE, n=n, density=n/cell.volume, cutoff=cutoff,
         noise_seed=1701+seed, bridge_seed=9071+seed, sampling='48 equal intervals per transit')
    if args.phase == 'prepare':
        pos = original.copy()
        neighbors = Neighbors(pos,box)
        for k in range(FRAMES+1):
            if k:
                brownian.frame([(cell,pos)])
            emit(test='warming', frame=k, energy_per_site=cell.read(pos)[0]/n, **neighbors.read(pos))
        emit(test='cooled_after_warming', **quench(cell,pos,original))
        np.savez(args.work/'prepared.npz', pos=pos, refine=args.refine, cold_sha=digest(args.cold),
                 rng=json.dumps(brownian.rng.bit_generator.state), bridge=json.dumps(brownian.bridge.bit_generator.state))
        return
    with source(args.prepared) as start:
        if int(start['refine']) != args.refine or str(start['cold_sha']) != digest(args.cold):
            raise ValueError('Temperature preparation does not match the cold source and step')
        pos = start['pos']
        brownian.rng.bit_generator.state = json.loads(str(start['rng']))
        brownian.bridge.bit_generator.state = json.loads(str(start['bridge']))
    mode = MODES[args.phase]
    a = np.diag(strain_basis()[mode if mode is not None else 1])
    amp = .00005 if args.phase=='half-shear' else .0001
    amplitudes = [0.] if mode is None else [amp, -amp]
    states, watches = [], []
    for value in amplitudes:
        stretch = np.exp(value*a)
        mat = Cell(lib, box*stretch,n,cutoff)
        p = pos*stretch
        states.append((mat,p))
        watches.append(Neighbors(p,mat.box))
    emit(test='strain', normalized_diagonal=a.tolist(), amplitudes=amplitudes, prepared_sha256=digest(args.prepared))
    slopes = []
    stressed = []
    for k in range(FRAMES+1):
        if k:
            brownian.frame(states)
        readings, diagnostics = [], []
        for (mat,p), watch in zip(states,watches):
            energy, _, diag = mat.read(p)
            readings.append(float((diag@a - n*TEMPERATURE*sum(a))/mat.volume))
            diagnostics.append(dict(energy_per_site=energy/n, **watch.read(p)))
        difference = None if mode is None else (readings[0]-readings[1])/(2*amp)
        if mode is not None:
            slopes.append(difference)
        stressed.append([p.copy() for _,p in states])
        emit(test='hold', frame=k, signed=readings, slope=difference, diagnostics=diagnostics)
    if slopes:
        emit(test='late', phase=args.phase, resolution=args.refine, initial=slopes[0],
             final_quarter=float(np.mean(slopes[-12:])), fraction=float(np.mean(slopes[-12:])/slopes[0]),
             last=slopes[-1], final_quarter_range=[min(slopes[-12:]),max(slopes[-12:])])
    if mode is None:
        emit(test='cooled_after_hold', **quench(cell,states[0][1],original))
    np.savez(args.work/'hold.npz', pos=np.array(stressed), boxes=np.array([m.box for m,_ in states]))


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('phase', choices=('prepare',*MODES))
    cli.add_argument('--cold', type=Path, required=True)
    cli.add_argument('--prepared', type=Path)
    cli.add_argument('--refine', type=int, choices=(1,2), required=True)
    cli.add_argument('--work', type=Path, required=True)
    args = cli.parse_args()
    if args.phase!='prepare' and args.prepared is None:
        cli.error('A hold needs the thermalized starting state.')
    run(args)
