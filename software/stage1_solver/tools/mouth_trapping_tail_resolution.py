"""Refine the wave mesh alone at the saved degree-twelve pole of each n=72 material."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn, spherical_yn

from mouth_exit_witness import inputs
from mouth_poly_shallow_release import Trial, emit
from mouth_shallow_puncture_fields import PORT
from mouth_trapping_angular_tail import integer_readout, recorded_seed
from mouth_trapping_mesh_check import project
from mouth_trapping_modes import Sector, pole


def wronskian_flux(op, seed):
    o, v = seed
    ph = op.port.project@v[op.port.ids]
    propagates = op.port.cutoff < o.real
    k = np.sqrt(o.real**2-op.port.cutoff[propagates]**2)
    z = PORT*k
    conductance = 1/(k*PORT**2*(spherical_jn(op.degree,z)**2+spherical_yn(op.degree,z)**2))
    power = .5*o.real*PORT**4*abs(ph[propagates])**2*conductance
    store = .25*(abs(o)**2*np.vdot(v,op.M@v).real+np.vdot(v,op.K@v).real
                  -2*PORT**3*np.vdot(ph,ph).real)
    return power/store


def run(radius):
    trial = Trial(radius,72,bonded_core=True)
    key = f'r{radius}n72'
    base = Path('_scratch/s11c/mouth_exit_release01')/key/'endpoint.npz'
    sh,c,_,digest,status,*_ = inputs(base,radius,72,trial)
    folder = Path('_scratch/s11c/mouth_trapping_tail01')/key
    guard = json.loads((folder/'outcome.json').read_bytes())
    stdout = (folder/'stdout').read_bytes()
    rows = list(map(json.loads,stdout.splitlines()))
    if (status!='balanced' or trial.mesh!=20 or rows[0]['endpoint_sha256']!=digest or
            guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']):
        raise ValueError('Invalid degree-twelve control')
    target, = [r for r in rows if r['kind']=='angular_mode' and r['degree']==12]
    previous = Sector(trial,sh,c,12)
    seed, _ = recorded_seed(previous,target)
    emit(kind='native',radius=radius,material_mesh=72,wave_mesh=20,endpoint_sha256=digest,
         tail_stdout_sha256=hashlib.sha256(stdout).hexdigest(),frequency=seed[0].real,
         loss=-2*seed[0].imag,Q=target['Q'])
    for mesh in (28,36,48):
        refined = copy.copy(trial)
        refined.mesh = mesh
        op = Sector(refined,sh,c,12)
        guess = (seed[0],project(previous,op,seed[1]))
        try:
            nxt, data = pole(op,guess)
            stats = integer_readout(op,nxt,data)
        except (ValueError, ArithmeticError, RuntimeError) as exc:
            emit(kind='unresolved_pole',radius=radius,wave_mesh=mesh,reason=str(exc))
            return
        backward = project(op,previous,nxt[1])
        reverse_overlap = abs(np.vdot(seed[1],previous.M@backward))**2/np.vdot(backward,previous.M@backward).real
        independent = wronskian_flux(op,nxt)
        difference = max(abs(independent-stats['physical_channel_power']))/sum(independent)
        emit(kind='wave_refinement',radius=radius,material_mesh=72,wave_mesh=mesh,
             active_wave_nodes=len(nxt[1]),reverse_projection_overlap=reverse_overlap,
             wronskian_power_per_energy=independent.tolist(),relative_channel_discrepancy=difference,**stats)
        if min(data['seed_overlap'],reverse_overlap)<.90 or data['detached'] or difference>1e-8:
            emit(kind='unresolved_correspondence_or_flux',radius=radius,wave_mesh=mesh)
            return
        previous,seed=op,nxt


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--radius',type=int,choices=(1,2),required=True)
    run(cli.parse_args().radius)
