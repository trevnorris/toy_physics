"""Measure layer velocities and traction during matched 0PN companion-flow runs.

Every observed displacement is relative to the affine expansion of its own
patch. Subtract matched isolated displacements before taking velocities, so
identical Brownian kicks cancel. All signed loads are bath forces on 12 sites.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_order_wave_pulse import digest, emit
from neutral_patch_companion_rim import CASES, run
from neutral_patch_conversion import LEVELS, SEEDS
from neutral_patch_sites import COUNT

ARCHIVE = Path('research/pde_ledger_v3/_measurements/neutral_patch_companion_rim.json')
LAYERS = ('brane_facing', 'wall_facing')
# Each row accumulates carrier displacement, the two material displacements,
# the bond impulse on the wall-facing sites and the outgoing grip impulse.
U, UPPER, LOWER, BOND, GRIP = range(5)


class LayerMeter:
    def __init__(self, histories):
        self.history_count = histories
        self.elapsed = 0.
        self.current = np.zeros((histories,len(CASES),5))
        self.windows = []
        self.times = []

    def __call__(self, dt, gate, carrier, motion, gradient, port):
        upper_motion = motion[:,COUNT:].mean(axis=1)
        lower_motion = motion[:,:COUNT].mean(axis=1)
        # Forces within either layer cancel after summing. The net upper energy
        # gradient equals the cross-layer force on the lower material.
        bond = gradient[:,COUNT:].sum(axis=1)
        increments = np.stack((dt*carrier,upper_motion,lower_motion,dt*bond,dt*port),axis=1)
        matched = increments.reshape(self.history_count,len(CASES),5)
        self.current += matched-matched[:,:1]
        self.elapsed += dt
        if gate is not None:
            if gate != len(self.windows)+1:
                raise ValueError('Missing conversion interval')
            self.windows.append(self.current.copy())
            self.times.append(self.elapsed)
            self.current.fill(0.)
            self.elapsed = 0.


def read_window(label, integrated, duration, exposure):
    # Shape (matched thermal history, field case, layer); retained until the
    # final per-history diagnostics to keep stochastic scatter visible.
    rate = integrated/duration
    flow = rate[:,:,U,None]
    velocity = rate[:,:,[UPPER,LOWER]]
    lag = flow-velocity
    mechanical = np.stack((rate[:,:,BOND],rate[:,:,GRIP]-rate[:,:,BOND]),axis=-1)
    viscous = COUNT*(flow*np.array([1.,exposure])-velocity)
    error = float(np.max(abs(mechanical-viscous)))
    mean_flow = flow.mean(axis=0)[:,0]
    nonzero = [1,2,3,4,5]
    rows = []
    for k in nonzero:
        rows.append(dict(case=CASES[k],carrier_speed=mean_flow[k],
                         layer_speed=velocity[:,k].mean(axis=0).tolist(),
                         lag_fraction=(lag[:,k]/flow[:,k]).mean(axis=0).tolist(),
                         drag=mechanical[:,k].mean(axis=0).tolist(),
                         seed_drag_range_width=np.ptp(mechanical[:,k],axis=0).tolist()))
    emit(test='layer_reading',interval=label,duration=duration,layer_order=LAYERS,
         max_layer_force_balance=error,cases=rows,
         seed_drag8=mechanical[:,1].tolist(),seed_lag_fraction8=(lag[:,1]/flow[:,1]).tolist())
    if error > 2e-9 or not np.all(np.isfinite(rate)):
        raise ArithmeticError('Measured velocity and mechanical force disagree')
    return rate.mean(axis=0)


def identify(stages):
    # Use only position/bond/grip meters. Fit constant coefficients using three
    # separated 8R windows. Validate against the other 8R windows and all 16R
    # windows with an independent source profile.
    for layer, pos in enumerate((UPPER,LOWER)):
        mechanical = stages[:,:,BOND] if layer==0 else stages[:,:,GRIP]-stages[:,:,BOND]
        x = stages[:,1,pos]/stages[:,1,U]
        y = mechanical[:,1]/(COUNT*stages[:,1,U])
        training = [0,2,4]
        design = np.column_stack((np.ones(len(x)),-x))
        fit,_,rank,_ = np.linalg.lstsq(design[training],y[training],rcond=None)
        if rank<2 or fit[1]<=0:
            raise ArithmeticError('Degenerate constant-mobility fit')
        def load(case):
            return COUNT*(fit[0]*stages[:,case,U]-fit[1]*stages[:,case,pos])
        holdout = np.array([1,3])
        e8 = np.max(abs(load(1)[holdout]-mechanical[holdout,1]))
        e16 = np.max(abs(load(3)-mechanical[:,3]))
        emit(test='layer_identification',layer=LAYERS[layer],training_8R_windows=training,
             speed_over_flow_span=[x.min(),x.max()],condition_number=np.linalg.cond(design[training]),
             inferred_drag_per_site=fit[1],inferred_exposure=fit[0]/fit[1],
             held_out_force_max_error_8R=e8,held_out_force_max_error_16R=e16)


def experiment(args):
    archive = json.loads(ARCHIVE.read_text())
    if digest(args.profile) != archive['output_profiles_sha256']['_scratch/s11c/companion_patch_fields_v1.npz']:
        raise ValueError('Changed companion velocity field')
    observer = LayerMeter(args.histories)
    emit(test='layer_protocol',seeds=SEEDS[:args.histories],law=args.law,exposure=args.exposure,
         frame='tangential with local affine expansion removed',layers=LAYERS,
         exposure_independently_observed=False,sites_per_layer=COUNT)
    run(args,observe_step=observer)
    if len(observer.windows)!=len(LEVELS)-1:
        raise ValueError('Incomplete layer histories')
    lower_exposure = float(args.exposure=='both')
    stages = [read_window([LEVELS[i],LEVELS[i+1]],window,time,lower_exposure)
              for i,(window,time) in enumerate(zip(observer.windows,observer.times,strict=True))]
    identify(np.array(stages))
    integrated = np.sum(observer.windows,axis=0)
    total_time = sum(observer.times)
    mean = read_window('transit',integrated,total_time,lower_exposure)
    # Independent prior-run impulse from the grip for every case and seed.
    suffix = 'both_01' if args.exposure=='both' else '01'
    reference = archive['runs'][f'companion_rim_{args.law}_{suffix}']['literal_stdout']
    final = next(a for a in map(json.loads,reference) if a['test']=='patch_gate' and a['index']==5)
    grip = np.array(final['grip_impulse'])[:args.histories]
    old_error = float(np.max(abs(integrated[:,:,GRIP]-(grip-grip[:,:1]))))
    f = np.stack([mean[:,BOND],mean[:,GRIP]-mean[:,BOND]],axis=-1)
    emit(test='layer_controls',prior_fine_grip_impulse_difference=old_error,
         twice_half_drive_over_full=(2*f[5]/f[1]).tolist(),
         reversed_sign_sum_over_full=((f[1]+f[2])/f[1]).tolist(),
         drag_16_over_8=(f[3]/f[1]).tolist())
    if args.refine==2 and args.pieces==160 and old_error>2e-10:
        raise ArithmeticError('Failed to reproduce prior grip impulse')


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile',type=Path,default=Path('_scratch/s11c/companion_patch_fields_v1.npz'))
    p.add_argument('--law',choices=('central','layered'),required=True)
    p.add_argument('--exposure',choices=('upper','both'),required=True)
    p.add_argument('--pieces',type=int,default=160)
    p.add_argument('--refine',type=int,choices=(1,2),default=2)
    p.set_defaults(histories=4,grip=1.,length_ratio=.01)
    experiment(p.parse_args())
