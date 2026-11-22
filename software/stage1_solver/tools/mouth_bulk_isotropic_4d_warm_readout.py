"""Compare signed 4D holds across strain amplitudes and matched Brownian time steps."""
import argparse
import json
from pathlib import Path

import numpy as np

from neutral_patch_sites import emit


def rows(folder):
    outcome = json.loads((folder/'outcome.json').read_text())
    if not outcome['limitsVerified'] or outcome['exitCode'] or outcome['childOutcome']['exitCode']:
        raise ValueError(f'Incomplete or unverified run {folder}')
    return list(map(json.loads,(folder/'stdout').read_text().splitlines()))


def slope_summary(s):
    return dict(initial=s[0], last=s[-1], final_quarter=float(np.mean(s[-12:])),
                retained=float(np.mean(s[-12:])/s[0]), final_quarter_range=[min(s[-12:]),max(s[-12:])])


def crystal(root):
    for index,name in enumerate(('bulk','shear')):
        data = rows(root/f'iso4_crystal{index}_01')
        setup, = [r for r in data if r['test']=='warm_setup']
        hold = [r for r in data if r['test']=='warm_hold']
        np.testing.assert_allclose([r['fraction'] for r in hold],np.arange(49)/48)
        signed = np.array([r['signed_stress'][1:] for r in hold])
        slopes = (signed[:,0]-signed[:,1])/(setup['amplitudes'][1]-setup['amplitudes'][2])
        np.testing.assert_allclose(slopes,[r['slope'] for r in hold],rtol=0.,atol=1e-13)
        emit(test='original_crystal',drive=name,dt=setup['dt'],n=setup['n'],T=setup['temperature'],
             max_original_lost=max(r['lost_original'] for r in hold),
             max_new_neighbors=max(r['new_neighbors'] for r in hold),
             terminal_rms=hold[-1]['rms_displacement'],**slope_summary(slopes))


def rms_frames(a,b):
    d = a-b
    d -= d.mean(axis=-2,keepdims=True)
    return np.sqrt(np.mean(np.sum(d*d,axis=-1),axis=-1))


def comparison(batch):
    data, zero = {},{}
    for refine,phases in [(1,['bulk','shear','half-shear']),(2,['shear','half-shear'])]:
        folder = batch/f'undriven_{refine}'
        undriven = rows(folder)
        with np.load(folder/'native/hold.npz') as saved:
            zero[refine] = saved['pos'][:,0]
            box = saved['boxes'][0]
        hold = [r['diagnostics'][0] for r in undriven if r['test']=='hold']
        emit(test='undriven',resolution=refine,starting_neighbors=hold[0]['initial_bonds'],
             most_original_lost=max(r['lost'] for r in hold),most_added=max(r['added'] for r in hold),
             terminal_rms=hold[-1]['rms_since_start'])
        for phase in phases:
            folder = batch/f'{phase}_{refine}'
            traj = rows(folder)
            setup, = [r for r in traj if r['test']=='strain']
            holds = [r for r in traj if r['test']=='hold']
            np.testing.assert_array_equal([r['frame'] for r in holds],np.arange(49))
            signed = np.array([r['signed'] for r in holds])
            slope = (signed[:,0]-signed[:,1])/(setup['amplitudes'][0]-setup['amplitudes'][1])
            np.testing.assert_array_equal(slope,[r['slope'] for r in holds])
            data[refine,phase] = slope
            with np.load(folder/'native/hold.npz') as saved:
                unstrained = saved['pos']*(box/saved['boxes'])[:,None,:]
            separation = rms_frames(unstrained,zero[refine][:,None])
            emit(test='signed_hold',resolution=refine,phase=phase,**slope_summary(slope),
                 max_signed_to_zero_nonaffine_rms=np.max(separation,axis=0).tolist(),
                 max_signed_to_signed_nonaffine_rms=max(rms_frames(unstrained[:,0],unstrained[:,1])))
        full,half = data[refine,'shear'],data[refine,'half-shear']
        emit(test='amplitude_comparison',resolution=refine,max_abs_full_minus_half=max(abs(full-half)),
             full_minus_half_final_quarter=float(np.mean((full-half)[-12:])))
    distances = rms_frames(zero[1],zero[2])
    for phase in ['shear','half-shear']:
        coarse,fine = data[1,phase],data[2,phase]
        emit(test='time_step_comparison',phase=phase,
             coarse_final_quarter=float(np.mean(coarse[-12:])),fine_final_quarter=float(np.mean(fine[-12:])),
             fine_minus_coarse_last_quarter=float(np.mean((fine-coarse)[-12:])),
             max_abs_pointwise=max(abs(coarse-fine)),unstrained_rms_at_probe_start=distances[0],
             unstrained_rms_at_probe_end=distances[-1],unstrained_rms_max=max(distances))


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('case',choices=('crystal','disordered'))
    cli.add_argument('root',type=Path)
    args = cli.parse_args()
    if args.case=='crystal':
        crystal(args.root)
    else:
        comparison(args.root)
