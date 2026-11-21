"""Compare the finite carrier runs to an independent fixed-bias virtual-work prediction."""
import argparse
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from mouth_empty_bend_stiffness import run as tangent
from mouth_flared_receiver import emit

ROOT=Path('_scratch/s11c')


def results(folder):
    status=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if (status['exitCode'] or status['childOutcome']['exitCode'] or not status['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Unverified source {folder}')
    return [json.loads(line) for line in (folder/'stdout').read_text().splitlines()]


def compare(a):
    def finite(label):
        records=results(ROOT/('bend_reaction_'+label))
        return {r['offset'][1]:r for r in records if r['test']=='native'}
    qrun,phase=finite(a.intake),finite(a.phase)
    # E*_delta=Q. Independently differentiate the bend Hessian at neutral
    # shape, before loading or using any of the finite-displacement currents.
    slopes=[]
    for step in (.005,.0025):
        states=[]
        for sign in (-1,1):
            args=SimpleNamespace(mesh=a.mesh,far=9.,roof=10.,profile='delayed',rate=.25+sign*step,
                                 source=None,extension='compact',save=None)
            buf=io.StringIO()
            with redirect_stdout(buf):
                tangent(args)
            prep,row=[json.loads(line) for line in buf.getvalue().splitlines()]
            states.append([prep['bias'],row['holding_stiffness'][1][1]])
            emit(test='straight_tangent',phase_bias=prep['bias'],intake=prep['bulk'],
                 deep_stiffness=states[-1][1],trace_value=row['surface_only'][1][1])
        (dm,km),(dp,kp)=states
        slopes.append((kp-km)/(dp-dm))
    dQdphase=2*step/(dp-dm)
    for size in (.01,.005):
        plus,minus=phase[size],phase[-size]
        deltaQ=(plus['bulk']+minus['bulk'])/2-phase[0.]['bulk']
        bias=(qrun[size]['phase_bias']+qrun[-size]['phase_bias'])/2-qrun[0.]['phase_bias']
        predicted=.5*slopes[-1]*size**2
        observed=qrun[size]['shell']
        # Raw interior current error includes the straight geometric quadrature error.
        # Relative changes in the control volume at fixed intake are evaluated
        # against that same straight volume; exterior boundary reactions define Q.
        center=qrun[0.]['shell']
        pm=lambda key: (np.array(observed[key])-np.array(qrun[-size]['shell'][key]))/2
        inc,out,side=map(pm,['low','high','side_out'])
        closure=out[2]-inc[2]+side[2]
        emit(test='conservative_control',mesh=a.mesh,offset=size,fixed_bias_intake=plus['bulk'],
             even_intake_change=deltaQ,virtual_work_prediction=predicted,
             ratio_to_prediction=deltaQ/predicted,fixed_rate_bias_change=bias,
             required_static_bias_prediction=-predicted/dQdphase,
             carrier_holder_odd=((np.array(qrun[size]['holding_force'])-np.array(qrun[-size]['holding_force']))/2).tolist(),
             inward_lower_x=inc[2],outward_upper_x=out[2],outward_wall_x=side[2],x_imbalance=closure,
             closure_fraction=abs(closure)/max(abs(inc[2]),abs(out[2]),abs(side[2])),
             lower_advected_x=inc[1],upper_advected_x=out[1],
             raw_number_closure_straight=center['closure'][0],
             number_closure_difference=observed['closure'][0]-center['closure'][0])
        if abs(deltaQ/predicted-1)>.01 or abs(bias*dQdphase/(-predicted)-1)>.01:
            raise ArithmeticError('Energy derivative failed the phase/flux comparison')
    emit(test='rate_halving',virtual_work_slopes=slopes,
         current_ratio=(phase[.01]['bulk']-phase[0.]['bulk'])/(phase[.005]['bulk']-phase[0.]['bulk']))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',required=True,type=int)
    p.add_argument('--phase',required=True)
    p.add_argument('--intake',required=True)
    compare(p.parse_args())
