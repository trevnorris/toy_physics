"""Local fixed-action balance on the archived bonded exits, measured in physical exit radius."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_exit_witness import inputs, validate_fields
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_trapping_angular_tail import recorded_seed
from mouth_trapping_modes import Sector, pole


def archived(radius,n):
    trial=Trial(radius,n,bonded_core=True)
    key=f'r{radius}n{n}'
    sh,chi,old,digest,status,*_=inputs(ROOT/f'_scratch/s11c/mouth_exit_release01/{key}/endpoint.npz',
                                        radius,n,trial)
    record=ROOT/'research/pde_ledger_v3/_measurements/mouth_wave_trapping.json'
    entry=json.loads(record.read_bytes())['supplementary']['tail_'+key]
    raw=entry['literal_stdout']
    if (status!='balanced' or entry['guard_exit'] or entry['child_exit'] or
            hashlib.sha256(raw.encode()).hexdigest()!=entry['stdout_sha256']):
        raise ValueError('Unqualified input preparation')
    rows=list(map(json.loads,raw.splitlines()))
    if rows[0]['endpoint_sha256']!=digest:
        raise ValueError('Order field does not match the modal scan')
    target,=[d for d in rows if d['kind']=='angular_mode' and d['degree']==12]
    return trial,sh,chi,old,digest,target


def measure(radius,n):
    trial,sh,chi,old,digest,target=archived(radius,n)
    op=Sector(trial,sh,chi,12)
    seed,info=recorded_seed(op,target)
    def readings(h):
        shifts=sh.copy()
        shifts[1]+=h
        obj=op if h==0. else Sector(trial,shifts,chi,12)
        if not np.array_equal(obj.live,op.live):
            raise ValueError('A shape move changed wave support')
        freq,stat=(seed,info) if h==0. else pole(obj,seed)
        if stat['seed_overlap']<.98 or stat['residual']>1e-7:
            raise ArithmeticError('Ambiguous continued pole')
        en,_,_,_,data,_=trial.evaluate(shifts,chi,old)
        validate_fields(data)
        terms=data['energy']
        nonwave=sum(v for k,v in terms.items() if k!='strained_wave')
        return dict(exit=obj.shape.exit,nonwave=nonwave,order4=terms['strained_wave'],
                    order12=freq[0].real,rate=-2*freq[0].imag,terms=terms,current=data['current'])
    center=readings(0.)
    side={h:readings(h) for h in (-.0025,-.00125,-.000625,.000625,.00125,.0025)}
    emit(kind='preparation',radius=radius,mesh=n,digest=digest,**center)
    for h in (.0025,.00125,.000625):
        low,high=side[-h],side[h]
        width=high['exit']-low['exit']
        slopes={k:(high[k]-low[k])/width for k in ('nonwave','order4','order12')}
        # Replace the previously prepared degree-four wave, retaining the
        # background density, inlet, material registration and order profile.
        action=-slopes['nonwave']/slopes['order12']
        if action<=0:
            raise ArithmeticError('No positive balancing action')
        a=center['exit']-low['exit']
        b=high['exit']-center['exit']
        E=lambda row: row['nonwave']+action*row['order12']
        curvature=2*((E(high)-E(center))/b-(E(center)-E(low))/a)/(a+b)
        emit(kind='exit_balance',log_step=h,physical_slopes=slopes,action=action,
             physical_closing_load=slopes['nonwave'],inward_slope_with_original_wave=slopes['nonwave']+slopes['order4'],
             cycle_energy=action*center['order12'],exit_curvature=curvature,
             threshold_work_minus_plus=[E(low)-E(center),E(high)-E(center)],
             at_half_action=slopes['nonwave']+.5*action*slopes['order12'],
             at_double_action=slopes['nonwave']+2*action*slopes['order12'])


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--radius',required=True,type=int,choices=(1,2))
    cli.add_argument('--mesh',required=True,type=int,choices=(48,72))
    a=cli.parse_args()
    measure(a.radius,a.mesh)
