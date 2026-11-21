"""Compare the energy left by the prescribed wall stroke and the two outgoing ports."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_order_wave_pulse import digest, emit


def read(folder):
    geom, stroke = map(json.loads,(folder/'stdout').read_text().splitlines())
    receipt = folder/'outcome.json'
    if receipt.exists():
        out = json.loads(receipt.read_text())
        if out['exitCode']:
            raise ValueError(f'{folder}: guard returned {out["exitCode"]}')
    else:
        # A separately documented terminal service may retain the native worker
        # result even though its launcher vanished. Keep this provenance visible.
        receipt = folder/'terminal_worker_verification.json'
        out = json.loads(receipt.read_text())
        if out['state']['ActiveState']!='inactive' or out['state']['MainPID']!='0':
            raise ValueError(f'{folder}: service not verified terminal')
    child = json.loads((folder/'child-outcome.json').read_text())
    limits = json.loads((folder/'limit-validation.json').read_text())
    if (not limits['verified'] or limits['actual']['memory.swap.max']!='0'
            or out['childOutcome']!=child or child['exitCode'] or not out['limitsVerified']):
        raise ValueError(f'{folder}: missing successful guarded receipt')
    trace = folder/'trace.npz'
    if geom['test']!='geometry' or stroke['test']!='stroke' or digest(trace)!=stroke['trace_sha256']:
        raise ValueError(f'{folder}: mismatched wave trace')
    with np.load(trace) as saved:
        data = dict(zip(saved['columns'],saved['readings'].T,strict=True))
    return receipt.name,stroke,data


def main(folders):
    controls = {}
    for folder in folders:
        receipt,stroke,data = read(folder)
        work = stroke['signed_wall_work']
        supplied = stroke['positive_work']
        windows = {}
        for lo,hi in ((4,8),(20,24),(44,48),(194,198)):
            t = data['time']
            if t[-1]<hi-.1:
                continue
            use = (t>=lo)&(t<hi)
            def percent(energy):
                return None if work==0 else 100*energy[use].mean()/work
            portions = dict(ordered=percent(data['ordered_throat']),weak=percent(data['weak_throat']),
                feed_net=percent(data['feeding_field']+data['remote_feed']),cut=percent(data['receiving_cut']),
                reached_feed_absorber=percent(data['remote_feed']))
            windows[f'{lo}-{hi}']=portions
        emit(case=folder.name, receipt=receipt, dt=stroke['dt'], net_wall_work=work,
             positive_work=supplied, returned_work=-stroke['negative_work'],
             returned_percent_of_positive=None if supplied==0 else -100*stroke['negative_work']/supplied,
             max_relative_budget_error=None if work==0 else stroke['max_balance_error']/work,
             peak_shear_stress=stroke['peak_stress'],window_percent_of_net_work=windows)
        controls[folder.name] = stroke,data
    full,df = controls['formation_wall_ax1_01']
    half,dh = controls['formation_wall_half_01']
    # The half stroke reverses at every time point. Both use the same mesh and dt.
    n=len(dh['time'])
    t_error=np.max(abs(dh['time']-df['time'][:n]))
    v_error=np.max(abs(2*dh['collar_stress']+df['collar_stress'][:n]))
    emit(test='opposite_half_stroke',time_alignment=t_error,
         full_work=full['signed_wall_work'],four_times_half=4*half['signed_wall_work'],
         max_signed_stress_sum=v_error,full_peak=np.max(abs(df['collar_stress'][:n])))
    if t_error>1e-10 or v_error>1e-12:
        raise ValueError('opposite half-stroke traces do not scale')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path,nargs='+')
    main(parser.parse_args().folder)
