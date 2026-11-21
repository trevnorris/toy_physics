"""Score completed blind middle-push forecasts against archived direct solves."""
import json
from pathlib import Path

from mouth_order_wave_pulse import ROOT, digest, emit


def read(path):
    return json.loads(path.read_text())


def rows(folder):
    outcome=read(folder/'outcome.json')
    limits=read(folder/'limit-validation.json')
    if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Unverified run {folder}')
    return [json.loads(row) for row in (folder/'stdout').read_text().splitlines()]


def archived(record,run,kind):
    r=record['runs'][run]
    if r['guard_exit'] or r['worker_exit'] or not r['verified_limits']:
        raise ValueError(f'Unverified reference {run}')
    return {d['t']:d['force_per_peak'] for s in r['literal_stdout'] if (d:=json.loads(s))['test']==kind}


def main():
    prior=Path('research/pde_ledger_v3/_measurements/mouth_middle_tag_grip.json')
    record=read(prior)
    for mesh in (26,34):
        for law in ('local','travelling'):
            label=f'n{mesh}_{law}'
            suffix='02' if label=='n34_travelling' else '01'
            hold=ROOT/f'mouth_wall_hold_{label}_{suffix}'
            calibration=rows(hold)[-1]
            prediction=rows(ROOT/f'mouth_hold_forecast_{label}_01')
            checksum=digest(hold/'hold.npz')
            if checksum!=calibration['sha256'] or checksum!=prediction[-1]['input_sha256']:
                raise ValueError('Frozen hold artifact does not match both processes')
            direct=archived(record,f'middle_replacement_{label}_01','middle_force')
            pred={r['time']:r['band_per_peak'] for r in prediction if r['test']=='prediction'}
            if not set(direct)<=set(pred) or len(direct)!=4:
                raise ValueError('Missing forecasts')
            error=max(abs(pred[t]-f) for t,f in direct.items())
            first=calibration['initial_midpoint'][2][0]
            late=calibration['last_midpoint'][2][0]
            emit(test='crosscheck',mesh=mesh,law=law,pulse_times=sorted(direct),
                 max_absolute_error_per_peak=error,prior_sha256=digest(prior),
                 held_entrance_band_initial=first,held_entrance_band_final=late,
                 first_midpoint_time=calibration['dt']/2,
                 last_midpoint_time=calibration['dt']*(calibration['hold_samples']-.5),
                 ratio_final_to_initial=late/first)
            if error>1e-8:
                raise ArithmeticError('Forecast disagrees with prior direct coupled calculation')


if __name__=='__main__':
    main()
