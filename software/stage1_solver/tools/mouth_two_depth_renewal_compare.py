"""Compare the independently guarded traction, carrier-subtraction and stability data."""
import json

import numpy as np

from mouth_bend_renewal import ROOT
from mouth_flared_receiver import emit
from mouth_two_depth_renewal import GAUGES, PORTS


def rows(name):
    folder=ROOT/name
    read=lambda filename: json.loads((folder/filename).read_text())
    guard=read('limit-validation.json')
    out=read('outcome.json')
    if (out['exitCode'] or out['childOutcome']['exitCode'] or not out['limitsVerified']
            or not guard['verified'] or guard['actual']['memory.swap.max']!='0'):
        raise ValueError(f'Unverified calculation {name}')
    return [json.loads(s) for s in (folder/'stdout').read_text().splitlines()]


def record(records,kind):
    item,=[r for r in records if r['test']==kind]
    return item


def pair(item,key):
    return np.array(item[key+'_real'])+1j*np.array(item[key+'_imag'])


def main():
    coarse=record(rows('two_depth_renewal_n26_local_01'),'driven')
    carrier=pair(coarse,'bare_carrier_matrix')
    for mesh in (26,34):
        for law in ('local','travelling'):
            main_rows=rows(f'two_depth_renewal_n{mesh}_{law}_01')
            direct=record(main_rows,'driven')
            if law=='local':
                check=rows(f'two_depth_local_stress_n{mesh}_01')
                alt=record(check,'local_stress_resolvent')
                pulse_type='local_stress_pulse'
                growth=record(check,'local_positive_mode')['rate']
            else:
                check=rows(f'two_depth_material_check_n{mesh}_{law}_01')
                alt=record(check,'material_resolvent')
                pulse_type='pulse_budget'
                growing,=[r for r in check if r['test']=='free_mode' and r['rate'][0]>.01]
                growth=growing['rate'][0]
            holding_error=np.max(abs(pair(alt,'bend_holding')-pair(direct,'bend_holding_response')))
            actual=pair(direct,'station_force')
            inferred=np.linalg.solve(GAUGES.T,carrier@pair(direct,'particular_q')-PORTS)
            pulse={r['t']:r for r in check if r['test']==pulse_type}
            prim={r['t']:r for r in main_rows if r['test']=='pulse_motion'}
            path_error=max(np.max(abs(np.array(r['gauge_per_peak'])-prim[t]['per_peak'])) for t,r in pulse.items())
            emit(test='comparison',mesh=mesh,law=law,holding_matrix_error=holding_error,
                 time_path_error_per_peak=float(path_error),
                 coarse_carrier_station_force_error_per_peak=abs(inferred-actual).tolist(),
                 maximum_replay_error=record(main_rows,'pulse_replay')['max_midpoint_force_error_per_peak'],
                 work_errors=[max(r['relative_work_error']) for r in pulse.values()],
                 clock=record(main_rows,'material')['clock'],end=max(pulse),
                 growth=growth,efold=1/growth)
            if holding_error>1e-7 or path_error>1e-4:
                raise ArithmeticError('Independent material solve disagrees with the free paths')


if __name__=='__main__':
    main()
