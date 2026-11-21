"""Read signed thickness/occupancy experiments and check the fixed-flux carrier limit."""
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from neutral_patch_sites import emit

ROOT = Path('_scratch/s11c')
WINDOWS = (1.,2.,4.,None)


def rows(folder):
    def read(name):
        return json.loads((folder/name).read_text())
    out = read('outcome.json')
    child = read('child-outcome.json')
    if not (out['exitCode']==child['exitCode']==0 and child==out['childOutcome']
            and child['guardReason'] is None and read('limit-validation.json')['verified']
            and read('duration-validation.json')['verified']):
        raise ValueError(str(folder))
    return list(map(json.loads,(folder/'stdout').read_text().splitlines()))


def window_overlap(sigma, period):
    def peaks(s):
        return sum(np.exp(-.5*(i*period/s)**2) for i in range(-6,7))
    if sigma is None:
        return 1.
    combined=np.hypot(sigma,1.)
    return sigma/combined*peaks(combined)/(peaks(1.)*peaks(sigma))


def patch(power):
    case = '01' if power==1 else '02'
    results = {}
    for center in ('main','shift'):
        data=rows(ROOT/f'neutral_width_stripe_{power}_{center}{case}')
        baseline=data[0]
        assert baseline['test']=='preparation'
        period=baseline['periods'][1]
        loads={(r['size'],r['drive'],r['load']):r for r in data if r['test']=='load'}
        assert len(loads)==28
        state={}
        for s in WINDOWS:
            for drive in (('pressure',) if s is None else ('pressure','affinity')):
                full,half=[],[]
                for h,vals in ((.001,full),(.0005,half)):
                    pair=[loads[s,drive,h*sign] for sign in (1,-1)]
                    for key in ('width','number','density'):
                        vals.append((pair[0]['center_log_changes'][key]-
                                     pair[1]['center_log_changes'][key])/(2*h))
                np.testing.assert_allclose(full,half,rtol=.0015,atol=2e-7)
                state[s,drive]=np.array(full)
        results[center]=state
        for s in WINDOWS:
            p=state[s,'pressure']
            overlap=window_overlap(s,period)
            emit(test='size_series', power=power, center=center, gaussian_width=s,
                 center_dlog_width_number_density_per_outward_pressure=p.tolist(),
                 ratio_to_uniform=p[0]/state[None,'pressure'][0],
                 local_window_null=overlap,
                 direction_reversed_half_step='checked',
                 site_bias_width_fraction=state[s,'affinity'][0]/state[s,'affinity'][1] if s else None)
    for key in results['main']:
        np.testing.assert_allclose(results['main'][key],results['shift'][key],rtol=.0015,atol=2e-7)


def active_annulus():
    # Unit far density and chemical potential of the historical calculation,
    # fixed intake Q=.25. This is the historical slowly varying exterior law.
    mp.mp.dps=42
    def state(r,width):
        current=mp.mpf('.25')/(4*mp.pi*r**2)
        def equation(n):
            return (n**4-1)/4+(current/(width*n))**2/2
        return mp.findroot(equation,1)
    for r in (8,16,32):
        n=state(r,1)
        v=mp.mpf('.25')/(4*mp.pi*r*r*n)
        slope=v*v/(n**4-v*v)
        h=mp.mpf('1e-6')
        fd=mp.log(state(r,mp.exp(h))/state(r,mp.exp(-h)))/(2*h)
        emit(test='annulus', radius=r, reference_density_deficit=float(1-n),
             dlog_density_per_log_width=float(slope),central_difference=float(fd),
             density_change_at_one_percent_narrower=float(state(r,mp.mpf('.99'))-n))
        if slope<=0 or abs(fd/slope-1)>mp.mpf('1e-11'):
            raise ArithmeticError('Carrier width derivative disagreement')


if __name__=='__main__':
    for power in (1,2):
        patch(power)
    active_annulus()
