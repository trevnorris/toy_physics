"""Read certified shear drives, partition work, and test for a repeatable ringdown."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks

from mouth_order_wave_pulse import digest,emit


def load(folder):
    outcome=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not limits['verified']:
        raise ValueError('Failed guard or worker '+str(folder))
    prep,drive,done=[json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    path=folder/'trace.npz'
    if digest(path)!=done['sha256'] or done['max_balance_over_final_work']>1e-8:
        raise ValueError('Trace digest or energy ledger failed')
    with np.load(path) as data:
        trace=dict(zip(data['columns'],data['trace'].T,strict=True))
    return prep,drive,trace


def fit_mode(t,y):
    """Fit signed velocity after switch-off, predict a later interval without refitting."""
    train=(t>=1)&(t<7)
    test=(t>=7)&(t<13)
    y=y/max(abs(y))
    def coefficients(rate,freq):
        phase=freq*(t-1)
        mat=np.column_stack((np.cos(phase),np.sin(phase)))*np.exp(-rate*(t-1))[:,None]
        c=np.linalg.lstsq(mat[train],y[train],rcond=None)[0]
        return mat@c
    def errors(params):
        return (coefficients(*params)-y)[train]
    scan=np.linspace(.8,32.,125)
    initial=min(scan,key=lambda f:np.linalg.norm(errors((.12,f))))
    opt=least_squares(errors,x0=(.12,initial),bounds=((0,.5),(2.5,33.)),ftol=1e-11)
    pred=coefficients(*opt.x)
    def relative(m):
        return np.linalg.norm((pred-y)[m])/max(1e-30,np.linalg.norm(y[m]))
    return dict(frequency=opt.x[1],envelope_rate=opt.x[0],train_error=relative(train),
                withheld_error=relative(test),train_rms=np.sqrt(np.mean(y[train]**2)),
                test_rms=np.sqrt(np.mean(y[test]**2)))


def report(path):
    prep,drive,trace=load(path)
    radius=prep['radius']
    t=trace['time']-drive['duration']
    e=trace['work'][-1]
    active=(trace['ordered']+trace['upper']+trace['lower'])/e
    exported=(trace['feed']+trace['loss'])/e
    below=trace['tail']/e
    receiving=trace['receiver']/e
    emit(test='input',folder=str(path),radius=radius,intake=prep['intake'],torque_frequency=drive['frequency'],
         drive_duration=drive['duration'],work=e,minimum_transit=prep['minimum_transit_time'])
    for lo,hi in ((0,1),(1,3),(3,6),(6,10),(10,14),(14,18)):
        window=(t>=lo)&(t<hi)
        if t[-1]<hi-.04:
            continue
        emit(test='after_switch_off',window=[lo,hi],active_mean=active[window].mean(),
             feed_mean=exported[window].mean(),tail_mean=below[window].mean(),
             receiving_mean=receiving[window].mean(),max_receiving_to_date=receiving[t<hi].max(),
             frozen_depth_work_per_input=trace['axial_force'][window].mean()/e)
    for name in ('near_velocity','far_velocity'):
        result=fit_mode(t,trace[name])
        emit(test='detector_fit',port=name,**result)
    # Peaks of the collocated post-drive signal identify a frequency to test
    # with a different waveform. The subsequent lifetime test uses new data.
    sample=t[1]-t[0]
    window=(t>=0)&(t<12)
    for name in ('near_velocity','far_velocity'):
        y=trace[name][window]
        spectrum=abs(np.fft.rfft(y*np.hanning(len(y)),n=4096))
        omega=2*np.pi*np.fft.rfftfreq(4096,sample)
        ids=find_peaks(spectrum)[0]
        ids=ids[(omega[ids]>=1)&(omega[ids]<=32)]
        ids=sorted(ids,key=lambda i:spectrum[i],reverse=True)[:3]
        emit(test='spectral_features',port=name,peaks=omega[ids].tolist(),
             peak_heights=spectrum[ids].tolist(),absolute_peak_velocity=float(max(abs(y))))
    for key in ('near_velocity','far_velocity'):
        if not np.any(np.isfinite(trace[key])) or not np.any(trace[key]):
            raise ArithmeticError('Missing detector')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folders',nargs='+',type=Path)
    for folder in p.parse_args().folders:
        report(folder)
