"""Independent prediction windows and force/energy ledgers for the six-cycle burst."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from mouth_order_wave_pulse import digest,emit


def read(folder):
    outcome=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not limits['verified']:
        raise ValueError('Guarded outcome failed')
    setup,done=[json.loads(s) for s in (folder/'stdout').read_text().splitlines()]
    p=folder/'trace.npz'
    if digest(p)!=done['sha256'] or done['max_balance_over_final_work']>1e-8:
        raise ValueError('Unverified time trace')
    with np.load(p) as f:
        columns=dict(zip(f['columns'],f['trace'].T,strict=True))
    return setup,columns


def fit(t,obs,start):
    # Two independent port amplitudes/phases, one angular frequency and rate.
    # Equal weighting in their training RMS. Report errors at each port.
    tr=(t>=start)&(t<start+4)
    guard=(t>=start+4)&(t<start+8)
    rms=np.sqrt(np.mean(obs[tr]**2,axis=0))
    y=obs/rms
    def predict(rate,frequency):
        z=t-start
        v=np.exp(-rate*z)[:,None]*np.column_stack((np.cos(frequency*z),np.sin(frequency*z)))
        amps=np.linalg.lstsq(v[tr],y[tr],rcond=None)[0]
        return v@amps
    def residual(q):
        return (predict(*q)-y)[tr].ravel()
    grid=np.linspace(1.,25.,97)
    best=min(grid,key=lambda w:np.linalg.norm(residual((.25,w))))
    ans=least_squares(residual,[.25,best],bounds=([0,.5],[4,30.]))
    pred=predict(*ans.x)
    def err(mask): return (np.linalg.norm((pred-y)[mask],axis=0)/np.linalg.norm(y[mask],axis=0)).tolist()
    return dict(train_start=start,holdout_start=start+4,omega=ans.x[1],amplitude_decay=ans.x[0],
                train_relative_errors=err(tr),holdout_relative_errors=err(guard),
                train_rms=rms.tolist(),holdout_rms=np.sqrt(np.mean(obs[guard]**2,axis=0)).tolist())


def measurements(key,setup,d):
    t=d['time']-setup['drive_duration']
    work=d['work'][-1]
    active=(d['ordered']+d['upper']+d['lower'])/work
    feed=(d['feed']+d['loss'])/work
    tail=d['tail']/work
    receiver=d['receiver']/work
    p=d['transport_momentum']/work
    frozen=d['axial_force']/work
    obs=np.column_stack((d['near_velocity'],d['far_velocity']))
    emit(**dict(setup,test='setup_reading',tag=key,work=work))
    for lo,hi in ((.5,4.5),(4.5,8.5),(8.5,12.5),(14,18)):
        i,j=np.searchsorted(t,[lo,hi])
        if j>=len(t):
            raise ValueError('Incomplete after-drive window')
        dt=t[j]-t[i]
        def mean(x): return np.trapz(x[i:j+1],t[i:j+1])/dt
        geom=mean(frozen)
        transport=(p[j]-p[i])/dt
        emit(test='window',tag=key,start=t[i],stop=t[j],active=mean(active),tail=mean(tail),
            feed=mean(feed),receiver=mean(receiver),receiver_max_through_window=receiver[:j+1].max(),
            active_at_start=active[i],active_at_stop=active[j],frozen_gradient=geom,
            momentum_term=transport,generalized_guide_force_per_work=geom+transport)
    for start in (.5,4.5,8.5):
        emit(test='two_probe_prediction',tag=key,**fit(t,obs,start))
    return t,obs


def main(a):
    results={}
    for folder in a.folders:
        setup,d=read(folder)
        name=folder.name
        results[name]=measurements(name,setup,d)
    for coarse,fine in a.pair:
        ta,oa=results[coarse]
        tb,ob=results[fine]
        interp=np.column_stack([np.interp(tb,ta,col) for col in oa.T])
        for start in (.5,4.5,8.5,14):
            use=(tb>=start)&(tb<start+4)
            mismatch=np.linalg.norm((interp-ob)[use],axis=0)/np.linalg.norm(ob[use],axis=0)
            emit(test='waveform_refinement',coarse=coarse,fine=fine,start=start,
                 port_relative_differences=mismatch.tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folders',type=Path,nargs='+')
    p.add_argument('--pair',action='append',nargs=2,default=[])
    main(p.parse_args())
