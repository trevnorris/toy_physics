"""Energy remaining after moving the six-cycle torque from the entrance into the throat."""
import argparse
from pathlib import Path

import numpy as np

from mouth_order_wave_pulse import emit
from mouth_shear_ringdown_compare import fit,read


def analyse(folder):
    s,d=read(folder)
    t=d['time']-s['drive_duration']
    w=d['work'][-1]
    active=d['ordered']+d['upper']+d['lower']
    tail=d['tail']
    feed=d['loss']+d['feed']
    cut=d.get('cut_export',np.zeros_like(t))
    receiver=d['receiver']
    ledger=(active+tail+feed+cut+receiver)/w
    err=np.max(abs(ledger[t>0]-1))
    emit(test='cavity_input',tag=folder.name,drive_depth=s.get('source_depth',.30),
         current=s['current'],residual_inertia=s['inertia_floor'],work=w,
         cut=s.get('cut'),max_error=err)
    if err>1e-8:
        raise ArithmeticError('Post-drive energy channels failed')
    for start in (0.,6.,14.,22.,34.,44.):
        stop=start+4
        m=(t>=start)&(t<=stop)
        if t[-1]<stop:
            continue
        def avg(x): return np.trapz(x[m],t[m])/(t[m][-1]-t[m][0])/w
        def net_rate(x):
            v=x[m]
            return (v[-1]-v[0])/(t[m][-1]-t[m][0])/w
        emit(test='cavity_budget',tag=folder.name,start=start,stop=stop,
             active=avg(active),tail=avg(tail),feed=avg(feed),cut_export=avg(cut),
             receiver=avg(receiver),net_feed_rate_per_input=net_rate(feed),
             active_loss_rate_per_input=-net_rate(active))
    obs=np.column_stack((d['near_velocity'],d['far_velocity']))
    for start in (7.,15.,27.,37.):
        if t[-1]>=start+8:
            emit(test='cavity_prediction',tag=folder.name,**fit(t,obs,start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folders',type=Path,nargs='+')
    for folder in p.parse_args().folders:
        analyse(folder)
