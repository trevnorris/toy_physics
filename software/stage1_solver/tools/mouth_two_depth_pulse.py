"""Independent smooth, single-pulse test of middle and distal loads on the unsupported mouth."""
import argparse
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_response import prepare
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response
from mouth_free_bend_release import pencil
from mouth_two_depth_load import PEAK_FORCE, observe_cuts


def main(args):
    dom, field=prepare(args)
    response=Response(dom,field,args.extension)
    energy,time=pencil(response,0.)
    read=np.tri(3)
    rhs=np.zeros((energy.shape[0],4))
    strengths=PEAK_FORCE*np.array([1.,1.,-.5,-.5])
    rhs[-3:]=read[[1,2,1,2]].T*strengths
    dt=args.dt
    factor=splu((energy*.5-time/dt).tocsc())
    state=np.zeros_like(rhs)
    x_history=[np.zeros((3,2))]
    work=np.zeros(4)
    abs_work=np.zeros(4)
    max_error=0.
    for k in range(round(16/dt)):
        t=(k+.5)*dt
        force=np.sin(np.pi*t/8)**2 if t<8 else 0.
        delta=factor.solve(force*rhs-energy@state)
        state+=delta
        dw=force*np.sum(rhs*delta,axis=0)
        work+=dw
        abs_work+=abs(dw)
        x=(read@state[-3:])/strengths
        np.testing.assert_allclose(x[:,:2],x[:,2:],rtol=1e-8,atol=1e-9)
        x_history.append(x[:,:2])
        t=(k+1)*dt
        if t in (4.,8.,16.):
            stored=.5*np.sum(state*(energy@state),axis=0)
            error=abs(stored-work)/abs_work
            max_error=max(max_error,max(error))
            emit(test='smooth_pulse',time=t,displacement_per_peak=x[:,:2].tolist(),
                 energy=stored.tolist(),work=work.tolist(),relative_error=error.tolist())
        if t==4.:
            # Real-time snapshot at peak force. The helper treats a real state
            # as a Fourier component with zero sine quadrature.
            observe_cuts(dom,field,response,state[:,:2]/PEAK_FORCE,2)
    if not np.isfinite(max_error) or max_error>2e-5:
        raise ArithmeticError('Short smooth-pulse energy balance failed')
    # An operational timing marker for a transient, independent of phase wrapping.
    series=np.array(x_history)[:round(8/dt)+1]
    times=np.arange(len(series))*dt
    crossings=np.zeros((3,2))
    for station in range(3):
        for drive in range(2):
            x=series[:,station,drive]
            scaled=x/x[-1]
            i=np.flatnonzero(scaled>=.5)[0]
            crossings[station,drive]=times[i-1]+dt*(.5-scaled[i-1])/(scaled[i]-scaled[i-1])
    emit(test='pulse_timing',width=8.,force_peak_time=4.,half_endpoint_time=crossings.tolist(),
         definition='First time reaching half the signed displacement at force shutoff')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    p.add_argument('--dt',type=float,choices=[.25,.125],required=True)
    main(p.parse_args())
