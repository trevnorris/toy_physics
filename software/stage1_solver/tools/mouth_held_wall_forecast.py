"""Measure two held-bend responses, then forecast a free middle push with the bare carrier.

The forecast process reads only the held-force histories and the original
carrier. Its comparison with an earlier free run is performed separately.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import splu

from mouth_order_wave_pulse import ROOT, digest, emit
from mouth_two_depth_load import DEPTHS, PEAK_FORCE


def measure(args, sub, dt, steps):
    from mouth_bend_renewal import volume_cost, wall_data
    from mouth_bend_static_shear import ANGULAR, archive_inputs
    order = ROOT/f'mouth_free_order_cubic_n{args.mesh}_d04/order_seed1.55.npz'
    geo, chi, _ = archive_inputs(order,None)
    with np.load(ROOT/f'mouth_flared_cubic_n{args.mesh}/field.npz') as saved:
        np.testing.assert_array_equal(saved['points'],geo.points)
        field = saved['psi']
    edges, live, bulk = volume_cost(geo,chi,sub,5,1e-8)
    contact,load,B,G,drift,tau,_ = wall_data(geo,chi,field,edges,live,5)
    if args.law == 'local':
        drift = np.zeros_like(drift)
    S = bulk+contact
    norm = 1/np.sqrt(S.diagonal())
    fac = linalg.cho_factor(norm[:,None]*S*norm)
    def solve(rhs):
        return norm[:,None]*linalg.cho_solve(fac,norm[:,None]*rhs)
    mb,mp = solve(load),solve(B)
    speed = drift+G@mp/tau
    advance = linalg.lu_factor(np.eye(len(speed))/dt+.5*speed)
    middle = DEPTHS[1]
    j = np.searchsorted(edges,middle)-1
    xi = (middle-edges[j])/(edges[j+1]-edges[j])
    clip = np.zeros(len(edges)); clip[j:j+2]=[1-xi,xi]
    cut = -ANGULAR*geo.radius(np.array([middle]))[0][0]*clip[live]
    read = np.vstack((ANGULAR*load.T,cut))
    # Full, reversed, half and undriven physical displacement copies, both shapes.
    scales = np.array([1.,-1.,.5,-.5,0.])*1e-5
    b = np.hstack([np.eye(2)*v for v in scales])
    p = np.zeros((len(speed),len(b.T)))
    source = G@mb@b/tau
    response = []
    err=0.
    for _ in range(steps):
        dp = linalg.lu_solve(advance,source-speed@p)
        readings = (read@(mb@b-mp@(p+dp/2))).reshape(3,5,2)
        full=(readings[:,0]-readings[:,1])/(2*scales[0])
        half=(readings[:,2]-readings[:,3])/scales[0]
        err=max(err,float(np.max(abs(full-half))),float(np.max(abs(readings[:,4]))))
        response.append(full)
        p+=dp
    if err>1e-9 or not np.isfinite(response).all():
        raise ArithmeticError('Hold response failed reversal/halving or finite-value checks')
    with args.file.open('xb') as fh:
        np.savez_compressed(fh,response=response,mesh=args.mesh,dt=dt,law=args.law)
    emit(test='hold',mesh=args.mesh,law=args.law,clock=tau,dt=dt,amplitudes=scales.tolist(),
         hold_samples=steps,hold_control_error=err,initial_midpoint=response[0].tolist(),
         last_midpoint=response[-1].tolist(),output=str(args.file),sha256=digest(args.file))


def forecast(args,extension,dt,steps):
    from mouth_carrier_lateral_response import prepare
    from mouth_free_bend_drive import Response
    from mouth_free_bend_release import pencil
    source=ROOT/f'mouth_flared_cubic_n{args.mesh}'
    dom,field=prepare(SimpleNamespace(source=source))
    energy,time=pencil(Response(dom,field,extension),0.)
    with np.load(args.file) as data:
        if (data['mesh']!=args.mesh or data['dt']!=dt or data['law']!=args.law):
            raise ValueError('Mismatched hold preparation')
        held=data['response']
    if held.shape!=(steps,3,2):
        raise ValueError('Incomplete held kernel')
    # The step was held at each midpoint. First differences give the response
    # to a single occupied input bin, including its immediate elastic term.
    kernel=np.diff(held,axis=0,prepend=np.zeros_like(held[:1]))
    direct=kernel[0,:2]
    padding=energy.shape[0]-2
    extra=sparse.block_diag((sparse.csc_matrix((padding,padding)),direct),format='csc')
    factor=splu(((energy+extra)/2-time/dt).tocsc())
    q=np.zeros(energy.shape[0]); history=np.zeros((steps,2)); band=0.
    for i in range(steps):
        mid=(i+.5)*dt
        drive=PEAK_FORCE*np.sin(np.pi*mid/8)**2 if mid<8 else 0.
        memory=np.einsum('sij,sj->i',kernel[1:i+1][::-1],history[:i])
        rhs=-(energy@q)
        rhs[-3:-1]+=drive  # point load at the middle: x_middle = X + first bend
        rhs[-2:]-=direct@q[-2:]+memory[:2]
        delta=factor.solve(rhs)
        history[i]=q[-2:]+delta[-2:]/2
        band_mid=kernel[0,2]@history[i]+memory[2]
        band=2*band_mid-band  # piecewise-linear states in the midpoint integrator
        q+=delta
        end=(i+1)*dt
        if end in (4.,8.,16.,18.5,steps*dt):
            emit(test='prediction',time=end,band_per_peak=band/PEAK_FORCE,
                 gauges_per_peak=(np.cumsum(q[-3:])/PEAK_FORCE).tolist())
    emit(test='forecast_source',law=args.law,mesh=args.mesh,dt=dt,input_sha256=digest(args.file))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('hold','forecast'))
    parser.add_argument('--law',choices=('local','travelling'),required=True)
    parser.add_argument('--mesh',type=int,choices=(26,34),required=True)
    parser.add_argument('--file',type=Path,required=True)
    args=parser.parse_args()
    extension,sub,dt,steps={26:('compact',2,.25,114),34:('wide',3,.125,227)}[args.mesh]
    if args.operation=='hold':
        measure(args,sub,dt,steps)
    else:
        forecast(args,extension,dt,steps)
