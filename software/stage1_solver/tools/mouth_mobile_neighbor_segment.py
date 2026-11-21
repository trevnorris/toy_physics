"""Use the separately measured mobile-neighbor response in a held 0PN bend.

Every material label carries one twelve-site coupon. The nine contacted
members can exchange force through that coupon during passage. Unit site
drag is included once in the coupon response. Translation and planar rotation
are free in the main preparation.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import digest,emit
from mouth_segment_contact_hold import AMPLITUDE,PITCH,Segment
from mouth_two_depth_load import DEPTHS


def material(args):
    manifest=json.loads((args.coupon.parent/'outcome.json').read_text())
    validation=json.loads((args.coupon.parent/'limit-validation.json').read_text())
    receipt=json.loads((args.coupon.parent/'stdout').read_text().splitlines()[-1])
    if (not manifest['limitsVerified'] or manifest['exitCode'] or manifest['childOutcome']['exitCode']
        or not validation['verified'] or int(validation['actual']['memory.swap.max'])
        or receipt['sha256']!=digest(args.coupon)):
        raise ValueError('Coupon artifact did not pass guarded calibration')
    with np.load(args.coupon) as saved:
        base='legacy_' if args.material=='scalar' else ''
        rate=saved[base+'rates']
        b=saved[base+'B']
        ids=saved['labels'].copy()
    live=rate>0 if args.frame=='fixed' else np.ones_like(rate,dtype=bool)
    if not np.any(rate==0): raise ArithmeticError('Missing free translation')
    directions,modes,n=b[:,live].shape
    inject=b[:,live].reshape(directions*modes,n)
    read=inject.T/directions
    if args.frame=='free':
        np.testing.assert_allclose(read@inject,np.eye(n),atol=1e-10)
    return np.tile(rate[live],directions),inject,read,directions,ids


def prepare_segment(args):
    s=Segment(args)
    trace,_=cubic_trace(args.mesh)
    w=.5 if args.pitch_at=='waist' else DEPTHS[1]
    a,slope=trace.geometry.radius(np.array([w]))
    period=PITCH/(trace.speed(w)*np.hypot(1.,slope[0]))
    s.age*=s.period/period
    s.period=period
    a0=trace.geometry.radius(np.array([.5]))[0]
    s.multiplicity*=float((a0/a)[0]**2)
    s.offset=int(np.ceil(max(s.age)))+3
    s.parcels=s.offset+int(np.ceil(s.transit/period))+5
    s.labels=np.arange(s.families)*s.parcels
    return s


def exponential_weights(rate,dt):
    decay=np.exp(-rate*dt)
    integral=np.full_like(rate,dt)
    np.divide(-np.expm1(-rate*dt),rate,out=integral,where=rate!=0)
    return decay[None,:,None],integral[None,:,None]


def run(args):
    rates,B,read,norient,ids=material(args)
    s=prepare_segment(args)
    if len(ids)!=s.families: raise ValueError('Wrong material family count')
    emit(test='transfer',kind=args.material,frame=args.frame,pitch_at=args.pitch_at,
         coupon_sha256=digest(args.coupon),directions=norient,zero_components=int((rates==0).sum()),
         tagged_sites=ids[:,1].tolist(),period=s.period,parcels=s.parcels,
         coverage=args.coverage,amplitude=args.amplitude,dt=args.dt)
    x=np.zeros((2,s.points,s.families))
    z=np.zeros((2,len(rates),s.parcels))
    def tip(z):
        return (read@z).reshape(2,-1)
    def forces(t,x,z):
        dx,external=s.rates(t,x,tip(z))
        return dx,B@external.reshape(2,s.families,s.parcels)
    def observe(t,x,z):
        y=tip(z)
        load=s.rates(t,x,y,True)
        _,external=s.rates(t,x,y)
        external=external.reshape(2,s.families,s.parcels)
        f=B@external
        elastic=rates[None,:,None]*z
        dz=f-elastic
        source=float(np.sum(external*(read@dz)))
        stored=float(np.sum(elastic*dz)/norient)
        dissipated=float(np.sum(dz*dz)/norient)
        denom=max(abs(source),abs(stored),abs(dissipated),1e-25)
        err=abs(source-stored-dissipated)/denom
        if not np.isfinite(err) or err>1e-7:
            raise ArithmeticError('Mobile-neighbor work budget failed')
        emit(test='coupon_power',t=t,source=source,energy_rate=stored,drag_power=dissipated,scaled_error=err)
        return load
    samples={0.:observe(0.,x,z)}
    start=0.
    for end in (.0625,1.,8.,28.3125,s.transit):
        steps=int(np.ceil((end-start)/args.dt-1e-12))
        dt=(end-start)/steps
        half,half_gain=exponential_weights(rates,dt/2)
        full,full_gain=exponential_weights(rates,dt)
        for step in range(steps):
            t=start+step*dt
            dx,f=forces(t,x,z)
            mid_x=x+dt/2*dx
            mid_z=half*z+half_gain*f
            vx,mid_f=forces(t+dt/2,mid_x,mid_z)
            x+=dt*vx
            z=full*z+full_gain*mid_f
        samples[end]=observe(end,x,z)
        start=end
    emit(test='predicted_retention',early=samples[.0625],late=samples[28.3125],full=samples[s.transit],
         matched=samples[28.3125]/samples[.0625],transit=samples[s.transit]/samples[.0625],
         from_initial=samples[s.transit]/samples[0.])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coupon',type=Path,required=True)
    p.add_argument('--pitch-at',choices=('waist','middle'),required=True)
    p.add_argument('--material',choices=('planar','scalar'),default='planar')
    p.add_argument('--frame',choices=('free','fixed'),default='free')
    p.add_argument('--dt',type=float,default=1/128)
    p.add_argument('--sub',type=int,default=5)
    p.add_argument('--coverage',type=float,default=1.)
    p.add_argument('--amplitude',type=float,default=AMPLITUDE)
    p.set_defaults(mesh=34,phase=0.,brane_grip=0.)
    run(p.parse_args())
