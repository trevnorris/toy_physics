"""Independently calibrate the brane end against held neighbors on the neutral patch."""
import argparse
import json

import numpy as np

from mouth_fixed_station_pair import NEUTRAL
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import digest,emit
from mouth_segment_contact_hold import AMPLITUDE,PITCH,Segment,coupon
from mouth_two_depth_load import DEPTHS
from neutral_patch_sites import COUNT, WITHIN, morse


def measure():
    ids,_,_=coupon()
    coords=np.array(json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']['rested_xy'])
    upper=coords[COUNT:]
    tangents=[]
    for _,index in ids:
        rel=np.delete(upper-upper[index-COUNT],index-COUNT,axis=0)
        vectors=np.column_stack((rel,np.zeros(len(rel))))
        radius=np.linalg.norm(vectors,axis=1)
        depth,alpha,rest=WITHIN
        exp=np.exp(-alpha*(radius-rest))
        _,first=morse(radius,*WITHIN)
        second=2*depth*alpha**2*exp*(2*exp-1)
        # Trace/3 under the same assumed isotropic orientation distribution
        # in the three surface tangents as used for the cross bonds.
        k=np.sum((second+2*first/radius)/3)
        def energy(axis,h):
            return np.sum(morse(np.linalg.norm(vectors-h*axis,axis=1),*WITHIN)[0])
        for h in (1e-4,5e-5):
            direct=np.mean([(energy(e,h)+energy(e,-h)-2*energy(e,0))/h**2 for e in np.eye(3)])
            np.testing.assert_allclose(direct,k,atol=3e-5,rtol=8e-7)
        tangents.append(float(k))
    if min(tangents)<=0: raise ValueError('Nonpositive brane cage stiffness in coupon')
    return ids,np.array(tangents)


def run(args):
    ids,cage=measure()
    segment=Segment(args)
    trace,_=cubic_trace(args.mesh)
    station=.5 if args.pitch_at=='waist' else DEPTHS[1]
    a,slope=trace.geometry.radius(np.array([station]))
    period=PITCH/(trace.speed(station)*np.hypot(1.,slope[0]))
    a0=trace.geometry.radius(np.array([.5]))[0]
    segment.age*=segment.period/period
    segment.period=period
    segment.multiplicity*=float((a0/a)[0]**2)
    segment.offset=int(np.ceil(max(segment.age)))+3
    segment.parcels=segment.offset+int(np.ceil(segment.transit/period))+5
    segment.labels=np.arange(segment.families)*segment.parcels
    segment.grip=np.repeat(cage,segment.parcels)
    emit(test='frozen_neighbor_calibration',archive=digest(NEUTRAL),ids=ids.tolist(),cage=cage.tolist(),
         pitch_station=station,period=period,coverage=args.coverage,wall_grip=1.,dt=args.dt)
    x=np.zeros((2,segment.points,segment.families))
    y=np.zeros((2,segment.parcels*segment.families))
    t=0.
    segment.rates(t,x,y,True)
    for end in (.0625,28.3125,segment.transit):
        steps=int(np.ceil((end-t)/args.dt-1e-12))
        dt=(end-t)/steps
        start=t
        for i in range(steps):
            now=start+i*dt
            dx,dy=segment.rates(now,x,y)
            vx,vy=segment.rates(now+dt/2,x+dt*dx/2,y+dt*dy/2)
            x+=dt*vx; y+=dt*vy
        t=end
        force=segment.rates(t,x,y,True)
        if t==.0625: initial=force
        else: emit(test='cage_retention',time=t,ratio=force/initial)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pitch-at',choices=('waist','middle'),required=True)
    p.add_argument('--dt',type=float,default=1/256)
    p.add_argument('--coverage',type=float,default=1.)
    p.set_defaults(mesh=34,sub=5,phase=0.,amplitude=AMPLITUDE,brane_grip=1.)
    run(p.parse_args())
