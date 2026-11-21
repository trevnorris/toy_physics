"""Refine the segment hold and test which station sets the pitch of arriving sites.

A row observed with a one-pitch spacing at the middle carries fewer labels
per unit time than one prepared with that spacing at the faster waist.
The carriers, wall-site coverage and microscopic force constants stay fixed.
"""
import argparse

import numpy as np

from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import emit
from mouth_segment_contact_hold import AMPLITUDE, PITCH, Segment
from mouth_two_depth_load import DEPTHS


def read_segment(args):
    segment=Segment(args)
    trace,_=cubic_trace(args.mesh)
    depth=.5 if args.pitch_at=='waist' else DEPTHS[1]
    radius,slope=trace.geometry.radius(np.array([depth]))
    physical_speed=trace.speed(depth)*np.hypot(1.,slope[0])
    period=PITCH/physical_speed
    # Age in time is measured from the waist in both preparations.
    segment.age*=segment.period/period
    segment.period=period
    waist_radius=trace.geometry.radius(np.array([.5]))[0][0]
    # A labeled incoming packet now represents one pitch of the wall
    # hypersurface area at the indicated depth in a given angular stream.
    segment.multiplicity*=(waist_radius/radius[0])**2
    segment.offset=int(np.ceil(max(segment.age)))+3
    segment.parcels=segment.offset+int(np.ceil(segment.transit/period))+5
    segment.labels=np.arange(segment.families)*segment.parcels
    emit(test='prepared_stream',pitch_at=depth,inlet_period=period,
         labels_per_model_time=1/period,parcel_area_factor=(radius[0]/waist_radius)**2,
         coverage_multiplier=args.coverage,brane_grip=args.brane_grip,dt=args.dt,sub=args.sub)
    x=np.zeros((2,segment.points,segment.families))
    y=np.zeros((2,segment.parcels*segment.families))
    early=None
    start=0.
    segment.rates(0.,x,y,True)
    for end in (.0625,28.3125,segment.transit):
        steps=int(np.ceil((end-start)/args.dt-1e-12))
        dt=(end-start)/steps
        for j in range(steps):
            now=start+j*dt
            # Explicit trapezoid corrector uses endpoint slopes; baseline
            # uses explicit midpoint. Do not interpolate the resulting forces.
            dx,dy=segment.rates(now,x,y)
            ex,ey=segment.rates(now+dt,x+dt*dx,y+dt*dy)
            x+=dt*(dx+ex)/2
            y+=dt*(dy+ey)/2
        force=segment.rates(end,x,y,True)
        if early is None:
            early=force
        else:
            emit(test='pitch_retention',pitch_at=depth,time=end,ratio=force/early)
        start=end


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pitch-at',choices=('waist','middle'),required=True)
    p.add_argument('--sub',choices=(3,4,5,6),type=int,default=5)
    p.add_argument('--dt',type=float,default=1/256)
    p.add_argument('--coverage',type=float,default=1.)
    p.set_defaults(mesh=34,amplitude=AMPLITUDE,phase=0.,brane_grip=1.)
    read_segment(p.parse_args())
