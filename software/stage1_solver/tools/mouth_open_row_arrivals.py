"""Count actual tag crossings in an axially mobile open Morse row on the 0PN carrier."""
import argparse

import numpy as np
from scipy.linalg import solve_banded

from mouth_connected_passage import MIDDLE, WAIST, arc_speed, bond_derivatives, integrals, path_grid
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import emit
from mouth_segment_contact_hold import PITCH


def experiment(args):
    trace,intake=cubic_trace(args.mesh)
    w=path_grid(trace,48)
    arc=np.r_[0., integrals(lambda v:np.hypot(1., trace.geometry.radius(v)[1]), w)]/PITCH
    speed=arc_speed(trace,w)/PITCH
    gates=[0.,arc[-1]]
    carrier_time=trace.transit(WAIST,MIDDLE,8)
    exit_lead=args.lead if args.exit_lead is None else args.exit_lead
    left,right=-args.lead,arc[-1]+exit_lead
    pos=np.arange(left,right,args.feed_pitch,dtype=float)
    ids=np.arange(len(pos))
    events=[{},{}]
    def readings(i,start,stop):
        a=sorted((t,g,tag) for tag,(t,g) in events[i].items() if start<=t<stop)
        if len(a)<16: raise ArithmeticError('Too few crossing events')
        ts,gap,_=zip(*a)
        return dict(count=len(a),rate=(len(a)-1)/(ts[-1]-ts[0]),mean_gap=float(np.mean(gap)),
                    min_gap=min(gap),max_gap=max(gap),period_scatter=float(np.std(np.diff(ts))))
    # Four carrier passages for startup, followed by two complete observation windows.
    window=carrier_time
    steps=int(np.ceil(6*window/args.dt))
    dt=6*window/steps
    lo_stretch=float('inf')
    for step in range(steps):
        now=step*dt
        dist=np.diff(pos)
        tension,k=bond_derivatives(dist)
        velocity=np.interp(pos,arc,speed)
        velocity[:-1]+=tension
        velocity[1:]-=tension
        mat=np.empty((3,len(pos)))
        mat[1]=1.
        mat[1,:-1]+=dt*k
        mat[1,1:]+=dt*k
        mat[0,0]=mat[2,-1]=0.
        mat[0,1:]=mat[2,:-1]=-dt*k
        new=pos+dt*solve_banded((1,1),mat,velocity,check_finite=False)
        if np.any(np.diff(new)<=0) or np.any(~np.isfinite(new)):
            raise ArithmeticError('Open row lost ordering')
        lo_stretch=min(lo_stretch,float(np.min(np.diff(new))))
        for gate,log in zip(gates,events):
            for j in np.flatnonzero((pos<gate)&(new>=gate)):
                if j+1==len(new) or ids[j] in log:
                    raise ArithmeticError('Crossing identity invalid')
                mix=(gate-pos[j])/(new[j]-pos[j])
                gap=dist[j]+mix*(new[j+1]-new[j]-dist[j])
                log[int(ids[j])]=(now+mix*dt,float(gap))
        pos=new
        if pos[-1]>=right:
            live=pos<right
            ids,pos=ids[live],pos[live]
        if pos[0]>=left+args.feed_pitch:
            # An unloaded reservoir supplies a label when enough space opens.
            # Its insertion time is an outcome of the preceding constituent's motion.
            ids=np.r_[ids[0]-1,ids]
            pos=np.r_[pos[0]-args.feed_pitch,pos]
    emit(test='open_row',intake=intake,mesh=args.mesh,actual_dt=dt,
         feed_pitch=args.feed_pitch,inlet_lead=args.lead,outlet_lead=exit_lead,model_time=6*window,
         minimum_spacing=lo_stretch,waist_to_middle_carrier_time=window)
    for j in (4,5):
        stats=[readings(i,j*window,(j+1)*window) for i in (0,1)]
        trips=[events[1][tag][0]-t for tag,(t,_) in events[0].items()
               if j*window<=t<(j+1)*window and tag in events[1]]
        if not trips: raise ArithmeticError('No complete tracked passage')
        emit(test='observed_crossings',window=j,waist=stats[0],middle=stats[1],
             paired_ids=len(trips),mean_tag_transit=float(np.mean(trips)),
             first_last_id_at_middle=[max(events[1]),min(events[1])])
    # Inventory of all labels initially and finally in the measured band,
    # including the exactly-on-waist initial label if present.
    initial=np.arange(left,right,args.feed_pitch,dtype=float)
    inventory=lambda x:int(np.count_nonzero((x>=0)&(x<gates[1])))
    error=inventory(pos)-inventory(initial)-len(events[0])+len(events[1])
    emit(test='number_balance',initial=inventory(initial),final=inventory(pos),
         waist_crossings=len(events[0]),middle_crossings=len(events[1]),error=error)
    if error: raise ArithmeticError('Unbalanced label flux')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--dt',type=float,default=.002)
    p.add_argument('--lead',type=float,default=8.)
    p.add_argument('--exit-lead',type=float)
    p.add_argument('--feed-pitch',type=float,default=1.)
    experiment(p.parse_args())
