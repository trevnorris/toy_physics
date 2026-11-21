"""Partition the middle-push wall force by the passage of a carrier cohort.

At pulse shutoff mark the carrier entering the waist. Integrate the actual
wall traction on either side of that front. Only the traveling-contact law
assigns these carrier labels to contact transport. The forecast uses the
original coupled mouth and applied drive with no observed force input.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from mouth_bend_renewal import volume_cost
from mouth_bend_static_shear import ANGULAR
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import emit
from mouth_two_depth_renewal import CoupledMouth, DEPTHS, PEAK_FORCE


def response(args):
    model=CoupledMouth(args)
    trace,_=cubic_trace(args.mesh)
    # Restore the full nodal order (the trajectory object uses only wall values).
    path=f'_scratch/s11c/mouth_free_order_cubic_n{args.mesh}_d04/order_seed1.55.npz'
    with np.load(path) as data:
        edges,keep,_=volume_cost(trace.geometry,data['chi'],args.sub,5,1e-8)
    waist, middle = edges[0], DEPTHS[1]
    crossing=trace.transit(waist,middle,8)
    roots,weights=leggauss(8)
    def side_force(left,right,stress):
        cuts=np.r_[left,edges[(edges>left)&(edges<right)],right]
        width=np.diff(cuts)
        w=cuts[:-1,None]+width[:,None]*(1+roots)/2
        m=np.interp(w,edges,stress)
        gr=np.diff(stress)/np.diff(edges)
        idx=np.clip(np.searchsorted(edges,w,side='right')-1,0,len(gr)-1)
        a,slope=trace.geometry.radius(w)
        return -ANGULAR*np.sum(width[:,None]*weights/2*(slope*m+a*gr[idx]))
    dt=args.dt
    advance=splu((model.energy/2-model.time/dt).tocsc())
    state=np.zeros(len(model.forces))
    for j in range(1,round(model.transit/dt)+1):
        t=j*dt
        mid=t-dt/2
        strength=np.sin(np.pi*mid/8)**2 if mid<8 else 0.
        force=PEAK_FORCE*strength*model.forces[:,0]
        state+=advance.solve(force-model.energy@state)
        if t not in (8.,16.,18.5,round(model.transit/dt)*dt):
            continue
        age=t-8.
        if age==0:
            front=waist
        elif age>=crossing:
            front=middle
        else:
            front=brentq(lambda w:trace.transit(waist,w,8)-age,waist,middle)
        stress=np.zeros(len(edges))
        stress[keep]=(model.dual_bend@state[model.nc-2:model.nc]-model.dual_ref@state[model.nc:])/PEAK_FORCE
        fresh=side_force(waist,front,stress)
        older=side_force(front,middle,stress)
        band=-ANGULAR*np.diff(trace.geometry.radius(np.array([waist,middle]))[0]
                              *np.interp([waist,middle],edges,stress))[0]
        if not np.isclose(band,fresh+older,atol=1e-9,rtol=1e-8):
            raise ArithmeticError('Area integral did not reproduce full wall resultant')
        emit(test='middle_force',law=args.reference,mesh=args.mesh,t=t,
             last_original_parcel_crosses_middle_at=8+crossing,
             front_clipped_to_middle=float(front),force_per_peak=band,
             force_on_waist_side_of_carrier_front=fresh,force_on_deep_side_of_carrier_front=older)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),required=True)
    p.add_argument('--reference',choices=('travelling','local'),required=True)
    args=p.parse_args()
    args.extension,args.sub,args.dt={26:('compact',2,.25),34:('wide',3,.125)}[args.mesh]
    response(args)
