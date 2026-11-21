"""Two mobile ends at a marked wall station, and the complete local support force.

A cold meridional row is a trial reduction of the neutral material. A row of
incoming brane constituents passes a single wall-side constituent. The axial
coordinates respond to drag, bond forces and material elasticity. The imposed
lateral row position and the wall fixture implement the displacement test.
"""
import argparse
import hashlib
import json

import numpy as np

from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import emit
from mouth_two_depth_load import DEPTHS
from neutral_patch_conversion import NEUTRAL
from neutral_patch_sites import BETWEEN, GAP, PITCH_UNIT, WITHIN, morse


class Station:
    def __init__(self,args):
        trace,intake=cubic_trace(34)
        w=DEPTHS[1]
        slope=trace.geometry.radius(w)[1]
        self.speed=float(trace.speed(w)*np.hypot(1.,slope)/.01)
        self.order=float(trace.sample(trace.chi,w)[0])
        self.depth=BETWEEN[0]*self.order
        frame=json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']
        lower,upper=frame['sites']
        xy=np.array(frame['rested_xy'])
        # Tap perpendicular to the saved tangential offset. The other
        # tangential separation and the interplane gap stay fixed in this test.
        self.gap2=GAP**2+sum((xy[upper]-xy[lower])**2)
        self.k_axial=args.axial_grip
        self.dt=args.dt
        count=args.padding+int(np.ceil(self.speed*12))
        self.tag=count
        self.initial=np.arange(-count,count+1,dtype=float)*PITCH_UNIT
        # Linear stretch modulus of the same-plane Morse bond at equilibrium.
        self.cohesion=2*WITHIN[0]*WITHIN[1]**2
        emit(test='local_preparation',mesh=34,station=float(w),intake=intake,order=self.order,
             carrier_in_pitches_per_time=self.speed,length_per_pitch_in_R=.01,
             gap_squared=self.gap2,central_ids=[lower,upper],axial_wall_grip=self.k_axial,
             lateral_wall_grip=1.,site_drag=1.,axial_row_modulus=self.cohesion,
             count=len(self.initial),dt=self.dt,temperature=0.)

    def rates(self,state,displacement):
        # x, z are displacements in the local tangential and meridional
        # directions, in patch units. Zero z marks the fixed wall station.
        wall_x=state[:,0]
        wall_z=state[:,1]
        upper_z=state[:,2:]
        sep=upper_z-wall_z[:,None]
        slip=displacement-wall_x
        length=np.sqrt(self.gap2+sep**2+slip[:,None]**2)
        energy,tension=morse(length,self.depth,BETWEEN[1],BETWEEN[2])
        stiffness=tension/length
        axial=stiffness*sep
        lateral=stiffness.sum(axis=1)*slip
        rates=np.empty_like(state)
        rates[:,0]=lateral-wall_x
        rates[:,1]=axial.sum(axis=1)-self.k_axial*wall_z
        rates[:,2:]=self.speed-axial
        stretch=self.cohesion*(np.diff(upper_z,axis=1)-PITCH_UNIT)
        rates[:,2:-1]+=stretch
        rates[:,3:]-=stretch
        return rates,lateral,energy

    def run(self,history,track):
        dt=self.dt
        time=np.arange(round(12/dt)+1)*dt
        amps=1e-4*np.array([0.,1.,-1.,.5,-.5])
        state=np.zeros((len(amps),len(self.initial)+2))
        state[:,2:]=self.initial
        position=np.empty((len(time),len(amps)))
        total=np.empty_like(position)
        first_loss=None
        loss_reason=None
        distinct=set()
        points=[]
        for i,t in enumerate(time):
            k1,force,energies=self.rates(state,amps*history(t))
            position[i]=state[:,0]
            total[i]=force
            if track:
                nearest=int(np.argmin(abs(state[0,2:]-state[0,1])))
                registered=-energies[0,nearest]/self.depth>.6
                if registered:
                    distinct.add(nearest)
                if first_loss is None and (nearest!=self.tag or not registered):
                    first_loss=[t,state[0,1],state[0,2+self.tag]]
                    loss_reason='different nearest' if nearest!=self.tag else 'well threshold'
                if i in {0,round(.05/dt),round(.1/dt),round(.25/dt),round(1/dt),round(4/dt),round(8/dt),len(time)-1}:
                    points.append([t,state[0,1],state[0,2+self.tag],nearest-self.tag,
                                   position[i,1],state[1,1],state[1,2+self.tag]])
            if i+1<len(time):
                k_mid=self.rates(state+.5*dt*k1,amps*history(t+.5*dt))[0]
                state+=dt*k_mid
        if first_loss is None and track:
            raise ValueError('Observation ended before loss of original partner')
        if not np.isfinite(total).all() or np.min(np.diff(state[:,2:],axis=1))<=0:
            raise ArithmeticError('Invalid material row')
        # Positive is the load on the wall. Recover all of it independently
        # from measured spring extension and wall-side velocity (unit drag).
        support=position+np.gradient(position,dt,axis=0,edge_order=2)
        def odd(values):
            return (values[:,1]-values[:,2])/(2*amps[1])
        reaction=odd(support)
        direct=odd(total)
        half_error=np.max(abs(direct-(total[:,3]-total[:,4])/amps[1]))
        support_error=np.max(abs(reaction-direct))
        if half_error>2e-5 or support_error>3e-3 or np.max(abs(total[:,0]))>1e-13:
            raise ArithmeticError(f'Local force controls failed: {half_error}, {support_error}')
        if track:
            emit(test='tracking',samples=points,columns=['t','unforced_wall_end_z','unforced_original_brane_end_z','nearest_index_from_original','held_wall_end_x','held_wall_end_z','held_original_brane_end_z'],
                 units='patch lengths from the fixed middle station',held_row_x=amps[1],
                 first_loss=first_loss,loss_reason=loss_reason,distinct_registered_partners=len(distinct))
        return time,odd(position),reaction,direct,support_error,half_error


def pulse(t):
    return np.sin(np.pi*t/8)**2 if t<8 else 0.


def experiment(args):
    station=Station(args)
    t,wall,reaction,_,ctrl,half=station.run(lambda _:1.,True)
    free_slip=1-wall
    if np.min(free_slip)<.02:
        raise ValueError('Calibrating slip too small')
    modulus=reaction/free_slip
    if min(modulus)<=0:
        raise ValueError('Nonpositive measured coupling')
    # The calibrated coefficient is frozen before the direct pulse trial.
    emit(test='held_calibration',coefficient_range=[min(modulus),max(modulus)],
         smallest_slip_per_applied=min(free_slip),peak_support_error=ctrl,
         halving_error=half,coefficient_sha256=hashlib.sha256(modulus.astype('<f8').tobytes()).hexdigest())
    x=0.
    prediction=np.empty(len(t))
    for i,now in enumerate(t):
        force=modulus[i]*(pulse(now)-x)
        prediction[i]=force
        if i+1<len(t):
            mid_x=x+args.dt/2*(force-x)
            mid_force=.5*(modulus[i]+modulus[i+1])*(pulse(now+args.dt/2)-mid_x)
            x+=args.dt*(mid_force-mid_x)
    _,_,_,direct,ctrl,half=station.run(pulse,False)
    err=max(abs(prediction-direct))/max(abs(direct))
    if not np.isfinite(err) or err>.002:
        raise ArithmeticError('Held response failed to transfer to local pulse')
    indices=np.rint(np.array([4.,8.,8.5,9.,12.])/args.dt).astype(int)
    emit(test='local_pulse',relative_peak_error=err,t=t[indices].tolist(),
         predicted_per_displacement=prediction[indices].tolist(),direct_per_displacement=direct[indices].tolist(),
         peak_support_error=ctrl,halving_error=half)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--axial-grip',type=float,choices=(0.,1.),required=True)
    p.add_argument('--dt',type=float,choices=(.001,.002),required=True)
    p.add_argument('--padding',type=int,choices=(48,80),required=True)
    experiment(p.parse_args())
