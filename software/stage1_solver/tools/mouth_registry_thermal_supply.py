"""Local flare traction for an independently specified thermal supply of facing sites.

The three inlet position distributions have the same transverse thermal speed
distribution. All follow the same Hamiltonian transverse motion downstream.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_renewing_registry_probe import PeriodicLining, path, read_inputs, cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, GAP, TEMPERATURE, emit, morse


def surface_energy(lining,xy):
    energy=np.zeros(xy.shape[:-1])
    for i in range(-9,11):
        for j in range(-9,11):
            delta=xy-np.array([i,j])@lining.basis
            r=np.sqrt((delta*delta).sum(axis=-1)+GAP**2)
            energy+=morse(r,*BETWEEN)[0]
    return energy


def main(args):
    patch,_,_=read_inputs(34)
    trace,_=cubic_trace(34)
    time,arc,chi=path(trace,args.steps)
    dt=time[1]-time[0]
    wall=PeriodicLining(patch['pitch'],96)
    e_mass=2*patch['shared_scale']
    speed_std=np.sqrt(TEMPERATURE*e_mass)
    grid=(np.arange(args.phases)+.5)/args.phases
    f=np.stack(np.meshgrid(grid,grid,indexing='ij'),axis=-1).reshape(-1,2)
    xy=f@wall.basis
    u=surface_energy(wall,xy)
    # Condition on an imposed speed along a primitive lattice row. The
    # perpendicular Maxwell distribution is integrated at fixed inlet flux.
    z,wz=leggauss(args.speeds)
    z=8*z
    weights=wz*np.exp(-.5*z*z)
    weights/=weights.sum()
    np.testing.assert_allclose(weights@(z*z),1.,atol=1e-12)
    v0=speed_std*z
    n=np.array([0.,1.])
    table=wall.projected_table(n)
    cross=n@wall.inverse/PITCH
    along=arc[:,None]*(np.array([1.,0.])@wall.inverse/PITCH)
    force0=wall.gradient(xy)[:,1]/PITCH
    y=np.zeros((len(z),len(f)))
    v=np.broadcast_to(v0[:,None],y.shape).copy()
    def acc(k):
        coord=f.T[:,None,:]+along[k,:,None,None]+cross[:,None,None]*y
        return -chi[k]*e_mass/PITCH*wall.interpolate(table,coord)
    a=acc(0)
    for j in range(1,len(time)):
        y+=dt*v+.5*dt**2*a
        anew=acc(j)
        v+=.5*dt*(a+anew)
        a=anew
    impulse=v-v0[:,None]
    if not np.isfinite(impulse).all(): raise ArithmeticError('Nonfinite passage')
    mass_flow=GAP*PITCH/time[-1]
    emit(test='thermal_setup',kBT_in_well_units=TEMPERATURE,transverse_speed_std=speed_std,
         phase_grid=args.phases,speed_nodes=args.speeds,integration_steps=args.steps,
         interval=[1.55,1.65],time=time[-1],inlet_order=chi[0],law='order weighted through the patch')
    def normalized(x):
        return x/x.sum()
    for name, strength in [('uniform',0.),('local_thermal',chi[0]),('ordered_inherited',1.)]:
        position=normalized(np.exp(-strength*(u-u.min())/TEMPERATURE))
        mean_vchange=weights@impulse
        def expectation(v_weights,x_weights):
            return (v_weights@impulse)@x_weights
        # In wall-relative variables a displacement translates the *inlet
        # probability density*. A speed increment shifts the Gaussian mean.
        dk=-strength/TEMPERATURE*(force0-position@force0)
        offset=mass_flow*((position*dk)@mean_vchange)
        speed=mass_flow*expectation(weights*z/speed_std,position)
        # Cell periodicity makes the static derivative vanish for uniform phases.
        # An independent centered change of both probability densities checks the slopes.
        q0=5.e-6
        checks=[]
        for factor in (1.,.5):
            q=factor*q0
            dv=q/time[-1]
            pos=[normalized(np.exp(-strength*(surface_energy(wall,xy+sgn*q*n/PITCH)-u.min())/TEMPERATURE))
                 for sgn in (-1,1)]
            vel=[normalized(weights*np.exp(sgn*z*dv/speed_std)) for sgn in (-1,1)]
            fk=mass_flow*((pos[1]-pos[0])@mean_vchange)/(2*q)
            fz=mass_flow*expectation(vel[1]-vel[0],position)/(2*dv)
            checks.append([fk,fz])
        np.testing.assert_allclose(checks[-1],[offset,speed],rtol=1.e-4,atol=2.e-10)
        emit(test='thermal_inlet',preparation=name,prep_strength=strength,
             peak_to_trough_probability=position.max()/position.min(),
             holding_force_offset_slope=offset,holding_force_inlet_speed_slope=speed,
             finite_difference_full_half=checks,unperturbed_mean_momentum_loading=mass_flow*(position@mean_vchange),
             rms_transverse_speed_change=np.sqrt((weights@(impulse*impulse))@position))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--steps',type=int,required=True)
    p.add_argument('--phases',type=int,required=True)
    p.add_argument('--speeds',type=int,required=True)
    main(p.parse_args())
