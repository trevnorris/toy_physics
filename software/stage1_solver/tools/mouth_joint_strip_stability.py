"""Test the sign of native bond work on the actual joint-strip transport.

Translate all members of each section together within one local tangent, with
the wall tips and exterior displacement held fixed. Negative curvature on this
admissible snapshot rules out a local energy minimum for the centered native
continuation, even after charging the perturbation its full attachment cost.
"""
import argparse
from pathlib import Path
import numpy as np

from mouth_joint_strip_hold import HeldWall
from mouth_middle_tag_grip import cubic_trace
from mouth_order_wave_pulse import emit
from neutral_patch_sites import COUNT, WITHIN, morse
from neutral_patch_turning_neighbor import ENTRY


def test_shape(system, trace, pos, center, width, wavelength):
    row,w=system.strip,system.wall
    j=np.arange(len(pos))-center
    # Smooth compact packet in material coordinates, with zero slope at ends.
    mask=abs(j)<width
    u=np.zeros_like(pos)
    u[mask]=np.cos(np.pi*j[mask]/(2*width))**2*np.cos(2*np.pi*j[mask]/wavelength)
    if center < width+2 or center+width+2>=len(pos):
        raise ArithmeticError('Perturbation touches open boundary')
    halfpos=(pos[1:]+pos[:-1])/2
    halfu=(u[1:]+u[:-1])/2
    axial=[]; hoop=[]; du=[]; orders=[]
    def bonds(s1,s2,b,v):
        mid=(s1+s2)/2
        use=v!=0
        mid=mid[use]
        depths=np.interp(mid,system.strip.arc,row.depths)
        axial.extend((s2-s1)[use])
        hoop.extend(b*trace.geometry.radius(depths)[0]/trace.geometry.radius(.5)[0])
        du.extend(v[use])
        orders.extend(np.interp(mid,row.arc,row.order)/ENTRY)
    # Two alternating hoop rows in the material unit cell. Each bond appears
    # once per two sites. Areas and the number of active ports match the wall model.
    bonds(pos[:-1],pos[1:],0.,np.diff(u))
    bonds(halfpos[:-1],halfpos[1:],0.,np.diff(halfu))
    for sign in (-1,1):
        b=sign*np.sqrt(3)/2
        bonds(pos[:-1],halfpos,b,halfu-u[:-1])
        bonds(pos[1:],halfpos,b,halfu-u[1:])
    ax,by,du,order=map(np.asarray,(axial,hoop,du,orders))
    r=np.hypot(ax,by)
    t=morse(r,*WITHIN)[1]
    depth,alpha,rest=WITHIN
    v=np.exp(-alpha*(r-rest))
    second=2*depth*alpha**2*v*(2*v-1)
    native=COUNT/2*np.sum(order*du**2*(t/r+(second-t/r)*(by/r)**2))
    def energy(h):
        return COUNT/2*(order @ morse(np.hypot(ax,by+h*du),*WITHIN)[0])
    checks=[]
    for h in (2e-4,1e-4):
        checks.append((energy(h)+energy(-h)-2*energy(0.))/h**2)
    np.testing.assert_allclose(checks,native,rtol=1e-5,atol=3e-4)
    idx=np.searchsorted(pos,system.stations)
    mix=(system.stations-pos[idx-1])/(pos[idx]-pos[idx-1])
    potential=w.multiplicity*w.k.sum(axis=1)
    wall=float(potential @ ((1-mix)*u[idx-1]**2+mix*u[idx]**2))
    # Independent quadratic form from each station's ordinary contact force.
    material_force=np.zeros(len(pos))
    np.add.at(material_force,idx-1,-potential*(1-mix)*u[idx-1])
    np.add.at(material_force,idx,-potential*mix*u[idx])
    np.testing.assert_allclose(wall,-u@material_force,rtol=1e-12)
    ref=COUNT*system.shear/ENTRY*np.sum(np.interp(halfpos,row.arc,row.order)*np.diff(u)**2)
    norm=COUNT*(u@u)
    emit(test='joint_trial_curvature',center_depth=float(np.interp(pos[center],row.arc,row.depths)),
         width_each_side=width,wavelength=wavelength,section_count=len(pos),native=native,
         native_centered_differences=checks,wall_if_held=wall,reference_shear=ref,
         total_native=native+wall,total_reference=ref+wall,
         native_per_unit_drag=(native+wall)/norm)


def run(args):
    from mouth_connected_passage import MIDDLE
    trace,_=cubic_trace(args.mesh)
    system=HeldWall(args)
    row=system.strip
    warm=4*system.wall.transit
    steps=int(np.ceil(warm/args.dt))
    for _ in range(steps):
        row.exchange(row.advance(warm/steps))
    tagged=None
    now=0.
    while True:
        new=row.advance(args.dt)
        if tagged is None:
            c=np.flatnonzero((row.pos<0)&(new>=0))
            if len(c): tagged=int(c[0]+row.id0)
        if tagged is not None:
            idx=tagged-row.id0
            if row.pos[idx]<row.gates[1]<=new[idx]:
                frac=(row.gates[1]-row.pos[idx])/(new[idx]-row.pos[idx])
                snapshot=row.pos+frac*(new-row.pos)
                emit(test='snapshot',id=tagged,time=now+frac*args.dt,depth=MIDDLE)
                for center_depth in (.65,.9,MIDDLE):
                    station=np.interp(center_depth,row.depths,row.arc)
                    center=int(np.argmin(abs(snapshot-station)))
                    for wave in (4.,8.,16.):
                        test_shape(system,trace,snapshot,center,12,wave)
                break
        row.exchange(new)
        now+=args.dt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,default=34)
    p.add_argument('--lead',type=float,default=8.)
    p.add_argument('--exit-lead',type=float,default=8.)
    p.add_argument('--dt',type=float,default=1/512)
    p.set_defaults(coupon=Path('_scratch/s11c/mobile_neighbor_probe_02/response.npz'),sub=5,
                   coverage=1.,brane_grip=0.,phase=0.,links=True,amplitude=1e-5)
    run(p.parse_args())
