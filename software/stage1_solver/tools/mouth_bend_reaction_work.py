"""Power, carried spring energy and momentum on the wall of the held two-reservoir bend.

The contact energy is for its mean reference only. A tag carrying that spring
energy and a carrier doing mechanical work on its reference are distinct ledgers.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve

from mouth_bend_renewal import ROOT, volume_cost, wall_data
from mouth_bend_static_shear import ANGULAR, archive_inputs
from mouth_carrier_lateral_inertia import taper
from mouth_curved_phase import quadratic
from mouth_empty_bend_stiffness import INTERVALS
from mouth_order_wave_pulse import digest, emit


def run(a):
    label=f'delayed_n{a.mesh}'
    folder=ROOT/f'mouth_free_order_{label}_d04'
    order_file=folder/'order_seed1.55.npz'
    geo,order,carrier=archive_inputs(order_file,ROOT/f'bend_static_shear_delcarrier{a.mesh}_01')
    meta=json.loads((folder/'stdout').read_text().splitlines()[0])
    path=Path(meta['source'])
    if digest(path)!=meta['source_sha256']:
        raise ValueError('Carrier and order mismatch')
    with np.load(path) as data:
        np.testing.assert_array_equal(geo.points,data['points'])
        psi=data['psi']
    edges,keep,bulk=volume_cost(geo,order,a.sub,a.quad,a.floor)
    contact,load,B,G,drift,clock,_=wall_data(geo,order,psi,edges,keep,a.quad)
    history=solve(drift*clock,G)
    operator=bulk+contact+B@history
    scale=1/np.sqrt((bulk+contact).diagonal())
    def inverse(rhs):
        return scale*solve(scale[:,None]*operator*scale,scale*rhs)
    m=inverse(load[:,1])
    inlet=np.zeros(len(drift))
    if a.inlet=='relaxed':
        inlet[len(inlet)//2:]=1.
        influence=inverse(B@inlet)
        ratio=m[0]/influence[0]
        m-=ratio*influence
        inlet*=ratio
    references=history@m+inlet
    dual=np.zeros(len(edges))
    dual[keep]=m
    last=np.flatnonzero(keep)[-1]+1
    slip_z,slip_a=[p.reshape(-1,2) for p in np.split(references,2)]
    po,vo,cells=geo.order_blocks[2]
    wall_nodes=cells[::len(vo)-1][:,[0,3,6]]
    def samples(z,index):
        i=np.clip(np.searchsorted(po,z,side='right')-1,0,len(po)-2)
        basis,deriv=quadratic(2*(z-po[i])/np.diff(po)[i]-1)
        def trace(data):
            return (basis*data[wall_nodes[i]]).sum(axis=-1)
        chi=np.clip(trace(order),0,1)
        grad=(deriv*psi[wall_nodes[i]]).sum(axis=-1)*2/np.diff(po)[i]
        radius,slope=geo.radius(z)
        arc=np.hypot(1.,slope)
        velocity=(grad/trace(psi)).imag/arc
        width=np.diff(edges)[index]
        u=(z-edges[index])/width
        mode=dual[index]*(1-u)+dual[index+1]*u
        mode_z=(dual[index+1]-dual[index])/width
        tm=mode/radius**2
        ta=mode_z/(2*radius*arc)
        return radius,slope,arc,velocity,chi,tm,ta,mode
    root,weight=leggauss(a.quad+2)
    def integrate(lo,hi):
        breaks=np.unique(np.r_[lo,edges[(edges>lo)&(edges<hi)],hi])
        dz=np.diff(breaks)[:,None]
        z=breaks[:-1,None]+dz*(root+1)/2
        idx=np.searchsorted(edges,z,side='right')-1
        r,rp,arc,V,chi,tm,ta,_=samples(z,idx)
        live=(tm!=0)|(ta!=0)
        inv=np.divide(1.,chi,out=np.zeros_like(chi),where=live)
        area=ANGULAR*dz*weight/2*r**2*arc
        loss=(tm*tm+2*ta*ta)*inv/clock
        # Meridional reaction on the carrier from virtual reference transport.
        # Only the component conjugate to its specified tangential velocity is used.
        drag=-loss/V
        dpz=(slip_z[idx,1]-slip_z[idx,0])/np.diff(edges)[idx]
        dpa=(slip_a[idx,1]-slip_a[idx,0])/np.diff(edges)[idx]
        gradient_work=(tm*dpz+2*ta*dpa)*V/arc
        tx=rp/arc*tm+2*ta
        shape=1-taper(z,*INTERVALS[1])[0]
        # The reference is represented in discontinuous linear elements. Its
        # derivative has interface jumps as well as the interior gradient.
        # The DG upwind inlet flux uses the incoming reference at the waist.
        k=np.flatnonzero((edges[:last]>=lo)&(edges[:last]<hi))
        rr,_,_,speed,_,tz,ta_edge,_=samples(edges[k],k)
        rz=np.r_[inlet[0],slip_z[:-1,1]]
        ra=np.r_[inlet[len(inlet)//2],slip_a[:-1,1]]
        jump=ANGULAR*np.sum(rr**2*speed*(tz*(slip_z[k,0]-rz[k])+2*ta_edge*(slip_a[k,0]-ra[k])))
        continuous_part=float(np.sum(area*gradient_work))
        return dict(loss=float(np.sum(area*loss)),work_from_travel_gradient=continuous_part+jump,
                    smooth_element_work=continuous_part,interface_jump_work=float(jump),
                    fourth_axis_force=float(np.sum(area*drag/arc)),angular_integrated_drag=float(np.sum(area*drag)),
                    tangential_holder_x=float(np.sum(area*tx)),deep_contact_force=float(np.sum(area*shape*tx)))
    def at_cut(z,side):
        idx=np.clip(np.searchsorted(edges,z,side=side)-1,0,len(edges)-2)
        r,rp,arc,V,chi,tm,ta,mode=samples(np.asarray(z),idx)
        incoming=ANGULAR*r**2*V*(tm**2+2*ta**2)/(2*chi) if chi>0 else 0.
        return dict(z=float(z),power_in_inward_tag_stream=float(incoming),shear_x_momentum_toward_bulk=float(ANGULAR*r*mode))
    whole=integrate(edges[0],edges[last])
    begin=at_cut(edges[0],'right')
    # Beyond the active numerical tail the mean-contact traction vanishes.
    # Use a one-sided inlet value and positivity of any outgoing energy for
    # the shortfall test, avoiding a division by vanishing order at the far cut.
    for lo,hi in ((1.173,2.067),(1.31,1.87)):
        low,high=at_cut(lo,'right'),at_cut(hi,'left')
        wall=integrate(lo,hi)
        closure=high['shear_x_momentum_toward_bulk']-low['shear_x_momentum_toward_bulk']-wall['tangential_holder_x']
        emit(test='shear_shell_per_unit_deep_offset',inlet=a.inlet,low=low,high=high,wall=wall,closure=closure)
        if abs(closure)>1e-7:
            raise ArithmeticError('Integrated shear traction and cut flux disagree')
    contact_coefficient=ANGULAR*load.T@m
    archived_power=ANGULAR*m.T@contact@m/clock
    emit(test='settings',inlet=a.inlet,mesh=a.mesh,sub=a.sub,clock=clock,active_tail=float(edges[last]),
         direct_loss=whole['loss'],compliance_loss=archived_power,force_from_wall=whole['deep_contact_force'],
         force_from_virtual_displacement=contact_coefficient[1],base_carrier_force=carrier[:,1].tolist())
    for b in (0.,.01,-.01,.005,-.005):
        loss=b*b*whole['loss']
        incoming=b*b*begin['power_in_inward_tag_stream']
        emit(test='load',inlet=a.inlet,bend=b,renewal_power=loss,
             head_if_carrier_funds_all=loss/.25,positive_mean_spring_energy_in=incoming,
             minimum_other_supply_if_tagged_energy_only=loss-incoming,
             power_from_reference_gradient=b*b*whole['work_from_travel_gradient'],
             component_along_w_on_carrier=b*b*whole['fourth_axis_force'],
             contact_holder=(b*contact_coefficient).tolist(),
             total_linear_holder=(b*(contact_coefficient+carrier[:,1])).tolist())
    if (abs(whole['loss']/archived_power-1)>1e-5 or
            abs(whole['deep_contact_force']-contact_coefficient[1])>2e-7 or
            abs(whole['work_from_travel_gradient']/archived_power-1)>0.025):
        raise ArithmeticError('Direct contact quadrature, transport or force check failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(34,48),default=48)
    p.add_argument('--sub',type=int,default=3)
    p.add_argument('--quad',type=int,default=5)
    p.add_argument('--floor',type=float,default=1e-8)
    p.add_argument('--inlet',choices=('straight','relaxed'),default='straight')
    run(p.parse_args())
