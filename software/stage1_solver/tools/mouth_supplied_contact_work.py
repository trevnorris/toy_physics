"""Contact-energy and reaction measurements on independently supplied straight order profiles.

The first nonzero power is quadratic in held displacement. Its coefficient uses
the axisymmetric driven profile and the elastic dipole at that supply. The
chemical-to-spring term is signed work on the spring along a material path.
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
from mouth_empty_bend_stiffness import INTERVALS
from mouth_order_wave_pulse import digest, emit


def verified(folder):
    outcome=json.loads((folder/'outcome.json').read_text())
    limits=json.loads((folder/'limit-validation.json').read_text())
    if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not outcome['limitsVerified']
            or not limits['verified'] or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Invalid guarded source: {folder}')
    return [json.loads(line) for line in (folder/'stdout').read_text().splitlines()]


def contact_readout(geo, order, psi, modes, keep, edges, clock):
    last=np.flatnonzero(keep)[-1]+1
    width=np.diff(edges)[:last,None]
    pos=edges[:last,None]
    pg,vg,cells=geo.order_blocks[2]
    parent=np.searchsorted(pg,(pos+width/2).ravel())-1
    nodes=cells[parent*(len(vg)-1)][:,[0,3,6]]

    def wall(data,z):
        h=np.diff(pg)[parent,None]
        t=(z-pg[parent,None])/h
        a,b,c=(v[:,None] for v in data[nodes].T)
        quad=2*(a+c-2*b)
        slope=4*b-3*a-c
        return a+t*(slope+t*quad),(slope+2*t*quad)/h

    dual=np.zeros(len(edges))
    dual[keep]=modes
    left=dual[:last,None]
    slope=np.diff(dual)[:last,None]/width

    def sample(z):
        chi,dchi=wall(order,z)
        field,dfield=wall(psi,z)
        radius,rp=geo.radius(z)
        arc=np.hypot(1.,rp)
        speed=(dfield/field).imag/arc
        mode=left+slope*(z-pos)
        tm=mode/radius**2
        ta=slope/(2*radius*arc)
        if np.any((chi<=0)&((tm!=0)|(ta!=0))):
            raise ValueError('Traction on unordered material')
        ext2=np.divide(tm**2+2*ta**2,chi**2,out=np.zeros_like(tm),where=chi>0)
        area=ANGULAR*radius**2*arc
        volume_flow=ANGULAR*radius**2*speed
        spring=chi*ext2/2
        return (area,volume_flow,spring,ext2,dchi,tm*rp/arc+2*ta)

    roots,wt=leggauss(10)
    z=pos+width*(1+roots)/2
    w=width*wt/2
    area,vol,spring,extension2,dchi,tx=sample(z)
    local=w*vol*extension2*dchi/2
    chemical_work=local.sum()
    renewal=np.sum(w*area*spring)*2/clock
    # Integrate the spatial change of the material surface flux, including
    # jumps in one-sided FE wall speed or strain at element boundaries.
    eps=width*1e-5
    dflow=(sample(z+eps)[1]-sample(z-eps)[1])/(2*eps)
    dil=np.sum(w*spring*dflow)
    _,fin,ein,*_=sample(pos)
    _,fout,eout,*_=sample(pos+width)
    dil+=np.sum((ein[1:]+eout[:-1])*(fin[1:]-fout[:-1])/2)
    incoming=float((fin*ein)[0,0])
    outgoing=float((fout*eout)[-1,0])
    slip=outgoing-incoming-dil+renewal-chemical_work
    hold=np.sum(w*area*tx*(1-taper(z,*INTERVALS[1])[0]))
    # Surface stress on the material = +t and on the wall = -t.
    # Integrating its x component on any segment gives ANGULAR*[a*m].
    def endpoint(z):
        return geo.radius(np.array([z]))[0][0]*np.interp(z,edges,dual)
    segment=ANGULAR*(endpoint(2.067)-endpoint(1.173))
    return dict(order_work_on_contacts=chemical_work,order_work_positive=np.maximum(local,0).sum(),
                order_work_negative=np.minimum(local,0).sum(),renewal=renewal,
                slip_work_from_material_balance=slip,area_flux_divergence=dil,
                incoming_spring_work=incoming,outgoing_spring_work=outgoing,
                generalized_contact_holder=hold,all_guide_net_x=np.sum(w*area*tx),
                segment_force_on_material_x=segment,segment_force_on_wall_x=-segment)


def main(a):
    records=verified(a.supply)
    setup=next(r for r in records if r['test']=='calibration')
    mesh=setup['mesh']
    geo,_,carrier=archive_inputs(Path(setup['order_source']),ROOT/f'bend_static_shear_delcarrier{mesh}_01')
    carrier_file=ROOT/f'mouth_flared_delayed_n{mesh}/field.npz'
    if digest(carrier_file)!=records[0]['sha256']:
        raise ValueError('Mismatched carrier')
    with np.load(carrier_file) as data:
        psi=data['psi']
    for record in records:
        if record['test']!='supplied_order' or (a.levels and record['multiple'] not in a.levels):
            continue
        p=Path(record['path'])
        if digest(p)!=record['sha256']:
            raise ValueError('Supplied field changed')
        with np.load(p) as data:
            order=data['chi']
        edges,keep,bulk=volume_cost(geo,order,a.sub,a.quad,1e-9)
        contact,load,B,G,drift,clock,_=wall_data(geo,order,psi,edges,keep,a.quad)
        matrix=bulk+contact+B@solve(clock*drift,G)
        norm=1/np.sqrt((bulk+contact).diagonal())
        entry=np.zeros(len(drift))
        entry[len(entry)//2:]=1.
        raw,correction=(norm[:,None]*solve(norm[:,None]*matrix*norm,
                             norm[:,None]*np.column_stack([load[:,1],B@entry]))).T
        for name,m in [('straight',raw),('traction_free',raw-correction*raw[0]/correction[0])]:
            r=contact_readout(geo,order,psi,m,keep,edges,clock)
            np.testing.assert_allclose(r['renewal'],ANGULAR*(m@contact@m)/clock,rtol=3e-5)
            np.testing.assert_allclose(r['generalized_contact_holder'],ANGULAR*(load[:,1]@m),atol=3e-8)
            if abs(r['all_guide_net_x'])>1e-7:
                raise ArithmeticError('Unweighted wall traction does not balance')
            emit(test='contact',mesh=mesh,sub=a.sub,multiple=record['multiple'],head=record['head'],inlet=name,
                 head_power=record['head_power'],passive_loss=record['passive_loss'],
                 signed_net_deep_holder_per_bend=r['generalized_contact_holder']+carrier[1,1],
                 convention='All powers per squared bend, all forces per bend, except head_power and passive_loss at zero bend',
                 **r)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--supply',type=Path,required=True)
    p.add_argument('--sub',type=int,default=3)
    p.add_argument('--quad',type=int,default=4)
    p.add_argument('--levels',nargs='+',type=float)
    main(p.parse_args())
