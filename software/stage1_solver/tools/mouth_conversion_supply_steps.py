"""Symmetric full/half supply perturbations at fixed carrier intake on a settled branch."""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq

from mouth_carrier_lateral_response import prepare
from mouth_conversion_supply import SuppliedOrder
from mouth_empty_bend_stiffness import quadrature
from mouth_free_order_diffuse import read, section_geometry
from mouth_supplied_contact_work import verified
from mouth_order_wave_pulse import digest,emit


def main(a):
    records=verified(a.source)
    setup=next(r for r in records if r['test']=='calibration')
    center=next(r for r in records if r.get('multiple')==a.center)
    dom,psi=prepare(SimpleNamespace(source=Path(records[0]['source']).parent))
    emit(**setup)
    if digest(Path(center['path']))!=center['sha256']:
        raise ValueError('Starting state changed')
    with np.load(center['path']) as data:
        seed=data['chi']
    model=SuppliedOrder(dom,psi,setup['width'],setup['mobility'])
    _,_,dr,dw=quadrature(dom)
    f=dom.N@psi
    ur,uw=(np.imag(g@(psi-1)/f) for g in (dr,dw))
    for factor in (-2,-1,0,1,2):
        delta=factor*a.step
        multiple=round(a.center+delta,10)
        h=setup['flat_head']*multiple
        if delta:
            result=model.equilibrium(seed,h)
            if result is None:
                raise ArithmeticError('Modulated layer did not settle')
            order=result[0]
            path=a.out/f'order_{multiple:.5f}.npz'
            np.savez_compressed(path,chi=order)
        else:
            order,path=seed,Path(center['path'])
        c=(dom.N@order)/model.quad_sum
        cr,cw=dr@order,dw@order
        material_rate=ur*cr+uw*cw
        fuel_per_chi=h*model.density*6*c*(1-c)
        mu_bare=-material_rate/setup['mobility']-fuel_per_chi
        supply=-fuel_per_chi*material_rate
        heat=material_rate**2/setup['mobility']
        axial_force=mu_bare*cw
        material_power=mu_bare*material_rate
        wt=dom.weight
        volume_supply=wt@supply
        depth=brentq(lambda z:read(dom,order,section_geometry(dom,z)[0],z)[0]-.5,.52,2.)
        emit(test='supplied_order',head=h,multiple=multiple,path=str(path),sha256=digest(path),
             head_power=h*setup['intake'],volume_supply=volume_supply,passive_loss=wt@heat,
             wall_half_order_depth=depth,bare_force_on_carrier_w=wt@axial_force,
             bare_force_work_on_carrier=wt@material_power,
             positive_local_carrier_work=wt@np.maximum(material_power,0),
             negative_local_carrier_work=wt@np.minimum(material_power,0))
        if abs(wt@(heat+material_power-supply))>5e-12 or abs(volume_supply/(h*setup['intake'])-1)>.001:
            raise ArithmeticError('Fuel, heat and carrier-work balance failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--center',type=float,default=.5)
    p.add_argument('--step',type=float,default=.01)
    main(p.parse_args())
