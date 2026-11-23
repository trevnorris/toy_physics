"""Low-frequency transfer of a prescribed weak breathing displacement into one interior order coordinate."""
import argparse
import hashlib
import json

import numpy as np
from scipy.special import expit, logit

from mouth_exit_release import Release
from mouth_exit_witness import inputs
from mouth_live_waist_fields import Carrier
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import BIAS, ROOT, Order, Trial, emit, source
from mouth_sealed_shape_relax import Meridian


def run(radius, mesh, directional):
    folder=ROOT/f'_scratch/s11c/mouth_breathing_matched01/r{radius}n{mesh}'
    if directional:
        folder=folder.with_name(folder.name+'_direction')
    guard=json.loads((folder/'outcome.json').read_text())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError('Invalid empty-background candidate')
    raw=(folder/'stdout').read_bytes()
    rows=list(map(json.loads,raw.splitlines()))
    refs,=[r for r in rows if r['kind']=='input']
    center,=[r for r in rows if r['kind']=='spectrum' and r['relative_exit']==0 and (r['shape_step']==.001 or directional)]
    waves=[r for r in rows if r['kind']=='signed' and r['exit_step']==.006 and r['outward_force_per_action']>0]
    trial=Trial(radius,mesh,bonded_core=True)
    basefile=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n{mesh}/endpoint.npz'
    sh,chi,_,digest,status,*_=inputs(basefile,radius,mesh,trial)
    if digest!=refs['material_sha256'] or status!='balanced':
        raise ValueError('Source preparation mismatch')
    emit(kind='input',radius=radius,mesh=mesh,input_sha256=digest,
         matched_sha256=hashlib.sha256(raw).hexdigest(),outward_candidates=len(waves))
    if not waves:
        return
    family=Release(radius,mesh)
    active=family.active
    logchi=logit(chi[active])
    chi_tangent=chi[active]*(1-chi[active])*family.rate
    depth=np.linspace(0,1,2*mesh+1)
    radial=np.linspace(0,1,mesh+1)
    inner=trial.co.copy()
    inner[0]+=sh[4]
    base=Meridian(radius,trial.x+sh[:4])
    b=base.exit
    rel=refs['original_exit_scale']/b

    def slope(config,order_coord):
        location=trial.x+sh[:4]
        location[0]+=np.log1p(config[0])
        location[1]+=np.log1p(rel*config[1])
        location[2]+=config[2]
        shape=Meridian(radius,location)
        c=chi.copy()
        c[active]=expit(logchi+order_coord*family.rate)
        tangent=c[active]*(1-c[active])*family.rate
        elastic=source(shape,inner,depth,radial,False)
        elastic+=core_source(shape,inner,*trial.core_reference,depth,radial)
        carrier=Carrier(shape,mesh,2*mesh,trial.mesh)
        carrier_grad=carrier.response(c,BIAS)[1]
        order_grad=Order(shape,mesh).energy_gradient(c)[1]
        return (order_grad+elastic+carrier_grad)[active] @ tangent

    zero=np.zeros(3)
    for wave in waves:
        idx=wave['index']
        mode=np.array(center['mass_normalized_vectors'][idx])
        amp=np.sqrt(2*wave['conditional_action']/wave['omega'])
        for h in (.0004,.0002):
            order_stiffness=(slope(zero,h)-slope(zero,-h))/(2*h)
            cross=(slope(h*mode,0.)-slope(-h*mode,0.))/(2*h)
            if order_stiffness<=0:
                raise ArithmeticError('The tested order direction is unstable')
            transfer=-cross/order_stiffness
            apparent=transfer*chi_tangent
            site=np.argmax(abs(apparent))
            ij=np.column_stack(np.nonzero(active))[site]
            emit(kind='quasistatic_order',index=idx,step=h,baseline_generalized_force=-slope(zero,0.),
                 order_stiffness=order_stiffness,mixed_stiffness=cross,
                 order_coordinate_per_unit_modal_displacement=transfer,
                 peak_material_order_per_unit_modal_displacement=abs(apparent[site]),
                 linear_extrapolation_to_conditional_action=amp*abs(apparent[site]),
                 location=ij.tolist(),local_chi=chi[tuple(ij)])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('radius',type=int,choices=(1,2))
    parser.add_argument('mesh',type=int,choices=(48,72))
    parser.add_argument('--directional',action='store_true')
    args=parser.parse_args()
    run(args.radius,args.mesh,args.directional)
