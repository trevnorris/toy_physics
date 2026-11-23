"""Check a proposed radius-two oscillation against exact finite geometry and held-shape stored energy."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import legval

from mouth_exit_witness import inputs
from mouth_live_waist_fields import Carrier
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import BIAS, CUT, ROOT, Order, Trial, emit, exterior, source
from mouth_sealed_shape_relax import Meridian


def check_source(path, digest):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError('Input readout digest mismatch')
    guard = json.loads((path.parent/'outcome.json').read_text())
    if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
        raise ValueError('Unverified input')


def finite_shapes(shape, exit_scale, registration, direction):
    # These gradients act on the current material. All signed increments are
    # affine in the three benchmark coordinates, including physical axial depth.
    t = np.linspace(0, 1, 2049)
    b, slope, _ = shape.profile(t)
    dz = direction[0]
    dr = exit_scale*t*direction[1] + 4*shape.neck*t*(1-t)*direction[2]
    beta = dr/b
    derivative = (exit_scale*direction[1]+4*shape.neck*(1-2*t)*direction[2])/shape.length-beta*slope
    outer = np.exp(np.exp(legval(2*t-1,registration))/2)
    measures = {}
    for name, x in [('wall',1.), ('outermost_registered',outer)]:
        shear = x*derivative
        deformation = np.zeros((len(t),2,2))
        deformation[:,0,0] = 1+beta
        deformation[:,0,1] = shear
        deformation[:,1,1] = 1+dz
        stretches = np.linalg.svd(deformation,compute_uv=False)
        strain = .5*(deformation+np.swapaxes(deformation,1,2))-np.eye(2)
        eig = np.linalg.eigvalsh(strain)
        jac = (1+beta)**3*(1+dz)
        inv = np.min(1+beta)<=0 or 1+dz<=0
        measures[name]=dict(min_jacobian=jac.min(), max_jacobian=jac.max(),
                            left_positive_shape_chart=bool(inv),
                            min_principal_stretch=min(stretches.min(), (1+beta).min()),
                            max_principal_stretch=max(stretches.max(), (1+beta).max()),
                            peak_absolute_infinitesimal_strain=max(abs(eig).max(),abs(beta).max()),
                            peak_displacement=np.hypot(x*dr,t*shape.length*dz).max())
    return measures


def run(spec):
    cfg=json.loads(spec.read_text())
    check_source(ROOT/cfg['source'],cfg['source_sha256'])
    radius, mesh = 2, 72
    trial=Trial(radius,mesh,bonded_core=True)
    endpoint=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n{mesh}/endpoint.npz'
    shifts,chi,_,digest,status,*_=inputs(endpoint,radius,mesh,trial)
    if digest != cfg['material_sha256'] or status!='balanced':
        raise ArithmeticError('Different strained shape')
    matched=ROOT/'research/pde_ledger_v3/_measurements/mouth_breathing_exit_support.json'
    old=json.loads(matched.read_bytes())['runs']['mouth_breathing_matched01/r2n72']
    check_source(ROOT/'_scratch/s11c/mouth_breathing_matched01/r2n72/stdout',old['stdout_sha256'])
    rows=list(map(json.loads,old['literal_stdout'].splitlines()))
    origin=rows[0]
    central, = [r for r in rows if r.get('relative_exit')==0 and r['shape_step']==.001]
    step = np.asarray(cfg['coordinate_peak'],dtype=float)
    if step.shape!=(3,) or not np.isfinite(step).all() or not step.any():
        raise ValueError('Expected a finite real three-coordinate displacement')
    c0 = trial.x+shifts[:4]
    fixed = Meridian(radius,c0)
    exit_scale=origin['original_exit_scale']
    order_coeff=trial.co.copy(); order_coeff[0]+=shifts[4]
    receiver_coeff=trial.cr.copy(); receiver_coeff[0]+=shifts[5]
    z=np.linspace(0,1,2*mesh+1)
    r=np.linspace(0,1,mesh+1)
    K=np.array(central['hessian'])
    grad=np.array(central['generalized_gradient'])
    emit(kind='input', material_sha256=digest, spec_sha256=hashlib.sha256(spec.read_bytes()).hexdigest(),
         source=cfg['source'], source_sha256=cfg['source_sha256'], spec=cfg,
         hessian_sha256=old['stdout_sha256'], centre=origin['center_material_energy'],
         predicted_unit_quadratic=.5*step@K@step)
    for s in (0.,-.125,.125,-.5,.5,-1.,1.):
        v=step*s
        metrics=finite_shapes(fixed,exit_scale,receiver_coeff,v)
        if metrics['wall']['left_positive_shape_chart']:
            emit(kind='chart_crossed',stroke_fraction=s,geometry=metrics)
            continue
        config=c0.copy();config[0]+=np.log1p(v[0]);config[1]+=np.log1p(exit_scale*v[1]/fixed.exit);config[2]+=v[2]
        shape=Meridian(radius,config)
        if shape.wall()['minimum_radius']<=0:
            raise ArithmeticError('Meridional shape crosses the axis')
        inner=source(shape,order_coeff,z,r,False)
        inner+=core_source(shape,order_coeff,*trial.core_reference,z,r)
        ext=exterior(shape,mesh,CUT,side_order_work=lambda w,x:source(shape,receiver_coeff,w/shape.length,x,True))
        carrier=Carrier(shape,mesh,2*mesh,trial.mesh)
        eflow,_,info=carrier.response(chi,BIAS)
        en=ext['energy']+ext['static_work']+np.sum(inner*chi)+Order(shape,mesh).energy_gradient(chi)[0]+eflow
        if abs(info['current']-info['terminal_current'])>1e-6*abs(info['current']):
            raise ArithmeticError('Stationary perturbation did not conserve current')
        if s==0 and abs(en-origin['center_material_energy'])>1e-8:
            raise ArithmeticError('Unmatched reference energy')
        work = en-origin['center_material_energy']-grad@v
        quad=.5*v@K@v
        emit(kind='stroke',stroke_fraction=s,geometry=metrics,
             harmonic_work=quad,finite_material_work=work,
             finite_over_harmonic=work/quad if s else None,
             intake=info['current'],terminal=info['terminal_current'],
             receiver_residual=ext['residual'],static_holding_work=-grad@v)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('spec',type=Path)
    run(parser.parse_args().spec)
