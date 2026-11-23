"""Independently integrate collar work from physical and reference Jacobian matrices."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import legder, legval
from scipy.interpolate import RegularGridInterpolator

from mouth_joint_collar_release import Release, core_source, Order


def check_input(path):
    parent = path.parent
    guard = json.loads((parent / "outcome.json").read_bytes())
    limits = json.loads((parent / "limit-validation.json").read_bytes())
    if (guard["exitCode"] or guard["childOutcome"]["exitCode"] or not guard["limitsVerified"]
            or not limits["verified"] or int(limits["actual"]["memory.swap.max"])):
        raise ValueError(f"Invalid guarded input {path}")


def material_work(pair, L, inner, outer):
    # Integrate in volume pulled back to the original *material* chart.
    # Resolve the two affine label derivatives separately on each triangle.
    nodes = L.reshape(-1, 2)[pair.id]
    e, t = pair.eta, pair.z
    cot = legval(2*t-1, pair.co)
    p = np.exp(cot)
    b, dbdw, _ = pair.shape.profile(t)
    radius = b*np.exp(-e*p)
    mat = np.stack([np.einsum('ki,kij->kj',f,nodes) for f in (pair.E,pair.T)], axis=-1)
    physical = np.zeros_like(mat)
    physical[:,0,0] = -p*radius
    physical[:,0,1] = radius*(dbdw*pair.shape.length/b - e*p*2*legval(2*t-1,legder(pair.co)))
    physical[:,1,1] = pair.shape.length
    S = pair.shape.neck+(pair.B*nodes[:,:,0]).sum(axis=1)
    # Increasing physical w corresponds to decreasing reference Z.
    material_det = np.linalg.det(mat)
    mat[:,1,:] *= -1
    F = physical @ np.linalg.inv(mat)
    h = radius/S
    det = np.linalg.det(F)*h*h
    if min(det) <= 0 or min(material_det)<=0:
        raise ValueError("Inverted material state")
    W = .5*((F*F).sum(axis=(1,2))+2*h*h-4)-np.log(det)+.25*(det-1)**2
    energy=[]
    for side,grid in [(e>0,inner),(e<0,outer)]:
        field = RegularGridInterpolator(grid[:2],grid[2],bounds_error=True)
        vals = field(np.column_stack((t[side],pair.coordinate[side])))
        energy.append(float((vals*W[side]*S[side]**2*material_det[side])@pair.qw[side]))
    return energy


def run(args):
    check_input(args.field)
    options=SimpleNamespace(state=args.state,order=args.order,flow=args.flow,polar=args.polar,
                            na=args.na,nt=args.nt)
    release=Release(options)
    c,v0=release.initial_chi,release.start
    with np.load(args.field) as z:
        if (int(z['radius'])!=4 or int(z['winding'])!=1 or z['mesh']!=args.flow
                or z['polar']!=args.polar):
            raise ValueError('Mismatched input for energy check')
        if (not np.allclose(z['shape_params'],release.model(v0).params) or
                z['shape_extra']!=release.model(v0).extra):
            raise ValueError('Energy fields and material prepared on different geometry')
        amp,phase=z['amplitude_ratio'],z['phase']
    center=release.model(v0)
    y=release.initial_labels
    if y is None:
        y=release.pair(center).identity
    a=np.array([.31,-.27,.17,-.11,.81])
    center_pair=release.pair(center)
    src=center_pair.evaluate(y,release.inner(c),release.outer(np.zeros_like(c)))['outer_source']
    outer0=release.outside(center,src,True)[1][2]
    offsets={}
    for step in [-.0004,-.0002,0.,.0002,.0004]:
        v=v0+step*a
        sh=release.model(v); pair=release.pair(sh)
        native=release.fixed_fields(v,c,y,amp,phase)
        src=pair.evaluate(y,release.inner(c),release.outer(outer0))['outer_source']
        ext,(_,_,outer)=release.outside(sh,src,True)
        parts=material_work(pair,y,release.inner(c),release.outer(outer))
        native_work=pair.evaluate(y,release.inner(c),release.outer(outer))['energy']
        core=(core_source(sh,release.co,release.ref,release.co,release.t,release.x)*c).sum()
        total=(sum(parts)+core+Order(sh,args.order).energy_gradient(c)[0]
               +ext['energy']+release.wound(sh,c).amplitude(amp,phase)[0])
        if abs(total-native)>2e-10*(1+abs(native)):
            raise ValueError('Independent partition mismatch')
        offsets[step]=[float(total),float(native)]
        print(json.dumps(dict(kind='partition',step=step,inner=parts[0],outer=parts[1],
                              native_collar=native_work,native_total=native,
                              independently_partitioned=total,receiver_order_residual=ext['residual'])),flush=True)
    print(json.dumps(dict(kind='energy_work',direction=a.tolist(),half=((np.array(offsets[-.0002])-offsets[.0002])/.0004).tolist(),
                          full=((np.array(offsets[-.0004])-offsets[.0004])/.0008).tolist(),radius=4,
                          input_sha256=[hashlib.sha256(x.read_bytes()).hexdigest() for x in [args.state,args.field]])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state',type=Path)
    parser.add_argument('field',type=Path)
    for k,v in [('order',16),('flow',13),('polar',10),('na',4),('nt',8)]:
        parser.add_argument('--'+k,type=int,default=v)
    run(parser.parse_args())
