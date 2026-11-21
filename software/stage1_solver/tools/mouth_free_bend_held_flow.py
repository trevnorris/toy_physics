"""Steady holding force on a fixed entrance with two freely adjusting bends.

Use virtual work in the carrier volume and independently sample reflecting-wall
pressure. The companion-like drive is one weak uniform brane stream at infinity
represented at the finite bath edge; the receiving bath has zero stream.
"""
import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_inertia import taper
from mouth_carrier_lateral_response import prepare, trace
from mouth_empty_bend_stiffness import INTERVALS, quadrature, shapes
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response


def wall_reaction(dom, field, perturbation):
    result=np.zeros(3)
    q, weights=leggauss(10)
    for block in (0,2):
        nodes=dom.blocks[block][0]
        half=np.diff(nodes)[:,None]/2
        p=(nodes[:-1,None]+half*(1+q)).ravel()
        r,w,_,_,wp,_,f,fr,fw=trace(dom,block,p,np.zeros_like(p),
                                      np.column_stack([field,perturbation]))
        shape=np.column_stack([np.ones_like(w)]+[1-taper(w,*span)[0] for span in INTERVALS])
        variation=(fr[:,0].conj()*fr[:,1]+fw[:,0].conj()*fw[:,1]
                   +(abs(f[:,0])**8-1)*f[:,0].conj()*f[:,1]/2).real
        result-=shape.T@((4*np.pi/3)*r*r*wp*(half*weights).ravel()*variation)
    return result


def volume_reaction(dom, field, perturbation, extension):
    r,w,gr,gw=quadrature(dom)
    _,sr,sw=shapes(r,w,extension,dom.far,dom.roof)
    radial,dr=taper(r,2.0 if extension=='compact' else 3.4,dom.far-.35)
    depth,dw=taper(w,3.3 if extension=='compact' else 4.3,dom.roof-.35)
    sr=np.column_stack((dr*depth,sr))
    sw=np.column_stack((radial*dw,sw))
    f, fr, fw=dom.N@field,gr@(field-1),gw@(field-1)
    h, hr, hw=dom.N@perturbation,gr@perturbation,gw@perturbation
    lag=(fr.conj()*hr+fw.conj()*hw+(abs(f)**8-1)*f.conj()*h/2).real
    dx=hr+2*h/r
    t_r=(dx.conj()*fr+fr.conj()*hr).real-lag
    t_w=(dx.conj()*fw+fr.conj()*hw).real
    return (sr.T@(dom.weight*t_r)+sw.T@(dom.weight*t_w))/3


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    args=p.parse_args()
    dom, field=prepare(args)
    response=Response(dom,field,args.extension)
    static=response.force_matrix(0.)
    for label, feed, bulk in [('brane',1.,0.),('bulk',0.,1.),('tow',-1.,-1.)]:
        perturbation, residual=response.dipole.steady(feed,bulk)
        volume=volume_reaction(dom,field,perturbation,args.extension)
        wall=wall_reaction(dom,field,perturbation)
        emit(test='incident_load',drive=label,volume=volume.tolist(),wall=wall.tolist(),
             field_residual=residual)
        if np.max(abs(volume-wall))>.005:
            raise ArithmeticError('Volume force and wall pressure disagree')
        for tension in (.02,.04):
            k=static+tension*response.area
            bend=np.linalg.solve(k[1:,1:],volume[1:])
            net=volume-k@np.r_[0.,bend]
            emit(test='relaxed_held',drive=label,tension=tension,force_per_speed=net[0],
                 holding_reaction_per_speed=-net[0],bend_per_speed=bend.tolist(),
                 internal_residual=max(abs(net[1:])))
            if max(abs(net[1:]))>1e-9 or not np.isfinite(net).all():
                raise ArithmeticError('Unbalanced free bend in the held-flow test')
