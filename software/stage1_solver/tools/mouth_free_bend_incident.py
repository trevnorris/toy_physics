"""Drive the brane reservoir and solve for the free entrance and two bend coordinates."""
import argparse
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_inertia import taper
from mouth_carrier_lateral_response import prepare, weighted_gram
from mouth_empty_bend_stiffness import quadrature, shapes
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response, real_pair
from mouth_free_bend_held_flow import volume_reaction
from mouth_free_bend_release import pencil


def boundary_momentum(dom, field, h, extension):
    r,w,gr,_=quadrature(dom)
    shape,sr,_=shapes(r,w,extension,dom.far,dom.roof)
    rg,dr=taper(r,2.0 if extension=='compact' else 3.4,dom.far-.35)
    zg,_=taper(w,3.3 if extension=='compact' else 4.3,dom.roof-.35)
    shape=np.column_stack((rg*zg,shape))
    sr=np.column_stack((dr*zg,sr))
    f,fr=dom.N@field,gr@(field-1)
    u,ur=dom.N@h,gr@h
    current=(u.conj()*fr+f.conj()*(ur+2*u/r)).imag
    return (shape.T@(dom.weight*current)+sr.T@(dom.weight*(f.conj()*u).imag))/3


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    args=p.parse_args()
    dom, field=prepare(args)
    response=Response(dom,field,args.extension)
    model=response.dipole
    mass=weighted_gram(dom.N,dom.weight)
    full_skew=sparse.bmat([[None,-mass],[mass,None]],format='csr')
    h=1j*dom.points[:,0]*field*dom.feed
    boundary=real_pair(h)
    b=volume_reaction(dom,field,h,args.extension)
    momentum=boundary_momentum(dom,field,h,args.extension)
    bulk_rhs=-(2/3)*(model.operator@boundary)[model.free]
    time_rhs=(2/3)*(full_skew@boundary)[model.free]
    for tension in (.02,.04):
        energy,time=pencil(response,tension)
        for omega in (0.,.02,.01,.005):
            rate=1j*omega
            rhs=np.r_[bulk_rhs+rate*time_rhs,b-rate*momentum]
            matrix=energy-rate*time
            if omega==0.:
                active=np.ones(len(rhs),bool)
                active[-3]=False  # entrance held; both bend forces vanish
                y=np.zeros(len(rhs),complex)
                y[active]=splu(matrix[active][:,active].tocsc()).solve(rhs[active])
                reaction=(matrix@y-rhs)[-3]
                emit(test='dc_holder',tension=tension,force_per_stream=float(reaction.real))
            else:
                y=splu(matrix.tocsc()).solve(rhs)
                error=float(np.max(abs(matrix@y-rhs)))
                accel=-omega**2*np.array([y[-3],sum(y[-3:])])
                emit(test='free_incident',tension=tension,omega=omega,
                     entrance_and_flare_real=accel.real.tolist(),
                     entrance_and_flare_imag=accel.imag.tolist(),residual=error)
                if not np.isfinite(accel).all() or error>3e-8:
                    raise ArithmeticError('Free field or coordinate balance failed')
