"""Isolate steady vortex work by turning off the axial reservoir drive on the same held mouth."""
import argparse
import numpy as np

from mouth_poly_shallow_release import emit
from mouth_sealed_shape_relax import Meridian
from mouth_wound_carrier import Wound, prepared


def run(radius,count,n,nt):
    config,order,digest=prepared(radius)
    b=Meridian(radius,config).exit
    emit(kind='input',radius=radius,winding=count,circulation=2*np.pi*count,mesh=n,polar=nt,
         reservoir_drop=0.,exit=b,material_sha256=digest)
    for h in (.006,.012):
        energies=[]
        for delta in (-h,h):
            shifted=config.copy()
            shifted[1]+=np.log1p(delta)
            shape=Meridian(radius,shifted)
            model=Wound(shape,order,count,n,nt,28.,22.)
            _,_,out=model.solve(0.)
            if abs(out['total_current'])>2e-12:
                raise ArithmeticError('Quiet axial drive carried net drain current')
            emit(kind='held',relative_exit=delta,excess_energy=out['excess_energy'],
                 meridional_intake=out['total_current'],work_terms=out['terms'],residual=out['residual'])
            energies.append(out['excess_energy'])
        force=(energies[0]-energies[1])/(2*h*b)
        emit(kind='vortex_only_force',relative_step=h,force=force)
        if count==0 and abs(force)>1e-9:
            raise ArithmeticError('Quiet unwound material generated a force')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('radius',choices=(1,2),type=int)
    p.add_argument('count',choices=(0,1,2),type=int)
    p.add_argument('mesh',type=int)
    p.add_argument('polar',type=int)
    args=p.parse_args()
    run(args.radius,args.count,args.mesh,args.polar)
