"""Fixed-depth, fixed-bias radius family with regional resistance and inertia checks.

Exterior mass regions in the zero-intake translation test are fixed by the
mesh blocks. At finite bias, phase-gradient energy divided by squared current
measures resistance in each region. The last identity is exact for stationary
continuity; measure the finite-element and quadrature error independently.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import spsolve

from mouth_carrier_lateral_inertia import volume
from mouth_carrier_lateral_response import Dipole, trace
from mouth_carrier_matched_bias import at_bias
from mouth_carrier_weak_radius import PhysicalDomain, translation
from mouth_flared_receiver import EXIT, emit

BIAS = .4510084282857874
WINDOWS = ((1.65, 2.55, EXIT+.4, EXIT+2.5), (1.9,3.25,EXIT+.7,EXIT+3.4))
REGIONS = ('feeding_slab_and_rounded_lip', 'conduit', 'receiver')


def regional_energy(dom, field, translation_field, rate):
    result=np.zeros((3,2))
    points, weights=leggauss(4)
    u,v=[d.ravel() for d in np.meshgrid(points,points)]
    fields=np.column_stack([field, translation_field])
    for block,(pg,vg,dofs) in dom.blocks.items():
        region = 0 if block < 2 else 1 if block==2 else 2
        for start in range(0,len(dofs),256):
            i,j=np.unravel_index(np.arange(start,min(start+256,len(dofs))),
                                 (len(pg)-1,len(vg)-1))
            hp,hv=np.diff(pg)[i,None]/2,np.diff(vg)[j,None]/2
            p,t=pg[i,None]+hp*(1+u),vg[j,None]+hv*(1+v)
            r,_,rp,rv,wp,wv,f,fr,fw=trace(dom,block,p,t,fields)
            weight=4*np.pi*r*r*hp*hv*abs(rp*wv-rv*wp)*np.outer(weights,weights).ravel()
            psi=f[...,0]
            phase_gradient2=(fr[...,0]/psi).imag**2+(fw[...,0]/psi).imag**2
            cost=abs(psi)**2*phase_gradient2/rate**2
            kinetic=(fr[...,1].real**2+fw[...,1].real**2+2*(f[...,1].real/r)**2)/3
            result[region] += (np.sum(weight*cost),np.sum(weight*kinetic))
    return result


class RefinedDipole(Dipole):
    """Use the same double-precision FE operator and solver, with a wider residual accumulator."""

    def solve(self, boundary, source=None):
        matrix = self.operator.astype(np.longdouble)
        values = np.r_[boundary.real, boundary.imag]
        rhs = np.zeros_like(values)
        if source is not None:
            rhs[:] = np.r_[source.real, source.imag]
        errors = []
        for step in range(4):
            residual = (matrix @ values.astype(np.longdouble) - rhs)[self.free]
            errors.append(float(np.max(np.abs(residual))))
            if step and errors[-1] < 3e-8:
                emit(test='refined_dipole_residual', corrections=step, errors=errors,
                     checked_vector_type=str(values.dtype))
                return values[:self.count]+1j*values[self.count:], errors[-1]
            correction = self.factor.solve(np.asarray(residual, dtype=float))
            values[self.free] -= correction
        raise ArithmeticError(f'Extended accumulation of physical equation residual: {errors}')


def steady_dipoles(dom,field,refine=False):
    model=(RefinedDipole if refine else Dipole)(dom,field)
    feed,_=model.steady(1.,0.)
    bulk,_=model.steady(0.,1.)
    tow,_=model.steady(-1.,-1.)
    source=dom.N.T@(dom.weight*(1j*(dom.N@tow)-model.rq*(dom.N@field)))
    acceleration,_=model.solve(np.zeros_like(field),source)
    return feed,bulk,tow,acceleration


def run(a):
    dom=PhysicalDomain(a.radius,1.,a.mesh,a.outer,a.roof,1.)
    phase=dom.bulk.astype(float)
    free=dom.free
    phase[free]=spsolve(dom.stiff[free][:,free],-(dom.stiff@phase)[free])
    z,rate,_=at_bias(dom,phase,BIAS,None)
    _,_,minus=at_bias(dom,phase,BIAS-.0005,z)
    _,_,plus=at_bias(dom,phase,BIAS+.0005,z)
    derivative=(plus-minus)/.001
    emit(test='pivot_preparation',radius=a.radius,mesh=a.mesh,outer=a.outer,roof=a.roof,
         axial_depth=EXIT,thickness=1.,bias=BIAS,material_intake=rate,
         action_derivative=derivative,nodes=len(z))
    if abs(derivative-rate)>2e-6*rate:
        raise ArithmeticError('Phase derivative of the stationary action and current differ')
    field=1+z
    zero_mass,work,phi=translation(dom,return_field=True)
    regions=regional_energy(dom,field,phi,rate)
    emit(test='regional_integrals',region_names=REGIONS,
         finite_flow_resistance=regions[:,0].tolist(),fixed_phase_drop_over_current=BIAS/rate,
         zero_intake_added_mass=regions[:,1].tolist(),zero_intake_energy=zero_mass,
         zero_intake_wall_work=work)
    if abs(regions[:,0].sum()*rate/BIAS-1)>5e-5:
        raise ArithmeticError('Phase energy and stationary throughflow fail regional balance')
    np.testing.assert_allclose([zero_mass,work],regions[:,1].sum(),rtol=3e-7)
    feed,bulk,tow,acceleration=steady_dipoles(dom,field,a.refine_dipole)
    for left,right,z0,z1 in WINDOWS:
        w=(left*a.radius,right*a.radius,z0,z1)
        if right*a.radius>=a.outer or z1>=a.roof:
            raise ValueError('Test window intersects driven bath')
        C,B,T=[volume(dom,np.c_[field,u],w,False)['virtual_wall_force'] for u in (feed,bulk,tow)]
        M=-volume(dom,np.c_[field,acceleration,tow],w,True)['virtual_wall_force']
        emit(test='finite_radius',window=w,radius=a.radius,radius_over_depth=a.radius/EXIT,
             intake=rate,C=C,B=B,tow=T,M=M,normalized_rate=rate/a.radius**3,
             normalized_mass=M/a.radius**3,C_over_D=C/rate,C_over_M=C/M,
             D_over_M=rate/M,weak_M_over_finite_M=zero_mass/M)
        if abs(C+B+T)>5e-8 or M<=0:
            raise ArithmeticError('Independent forcing columns or measured inertia failed')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--outer',type=float,required=True)
    p.add_argument('--roof',type=float,required=True)
    p.add_argument('--refine-dipole', action='store_true')
    run(p.parse_args())
