"""Independent weak-intake test using hydraulic conductance and translation energy.

Two scalar variational problems, solved with physical length coordinates. A
small common axial phase bias drives the monopole. A transverse wall velocity
drives a laboratory dipolar velocity potential, with the baths at rest.
"""
import argparse
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import spsolve
from scipy.integrate import quad

from mouth_carrier_lateral_response import quadrature_radii, weighted_gram
from mouth_carrier_radius_family import RadiusFamily
from mouth_curved_phase import quadratic, tensor
from mouth_flared_receiver import EXIT, ROUNDING, FLARE_START, FlaredDomain, emit


class PhysicalDomain(RadiusFamily):
    def __init__(self,radius,depth,mesh,outer,roof,thickness):
        self.length_unit=radius
        self.depth_stretch=depth/radius
        self.feed_depth=thickness/depth
        FlaredDomain.__init__(self,mesh,'cubic',round(roof/depth,11),far=round(outer/radius,11))
        self.points=self.stretched_points
        del self.stretched_points
        r,w=self.points.T
        side=np.isclose(r,outer,atol=2e-9,rtol=0.)
        np.testing.assert_array_equal(self.feed, side & (w<=1e-10))
        np.testing.assert_array_equal(self.bulk,(side & (w>=EXIT*depth-1e-10)) |
                                      np.isclose(w,roof,atol=2e-9,rtol=0.))
        np.testing.assert_allclose(w.min(),-thickness,atol=2e-9)
    def assemble(self,grids):
        super().assemble(grids)
        # Only chart coordinates enter the parent's boundary classification.
        self.points=np.round(self.stretched_points/[self.length_unit,
                              self.length_unit*self.depth_stretch],11)


def conductance(dom):
    k=dom.stiff
    phase=dom.bulk.astype(float)
    mask=dom.free
    phase[mask]=spsolve(k[mask][:,mask],-(k@phase)[mask])
    reaction=k@phase
    return float(phase@reaction),float(reaction[dom.bulk].sum()),float(-reaction[dom.feed].sum())


def translation(dom):
    k=dom.stiff+weighted_gram(dom.N,2*dom.weight/quadrature_radii(dom)**2)
    # Lab potential has zero outer value and normal velocity n_x on the wall.
    load=np.zeros(len(dom.points))
    root,weights=leggauss(8)
    for block in (0,2):
        pg,vg,dofs=dom.blocks[block]
        half=np.diff(pg)[:,None]/2
        p=pg[:-1,None]+half*(1+root)
        r,_,_,_,wp,_=dom.chart(block,p,np.zeros_like(p))
        n=tensor(quadratic(root)[0],quadratic(np.full_like(root,-1.))[0])
        assembled=(4*np.pi*r*r*wp*half*weights)@n
        ids=dofs[np.arange(len(pg)-1)*(len(vg)-1)]
        np.add.at(load,ids.ravel(),assembled.ravel())
    free=dom.free & (dom.points[:,0]>0)
    phi=np.zeros_like(load)
    phi[free]=spsolve(k[free][:,free],load[free])
    return float(phi@(k@phi))/3, float(load@phi)/3


def neck_integrals(dom):
    def radius(w):
        if w<ROUNDING:
            return dom.rc-np.sqrt(ROUNDING**2-(w-ROUNDING)**2)
        return dom.radius(w)[0]
    def mean(power):
        return sum(quad(lambda w:radius(w)**power,lo,hi,epsabs=1e-11)[0]
                   for lo,hi in ((0.,ROUNDING),(ROUNDING,FLARE_START),(FLARE_START,EXIT)))
    return mean(-3),mean(3)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=float,required=True)
    p.add_argument('--depth',type=float,required=True)
    p.add_argument('--mesh',type=int,required=True)
    p.add_argument('--outer',type=float,required=True)
    p.add_argument('--roof',type=float,required=True)
    p.add_argument('--reference-bias',type=float,default=.4510084282857874)
    a=p.parse_args()
    dom=PhysicalDomain(a.radius,a.depth,a.mesh,a.outer,a.roof,1.)
    c,flux,feed=conductance(dom)
    mass,work=translation(dom)
    inv,vol=neck_integrals(dom)
    conduit_c=(4*np.pi/3)*a.radius**3/(a.depth*inv)
    conduit_mass=(4*np.pi/3)*a.radius**3*a.depth*vol
    emit(test='weak_geometry',radius=a.radius,depth=a.depth,mesh=a.mesh,outer=a.outer,
         roof=a.roof,slab_thickness=1.,conductance_energy=c,bulk_flux=flux,feed_flux=feed,
         mass_from_energy=mass,mass_from_wall_work=work,conductance_over_mass=c/mass,
         weak_ratio_at_reference_bias=a.reference_bias*c/mass,
         slender_conductance=conduit_c,slender_mass=conduit_mass,slender_ratio=conduit_c/conduit_mass)
    np.testing.assert_allclose([c,c,work],[feed,flux,mass],rtol=3e-7,atol=2e-10)
    if min(c,mass)<=0:
        raise ArithmeticError('Unphysical kinetic norm or conductance')
