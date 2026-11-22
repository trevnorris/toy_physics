"""Signed depth and mid-flare sensitivities of the fixed-phase rigid comparator.

The partial derivatives measure properties of the carrier field on controlled
geometries. No material energy is tuned to the gravity measurements.
"""
import argparse

import numpy as np
from scipy.integrate import quad
from scipy.sparse.linalg import spsolve

from mouth_carrier_lateral_inertia import volume
from mouth_carrier_matched_bias import at_bias
from mouth_carrier_weak_radius import PhysicalDomain
from mouth_flared_receiver import EXIT, FLARE_START, FLARE_LENGTH, ROUNDING, emit
from mouth_radius_pivot import BIAS, WINDOWS, steady_dipoles


class TrialDomain(PhysicalDomain):
    def __init__(self, radius, depth, bulge, mesh):
        self.bulge = bulge
        super().__init__(radius, depth, mesh, 28., 22., 1.)

    def assemble(self, grids):
        super().assemble(grids)
        # The parent class classifies bath nodes in native coordinates before
        # restoring physical coordinates and validating both boundary sets.
        # Snap the known receiver seam in this classification view. Roundoff
        # in the physical-to-native mapping can put its corner just below EXIT.
        seam = np.isclose(self.points[:,1], EXIT, rtol=0., atol=2e-10)
        self.points[seam,1] = round(EXIT,11)

    def radius(self, z):
        base, slope = super().radius(z)
        x = np.clip((z-FLARE_START)/FLARE_LENGTH, 0., 1.)
        y = x*(1-x)
        envelope = 64*y**3
        derivative = 192*y**2*(1-2*x)/FLARE_LENGTH
        factor = 1+self.bulge*envelope
        return base*factor, slope*factor+base*self.bulge*derivative


def material_area(dom, radius, depth):
    # Integrate the 3D area of the wall in four spatial dimensions, including
    # the rounded lip. Use its angle to avoid an infinite radial slope at w=0.
    def rounded(angle):
        r = radius*(dom.rc-ROUNDING*np.sin(angle))
        ds = ROUNDING*np.hypot(radius*np.cos(angle), depth*np.sin(angle))
        return 4*np.pi*r*r*ds

    def wall(w):
        shape, slope = dom.radius(w)
        return 4*np.pi*(radius*shape)**2*np.hypot(depth, radius*slope)

    return (quad(rounded, 0., np.pi/2, epsabs=1e-10)[0]
            + quad(wall, ROUNDING, FLARE_START, epsabs=1e-10)[0]
            + quad(wall, FLARE_START, EXIT, epsabs=1e-10)[0])


def solve(radius, mode, step, mesh):
    # Choose the native axial coordinate of the roof before mapping it to the
    # exact physical height 22. This gives unique rounded bath membership.
    depth = 22./round(22./(1+step), 8) if mode=='depth' else 1.
    bulge = step if mode=='bulge' else 0.
    dom = TrialDomain(radius, depth, bulge, mesh)
    phase = dom.bulk.astype(float)
    free = dom.free
    phase[free] = spsolve(dom.stiff[free][:,free], -(dom.stiff@phase)[free])
    z, intake, energy = at_bias(dom, phase, BIAS, None)
    field = 1+z
    feed, bulk, tow, acc = steady_dipoles(dom,field,refine=True)
    emit(test='shape', radius=radius, mode=mode, step=step, depth_scale=depth,
         physical_depth=depth*EXIT, bulge=bulge, phase_bias=BIAS, mesh=mesh,
         outer=28., roof=22., slab_thickness=1., D=intake, finite_phase_energy=energy,
         material_wall_area=material_area(dom,radius,depth))
    for left,right,z0,z1 in WINDOWS:
        window = (left*radius, right*radius, z0+(depth-1)*EXIT, z1+(depth-1)*EXIT)
        C, B, T = [volume(dom,np.c_[field,h],window,False)['virtual_wall_force']
                   for h in (feed,bulk,tow)]
        mass = -volume(dom,np.c_[field,acc,tow],window,True)['virtual_wall_force']
        emit(test='shape_response', mode=mode, step=step, radius=radius, window=window,
             D=intake, M=mass, C=C, B=B, T=T, D_over_M=intake/mass, C_over_M=C/mass,
             D_over_R3=intake/radius**3, M_over_R3=mass/radius**3)
        if abs(C+B+T)>5e-8 or min(intake,mass)<=0:
            raise ArithmeticError('Nonpositive current or inertia, or failed force superposition')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',type=int,choices=(1,2),required=True)
    p.add_argument('--mode',choices=('depth','bulge'),required=True)
    p.add_argument('--step',type=float,choices=(-.05,-.025,.025,.05),required=True)
    p.add_argument('--mesh',type=int,default=30)
    args=p.parse_args()
    solve(args.radius,args.mode,args.step,args.mesh)
