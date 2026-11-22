"""Prepared bent, intact or perforated order surface in one coherent 4D carrier.

Opposite spatial phase steps prepare uniform bulkward or brane-ward velocities.
Measure instantaneous order kinetics and the external distributed order force
required to keep the given configuration still. The puncture removes order
through the whole active depth; two smooth rim tapers control its edge cost.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

from mouth_flat_phase_control import saved_bath
from mouth_free_front_wave_head import emit


RADIUS = 100.
BOX = 16.
HOLE = 8.
PROBE = 14.
LENGTH = 32.
MOBILITY = .5
KAPPA = 2.


def gauss(start, stop, count):
    x, weight = leggauss(count)
    return (start+stop+(stop-start)*x)/2, (stop-start)*weight/2


def order(r, w, width):
    """Return the value, two cylindrical gradients, bare 4D affinity and free energy."""
    distance = np.hypot(r, RADIUS+w)
    normal_r = r/distance
    normal_w = (RADIUS+w)/distance
    q = expit(RADIUS-distance)
    first = -q*(1.-q)
    second = q*(1.-q)*(1.-2*q)
    if width:
        t = np.clip((r-HOLE)/width, 0., 1.)
        taper = t**3*(10.-15.*t+6.*t*t)
        dt = 30.*t*t*(1.-t)**2 / width
        ddt = 60.*t*(1.-t)*(1.-2*t) / width**2
    else:
        taper, dt, ddt = 1., 0., 0.
    chi = taper*q
    cr = taper*first*normal_r + dt*q
    cw = taper*first*normal_w
    laplacian = (taper*(second+3*first/distance) + 2*dt*first*normal_r
                 + q*(ddt+2*dt/r))
    potential = (chi*(1.-chi))**2
    mu = 2.*chi*(1.-chi)*(1.-2.*chi) - KAPPA*laplacian
    return chi, cr, cw, mu, potential


def z_stress(r, w, width):
    _, cr, cw, _, potential = order(r,w,width)
    return (KAPPA*cw*cw - potential - KAPPA*(cr*cr+cw*cw)/2, KAPPA*cr*cw)


def measure(mesh, width, phase):
    regions = [0.,HOLE,BOX] if not width else [0.,HOLE,HOLE+width,BOX]
    parts = [gauss(left,right,mesh) for left,right in zip(regions,regions[1:])]
    r, wr = (np.concatenate(v) for v in zip(*parts,strict=True))
    w, ww = gauss(-LENGTH/2,LENGTH/2,2*mesh)
    darea = 4*np.pi*r*r*wr
    weight = darea[:,None]*ww
    chi, cr, cw, mu, _ = order(r[:,None],w,width)
    u = phase/LENGTH
    hold_static = np.sum(weight*(-mu*cw))
    hold_increment = -np.sum(weight*cw*cw)*u/MOBILITY
    # Measured independently by evaluating pressure and order stresses at all boundaries.
    top, _ = z_stress(r,LENGTH/2,width)
    bottom, _ = z_stress(r,-LENGTH/2,width)
    _, side = z_stress(BOX,w,width)
    surface_stress = np.dot(darea,top-bottom) + 4*np.pi*BOX**2*np.dot(ww,side)
    held_mu = -u*cw/MOBILITY
    distributed_holder = np.sum(weight*(held_mu-mu)*cw)
    loss = np.sum(weight*MOBILITY*held_mu**2)

    # Same graph in both cases, continued as a virtual flux cut where order is absent.
    nz = np.sqrt(1-(r/RADIUS)**2)
    projected = darea
    cut = r<HOLE
    occupied = (r>=HOLE+width) if width else np.ones_like(r,dtype=bool)
    central_rate = np.dot((projected/nz)[cut], (u*nz)[cut])
    common = r>=HOLE
    annular_rate = np.dot((projected/nz)[common], (u*nz)[common])
    probe_w = np.sqrt(RADIUS**2-PROBE**2)-RADIUS
    probe_chi, rr, zz, bare_mu, _ = order(PROBE,probe_w,width)
    mag = np.hypot(rr,zz)
    normal_speed = -u*zz/mag
    unheld = normal_speed - MOBILITY*bare_mu/mag
    midpoint_w = np.sqrt(RADIUS**2-1.)-RADIUS
    aperture_chi = order(1.,midpoint_w,width)[0]
    emit(test='matched_bend', mesh=mesh, rim_taper=width, phase=phase,
         phase_length=LENGTH, density=1., static_pressure=.2, carrier_speed=u,
         projected_hole_radius=HOLE, outer_radius=BOX, sphere_radius=RADIUS,
         principal_curvature_each=1/RADIUS, retained_probe_radius=PROBE,
         retained_half_order=probe_chi, central_reference_order=aperture_chi,
         center_contains_order_face=not bool(width),
         jn_retained_probe=normal_speed, jn_virtual_cut_center=u,
         retained_velocity_with_static_holder=normal_speed,
         retained_velocity_all_support_removed=unheld,
         retained_velocity_with_full_distributed_holder=0.,
         center_order_velocity_with_static_holder=u if not width else None,
         center_order_velocity_without_support=u-3*MOBILITY*KAPPA/RADIUS if not width else None,
         virtual_central_flux=central_rate, flux_outer_annulus=annular_rate,
         carrier_normal_stress_retained=.2+normal_speed**2,
         two_sided_carrier_stress_jump=0.,
         axial_holder_whole_order_profile=distributed_holder,
         axial_holder_zero_flow=hold_static, axial_holder_increment=hold_increment,
         static_support_from_boundary_stresses=surface_stress,
         power_for_sustained_conversion_at_frozen_fields=loss,
         occupied_curved_area_away_from_taper=sum((projected/nz)[occupied]))
    if (abs(hold_static-surface_stress)>2e-7 or abs(hold_static+hold_increment-distributed_holder)>2e-8
            or abs(loss+u*hold_increment)>1e-10 or abs(probe_chi-.5)>1e-12
            or abs(unheld - (normal_speed-3*MOBILITY*KAPPA/RADIUS))>1e-12
            or abs(central_rate+annular_rate-u*4*np.pi*BOX**3/3)>1e-9):
        raise ArithmeticError('Order traction, kinematics, conversion power or cut flux imbalance')


def local_load(mesh):
    # This column follows the common outward normal and never intersects
    # the widest puncture rim or the outer wall. Force is per midsurface area.
    nr = PROBE/RADIUS
    nw = np.sqrt(1.-nr**2)
    for halfspan in (8., 10.):
        s, weight = gauss(-halfspan, halfspan, mesh)
        r, w = nr*(RADIUS+s), nw*(RADIUS+s)-RADIUS
        if min(r)<=HOLE+4. or max(r)>=BOX:
            raise ValueError('Normal sampling interval enters a rim')
        for width in (0.,2.,4.):
            _, cr, cw, mu, _ = order(r,w,width)
            slope = nr*cr+nw*cw
            jac = (1.+s/RADIUS)**3
            for phase in (0., saved_bath(), -saved_bath(), saved_bath()/2, -saved_bath()/2):
                u = phase/LENGTH
                static = weight@(-mu*slope*jac)
                total = weight@((-mu-u*cw/MOBILITY)*slope*jac)
                analytic_increment = -u*nw/MOBILITY*(weight@(jac*slope**2))
                emit(test='retained_normal_column',rim_taper=width,phase=phase,
                     halfspan=halfspan, probe_radius=PROBE, radial_domain=[min(r),max(r)],
                     static_traction=static,total_traction=total,
                     traction_change=total-static,projected_prediction=analytic_increment)
                if abs(total-static-analytic_increment)>1e-13:
                    raise ArithmeticError('Local normal projection disagreement')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=int,choices=(48,96),required=True)
    parser.add_argument('--local-column', action='store_true')
    args=parser.parse_args()
    if args.local_column:
        local_load(args.mesh)
    else:
        for width in [0.,2.,4.]:
            for drive in [0.,1.,-1.,.5,-.5]:
                measure(args.mesh,width,drive*saved_bath())
