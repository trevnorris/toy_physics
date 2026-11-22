"""Normally retained, freely sliding spherical guide ending inside the open outlet.

This extra contact law projects the incident radial velocity at entry. Every tag
then remains on r=R(w)-offset and feels only normal reaction until the true outlet.
No previous Morse-well stage is added. In the receiver each tag coasts in four
coordinates. All force results are normalized by the assumed tagged mass flux.
"""
import argparse

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

from mouth_calibrated_flare_coating import read_inputs
from mouth_flared_receiver import EXIT
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import TEMPERATURE, emit


def run(mesh, subdivisions, offset):
    trace, intake = cubic_trace(mesh)
    patch, _, _ = read_inputs(mesh)
    spread = np.sqrt(2*patch['shared_scale']*TEMPERATURE)
    entry = 1.175
    targets = np.array([trace.crossing(.5), trace.crossing(.01), 1.65, 2.4, EXIT])
    edges = np.unique(np.r_[entry, targets, trace.edges[(trace.edges>entry) & (trace.edges<EXIT)]])
    mu, weight = leggauss(24)
    weight /= 2
    angular = np.sqrt(1-mu*mu)
    speed = np.array([0., 2.e-5, -2.e-5, 1.e-5, -1.e-5])
    U = speed[:, None]
    radial, radial_weight = hermgauss(12)
    radial = radial*np.sqrt(2)*spread
    radial_weight /= np.sqrt(np.pi)
    a0 = trace.geometry.radius(np.array([entry]))[0][0]-offset
    angular_momentum = a0*U*angular
    # Projection onto moving surface: discard relative normal speed but transmit
    # the tangential part. Relative lost kinetic energy includes the undriven heat.
    relative_radial = radial[None,:,None]+U[:,None,:]*mu
    joint = radial_weight[:,None]*weight
    capture = np.sum(joint*relative_radial*mu, axis=(1,2))
    dissipated = .5*np.sum(joint*relative_radial**2, axis=(1,2))
    np.testing.assert_allclose(capture,speed/3,atol=1.e-16)
    np.testing.assert_allclose(dissipated,.5*(spread**2+speed**2/3),rtol=1.e-13)
    # Parameter phi advances around each initial great circle. Its unit radial
    # and tangent vectors have x projections nx(phi) and tx(phi).
    def direction(tau):
        phase = tau[...,None,None]*angular_momentum
        return (mu*np.cos(phase)+angular*np.sin(phase),
                -mu*np.sin(phase)+angular*np.cos(phase))

    ratio = lambda x: [(x[1]-x[2])/(2*speed[1]),(x[3]-x[4])/(2*speed[3])]
    emit(test='normal_guide_setup', mesh=mesh, radius=1., intake=intake, entry=entry, outlet=EXIT,
         offset=offset, nodes_per_element=subdivisions, entering_radius=a0, speeds=speed.tolist(),
         radial_sigma=spread, absorbed_relative_energy_per_tag_mass=dissipated.tolist(),
         additional_absorbed_energy_per_incident_speed_squared=1/6,
         inlet_accommodation_momentum=capture.tolist(), accommodation_fraction=ratio(capture))
    total_tau = 0.
    pressure_on_tag = np.zeros_like(speed)
    last_radial_speed = None
    roots, weights = leggauss(4)
    for lo,hi in zip(edges[:-1],edges[1:],strict=True):
        nodes = np.linspace(lo,hi,subdivisions+1)
        widths = np.diff(nodes)
        gauss = nodes[:-1,None] + .5*widths[:,None]*(1+roots)
        ag = trace.geometry.radius(gauss)[0]-offset
        ug = trace.speed(gauss)
        dt = .5*widths*((1/ug)@weights)
        step_tau = .5*widths*((1/(ag**2*ug))@weights)
        times = np.r_[total_tau, total_tau+np.cumsum(step_tau)]
        nx, tx = direction(times)
        sample = nodes.copy()
        sample[0],sample[-1] = np.nextafter(lo,hi),np.nextafter(hi,lo)
        radii, slopes = trace.geometry.radius(sample)
        radii -= offset
        uz = trace.speed(sample)
        rad_speed = slopes*uz
        if last_radial_speed is not None:
            pressure_on_tag += (rad_speed[0]-last_radial_speed)*(nx[0]@weight)
        # d(v_radial) n_x and the centripetal force −L² n_x/a³ dt are
        # measured separately from the outgoing particle momentum.
        pressure_on_tag += np.diff(rad_speed)@((nx[:-1]+nx[1:])*.5@weight)
        nx_mid,_ = direction((times[:-1]+times[1:])*.5)
        middle_radius = trace.geometry.radius((nodes[:-1]+nodes[1:])*.5)[0]-offset
        pressure_on_tag -= (dt/(middle_radius**3))@((nx_mid*angular_momentum**2)@weight)
        total_tau = times[-1]
        last_radial_speed = rad_speed[-1]
        if hi in targets:
            vx = last_radial_speed*nx[-1]+angular_momentum/radii[-1]*tx[-1]
            mean_exit = vx@weight
            wall = capture-pressure_on_tag
            exported = speed-mean_exit
            emit(test='guide_cut', depth=hi, radius=radii[-1],
                 wall_momentum_per_flux=wall.tolist(), tag_velocity_budget=exported.tolist(),
                 pressure_transfer_fraction=ratio(wall), momentum_transfer_fraction=ratio(exported),
                 max_per_flux_discrepancy=float(max(abs(wall-exported))))
    # At the outlet the radius has zero slope. The exact weak-drive response
    # to an inlet angular speed follows from constant angular momentum.
    linear_exit = 1-2*a0/(3*radii[-1])
    np.testing.assert_allclose(ratio(exported),linear_exit,atol=3e-6)
    np.testing.assert_allclose(ratio(wall),ratio(exported),atol=3e-5)
    if offset<=0: raise ValueError('Require finite clearance from the reentrant exit edge')
    if np.min(uz)<=0: raise ArithmeticError('Invalid exit velocity')
    for plane in (EXIT+.10,EXIT+1.,7.):
        flight=(plane-EXIT)/uz[-1]
        # Finite-U free flight to a full receiver plane; all tags retain vz.
        radii_at_plane=np.hypot(radii[-1]+flight*last_radial_speed,flight*angular_momentum/radii[-1])
        if np.max(radii_at_plane)>=trace.geometry.far:
            raise ArithmeticError('Missing a lateral crossing before the plane')
        emit(test='open_receiver', depth=plane, flight_time=flight, clearance_at_opening=offset,
             crossed_tag_fraction=float(np.sum(weight)),
             radial_extent=[float(radii_at_plane.min()),float(radii_at_plane.max())],
             pressure_transfer_fraction=ratio(wall), momentum_transfer_fraction=ratio(exported),
             linear_endpoint_formula=linear_exit, receiver_force=0.)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,required=True,choices=[26,34])
    p.add_argument('--subdivisions',type=int,required=True)
    p.add_argument('--offset-pitches',type=float,required=True)
    args=p.parse_args()
    run(args.mesh,args.subdivisions,args.offset_pitches*PITCH)
