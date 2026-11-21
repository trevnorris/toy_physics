"""Excite a local toroidal stress with finite duration; read the surviving spatial pattern."""
import argparse
from pathlib import Path

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

from mouth_shear_pattern_release import angular_stiffness, assemble, digest, emit, impedance, momentum_and_mass, prepared, truncate


def run(a):
    dom, chi, intake, src = prepared(a)
    dom, chi, cut = truncate(dom, chi, a.cut)
    points, blocks, _, parts, remote, _ = assemble(dom, chi, a.sub, 1, a.feed_spacing)
    masses, _ = momentum_and_mass(dom, blocks, chi, len(points), 0.)
    cut_drag, _, _ = impedance(dom, blocks, points, chi, cut, 0.)
    if a.degree > 1:
        angle = angular_stiffness(dom, blocks, chi, len(points))
        parts = [s+(a.degree*(a.degree+1)-2)*b for s, b in zip(parts, angle, strict=True)]
    full_mass = masses.sum(axis=0)
    live = (full_mass > 0) & ((points[:, 0] > 0) | (a.degree == 1))
    mass, (r, w) = full_mass[live], points[live].T
    remote, cut_drag = remote[live], cut_drag[live]
    k = sum(parts).tocsr()[live][:, live]
    feed_stiff = parts[0][live][:, live]
    feed_mass = masses[0, live]
    # A clamped reference is used only as a spatial readout, never as the drive.
    behind = w > .5000000001
    scale = diags(1/np.sqrt(mass[behind]))
    op = scale@k[behind][:, behind]@scale
    eig, vec = eigsh(op, k=1, sigma=0, v0=np.ones(op.shape[0]), tol=2e-12)
    omega = np.sqrt(eig[0])
    reference = np.zeros_like(r)
    reference[behind] = scale@vec[:, 0]
    mr = mass*reference
    # Local torque density factors into a specified angular harmonic and this
    # smooth interior profile; regular at the axis. Its overall scale is arbitrary.
    g = r**(a.degree-1)*np.exp(-(r/.7)**2)*np.maximum(0., 1-((w-.95)/.18)**2)**3
    force = mass*g
    force /= np.sqrt(force@(force/mass))
    forcing = np.flatnonzero(force)
    edge, far = np.flatnonzero(cut_drag), np.flatnonzero(remote)
    ro = np.flatnonzero(reference)
    # An independent ordered patch on the flare side of the imposed torque.
    sensor = mass*np.maximum(0., 1-((w-1.19)/.055)**2)*r**(a.degree-1)*np.exp(-(r/.7)**2)
    sensor /= sensor.sum()
    so = np.flatnonzero(sensor)
    dt_sample = .08
    bound = np.max(np.asarray(abs(k).sum(axis=1)).ravel()/mass)
    stride = int(np.ceil(dt_sample*np.sqrt(bound)/1.4))
    dt = dt_sample/stride
    drive = a.drive
    duration = 12*np.pi/drive
    damping = remote+cut_drag
    update = dt/(mass+.5*dt*damping)
    damp_ratio = (mass-.5*dt*damping)/(mass+.5*dt*damping)
    emit(test='local_source', degree=a.degree, radius=a.radius, sub=a.sub, source_depth=.95,
         readout_depth=1.19, end=a.end, outer=a.outer, feed_spacing=a.feed_spacing, cut=cut,
         order_file=str(src), order_sha256=digest(src), intake=intake, drive_frequency=drive,
         drive_duration=duration, clamp_frequency=omega, initial_overlap_squared=(reference@force)**2,
         trace_spacing=dt_sample, dt=dt)
    u = np.zeros_like(mass)
    vel = np.zeros_like(mass)
    work = feed_loss = cut_loss = 0.
    trace = []
    for i in range(round(a.end/dt_sample)*stride):
        t = i*dt
        residual = k@u
        amp = 1e-4*np.sin(np.pi*t/duration)**2*np.cos(drive*(t-duration/2)) if t < duration else 0.
        if amp:
            residual[forcing] -= amp*force[forcing]
        vn = damp_ratio*vel-update*residual
        un = u+dt*vn
        if amp:
            work += .5*dt*amp*(force[forcing]@(vn[forcing]+vel[forcing]))
        feed_loss += .25*dt*(remote[far]@(vn[far]+vel[far])**2)
        cut_loss += .25*dt*(cut_drag[edge]@(vn[edge]+vel[edge])**2)
        if (i+1)%stride == 0:
            energy = .5*(mass@vn**2+u@(k@un))
            feed = .5*(feed_mass@vn**2+u@(feed_stiff@un))+feed_loss
            trace.append([(i+.5)*dt, work, energy, feed, feed_loss, cut_loss, mr[ro]@vn[ro],
                          mr[ro]@(.5*(u[ro]+un[ro])), sensor[so]@vn[so],
                          mass[behind]@(vn[behind]**2), mass@(vn**2), energy+feed_loss+cut_loss-work])
        u, vel = un, vn
    names = ['time','work','energy','feed_net','remote_loss','cut_loss','projection_velocity',
             'projection_position','patch_velocity','behind_twice_kinetic','twice_kinetic','closure']
    trace = np.array(trace)
    err = np.max(abs(trace[:, -1]))/work
    if work <= 0 or err > 1e-8 or not np.isfinite(trace).all():
        raise ArithmeticError('Driven energy ledger failed')
    file = a.outdir/'probe.npz'
    np.savez_compressed(file, columns=names, trace=trace)
    emit(test='local_release', input_work=work, max_energy_error=err, trace_sha256=digest(file))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.set_defaults(radius=1., mesh=34, state=None)
    p.add_argument('--degree', type=int, required=True)
    p.add_argument('--drive', type=float, required=True)
    p.add_argument('--cut', type=float, default=1.85)
    p.add_argument('--sub', type=int, default=1)
    p.add_argument('--outer', type=float, default=35.)
    p.add_argument('--feed-spacing', type=float, default=.08)
    p.add_argument('--end', type=float, default=60.)
    p.add_argument('--outdir', type=Path, required=True)
    a = p.parse_args()
    if a.degree < 1 or a.drive <= 0 or a.outer < 10 or a.end <= 12*np.pi/a.drive:
        p.error('Choose an angular degree, open feed and positive drive with a release interval.')
    run(a)
