"""Closed-loop velocity on the shifted material graph of a held R=1 mouth."""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from scipy.sparse.linalg import spsolve

from lateral_source_rate import verified_packet
from mouth_advected_rim import EVIDENCE, RATE, Sphere, archive_drive, grid


K = .4
LOOPS = ((.25, .8, .75), (.35, .9, .55))


def flat_fields(path, pins, slab):
    ix = abs(slab.x) < 1 + 5*slab.dx
    ir = slab.r < 1 + 5*slab.dx
    with verified_packet(path, pins) as packet:
        face = packet['layers' if 'layers' in packet else 'field'][-1][ix][:, ir]
    x, r = slab.x[ix], slab.r[ir]
    density = abs(face)**2
    phase = np.angle(face * face[len(x)//2, 0].conjugate())
    if density.min() <= 0 or abs(phase).max() >= np.pi/2:
        raise ValueError('A single smooth phase branch has not been verified')
    even_phase = np.concatenate((phase[:, ::-1], phase), axis=1)
    spline = RectBivariateSpline(x, np.r_[-r[::-1], r], even_phase)
    density_itp = RegularGridInterpolator((x, np.r_[0, r]),
                         np.concatenate((density[:, :1], density), axis=1))
    return density_itp, spline


def lifted_shape(sphere, data):
    # The transport=False system of Sphere.solve, retaining the solution coefficients.
    active = sphere.free.ravel()
    matrix = K*data.diffusion
    result = np.zeros(active.size)
    result[active] = spsolve(matrix[active][:, active], data.rhs[active])
    return result.reshape(sphere.indices.shape)


def points(a, b, pieces):
    roots, weights = leggauss(5)
    width = (b-a)/pieces
    nodes = (a + width*(np.arange(pieces)[:, None]+(roots+1)/2)).ravel()
    return nodes, np.tile(weights*width/2, pieces)


def circulation(sphere, coeffs, sample, normalization, fields, loop, pieces):
    lo, hi, spread = loop
    radius, wr = points(lo, hi, pieces)
    mu, wm = points(-spread, spread, pieces)
    # Increasing radius toward the companion, then across on the outer arc.
    rr = np.r_[radius, hi+0*mu, radius, lo+0*mu]
    mm = np.r_[spread+0*radius, mu, -spread+0*radius, mu]
    dr = np.r_[wr, 0*wm, -wr, 0*wm]
    dm = np.r_[0*wr, -wm, 0*wr, wm]
    p = np.column_stack([Legendre.basis(i)(mm) for i in range(sphere.modes)])
    dp = np.column_stack([Legendre.basis(i).deriv()(mm) for i in range(sphere.modes)])
    hr = np.einsum('ij,ij->i', sphere.basis(rr, 1)@coeffs, p)
    hm = np.einsum('ij,ij->i', sphere.basis(rr)@coeffs, dp)
    dh = hr*dr + hm*dm
    sine = np.sqrt(1-mm**2)
    x, transverse = -rr*mm, rr*sine
    density_itp, phase = fields
    n = density_itp(np.column_stack((x, transverse)))
    q0, ja, jt = sample(rr, mm)
    vz = normalization*q0/n
    phase_a = -phase.ev(x, transverse, dx=1)
    phase_t = phase.ev(x, transverse, dy=1)
    da, dt = mm*dr+rr*dm, sine*dr-rr*mm/sine*dm
    grad_a = mm*hr + sine**2*hm/rr
    grad_t = sine*hr - mm*sine*hm/rr
    # For the gradient reconstruction the flat line element is exactly d(phi).
    # Also evaluate it by line quadrature as a control on the closed-loop geometry.
    normal_part = np.dot(vz, dh)
    from_phase = np.dot(phase_a*grad_a + phase_t*grad_t, dh)
    from_currents = np.dot((ja*grad_a + jt*grad_t)/n, dh)
    return dict(normal=float(normal_part), phase_tilt=float(from_phase),
                current_tilt=float(from_currents), lifted=float(normal_part+from_phase),
                lifted_fd=float(normal_part+from_currents), closure_h=float(dh.sum()),
                closure_flat_phase=float(np.dot(phase_a, da)+np.dot(phase_t, dt)))


def main():
    evidence = json.loads(EVIDENCE.read_text())['runs']['mouth_front_pressure_dipole02']
    assert evidence['limits_verified'] and evidence['guard_exit'] == evidence['child_exit'] == 0
    for raw in evidence['stdout_lines']:
        record = json.loads(raw)
        if record['grid'][0] != 96:
            continue
        slab = grid(*record['grid'])
        pins = {Path(name): digest for name, digest in record['inputs'].items()}
        samples = archive_drive(record)
        for path, sample, label in zip(pins, samples, ('single', 'pair'), strict=True):
            if label == 'single' and record['d'] != 8:
                continue
            fields = flat_fields(path, pins, slab)
            for elements, modes, angles in ((16, 5, 48), (32, 8, 96)):
                sphere = Sphere(elements, modes, angles)
                data = sphere.prepare(sample)
                coeffs = lifted_shape(sphere, data)
                for loop in LOOPS:
                    for pieces in (32, 64):
                        out = circulation(sphere, coeffs, sample, RATE/data.projected_rate,
                                          fields, loop, pieces)
                        print(json.dumps(dict(grid=record['grid'], d=record['d'] if label=='pair' else 0,
                              elements=elements, modes=modes, angles=angles, loop=loop,
                              pieces=pieces, K=K, **out), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
