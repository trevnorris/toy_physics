"""Map a measured paired-minus-alone aperture field to the radius-one curved carrier.

The two solvers meet by a trial Dirichlet condition on equal fractional radii
of their respective openings. The held receiver retains its isolated bath.
Only odd angular modes have tangential velocity at the equator at first order.
"""
import argparse
import json
from pathlib import Path
import shlex
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.interpolate import RegularGridInterpolator

from incident_flow_falloff import check_current, checked_grid
from lateral_source_rate import verified_packet
from mouth_carrier_lateral_response import Dipole, prepare
from mouth_contact_tracking import WallTrace
from mouth_flared_receiver import ROUNDING
from mouth_front_pressure_dipole import MANIFEST, grid
from mouth_order_wave_pulse import digest, emit


DEGREES = (1, 3, 5, 7, 9)


def face(path, pins, g, centers):
    with verified_packet(path, pins) as p:
        checked_grid(g, p)
        if 'centers' in p:
            np.testing.assert_array_equal(p['centers'], centers)
            np.testing.assert_array_equal(p['rates'], np.full(len(centers), .25))
        elif int(p['d']) != -centers[0]:
            raise ValueError('Archive separation differs')
        currents = p['currents']
        if len(currents) != len(centers):
            raise ValueError('Archive mouth count differs')
        rates = []
        for i, center in enumerate(centers):
            current = check_current(g, center, currents[i:i+1])
            rates.append(float(np.sum(g.dv * current)))
        ix = abs(g.x) < 1 + 2*g.dx
        ir = g.r < 1 + 2*g.dx
        field = p['layers' if 'layers' in p else 'field'][-1, ix, :][:, ir]
    if not np.isfinite(field).all() or np.min(field.real) <= 0:
        raise ValueError('Invalid archive phase reference')
    return g.x[ix], g.r[ir], field, rates


def project(x, radial, log_ratio, scaled_r, count):
    """Integrate in direction cosine about the companion-facing (-x) axis."""
    interp = RegularGridInterpolator(
        (x, np.r_[0, radial]), np.concatenate((log_ratio[:, :1], log_ratio), axis=1))
    root, weight = leggauss(count)
    coefficients = np.zeros((len(DEGREES), len(scaled_r)), complex)
    for i, r in enumerate(scaled_r):
        if r == 0:
            continue
        # At each interpolant cell boundary restart the Gauss panel.
        turns = np.sqrt(1-(radial[radial < r]/r)**2)
        edges = np.unique(np.r_[-1, 1, -x[abs(x) < r]/r, turns, -turns])
        half = np.diff(edges)/2
        mu = ((edges[1:] + edges[:-1])[:, None]/2 + half[:, None]*root).ravel()
        weights = (half[:, None]*weight).ravel()
        z = interp(np.c_[-r*mu, r*np.sqrt(1-mu**2)])
        for k, l in enumerate(DEGREES):
            coefficients[k, i] = (2*l+1)/2 * (weights*Legendre.basis(l)(mu)) @ z
    return coefficients


def inputs(scaled_r):
    coeff = {}
    records = json.loads(MANIFEST.read_text())['runs']
    for d in (8, 16):
        entry, = [a for a in records if a['name'] == f'L96C6_d{d}']
        if entry['guard_exit'] or entry['worker_exit'] or not entry['limits_verified']:
            raise ValueError('Guarded phase archive failed')
        reference = json.loads(entry['stdout'][-1])
        if reference['d'] != d:
            raise ValueError('Archive run distance mismatch')
        g = grid(*reference['grid'])
        pins = {Path(k): h for k,h in (entry['input_sha256']|entry['output_fields_sha256']).items()}
        command = shlex.split(entry['command'])
        baseline, paired = [Path(command[command.index(flag)+1]) for flag in ('--baseline', '--output')]
        x, radial, isolated, rate0 = face(baseline, pins, g, [0])
        xp, rp, other, rates = face(paired, pins, g, [-d, 0])
        np.testing.assert_array_equal(x, xp)
        np.testing.assert_array_equal(radial, rp)
        # Real part is delta log amplitude, imaginary part is phase difference.
        change = np.log(other/isolated)
        coeff[d] = project(x, radial, change, scaled_r, 16)
        coarse = project(x, radial, change, scaled_r, 8)
        error = np.max(abs(coeff[d]-coarse))
        emit(test='feeding', d=d, grid=reference['grid'], pins={str(p):pins[p] for p in (baseline,paired)},
             alone_rate=rate0, paired_rates=rates, phase_coeff_at_wall=coeff[d][:,-1].imag.tolist(),
             input_quadrature_max_delta=float(error), max_abs_log_ratio=float(np.max(abs(change))))
        if error > 1e-9 or np.max(abs(change)) > .1:
            raise ArithmeticError('Unresolved phase projection or large perturbation')
    return coeff


def main(args):
    trace = WallTrace(48, .5)
    dom, field = prepare(SimpleNamespace(source=Path('_scratch/s11c/mouth_flared_delayed_n48')))
    np.testing.assert_array_equal(dom.points, trace.geometry.points)
    np.testing.assert_array_equal(trace.psi, field[trace.nodes])
    radius, w = dom.points.T
    waist = np.flatnonzero(w == ROUNDING)
    waist = waist[np.argsort(radius[waist])]
    scaled = radius[waist]/radius[waist[-1]]
    np.testing.assert_allclose([scaled[0], radius[waist[-1]]], [0, dom.neck], atol=1e-11)
    wg = w[trace.nodes]
    np.testing.assert_allclose(wg, np.column_stack((trace.edges[:-1], (trace.edges[1:]+trace.edges[:-1])/2, trace.edges[1:])), atol=1e-10, rtol=0)
    coeffs = inputs(scaled)
    result = dict(w_nodes=wg, wall_edges=trace.edges, degrees=DEGREES)
    for d in coeffs:
        result[f'modes{d}'] = []
        result[f'phase_modes{d}'] = []
    # A reflecting wall is the natural boundary of the parent Galerkin operator.
    for j, l in enumerate(DEGREES):
        solver = Dipole(dom, field, degree=l, active=w > ROUNDING)
        equatorial = Legendre.basis(l).deriv()(0)
        for d, values in coeffs.items():
            errors = []
            for mode in ('modes', 'phase_modes'):
                bnd = np.zeros_like(field)
                z = values[j] if mode == 'modes' else 1j*values[j].imag
                bnd[waist] = field[waist]*z
                sol, err = solver.solve(bnd)
                np.testing.assert_array_equal(sol[waist], bnd[waist])
                np.testing.assert_array_equal(sol[dom.bulk], 0.)
                phase = (sol[trace.nodes]/trace.psi).imag
                v = equatorial*phase/radius[trace.nodes]
                result[f'{mode}{d}'].append(v)
                errors.append(err)
            emit(test='mode', d=d, degree=l, max_equation_residual=max(errors))
        del solver
    levels = (.95,.75,.5,.3,.1,.01)
    stations = np.array([trace.crossing(c) for c in levels])
    sample_w = np.linspace(stations[0],stations[-1],401)
    for d in coeffs:
        full = np.sum(result[f'modes{d}'],axis=0)
        phase = np.sum(result[f'phase_modes{d}'],axis=0)
        result[f'velocity{d}'] = full
        result[f'phase_velocity{d}'] = phase
        v = trace.sample(full,sample_w)[0]
        scale = np.max(abs(v))
        def relative(a):
            return float(np.max(abs(trace.sample(a,sample_w)[0]))/scale)
        emit(test='wall_reading', d=d, order=levels, waist_radius=dom.neck,
             velocity_at_waist=float(full[0,0]), velocity_at_order=trace.sample(full,stations)[0].tolist(),
             velocity_min=float(v.min()), velocity_max=float(v.max()),
             phase_only_difference_relative=relative(full-phase),
             dipole_only_difference_relative=relative(full-result[f'modes{d}'][0]),
             last_mode_relative=relative(result[f'modes{d}'][-1]),
             native_carrier_pins=trace.hashes)
        if not np.all(np.isfinite(v)) or scale <= 0:
            raise ArithmeticError('Invalid predicted transport')
    np.savez_compressed(args.output, **result)
    emit(test='output', file=str(args.output), sha256=digest(args.output))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    main(p.parse_args())
