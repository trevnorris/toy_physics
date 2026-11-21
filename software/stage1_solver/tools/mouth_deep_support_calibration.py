"""Phase 1: test local radial equilibrium and forecast a differential sideways push.

All four radial bands are inside the deeper lateral bend. Fit only the first and
last bands. No lateral force, free acceleration or release result enters the fit.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_carrier_lateral_response import prepare, trace
from mouth_flared_receiver import emit

BANDS = ((1.22, 1.40), (1.42, 1.59), (1.62, 1.79), (1.81, 1.99))
ARCHIVE = Path('research/pde_ledger_v3/_measurements/mouth_free_bend_response.json')


def quadrature(edges, lower, upper):
    cuts = np.r_[lower, edges[(edges > lower) & (edges < upper)], upper]
    half = np.diff(cuts)[:, None] / 2
    nodes, weights = leggauss(8)
    return (cuts[:-1, None] + half*(nodes+1)).ravel(), (half*weights).ravel()


def stationary_energy(field, radial, axial):
    density = abs(field)**2
    potential = (density-1)**2*(density**3 + 2*density**2 + 3*density + 4)/20
    return (abs(radial)**2+abs(axial)**2)/2 + potential


def band_work(dom, field, lo, hi):
    axial_grid, radial_grid, _ = dom.blocks[2]
    w, dw = quadrature(axial_grid, lo, hi)
    t = (w-lo)/(hi-lo)
    bump = 16*t**2*(1-t)**2
    derivative = 32*t*(1-t)*(1-2*t)/(hi-lo)
    radius, slope = dom.radius(w)
    root = np.sqrt(1+slope**2)
    area_work = 4*np.pi * dw @ (2*radius*bump*root + radius**2*slope*derivative/root)
    volume_work = 4*np.pi * dw @ (radius**2*bump)

    def area_difference(step):
        def area(sign):
            return 4*np.pi * dw @ ((radius+sign*step*bump)**2 *
                                   np.sqrt(1+(slope+sign*step*derivative)**2))
        return (area(1)-area(-1))/(2*step)

    np.testing.assert_allclose(area_difference(1e-6), area_work, rtol=1e-7)
    # On a reflecting wall the continuum normal derivative vanishes. The
    # radial surface variation then measures the outward pressure -L.
    wall = trace(dom, 2, w, np.zeros_like(w), field[:, None])
    f, fr, fw = (a[:, 0] for a in wall[-3:])
    pressure = -stationary_energy(f, fr, fw)
    surface_work = 4*np.pi * dw @ (radius**2*bump*pressure)

    # Independent bulk shape derivative of E-Q*delta at fixed intake.
    # Radially extend the compact displacement as eta_r = r*bump/a.
    v, dv = quadrature(radial_grid, 0., 1.)
    z, u = np.meshgrid(w, v, indexing='ij')
    nodes = trace(dom, 2, z.ravel(), u.ravel(), field[:, None])
    r = nodes[0].reshape(z.shape)
    f, fr, fw = (a[:, 0].reshape(z.shape) for a in nodes[-3:])
    scale = (bump/radius)[:, None]
    axial_slope = r*((derivative-slope*bump/radius)/radius)[:, None]
    shape_density = ((3*stationary_energy(f, fr, fw)-abs(fr)**2)*scale
                     - (fr.conj()*fw).real*axial_slope)
    measure = 4*np.pi*r**2 * (dw*radius)[:, None] * dv
    bulk_work = -np.sum(measure*shape_density)
    return area_work/volume_work, surface_work/volume_work, bulk_work/volume_work


def radial_test(mesh):
    source = Path(f'_scratch/s11c/mouth_flared_cubic_n{mesh}')
    dom, field = prepare(SimpleNamespace(source=source))
    assert dom.profile == 'cubic'
    curvature, surface, bulk = np.array([band_work(dom, field, *b) for b in BANDS]).T
    sample = [0, 3]
    design = np.column_stack((np.ones(len(BANDS)), curvature))
    for method, pressure in [('surface', surface), ('volume', bulk)]:
        backing, tension = np.linalg.solve(design[sample], pressure[sample])
        predicted = backing+tension*curvature
        emit(test='radial_equilibrium', mesh=mesh, method=method,
             bands=BANDS, fit_indices=sample, curvature=curvature.tolist(),
             actual_pressure=pressure.tolist(), backing=backing, tension=tension,
             predicted_pressure=predicted.tolist(), residual=(pressure-predicted).tolist())


def lateral_test(mesh, rows):
    def unique(test, **criteria):
        row, = [r for r in rows if r['test']==test and all(r[k]==v for k,v in criteria.items())]
        return row
    area = np.array(unique('static')['area'])
    data = unique('impedance', omega=.005)
    carrier = np.array(data['real']) + 1j*np.array(data['imag'])
    # q=(entrance translation, middle-entrance, flare-middle). Unit physical
    # forces on adjacent stations with opposite signs do work on one offset.
    read = np.tri(3)
    drives = np.diff(np.eye(3), axis=0).T
    loads = read.T @ drives
    for tension in (.02, .04):
        z = carrier + np.diag(area*np.array([0., .02, tension]))
        q = np.linalg.solve(z, loads)
        gauges = read @ q
        np.testing.assert_allclose(z@q, loads, atol=1e-10)
        np.testing.assert_allclose(z@np.linalg.solve(z, -.5*loads), -.5*loads, atol=1e-10)
        imp = np.linalg.solve(q[1:].T, (loads-carrier@q)[1:].T).T
        np.testing.assert_allclose(imp, np.diag(area[1:]*[.02,tension]), atol=1e-10)
        emit(test='balanced_push_forecast', mesh=mesh, first=.02, deeper=tension,
             omega=.005, peak_force=1e-9, applied_station_forces_per_peak=drives.tolist(),
             positions_per_force_real=gauges.real.tolist(), positions_per_force_imag=gauges.imag.tolist(),
             signed_phase_lag_degrees=(-np.angle(gauges, deg=True)).tolist(),
             deep_offset_per_flare_pair=[q[2,1].real,q[2,1].imag],
             inferred_from_own_generated_motion=imp.real.tolist())


if __name__ == '__main__':
    for n in (26, 34, 48):
        radial_test(n)
    archived = json.loads(ARCHIVE.read_text())['runs']
    for n, extension in ((26,'compact'),(34,'wide')):
        run = archived[f'free_center_bend_n{n}_{extension}_01']
        assert run['guarded_exit']==run['worker_exit']==0 and run['limits_verified']
        lateral_test(n, [json.loads(line) for line in run['literal_stdout']])
