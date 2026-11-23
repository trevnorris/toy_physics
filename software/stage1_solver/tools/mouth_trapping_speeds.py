"""Map the local transverse principal symbol of a saved strained crystal and its occupied wave."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import legder, leggauss, legval
from scipy.interpolate import RegularGridInterpolator

from mouth_exit_witness import inputs
from mouth_finite_bore_support import graph_values
from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_wave import OccupiedWave
from mouth_poly_bonded_core import core_stretches
from mouth_poly_shallow_release import Trial
from mouth_polyconvex_neutral import stress
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT
from mouth_strained_wave import install


def emit(**items):
    print(json.dumps(items, allow_nan=False), flush=True)


def shear_check(F, chi, normal):
    # Incremental stress on a transverse plane wave in the current configuration.
    # Compare work per reference mass with the speed from the optical pencil.
    a = np.eye(4)[np.argmin(abs(normal))]
    a -= normal @ a * normal
    a /= np.linalg.norm(a)
    H = np.outer(a, normal @ F)
    step = 2e-5
    change = stress(F+step*H, chi)-stress(F-step*H, chi)
    measured = np.sum(change*H)/(2*step*chi)
    return float(measured), abs(measured/(np.linalg.norm(normal @ F)**2)-1)


def sample(trial, shifts, chi, seed, quadrature):
    shape = Meridian(trial.radius, trial.x+shifts[:4])
    coeff = trial.co.copy()
    coeff[0] += shifts[4]
    wave = OccupiedWave(shape, trial.n, 2*trial.n, trial.mesh, PORT, seed)
    install(wave, shape, coeff, core_reference=trial.core_reference)
    wave = ActiveWave.install(SimpleNamespace(wave=wave)).wave
    cell_chi = np.asarray(wave.project @ chi.T.ravel())
    vector = seed[1]
    mode = np.r_[vector, 0j][wave.ids]
    material = RegularGridInterpolator((np.linspace(0, 1, chi.shape[0]),
                                        np.linspace(0, 1, chi.shape[1])), chi)
    mesh = wave.mesh
    nt = len(cell_chi)//mesh
    i, j = np.indices((mesh, nt)).reshape(2, -1)
    x0, x1 = i/mesh, (i+1)/mesh
    points, weights = leggauss(quadrature)
    points, weights = (points+1)/2, weights/2
    chunks = []
    checked = 0.
    flat_error = 0.
    delta_mass = shape.length/nt
    for z, wz in zip(points, weights):
        t = (j+z)/nt
        b, db, _ = shape.profile(t)
        lp = legval(2*t-1, coeff)
        dlp = 2*legval(2*t-1, legder(coeff))
        split = np.clip(np.exp(-.5*np.exp(lp)), x0, x1)
        for collar, lo, hi in ((False, x0, split), (True, split, x1)):
            keep = hi > lo
            for y, wy in zip(points, weights):
                x = (lo + (hi-lo)*y)[keep]
                r, tk, widths = x*b[keep], t[keep], b[keep]
                if collar:
                    a = -np.log(x)/np.exp(lp[keep])
                    radius, rs, rz, hoop, jac, _ = graph_values(a, tk, shape, lp[keep], dlp[keep])
                    axial = np.full(r.shape, -shape.length)
                    if np.max(abs(radius-r)) > 1e-11:
                        raise ArithmeticError('Material inverse missed the collar')
                else:
                    rs, rz, axial, jac = core_stretches(shape, coeff, *trial.core_reference, tk, r)
                    hoop = rs
                    axial = np.full(r.shape, axial)
                brr, bww = rs**2+rz**2, axial**2
                brw = rz*axial
                gram = np.empty((len(r), 2, 2))
                gram[:, 0, 0], gram[:, 0, 1] = brr, brw
                gram[:, 1, 0], gram[:, 1, 1] = brw, bww
                eigen, axes = np.linalg.eigh(gram)
                cmin, cmax = np.sqrt(eigen[:, 0]), np.sqrt(eigen[:, 1])
                u = (x-x0[keep])*mesh
                N = np.stack(((1-u)*(1-z), u*(1-z), (1-u)*z, u*z), axis=1)
                volume = wy*wz*delta_mass*widths*(hi-lo)[keep]*r**2
                kinetic = volume*cell_chi[keep]/jac*r**2*abs(np.sum(N*mode[keep], axis=1))**2
                local_chi = material(np.column_stack((tk, x)))
                direction = axes[:, :, 0]  # sign arbitrary; angle is relative to increasing r.
                cslow = np.minimum(cmin, hoop)
                cfast = np.maximum(cmax, hoop)
                chunks.append(dict(radius=r, depth=shape.length*tk, t=tk, x=x,
                    collar=np.full(r.shape, collar), cell=i[keep]*nt+j[keep],
                    chi=local_chi, wave_cell_chi=cell_chi[keep], jacobian=jac,
                    current_volume_per_solid_angle=volume, cavity_kinetic_weight=kinetic,
                    c_radial=np.sqrt(brr), c_depth=np.sqrt(bww), c_hoop=hoop,
                    c_meridional_min=cmin, c_meridional_max=cmax,
                    c_min=cslow, c_max=cfast,
                    slow_meridional_nr=direction[:, 0], slow_meridional_nw=direction[:, 1]))
                if np.min(jac) <= 0 or np.min(cslow) <= 0:
                    raise ArithmeticError('Invalid local elastic deformation')
                # Select a nonzero-order sample in each material and quadrature panel.
                pos = np.flatnonzero(cell_chi[keep] > 0)
                if len(pos):
                    k = pos[len(pos)//2]
                    F = np.diag([rs[k], axial[k], hoop[k], hoop[k]])
                    F[0, 1] = rz[k]
                    for normal in (np.array([0., 0., 1., 0.]),
                                   np.r_[direction[k], 0., 0.], np.r_[axes[k, :, 1], 0., 0.]):
                        _, err = shear_check(F, cell_chi[keep][k], normal)
                        checked = max(checked, err)
                        flat, _ = shear_check(np.eye(4), cell_chi[keep][k], normal)
                        flat_error = max(flat_error, abs(flat-1))
    grid = {key: np.concatenate([row[key] for row in chunks]) for key in chunks[0]}
    cavity = wave.matrix(cell_chi[:, None, None]*wave.mc)
    native = np.vdot(vector, cavity@vector).real
    frac = native/np.vdot(vector, (wave.cut_m+cavity)@vector).real
    error = np.sum(grid['cavity_kinetic_weight'])/native-1
    emit(kind='verification', radius=trial.radius, order_mesh=trial.n, field_mesh=mesh,
         quadrature=quadrature, matched_flat_speed=1., matched_flat_speed_squared_error=flat_error,
         max_stress_speed_squared_relative_error=checked,
         quadrature_kinetic_relative_error=error, native_cavity_fraction=frac)
    if abs(error) > 1e-10 or max(checked, flat_error) > 1e-8:
        raise ArithmeticError('Archived wave mass or constitutive speed check failed')
    return grid, frac


def stats(grid, mask, region):
    kinetic = grid['cavity_kinetic_weight']
    volume = grid['current_volume_per_solid_angle']
    active = mask & (grid['wave_cell_chi'] > 0)
    total_mass = np.sum(kinetic)
    if not np.any(active):
        return dict(region=region, active=False, kinetic_fraction_of_cavity=0.)
    v, k = volume[active], kinetic[active]
    summary = dict(region=region, active=True, kinetic_fraction_of_cavity=float(k.sum()/total_mass),
                   active_volume_fraction_in_region=float(v.sum()/volume[mask].sum()))
    for name in ('c_radial', 'c_depth', 'c_hoop', 'c_min', 'c_max'):
        speed = grid[name][active]
        ix = np.argsort(speed)
        q = np.interp([.1, .5, .9], (np.cumsum(v[ix])-.5*v[ix])/v.sum(), speed[ix])
        slow = speed < 1.-1e-12
        summary[name] = dict(range=[float(min(speed)), float(max(speed))],
                             active_volume_quantiles=q.tolist(),
                             slower_active_volume_fraction=float(v[slow].sum()/v.sum()),
                             slower_cavity_kinetic_fraction=float(k[slow].sum()/total_mass))
    return summary


def run(radius, n, output):
    trial = Trial(radius, n, bonded_core=True)
    path = Path(f'_scratch/s11c/mouth_exit_release01/r{radius}n{n}/endpoint.npz')
    shifts, chi, seed, digest, status, _, _ = inputs(path, radius, n, trial)
    if status != 'balanced':
        raise ValueError('Expected saved four-coordinate balanced state')
    endpoint, = [json.loads(row) for row in (path.parent/'stdout').read_text().splitlines()
                 if json.loads(row)['kind'] == 'endpoint']
    grid, fraction = sample(trial, shifts, chi, seed, 3)
    if abs(fraction-endpoint['cavity_fraction']) > 1e-8:
        raise ArithmeticError('Mismatch with archived cavity kinetic fraction')
    out = output/f'r{radius}n{n}.npz'
    np.savez_compressed(out, **grid, endpoint_sha256=digest)
    emit(kind='map', radius=radius, order_mesh=n, endpoint_sha256=digest, path=str(out),
         sha256=hashlib.sha256(out.read_bytes()).hexdigest(), points=len(grid['radius']),
         wave_frequency=seed[0].real, wave_amplitude_loss=-seed[0].imag,
         flat_reference='same positive order, F=identity, c=1 in every direction',
         zero_order='no shear mode in zero-weight cells; speeds there are latent F factors')
    t, order = grid['t'], grid['wave_cell_chi']
    regions = dict(cavity=np.ones(t.shape, bool), core=~grid['collar'], collar=grid['collar'],
                   entry_half=t<.5, exit_half=t>=.5, at_most_half_order=order<=.5,
                   above_half_order=order>.5)
    for region, mask in regions.items():
        emit(kind='speeds', radius=radius, order_mesh=n, **stats(grid, mask, region))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for case in ((1, 48), (1, 72), (2, 48), (2, 72)):
        run(*case, args.output)
