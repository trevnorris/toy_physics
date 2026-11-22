"""Signed shape work on a transported, partially shear-free field using closed P1 triangles."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import eigsh
from closed_wave import Shape, arrays, matrix
from shape_work import energy, duct
from mouth_live_waist_wave import OccupiedWave
from mouth_sealed_shape_relax import Meridian


def wave_forms(domain, chi):
    xy, ids, tri, gradient, quadrature = domain
    count = len(xy)
    shape = (len(tri), 3, 3)
    stiffness, mass, axial = (np.zeros(shape) for _ in range(3))
    for interp, radius, volume in quadrature:
        angular = 20 * np.outer(interp, interp)
        curl = 2 * interp + radius[:, None] * gradient[:, :, 0]
        dz = radius[:, None] * gradient[:, :, 1]
        zz = dz[:, :, None] * dz[:, None, :]
        weight = (volume * (chi[tri] @ interp))[:, None, None]
        stiffness += weight * (angular + curl[:, :, None] * curl[:, None, :] + zz)
        mass += weight * radius[:, None, None]**2 * np.outer(interp, interp)
        axial += weight * zz
    return tuple(matrix(a, tri, count) for a in (stiffness, mass, axial))


def sampled_order(source, n):
    with np.load(source) as data:
        raw, target = data['order'], float(data['frequency'])
    grid = tuple(np.linspace(0, 1, s) for s in raw.shape)
    z, r = np.meshgrid(np.linspace(0, 1, 2*n+1), np.linspace(0, 1, n+1), indexing='ij')
    sample = RegularGridInterpolator(grid, raw)(np.c_[z.ravel(), r.ravel()])
    return sample, target


def measure(source, state, n):
    archive = Path('research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json')
    training = json.loads(archive.read_text())['training']['1']
    center = training['extension_states'][1 if state == 'mid' else -1]
    ref = np.asarray(center['coordinates'])
    chi, target = sampled_order(source, n)
    initial = arrays(Shape(1, ref), n)
    xy = initial[0]
    stiff, mass, axial = wave_forms(initial, chi)
    boundary_free = (xy[:, 0] > 0) & (xy[:, 1] > 0)
    supported = boundary_free & (mass.diagonal() > 0)
    reduce = lambda M: M[supported][:, supported]
    vals, modes = eigsh(reduce(stiff), M=reduce(mass), k=4, sigma=target**2,
                        v0=np.ones(int(supported.sum())), tol=1e-11)
    # The physical P1 eigenvectors and eigenvalues are independent. Read the observed
    # pole amplitude to match the same occupation across different discretizations.
    with np.load(source) as data:
        nas, nrs = np.asarray(data['order'].shape)-1
        fieldmesh = {48: 14, 72: 20}[nrs]
        model = OccupiedWave(Meridian(1, ref), nrs, nas, fieldmesh, 4.2)
        cells = np.r_[data['vector'], 0j][model.ids].reshape(fieldmesh, -1, 4)
    numz = cells.shape[1]
    wav = np.zeros((fieldmesh+1, numz+1), complex)
    wav[:-1, :-1],wav[1:, :-1] = cells[:, :, 0],cells[:, :, 1]
    wav[:-1, 1:],wav[1:, 1:] = cells[:, :, 2],cells[:, :, 3]
    grid = (np.linspace(0, 1, numz+1), np.linspace(0, 1, fieldmesh+1))
    z, r = np.meshgrid(np.linspace(0, 1, 2*n+1), np.linspace(0, 1, n+1), indexing='ij')
    srcmode = RegularGridInterpolator(grid, wav.T)(np.c_[z.ravel(), r.ravel()])[supported]
    metric = reduce(mass)
    overlap = abs(modes.T @ (metric @ srcmode))**2 / np.vdot(srcmode, metric @ srcmode).real
    pick = int(np.argmax(overlap))
    freq = float(np.sqrt(vals[pick]))
    full = np.zeros(len(xy))
    full[supported] = modes[:, pick]
    def expectation(M):
        return float(full @ (M @ full))
    exact_axial = -freq * expectation(axial) / expectation(stiff)
    print(json.dumps(dict(source=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                          target=target, closed_nearby=np.sqrt(vals).tolist(), selected=freq,
                          overlap=overlap.tolist(), removed_empty_nodes=int(sum(boundary_free & ~supported)))), flush=True)
    for name, index in (('radial', 2), ('axial', 0)):
        energies = []
        for sign in (-1., 1.):
            at = ref.copy()
            at[index] += sign*.002
            shape = Shape(1, at)
            domain = arrays(shape, n)
            shifted_k, shifted_m, _ = wave_forms(domain, chi)
            omega = np.sqrt(expectation(shifted_k) / expectation(shifted_m))
            energies.append([omega, energy(domain, chi), duct.solve(shape, n//2)])
        dw, do, dc = (np.asarray(energies[1]) - energies[0]) / .004
        output = dict(coord=name, wave=dw, interior_order=do, weak_carrier=dc, total=dw+do+dc)
        if name == 'axial':
            output.update(exact_axial_wave=exact_axial, analytic_residual=dw-exact_axial)
            assert abs(exact_axial-dw) < 4e-6
        print(json.dumps(output), flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('state', choices=('mid', 'tight'))
    p.add_argument('n', type=int)
    a=p.parse_args()
    measure(a.source, a.state, a.n)
