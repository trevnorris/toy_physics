"""Release a period-24 wrinkle and measure shear in the full two-layer patch."""
import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from neutral_patch_open_w import OpenCell
from neutral_patch_sites import emit
from neutral_patch_structural_endpoints import PITCH


class WrinkleCell(OpenCell):
    """Reuse the small-cell pair law with a geometric neighbor list."""
    def __init__(self, power, cutoff=6.):
        self.side = 24
        self.n = 2*self.side**2
        self.power = power
        self.box = self.side*PITCH*np.array([1., np.sqrt(3)/2])
        self.area = np.prod(self.box)
        self.cutoff = cutoff
        self.reference = None
        self.shifts = np.zeros(2)

    def evaluate(self, pos, shear=0., tangent=False):
        if abs(shear) > .001:
            raise ValueError('Neighbor radius covers only these small imposed strains')
        xy = pos[:, :2]
        if self.reference is None or np.linalg.norm(xy-self.reference, axis=1).max() > .25:
            self.reference = xy.copy()
            # In this box the cutoff is strictly inside half the shortest period.
            # Include .5 for relative pair motion and .001*cutoff for shear.
            radius = self.cutoff*1.001 + .51
            wrapped = xy % self.box
            wrapped = np.where(wrapped >= self.box, 0., wrapped)
            tree = cKDTree(wrapped, boxsize=self.box)
            self.i, self.j = tree.query_pairs(radius, output_type='ndarray').T
        return super().evaluate(pos, shear, tangent)

    def relax_sites(self, initial, shear=0.):
        def objective(flat):
            energy, gradient, _ = self.evaluate(flat.reshape(-1, 3), shear)
            return energy, gradient.ravel()
        fit = minimize(objective, initial.ravel(), jac=True, method='L-BFGS-B',
                       options=dict(ftol=3e-16, gtol=2e-10, maxiter=10000, maxls=70, maxcor=30))
        return fit.x.reshape(-1, 3)


def hollow_stack(cell):
    i, j = np.indices((cell.side, cell.side))
    xy = PITCH*np.column_stack(((i+.5*(j % 2)).ravel(), (j*np.sqrt(3)/2).ravel()))
    offset = PITCH*np.array([.5, 1/(2*np.sqrt(3))])
    # Independent zero-normal-load values for this flat hollow crystal.
    gap = {1: .7481782361436098, 2: .752300839061282}[cell.power]
    return np.column_stack((np.vstack((xy, xy+offset)), np.repeat([-.5*gap, .5*gap], len(xy))))


def initial_mode(cell, base, axis):
    h, _, _ = cell.evaluate(base, tangent=True)
    phase = 2*np.pi*base[:cell.n//2, axis]/cell.box[axis]
    cols = []
    for layer in range(2):
        for coordinate in range(3):
            for wave in (np.cos(phase), np.sin(phase)):
                v = np.zeros_like(base)
                v.reshape(2, -1, 3)[layer, :, coordinate] = wave/np.linalg.norm(wave)
                cols.append(v.ravel())
    basis = np.column_stack(cols)
    eigen, vec = eigh(basis.T@h@basis)
    lowest = basis@vec[:, eigen < eigen[0]+1e-7]
    # Choose the lowest mode with the crest initially at the coordinate origin.
    trial = np.zeros_like(base)
    trial[:, 2] = np.tile(np.cos(phase), 2)
    mode = lowest@(lowest.T@trial.ravel())
    mode /= np.linalg.norm(mode)
    np.testing.assert_allclose(h@mode, eigen[0]*mode, atol=1e-9)
    return eigen[0], mode.reshape(-1, 3)


def spectrum(cell, pos):
    h, mixed, born = cell.evaluate(pos, tangent=True)
    translations = np.tile(np.eye(3), (cell.n, 1))/np.sqrt(cell.n)
    np.testing.assert_allclose(h@translations, 0, atol=1e-9)
    np.testing.assert_allclose(translations.T@mixed, 0, atol=1e-8)
    # Shift just the three known translations above the internal spectrum.
    for dim in range(3):
        h[dim::3, dim::3] += 2000./cell.n
    eigen, modes = eigh(h, overwrite_a=True, driver='evd')
    return eigen, modes, mixed, born


def downhill_step(cell, pos, mode, energy):
    mode = mode/np.max(abs(mode))
    choices = [pos+s*a*mode for a in (.00003, .0001, .0003, .001, .003, .02, .1, .24)
               for s in (-1, 1)]
    # A shift of the wrinkle costs little even when moving sites along a
    # straight tangent to that shift stretches the bonds. Trial translations
    # of the periodic displacement profile follow its curved path. They only
    # seed the next unrestricted minimization if the actual pair energy drops.
    flat = hollow_stack(cell)
    disp = (pos-flat).reshape(2, cell.side, cell.side//2, 2, 3)
    # Each row parity has uniform reference x and y spacings.
    for axis in (1, 2):
        size = disp.shape[axis]
        coeff = np.fft.rfft(disp, axis=axis)
        shape = [1]*disp.ndim
        shape[axis] = size//2+1
        frequency = np.fft.rfftfreq(size).reshape(shape)
        for shift in (-.5, -.25, -.125, .125, .25, .5):
            moved = np.fft.irfft(coeff*np.exp(2j*np.pi*shift*frequency), n=size, axis=axis)
            choices.append(flat+moved.reshape(pos.shape))
    trial = min(choices, key=lambda q: cell.evaluate(q)[0])
    change = cell.evaluate(trial)[0]-energy
    emit(test='saddle_step', energy_change=change, peak_displacement=np.max(abs(trial-pos)))
    if change >= -8*abs(np.spacing(energy)):
        raise ArithmeticError('No resolved energy decrease out of the apparent saddle')
    return trial


def settle(cell, pos, folder):
    for attempt in range(10):
        pos = cell.relax_sites(pos)
        for polish in range(4):
            energy, grad, _ = cell.evaluate(pos)
            eigen, modes, mixed, born = spectrum(cell, pos)
            emit(test='stationarity', attempt=attempt, polish=polish, energy_per_site=energy/cell.n,
                 residual=np.max(abs(grad)), least=eigen[:5].tolist())
            np.savez(folder/'latest.npz', position=pos, power=cell.power)
            if eigen[0] < -2e-8:
                pos = downhill_step(cell, pos, modes[:, 0].reshape(-1, 3), energy)
                break
            if np.max(abs(grad)) < 2e-9:
                return pos, eigen, modes, mixed, born
            inverse = np.zeros_like(eigen)
            stable = eigen > 2e-8
            inverse[stable] = 1/eigen[stable]
            step = -(modes@(inverse*(modes.T@grad.ravel()))).reshape(-1, 3)
            if np.max(abs(step)) > .03:
                step *= .03/np.max(abs(step))
            pos = pos+step
        else:
            residual = np.max(abs(cell.evaluate(pos)[1]))
            if residual > 2e-7:
                raise ArithmeticError(f'Force persists after local refinement: {residual}')
    raise ArithmeticError('No stationary stable branch identified')


def run(power, folder, axis, checkpoint=None):
    cell = WrinkleCell(power)
    flat = hollow_stack(cell)
    e0, g0, _ = cell.evaluate(flat)
    value, wave = initial_mode(cell, flat, axis)
    peak = abs(wave).max()
    step = 1e-4/peak
    curv = sum(cell.evaluate(flat+s*step*wave)[0]-e0 for s in (-1, 1))/step**2
    emit(test='flat', power=power, axis=axis, sites=cell.n, force=np.max(abs(g0)),
         energy_per_site=e0/cell.n, wrinkle_curvature=value, direct_curvature=curv,
         normal_fraction=np.sum(wave[:, 2]**2))
    if value >= 0 or abs(curv-value) > .0002 or np.max(abs(g0)) > 1e-8:
        raise ArithmeticError('Initial state or unstable mode inconsistent')
    noise = np.random.default_rng(319).normal(size=flat.shape)
    start = flat+.08*wave/peak+.0005*(noise-noise.mean(axis=0))
    if checkpoint:
        with np.load(checkpoint) as saved:
            if int(saved['power']) != power or saved['position'].shape != flat.shape:
                raise ValueError('Checkpoint material or cell differs')
            start = saved['position'].copy()
        emit(test='continuation', input=str(checkpoint), energy_per_site=cell.evaluate(start)[0]/cell.n)
    pos, eigen, modes, mixed, born = settle(cell, start, folder)
    energy, gradient, preload = cell.evaluate(pos)
    stable = eigen > 2e-8
    projections = modes.T@mixed
    soft_force = np.linalg.norm(projections[~stable])
    if soft_force > 1e-6:
        raise ArithmeticError(f'Unresolved shear coupling to a soft mode: {soft_force}')
    inv = np.zeros_like(eigen)
    inv[stable] = 1/eigen[stable]
    response = -(modes@(inv*projections)).reshape(-1, 3)
    predicted = (born+mixed@response.ravel())/cell.area
    mid = (pos[:cell.n//2, 2]+pos[cell.n//2:, 2])/2
    np.savez(folder/'state.npz', position=pos, flat=flat, response=response, box=cell.box,
             power=power, axis=axis, soft_eigenvalues=eigen[:12])
    emit(test='wrinkled', power=power, reference_pair_midplane_span=np.ptp(mid),
         reference_pair_midplane_rms=np.std(mid), energy_change_per_site=(energy-e0)/cell.n,
         force_residual=np.max(abs(gradient)), baseline_shear=preload, least=eigen[:6].tolist(),
         soft_count=int(sum(~stable)), soft_shear_coupling=soft_force,
         affine=born/cell.area, relaxed=predicted)
    checks = []
    for eps in (1e-4, 5e-5):
        sides = []
        for sign in (-1, 1):
            shear = sign*eps
            q = cell.relax_sites(pos+shear*response, shear)
            # Use the unsheared tangent as a local preconditioner. Evaluate every
            # correction and residual with the actual strained nonlinear forces.
            for _ in range(6):
                es, gs, stress = cell.evaluate(q, shear)
                if np.max(abs(gs)) < 4e-9:
                    break
                dq = (modes@(inv*(modes.T@gs.ravel()))).reshape(-1, 3)
                if np.max(abs(dq)) > .02:
                    break
                q -= dq
            es, gs, stress = cell.evaluate(q, shear)
            sides.append((es, stress, np.max(abs(gs))))
        (em, sm, gm), (ep, sp, gp) = sides
        row = dict(strain=eps, energy_curvature=(ep+em-2*energy)/(eps**2*cell.area),
                   stress_slope=(sp-sm)/(2*eps), residual=max(gm,gp))
        checks.append(row)
        emit(test='direct_shear', **row)
    extended = WrinkleCell(power, 7.)
    ex, gx, sx = extended.evaluate(pos)
    emit(test='cutoff', energy_difference_per_site=(ex-energy)/cell.n, shear_difference=sx-preload,
         force_residual=np.max(abs(gx)))
    if energy >= e0 or np.max(abs(gx)) > 2e-7 or abs(ex-energy) > 2e-7 or abs(sx-preload) > 1e-7:
        raise ArithmeticError('Endpoint or cutoff check unresolved')
    for row in checks:
        np.testing.assert_allclose([row['energy_curvature'], row['stress_slope']], predicted,
                                   rtol=.002, atol=.005)
        if row['residual'] > 2e-7:
            raise ArithmeticError('Finite-strain equilibrium unresolved')


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('power', type=int, choices=[1, 2])
    cli.add_argument('folder', type=Path)
    cli.add_argument('--axis', choices=['x','y'], help='Direction along which the seed varies')
    cli.add_argument('--checkpoint', type=Path, help='Continue the recorded positions of this preparation')
    args = cli.parse_args()
    axis = args.power-1 if args.axis is None else ('x','y').index(args.axis)
    run(args.power, args.folder, axis, args.checkpoint)
