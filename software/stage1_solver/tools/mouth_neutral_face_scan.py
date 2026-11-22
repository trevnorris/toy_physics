"""Scan a weak clamped neutral probe above a uniformly braced OpenCell bilayer.

The test constituent samples the existing unscaled OpenCell pair gradient.
The material order acts on the existing external bracing during preparation.
The scan reports upward holding force per unit probe coupling in patch units.
"""
import argparse
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from mouth_brace_grip import Material, evolve_interval
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_open_w import mix
from neutral_patch_sites import BETWEEN, GAP, TEMPERATURE, WITHIN, emit, morse
from neutral_patch_structural_endpoints import PITCH


HEIGHT = 1.1
READS_PER_TRANSIT = 16
BURN_TRANSITS = 4
TRANSITS_PER_BLOCK = 2
FIELDS = ('upper', 'lower', 'total', 'flat_upper', 'flat_lower', 'flat_total',
          'upper_rms_height', 'lower_rms_height', 'mean_face_gap')


def radial_integrals(a, cutoff, well):
    """Return integral of M'(r) dr and integral of r M(r) dr from a to cutoff."""
    depth, alpha, length = well

    def tails(r):
        w = np.exp(alpha * (length-r))
        energy = depth*w*(w-2)
        integral = depth * (w*w*(r/(2*alpha) + 1/(2*alpha)**2)
                            - 2*w*(r/alpha + 1/alpha**2))
        return energy, integral

    ma, ia = tails(a)
    mc, ic = tails(cutoff)
    return mc-ma, ia-ic


def area_force(gap, cell):
    """Uniform XY area average per site from the polar area element 2 pi r dr.

    dz is fixed during this integral. The native truncated force has no surface
    delta function at cutoff. Its radial integral retains M(cutoff).
    """
    a = np.minimum(np.abs(gap), cell.cutoff)
    da, ia = radial_integrals(a, cell.cutoff, WITHIN)
    db, ib = radial_integrals(a, cell.cutoff, BETWEEN)
    f, df, _ = mix(gap, cell.power)
    return 2*np.pi/cell.area * (gap*((1-f)*da + f*db) + df*(ib-ia))


def observe(pos, mat):
    lower, upper = np.split(pos[:, 2], 2)
    zprobe = upper.mean() + HEIGHT
    gap = zprobe-pos[:, 2]
    if gap.min() <= 0 or upper.min() <= lower.max():
        raise ArithmeticError('Faces crossed or probe plane intersects slab')
    f_low, f_up = (f.sum() for f in np.split(area_force(gap, mat.cell), 2))
    flat_low, flat_up = mat.half*area_force(np.array([zprobe-lower.mean(), HEIGHT]), mat.cell)
    return np.array([f_up, f_low, f_up+f_low, flat_up, flat_low, flat_up+flat_low,
                     upper.std(), lower.std(), upper.mean()-lower.mean()])


def direct_grid(pos, mat, count):
    """Frozen point-by-point sum over the periodic cell using native pair gradients."""
    cell = mat.cell
    gap = pos[mat.half:, 2].mean()+HEIGHT-pos[:, 2]
    f, df, _ = mix(gap, cell.power)
    # With nearest-image lateral displacements, this bound includes every tile
    # with any point within the spherical cutoff.
    extent = np.ceil(cell.cutoff/cell.box + .5).astype(int)
    tiles = [np.asarray(tile)*cell.box for tile in
             product(*(range(-s, s+1) for s in extent))
             if np.linalg.norm(np.maximum(np.abs(tile)-.5, 0)*cell.box) < cell.cutoff]
    xy = (np.indices((count, count)).reshape(2, -1).T+.5)*cell.box/count
    readings = []
    for start in range(0, len(xy), 96):
        nearest = xy[start:start+96, None, :] - pos[:, :2]
        nearest -= cell.box*np.rint(nearest/cell.box)
        forces = np.zeros(nearest.shape[:2])
        for tile in tiles:
            delta = nearest+tile
            rr = (delta*delta).sum(axis=2) + gap**2
            probe, site = np.nonzero(rr < cell.cutoff**2)
            r = np.sqrt(rr[probe, site])
            ea, ga = morse(r, *WITHIN)
            eb, gb = morse(r, *BETWEEN)
            grad = ((1-f[site])*ga + f[site]*gb)*gap[site]/r + df[site]*(eb-ea)
            forces[probe, site] += grad
        readings.append(np.column_stack([forces[:, mat.half:].sum(axis=1),
                                         forces[:, :mat.half].sum(axis=1)]))
    grid = np.vstack(readings)
    return grid.mean(axis=0), grid.sum(axis=1).std()


def check_grid(pos, mat, count, label, **fields):
    expected = observe(pos, mat)[:2]
    errors = []
    for n in (count//2, count):
        direct, spread = direct_grid(pos, mat, n)
        errors.append((direct-expected).tolist())
    emit(test='frozen_area_scan', state=label, grids=[count//2, count],
         upper_lower_area=expected.tolist(), fine_direct=direct.tolist(),
         direct_minus_area=errors, fine_lateral_force_sd=spread, **fields)
    np.testing.assert_allclose(direct, expected, atol=2.e-6, rtol=1.e-6)


def cold_registered(mat, chi):
    """Set the open-normal registered gap by normal force; check all site modes."""
    def state(gap):
        pos = mat.flat.copy()
        pos[:, 2] = mat.sign*gap
        return pos

    def gradient(pos):
        grad = mat.cell.evaluate(pos)[1]
        grad[:, 2] += chi*(mat.brace*(pos[:, 2]-mat.height)-mat.load)
        return grad

    gap = brentq(lambda d: gradient(state(d))[mat.half:, 2].mean(), .70, 1.20,
                 xtol=2.e-14)
    pos = state(gap)
    err = np.max(np.abs(gradient(pos)))
    hess = mat.cell.evaluate(pos, tangent=True)[0]
    hess[2::3, 2::3] += chi*mat.brace*np.eye(mat.n)
    eig = np.linalg.eigvalsh(hess)
    stable = bool(np.max(abs(eig[:2])) < 1.e-7 and eig[2] > 1.e-7)
    out = observe(pos, mat)
    emit(test='cold_registered', order=chi, max_force=err, locally_stable=stable,
         lowest_modes=eig[:4].tolist(), **dict(zip(FIELDS, out)))
    if err > 1.e-8:
        raise ArithmeticError('Cold branch has unresolved site forces')
    return pos


def warm_history(mat, args, seed, name, initial, orders, transit):
    rng = np.random.default_rng(seed)
    bridge = np.random.default_rng(seed+7013)
    states = [p.copy() for p in initial]
    schedules = [lambda times, c=c: np.full_like(times, c) for c in orders]
    burn = BURN_TRANSITS*transit
    evolve_interval(mat, states, 0., burn, schedules, rng, bridge, args.refine)
    old = burn
    blocks = []
    for block in range(2):
        readings = []
        for _ in range(TRANSITS_PER_BLOCK*READS_PER_TRANSIT):
            new = old+transit/READS_PER_TRANSIT
            evolve_interval(mat, states, old, new, schedules, rng, bridge, args.refine)
            readings.append([observe(p, mat) for p in states])
            old = new
        avg = np.mean(readings, axis=0)
        blocks.append(avg)
        emit(test='warm_block', seed=seed, history=name, block=block+1, end_time=old,
             order_values=orders, fields=FIELDS, means=avg.tolist(),
             half_to_entry=avg[1, 2]/avg[0, 2])
    for chi, pos in zip(orders, states, strict=True):
        check_grid(pos, mat, args.grid, 'last_warm_snapshot', seed=seed, history=name, order=chi)
    return np.mean(blocks, axis=0)


def group_summary(data):
    """Ratio of means; paired leave-one-independent-seed-out standard error."""
    data = np.asarray(data)
    means = data.mean(axis=0)

    def force_ratios(x):
        return x[1, :6]/x[0, :6]

    ratios = force_ratios(means)
    jack_se = None
    if len(data) > 1:
        jack = np.array([force_ratios((data.sum(axis=0)-v)/(len(data)-1)) for v in data])
        jack_se = np.sqrt((len(data)-1)*np.mean((jack-jack.mean(axis=0))**2, axis=0)).tolist()
    return dict(fields=FIELDS, means=means.tolist(), ratio_fields=FIELDS[:6],
                half_to_entry=ratios.tolist(), paired_seed_jackknife_se=jack_se)


def run(args):
    args.work.mkdir(parents=True, exist_ok=True)
    args.brace, args.temperature, args.mobility = 'elastic', TEMPERATURE, 1.
    trace, rate = cubic_trace(34)
    start, stop = 1.175, trace.crossing(.5)
    transit = trace.transit(start, stop, 8)
    orders = [float(trace.sample(trace.chi, np.array([start]))[0][0]), .5]
    np.testing.assert_allclose(orders[0], .875299070587, atol=1.e-11, rtol=0.)
    mat = Material(args)
    emit(test='setup', order_values=orders, order_action='external preload and ties only',
         probe='unscaled OpenCell gradient, zero back action', sign='upward holding',
         force_units='within-well energy / native length per unit weak coupling',
         height=HEIGHT, height_reference='instantaneous mean upper original face', side=args.side,
         cells_area=mat.cell.area, sites_per_face=mat.half, per_face_density=mat.half/mat.cell.area,
         wells=[WITHIN, BETWEEN], mix_power=args.power, cutoff=args.cutoff, mobility=mat.mobility,
         reference_gap=GAP, brace_stiffness=mat.brace, temperature=mat.temperature, seeds=args.seeds,
         carrier_intake_for_clock=rate, clock_depths=[start, stop], transit=transit,
         burn_transits=BURN_TRANSITS, read_block_transits=TRANSITS_PER_BLOCK, read_blocks=2,
         reads_per_transit=READS_PER_TRANSIT, dt_max=.001/args.refine)
    # A corrugation reaches the mixing interval; the signed derivative of the
    # mixing weight would be invisible in scans of an exactly flat upper face.
    calibration = mat.flat.copy()
    calibration[mat.half:, 2] += .30*np.cos(np.arange(mat.half)*2*np.pi/mat.half)
    check_grid(calibration, mat, args.grid, 'manufactured_corrugation')
    cold = [cold_registered(mat, chi) for chi in orders]
    cold_means = [observe(p, mat) for p in cold]
    emit(test='cold_ratio', **group_summary([cold_means]))
    samples = {'registered_cold': [], 'shifted_rough': []}
    for seed in args.seeds:
        rough = mat.flat.copy()
        rough[mat.half:, 0] += .25*PITCH
        rough[:, 2] *= .90
        rough += .025*np.random.default_rng(seed+23).normal(size=rough.shape)
        for name, initial in [('registered_cold', cold), ('shifted_rough', [rough,rough])]:
            samples[name].append(warm_history(mat, args, seed, name, initial, orders, transit))
        emit(test='history_read_difference', seed=seed, fields=FIELDS,
             shifted_minus_registered=(samples['shifted_rough'][-1]-samples['registered_cold'][-1]).tolist())
    for name, values in samples.items():
        emit(test='warm_summary', history=name, seeds=args.seeds, **group_summary(values))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--seeds', type=int, nargs='+', default=[1931, 2953, 3989, 5021])
    parser.add_argument('--power', type=int, choices=(1, 2), default=1)
    parser.add_argument('--side', type=int, choices=(4, 6), default=4)
    parser.add_argument('--refine', type=int, choices=(1, 2), default=1)
    parser.add_argument('--cutoff', type=float, default=4.5)
    parser.add_argument('--grid', type=int, default=48)
    args = parser.parse_args()
    if len(args.seeds) != len(set(args.seeds)) or args.grid < 24 or args.cutoff <= HEIGHT:
        parser.error('Use distinct seeds, grid >= 24 and cutoff above the probe height')
    run(args)
