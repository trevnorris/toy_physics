"""Force on a held staggered row feature, per unit steady flux of fresh mass.

A cold uniform supply first crosses the stationary 1.175--1.55 strip. At the
gauge, copy each actual (unwrapped transverse position, velocity) into a local
1.55--1.65 continuation. Every row within interaction range of that continuation
belongs to its fixture, including rows on the upstream side of the section.
This is a two-stage injection experiment. No step, kick, or upstream force is
charged to the local holding load. Row phases are read independently of tags.
"""
import argparse
from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_mobile_lining import ROW_SPACING, SUPPORT, row_formation, row_potential
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, GAP, emit, morse


START, GAUGE, EXIT = 1.175, 1.55, 1.65


@dataclass
class Track:
    t: np.ndarray
    depth: np.ndarray
    x: np.ndarray
    chi: np.ndarray


def tracks(trace, pitch, stride):
    depth = np.unique(np.r_[np.linspace(START, EXIT, 8193), GAUGE,
                             trace.edges[(trace.edges > START) & (trace.edges < EXIT)]])
    clock = cumulative_trapezoid(1/trace.speed(depth), depth, initial=0.)
    arc = cumulative_trapezoid(np.hypot(1., trace.geometry.radius(depth)[1]), depth, initial=0.)/pitch
    rows = int(np.floor(arc[-1]/ROW_SPACING)) + 1
    dt = clock[-1]/(rows*stride)
    time_to_depth = PchipInterpolator(clock, depth)
    time_to_x = PchipInterpolator(clock, arc)
    times = clock[np.searchsorted(depth, [START, GAUGE, EXIT])]
    parts = []
    for a, b, lo, hi in zip(times[:-1], times[1:], [START, GAUGE], [GAUGE, EXIT], strict=True):
        t = np.linspace(a, b, int(np.ceil((b-a)/dt))+1)
        z = time_to_depth(t)
        z[[0, -1]] = [lo, hi]
        chi = trace.sample(trace.chi, z)[0]
        np.testing.assert_allclose(b-a, trace.transit(lo, hi, 12), rtol=2e-6)
        if np.any((chi < 0) | (chi > 1)):
            raise ValueError('Unphysical shear-order profile')
        parts.append(Track(t-a, z, time_to_x(t), chi))
    return rows, parts


def pair_tables(p, grid):
    """Force derivatives for one periodic row, per transverse lattice pitch."""
    u = np.arange(grid+1)[None, :]/grid
    x = np.arange(SUPPORT*(grid//2)+1)[:, None]/(grid//2)
    tangent = row_potential(x, u, p, True)[1]
    # The positive gap coordinate points from the rows to the stream. Multiply
    # the native dU/dgap by p for the derivative per normal lattice pitch.
    normal = np.zeros_like(tangent)
    for image in range(-10, 11):
        r = np.sqrt(p*p*(x*x+(u-image)**2)+GAP**2)
        normal += morse(r, *BETWEEN)[1]*(p*GAP/r)
    window = np.clip(SUPPORT-x, 0, 1)
    taper = window**2*(3-2*window)
    np.testing.assert_allclose(tangent, -tangent[:, ::-1], atol=3e-13)
    np.testing.assert_allclose(normal, normal[:, ::-1], atol=3e-13)
    return tangent*taper, normal*taper


def row_origin(wall, indices):
    """Phase and magnitude of the wall's staggered Fourier density component."""
    component = np.exp(2j*np.pi*(wall[indices]-.5*(indices % 2))).mean()
    return np.angle(component)/(2*np.pi), abs(component)


def active_rows(x, n):
    return np.flatnonzero(np.abs(x-np.arange(n)*ROW_SPACING) < SUPPORT)


def field(track, n, tables, scale, resolve_rows):
    """Propagator force and, locally, separate wall row pair-force profiles.

    Table coordinates are relative to the staggered origin. Row j is displaced
    by j%2 half-periods. Bilinear interpolation commutes with summing the
    periodic force on this rigid branch.
    """
    tangent, normal = tables
    grid = tangent.shape[1]-1
    ny = grid//2
    total = np.empty((len(track.t), grid+1))
    samples = []
    for k, (x, chi) in enumerate(zip(track.x, track.chi, strict=True)):
        ids = active_rows(x, n)
        axial = np.abs(x-ids*ROW_SPACING)*ny
        i = np.minimum(axial.astype(np.intp), tangent.shape[0]-2)
        frac = (axial-i)[:, None]

        def row_forces(table):
            native = table[i, :-1] + frac*(table[i+1, :-1]-table[i, :-1])
            odd = (ids % 2 == 1)
            native[odd] = np.roll(native[odd], ny, axis=1)
            return scale*chi*native

        tan = row_forces(tangent)
        total[k, :-1] = -tan.sum(axis=0)
        total[k, -1] = total[k, 0]
        if resolve_rows:
            samples.append((ids, tan, row_forces(normal)))
    return total, samples


def shape(y, origin, grid):
    u = (y-origin) % 1.*grid
    i = np.minimum(u.astype(np.intp), grid-1)
    return i, u-i


def gather(force, i, f):
    return force[i] + f*(force[i+1]-force[i])


def developed_inlet(track, force, count, offset):
    """Follow half the cold midpoint pairs, then reflect each entire trajectory."""
    y = (np.arange(count//2)+offset)/count
    v = np.zeros_like(y)
    dt = track.t[1]
    kick = gather(force[0], *shape(y, 0., force.shape[1]-1))
    for grad in force[1:]:
        v += .5*dt*kick
        y += dt*v
        kick = gather(grad, *shape(y, 0., force.shape[1]-1))
        v += .5*dt*kick
    return np.r_[y, -y], np.r_[v, -v]


def registry(y, v, wall_y):
    angle = 2*np.pi*(y-wall_y)
    wrapped = (y-wall_y+.5) % 1.-.5
    return dict(phase_cos=np.mean(np.cos(angle)), phase_sin=np.mean(np.sin(angle)),
                fraction_within_eighth=np.mean(np.abs(wrapped) <= .125),
                mean_slip=np.mean(v), rms_slip=np.sqrt(np.mean(v*v)),
                mean_slip_times_phase_sin=np.mean(v*np.sin(angle)))


def passage(track, force, row_samples, wall, inlet_y, inlet_v):
    """Holding impulse by row-force quadrature, and a separate tag momentum sum."""
    y, v = inlet_y.copy(), inlet_v.copy()
    grid = force.shape[1]-1
    used = np.unique(np.concatenate([s[0] for s in row_samples]))
    origin, coherence = row_origin(wall, used)
    ref = int(np.rint(track.x[0]/ROW_SPACING))
    last = int(np.rint(track.x[-1]/ROW_SPACING))
    last = min(last, len(wall)-1)
    inlet = registry(y, v, wall[ref])
    dt = track.t[1]
    wall_impulse = np.zeros((2, len(wall)))
    normal_t = np.empty(len(track.t))
    for k, (ids, tangent, normal) in enumerate(row_samples):
        i, f = shape(y, origin, grid)
        kick = gather(force[k], i, f)
        # Scatter the ensemble's shape weights to each phase node. Contract with
        # the raw forces on individual rows before summing. This row-side readout
        # uses no particle velocity or momentum difference.
        left = np.bincount(i, weights=1.-f, minlength=grid)
        right = np.bincount(i, weights=f, minlength=grid)
        occupancy = (left+np.roll(right, 1))/len(v)
        row_tan = tangent @ occupancy
        row_normal = normal @ occupancy
        normal_t[k] = row_normal.sum()
        weight = .5*dt if k in (0, len(track.t)-1) else dt
        wall_impulse[0, ids] += weight*row_tan
        wall_impulse[1, ids] += weight*row_normal
        if k:
            v += .5*dt*kick
        if k < len(track.t)-1:
            v += .5*dt*kick
            y += dt*v

    holder = -wall_impulse[0].sum()
    particle_gain = np.mean(v-inlet_v)
    result = dict(inlet_count=len(inlet_v), outlet_count=len(v), inlet=inlet,
                  outlet=registry(y, v, wall[last]), wall_stagger_phase=origin,
                  wall_stagger_magnitude=coherence, feature_row=ref, feature_row_phase=wall[ref] % 1.,
                  inlet_mean_momentum=np.mean(inlet_v), outlet_mean_momentum=np.mean(v),
                  rms_unwrapped_slip_distance=np.sqrt(np.mean((y-inlet_y)**2)),
                  holding_tangent=holder, incident_force_on_rows=-holder,
                  particle_momentum_gain=particle_gain, closure=holder-particle_gain,
                  wall_gap_normal=wall_impulse[1].sum(), active_rows=used.tolist(),
                  incident_tangent_by_row=wall_impulse[0, used].tolist(),
                  gap_normal_by_row=wall_impulse[1, used].tolist())
    return result, normal_t


def controls(results, q, speed, translated):
    """Emit paired operands, raw zero levels and control residuals before checks."""
    zero, normal0 = results['shared', 0., 0.]
    case, common = translated
    original, orig_n = results['fixed', q, -speed]
    emit(test='global_translation', gauge_translation_pitches=case['global_shift'],
         original_holding=original['holding_tangent'], translated_holding=case['holding_tangent'],
         holding_difference=case['holding_tangent']-original['holding_tangent'],
         original_normal=original['wall_gap_normal'], translated_normal=case['wall_gap_normal'],
         max_normal_history_difference=np.max(abs(common-orig_n)))
    pairs = []
    for (supply, offset, drift), (row, normal) in results.items():
        if offset < 0 or (offset == 0 and drift <= 0):
            continue
        mirror, mirrored_normal = results[supply, -offset, -drift]
        normal_diff = np.max(np.abs(normal-mirrored_normal))
        pairs.append((row['holding_tangent']+mirror['holding_tangent'], normal_diff))
        emit(test='sign_pair', supply=supply, offset=offset, inlet_drift=drift,
             reflected_holdings=[row['holding_tangent'], mirror['holding_tangent']],
             holding_sum=pairs[-1][0], reflected_normals=[row['wall_gap_normal'], mirror['wall_gap_normal']],
             max_normal_history_difference=normal_diff,
             max_normal_history_change_from_zero=np.max(np.abs(normal-normal0)),
             integrated_normal_change_from_zero=row['wall_gap_normal']-zero['wall_gap_normal'])
    for amplitude in (q, q/2):
        for supply in ('fixed', 'translated'):
            plus = results[supply, amplitude, 0.][0]['holding_tangent']
            minus = results[supply, -amplitude, 0.][0]['holding_tangent']
            emit(test='offset_response', supply=supply, offset_pitches=amplitude,
                 holdings=[plus, zero['holding_tangent'], minus], slope=(plus-minus)/(2*amplitude))
    for v in (speed, speed/2):
        plus = results['shared', 0., v][0]['holding_tangent']
        minus = results['shared', 0., -v][0]['holding_tangent']
        emit(test='drift_response', inlet_drift=v, holdings=[plus, zero['holding_tangent'], minus],
             slope=(plus-minus)/(2*v))
    tolerance = 2e-10
    all_rows = [r for r, _ in results.values()]+[case]
    if not np.all(np.isfinite([r['closure'] for r in all_rows])):
        raise ArithmeticError('Nonfinite reaction or momentum')
    if max(abs(r['closure']) for r in all_rows) > tolerance:
        raise ArithmeticError('Row-side reaction and particle momentum ledger disagree')
    if max(max(abs(h), n) for h, n in pairs) > tolerance or abs(zero['holding_tangent']) > tolerance:
        raise ArithmeticError('Reflected force or zero control failed')
    if (abs(case['holding_tangent']-original['holding_tangent']) > tolerance
            or np.max(np.abs(common-orig_n)) > tolerance):
        raise ArithmeticError('Common translation of the local continuation changed the reaction')


def run(args):
    if args.grid < 8 or args.grid % 2 or args.phases < 4 or args.phases % 2 or args.stride < 1:
        raise ValueError('Use an even grid and phase count, and a positive stride')
    if not (0 < args.q < .25 and 0 < args.phase_offset < 1):
        raise ValueError('Offset drive must be in (0, 1/4) pitches and quadrature offset in (0, 1)')
    tr, intake = cubic_trace(args.mesh)
    patch, _, _ = read_inputs(args.mesh)
    p = patch['pitch']
    pitch = p*PITCH
    scale = 2*patch['shared_scale']/pitch**2
    n, (prep, local) = tracks(tr, pitch, args.stride)
    relaxed = row_formation(n, p, 191)
    stagger = .5*(np.arange(n) % 2)
    origin = np.mean(relaxed-stagger)
    error = np.max(np.abs(relaxed-stagger-origin))
    emit(test='wall_relaxation_snap', phase=origin, max_snap_pitches=error)
    if error > 2e-5:
        raise ValueError('Staggered branch did not relax')
    wall = origin+stagger
    tables = pair_tables(p, args.grid)
    upstream, _ = field(prep, n, tables, scale, False)
    relative_y, inlet_v = developed_inlet(prep, upstream, args.phases, args.phase_offset)
    del upstream
    local_force, row_samples = field(local, n, tables, scale, True)
    q = args.q
    speed = q/local.t[-1] if args.drift is None else args.drift
    if not np.isfinite(speed) or speed <= 0:
        raise ValueError('Use a positive finite inlet drift')
    emit(test='setup', mesh=args.mesh, intake=intake, grid=args.grid, axial_grid=args.grid//2,
         stride=args.stride, phases=args.phases, reflected_half_cell_offset=args.phase_offset,
         upstream_steps=len(prep.t)-1, downstream_steps=len(local.t)-1,
         upstream_time=prep.t[-1], downstream_time=local.t[-1], depth=[START, GAUGE, EXIT],
         chi_at_sections=[prep.chi[0], local.chi[0], local.chi[-1]], gap_native=GAP,
         pitch=pitch, rows=n, q_pitches=q, drift_pitches_per_time=speed,
         represented_mass_flux=1., weight_per_inlet_phase=1./args.phases,
         represented_mass_inventory=local.t[-1],
         normalization='holding force per unit steady fresh mass flux, in pitches per carrier time',
         normal_sign='positive gap-normal points from the wall toward passing material',
         feature='all rigid rows in reach after 1.55 including overlap across injection section')
    results = {}
    def measure(supply, offset, drift, global_shift=0.):
        shift = offset if supply == 'translated' else 0.
        inlet_y = (relative_y+origin)+global_shift+shift
        moved_wall = wall+global_shift+offset
        row, normal = passage(local, local_force, row_samples, moved_wall, inlet_y, inlet_v+drift)
        row.update(supply=supply, offset=offset, inlet_drift=drift, global_shift=global_shift)
        emit(test='local_feature_traction', **row)
        return row, normal

    for drift in (0., speed, -speed, speed/2, -speed/2):
        results['shared', 0., drift] = measure('shared', 0., drift)
    for supply in ('fixed', 'translated'):
        for offset in (q, -q, q/2, -q/2):
            results[supply, offset, 0.] = measure(supply, offset, 0.)
        for offset, drift in ((q, -speed), (-q, speed)):
            results[supply, offset, drift] = measure(supply, offset, drift)
    translated = measure('fixed', q, -speed, .273)
    controls(results, q, speed, translated)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, choices=(26, 34), required=True)
    parser.add_argument('--phases', type=int, required=True)
    parser.add_argument('--stride', type=int, required=True, help='time steps per row on the full strip')
    parser.add_argument('--grid', type=int, required=True, help='phase nodes per pitch, twice the axial density')
    parser.add_argument('--q', type=float, default=.002, help='static fixture displacement in pitches')
    parser.add_argument('--drift', type=float, help='inlet velocity step in pitches/time, default q/deep transit')
    parser.add_argument('--phase-offset', type=float, default=.5,
                        help='node fraction of each positive-phase cell, then reflected')
    run(parser.parse_args())
