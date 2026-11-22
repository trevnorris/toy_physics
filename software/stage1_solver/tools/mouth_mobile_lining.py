"""Local mobile lining and fresh constituent crossings on the radius-one cubic mouth."""
import argparse
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from mouth_segment_contact_hold import PITCH
from neutral_patch_sites import BETWEEN, GAP, WITHIN, emit, morse


ROW_SPACING = np.sqrt(3) / 2
SUPPORT = 5  # axial cutoff in measured nearest-neighbor pitches
NY = 256
NX_PER_PITCH = 96


def row_potential(x, y, pitch, cross):
    """Energy per site and derivative w.r.t. transverse separation in pitches."""
    energy = np.zeros(np.broadcast_shapes(np.shape(x), np.shape(y)))
    derivative = np.zeros_like(energy)
    for image in range(-10, 11):
        dy = y - image
        r = np.sqrt(pitch**2 * (x*x + dy*dy) + (GAP**2 if cross else 0))
        term, radial = morse(r, *(BETWEEN if cross else WITHIN))
        energy += term
        derivative += radial * pitch**2 * dy / r
    return energy, derivative


def row_formation(n, pitch, seed):
    # Only relative row phases relax. Longitudinal rows and their transverse
    # repeat period define the candidate material and its crystallographic axis.
    initial = np.arange(n) % 2 / 2 + np.random.default_rng(seed).uniform(-.18, .18, n)
    i, j = np.triu_indices(n, k=1)
    keep = (j-i)*ROW_SPACING < SUPPORT
    i, j = i[keep], j[keep]
    axial = (j-i)*ROW_SPACING

    def energy(y):
        dy = y[j] - y[i]
        terms, dy_grad = row_potential(axial, dy - np.rint(dy), pitch, False)
        grad = np.bincount(j, weights=dy_grad, minlength=n) - np.bincount(i, weights=dy_grad, minlength=n)
        return terms.sum(), grad

    initial_energy = energy(initial)[0]
    sol = minimize(energy, initial, jac=True, method='L-BFGS-B',
                   options={'ftol': 1e-15, 'gtol': 1e-10, 'maxiter': 5000, 'maxcor': 30})
    y = sol.x
    coherence = lambda v: abs(np.exp(2j*np.pi*(v - np.arange(n)%2/2)).mean())
    residual = np.max(abs(energy(y)[1]))
    emit(test='row_formation', seed=seed, initial_coherence=coherence(initial),
         relaxed_coherence=coherence(y), energy_decrease_per_row=(initial_energy-sol.fun)/n,
         max_native_phase_force=residual, optimizer=sol.message)
    if not sol.success or residual > 1e-5:
        raise ValueError('Row phases did not relax')
    return y


def run(args):
    trace, intake = cubic_trace(args.mesh)
    patch, _, _ = read_inputs(args.mesh)
    pitch_native = patch['pitch']
    pitch = pitch_native*PITCH
    lo, gauge, hi = 1.175, 1.55, 1.65
    w = np.unique(np.r_[np.linspace(lo, hi, 8193), gauge, trace.edges[(trace.edges>lo)&(trace.edges<hi)]])
    clock = cumulative_trapezoid(1/trace.speed(w), w, initial=0.)
    arc = cumulative_trapezoid(np.hypot(1., trace.geometry.radius(w)[1]), w, initial=0.)
    np.testing.assert_allclose(clock[-1], trace.transit(lo, hi, 12), rtol=1e-6)
    n = int(np.floor(arc[-1]/pitch/ROW_SPACING))+1
    steps = n*args.stride
    times = np.linspace(0, clock[-1], steps+1)
    depths = PchipInterpolator(clock, w)(times)
    x = PchipInterpolator(clock, arc)(times)/pitch
    chi = trace.sample(trace.chi, depths)[0]
    gauge_index = int(np.argmin(abs(depths-gauge)))
    gauge_row = int(np.rint(x[gauge_index]/ROW_SPACING))
    np.testing.assert_array_less(-1e-10, chi)
    np.testing.assert_array_less(chi, 1.00000001)
    y = row_formation(n, pitch_native, 191)
    # A cold relaxed row is a controllable symmetry branch. Keep its arbitrary
    # mean phase, and remove only sub-tolerance optimizer errors for mirrored runs.
    stagger = np.arange(n)%2/2
    phase = (y-stagger).mean()
    snap = np.max(abs(y-stagger-phase))
    if snap > 2e-5:
        raise ValueError('Independent relaxation did not select the staggered branch')
    xgrid = np.arange(SUPPORT*NX_PER_PITCH+1)[:,None]/NX_PER_PITCH
    ygrid = np.arange(NY+1)[None,:]/NY
    dist = np.arange(1, int(SUPPORT/ROW_SPACING)+1)[:,None]*ROW_SPACING
    cross_table = row_potential(xgrid, ygrid, pitch_native, True)[1]
    same_table = row_potential(dist, ygrid, pitch_native, False)[1]
    for table in (cross_table, same_table):
        np.testing.assert_allclose(table, -table[:,::-1], atol=3e-13)
    # Window the last axial interval in the passing-site force continuously.
    # The window depends only on prescribed axial position, preserving reciprocity in y.
    cutoff = np.clip(SUPPORT - xgrid, 0, 1)
    cross_table *= cutoff**2*(3-2*cutoff)
    emit(test='strip_setup', mesh=args.mesh, intake=intake, transit=clock[-1], rows=n,
         length_pitches=x[-1], dt=times[1], phases=args.phases, stride=args.stride,
         gauge_depth=depths[gauge_index], gauge_row=gauge_row, gauge_order=chi[gauge_index],
         inlet_order=chi[0], exit_order=chi[-1], pitch=pitch, mass_scale=2*patch['shared_scale'],
         max_snap_pitches=snap, native_pitch=pitch_native, relaxed_global_phase=phase)
    arrays = b''.join(np.asarray(a, dtype=np.float64).tobytes() for a in (x, chi, cross_table, same_table))
    source = Path(__file__).with_suffix('.cpp')
    with tempfile.TemporaryDirectory() as directory:
        binary = Path(directory)/'lining'
        subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra',str(source),'-o',str(binary)], check=True)
        for case in args.cases:
            cmd = [str(binary), n, args.stride, args.phases, NY, NX_PER_PITCH, SUPPORT,
                   len(dist), times[1], 2*patch['shared_scale']/pitch**2, gauge_index, phase, case]
            subprocess.run(list(map(str,cmd)),input=arrays,check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, required=True, choices=[26,34])
    parser.add_argument('--stride', type=int, required=True)
    parser.add_argument('--phases', type=int, required=True)
    parser.add_argument('--cases', nargs='+', required=True,
                        choices=['free','edge','reverse','half','rest','smooth','shift','boost','edge-half','edge-reverse'])
    run(parser.parse_args())
