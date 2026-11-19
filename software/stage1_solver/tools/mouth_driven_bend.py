"""Drive and release the previously proposed antisymmetric side coordinate.

The carrier control fixes its stationary inlet current. The odd-sector test
holds the reference depth with an independent static support and postulates
the two earlier bend responses at that same current. No wall damping is
inferred from the stationary flow. All frequencies are angular frequencies.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from mouth_dipole_bow_decay import generator
from mouth_trapped_wave_ports import ARCHIVE

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from stage1_solver.config import TensorGridSpec
from stage1_solver.open_sleeve import OpenSleeve


def saved_row(manifest, run, prefix):
    lines = manifest["runs"][run]["stdout_lines"]
    line, = [line for line in lines if line.startswith(prefix)]
    return dict(pair.split("=", 1) for pair in line.split() if "=" in pair)


def steady_current(archive, radius, length, current):
    p = archive["parameters"]
    for cells in (8, 16):
        # One open outlet. A section translated sideways preserves its area.
        # The grid supplies the axisymmetric carrier control only.
        grid = TensorGridSpec(2*radius, 2*cells, 0., length, 2*cells)
        solver = OpenSleeve(grid, cells, cells, hbar=p["hbar"], constituent_mass=p["m"],
                            eos_k=p["K_eos"], bath_density=p["rho0"]/p["m"])
        psi, phase, residual = solver.solve(current/p["m"])
        inlet, outlet, far, error = solver.budget(psi, phase)
        print(f"carrier cells={cells} mass_target={current:.12g} mass_in={p['m']*inlet:.12g} "
              f"mass_through={p['m']*outlet:.12g} mass_out={p['m']*far:.12g} "
              f"phase={phase:.12g} cell_error={error:.3g} residual={residual:.3g}", flush=True)
        assert np.allclose(p["m"]*np.array([inlet,outlet,far]), current, rtol=0, atol=1e-10)
        assert error < 1e-8 and residual < 1e-11


def integrate(rhs, period, init, dense=False):
    sol = solve_ivp(rhs, (0, period), init, method="DOP853", rtol=3e-13, atol=1e-15,
                    dense_output=dense)
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol


def drive_and_release(ws, wd, coupling, damping, frequency, force):
    """Use the actual periodic steady response, switch force off and sample b.

    A complex drive represents two real runs in quadrature. At phase zero,
    positive lag is -arg(b/F) for F(t)=Re[F exp(+i frequency t)].
    """
    mass, stiff = 12., 6.
    period = 2*np.pi/ws
    def matrix(t):
        return generator(t, ws, wd, coupling, mass, stiff, damping)
    def homogeneous(t, y):
        return (matrix(t) @ y.reshape(4,4)).ravel()
    monodromy = integrate(homogeneous, period, np.eye(4).ravel()).y[:,-1].reshape(4,4)
    def forced(t, y):
        dy = matrix(t) @ y
        dy[3] += force*np.exp(1j*frequency*t)/mass
        return dy
    period_phase = np.exp(1j*frequency*period)
    offset = integrate(forced, period, np.zeros(4,complex)).y[:,-1]
    start = np.linalg.solve(period_phase*np.eye(4)-monodromy, offset)
    sol = integrate(forced, period, start, dense=True)
    assert np.max(abs(sol.y[:,-1]-period_phase*start)) < 1e-10*force
    time = period*np.arange(128)/128
    envelope = sol.sol(time)[2] * np.exp(-1j*frequency*time)
    response = np.mean(envelope)/force
    bare = 1/(stiff-mass*frequency**2+1j*damping*frequency)
    print(f"drive C={damping:g} omega={frequency:.12g} force={force:g} "
          f"same_tone_displacement={force*abs(response):.12g} lag_deg={-np.angle(response,deg=True):.12g} "
          f"abs_chi={abs(response):.12g} bare_abs_chi={abs(bare):.12g} "
          f"bare_lag_deg={-np.angle(bare,deg=True):.12g}", flush=True)
    # Sample the switch-off wall trace at the same support phase. Extract the
    # four damped roots from the signal, with no supplied loss or resonance.
    state = start.copy()
    signal = []
    for _ in range(180):
        signal.append(state[2]/force)
        state = monodromy @ state
    signal = np.asarray(signal)
    def poles(data):
        history = np.lib.stride_tricks.sliding_window_view(data, 5)
        coeff, _, rank, _ = np.linalg.lstsq(history[:,:4], history[:,4], rcond=None)
        assert rank == 4
        return np.roots(np.r_[1., -coeff[::-1]])
    targets = []
    for count in (72, 180):
        roots = poles(signal[:count])
        # Unfold within half a pump frequency of the archived dipole.
        phase = np.angle(roots)/period
        unfolded = phase + ws*np.rint((wd-phase)/ws)
        i = np.argmin(abs(unfolded-wd))
        targets.append((-2*np.log(abs(roots[i]))/period, unfolded[i]))
    np.testing.assert_allclose(targets[0], targets[1], atol=2e-9, rtol=1e-8)
    gamma, resonance = targets[-1]
    assert gamma>0 and abs(resonance-wd)<.02
    print(f"switch_off C={damping:g} drive={frequency:.12g} energy_rate={gamma:.12g} "
          f"dipole_frequency={resonance:.12g}", flush=True)
    return response, gamma, resonance


def line_test(first, second, gamma, shift):
    """Conditional selective probes of adjacent coherences, equal integrated areas.

    Impose one cold Markov jump sqrt(gamma)*a, no pure phase noise and a common
    linear frequency shift. In the three-state basis D[a] gives diagonal
    coherence rates gamma/2 and 3*gamma/2. The 12 coherence also feeds the 01
    coherence; selective detection of 12 excludes that lower-line contribution.
    """
    centers = np.array([first, second])+shift
    halfwidths = .5*gamma*np.array([1.,3.])
    gap = first-second
    counts = []
    for points in (8001, 16001):
        span = max(gap, max(halfwidths))
        x = np.linspace(min(centers)-4*span, max(centers)+4*span, points)
        spectrum = np.sum(halfwidths/((x[:,None]-centers)**2+halfwidths**2), axis=1)
        peaks, _ = find_peaks(spectrum)
        counts.append(len(peaks))
    assert counts[0]==counts[1]
    print(f"selective_lines conditional_common_shift={shift:.12g} first_center={centers[0]:.12g} "
          f"second_center={centers[1]:.12g} first_FWHM={gamma:.12g} second_FWHM={3*gamma:.12g} "
          f"equal_area_peak_count={counts[-1]} peak_omegas={','.join(f'{v:.10g}' for v in x[peaks])}")


def main(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    archive = yaml.safe_load(ARCHIVE.read_text())
    point = saved_row(manifest, "smooth_dipole_l015_01", "point n=1 ")
    levels = saved_row(manifest, "smooth_dipole_l015_01", "levels ")
    bend = saved_row(manifest, "mouth_wave_linewidth_budget01", "bow h=0.015 ")
    current = saved_row(manifest, "mouth_wave_depth_response01", "archive ")
    ws, wd = float(point["omega_s"]), float(point["omega_d"])
    p = archive["parameters"]
    steady_current(archive, archive["phase_A"]["q_star"][0],float(point["L"]),
                   float(current["reference_inward_current"]))
    coupling = float(bend["overlap"])*np.sqrt(2*p["I_wave"]/ws)
    for damping in (.1, 6.5):
        fits = []
        for freq in (float(bend["difference"]), float(bend["sum"])):
            response, gamma, resonance = drive_and_release(ws,wd,coupling,damping,freq,.001)
            half = drive_and_release(ws,wd,coupling,damping,freq,.0005)
            np.testing.assert_allclose((response,gamma,resonance), half, atol=3e-9, rtol=2e-7)
            fits.append([gamma,resonance])
        np.testing.assert_allclose(fits[0],fits[1],atol=3e-9,rtol=2e-7)
        line_test(float(levels["relaxed_first"])/p["hbar"],
                  float(levels["relaxed_second"])/p["hbar"],fits[0][0],fits[0][1]-wd)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spectrum_manifest",type=Path)
    main(parser.parse_args().spectrum_manifest)
