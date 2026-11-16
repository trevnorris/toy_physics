"""Sideways carrier momentum across a stationary, phase-continuous mouth.

The thin real surface energy is a conditional example. A single feeding-side
trace is used for different receiving derivatives, without solving a receiver.
"""

import numpy as np
import sympy as sp

from held_crossflow_budget import archived_packet


def trace_balance():
    c, amp = sp.symbols("c u", positive=True)  # c = hbar^2/(2m)
    left, right, phase_left, phase_right = sp.symbols("uL uR qL qR", real=True)
    y, n, nx = sp.symbols("y n nx", real=True)
    energy = sp.Function("e")(y, n)
    def dx(value):
        return sp.diff(value, y) + nx * sp.diff(value, n)

    # Cancel the common phase in c*[d_N psi] = e_n*psi.
    derivative_jump = right - left + sp.I * amp * (phase_right - phase_left)
    condition = c * derivative_jump - sp.diff(energy, n) * amp
    energy_n = sp.Symbol("energy_n", real=True)
    real, imag = condition.subs(sp.diff(energy, n), energy_n).as_real_imag()
    match = sp.solve([real, imag], [right, phase_right], dict=True)[0]
    log_jump = sp.simplify((2 * (right - left) / amp).subs(match))
    log_jump = log_jump.subs(energy_n, sp.diff(energy, n))
    current_jump = sp.simplify((amp**2 * (phase_right - phase_left)).subs(match))
    # A common tangential phase derivative multiplies the normal current.
    # The remaining feeding-minus-receiving flux is the mixed density stress.
    drive = sp.simplify(c * n * dx(log_jump) / 2)
    edge_flux = n * sp.diff(energy, n) - energy
    translated = sp.simplify(drive - dx(edge_flux))
    # Build the full stress on each side from n, n_N and their tangential
    # derivatives. Vary both density and coupling along a trace inside its support.
    u_trace = 1 + y/7
    n_trace = u_trace**2
    trial_energy = (1-y**2)**3 * (n**2/3 + n**3/7 + y*n**2/9)
    trial_en = sp.diff(trial_energy, n).subs(n, n_trace)
    left_slope = y/11
    right_slope = left_slope + u_trace * trial_en/c
    carry = n_trace * (1+y/5) * (1-y/3)
    def normal_flux(slope):
        nz = 2*u_trace*slope
        return carry + c*(sp.diff(n_trace, y)*nz/n_trace - sp.diff(nz, y))/2

    stress_difference = normal_flux(left_slope) - normal_flux(right_slope)
    trial_drive = drive.subs(energy, trial_energy).doit().subs(
        {nx: sp.diff(n_trace, y), n: n_trace}, simultaneous=True)
    stress_residual = sp.factor(stress_difference - trial_drive)
    lam, g = sp.symbols("lambda g", real=True)
    restoring = lam * c * g * (sp.sqrt(n) - 1)**2
    amp_jump = sp.simplify(match[right] - left).subs(energy_n, sp.diff(restoring, n)).subs(n, amp**2)
    receiving = sp.factor(-g * (amp-1) + amp_jump)
    print("outgoing minus incoming normal number current, in units hbar/m =", current_jump)
    print("outgoing minus incoming normal log-density slope =", log_jump)
    print("feeding minus receiving x momentum flux =", drive)
    print("after subtracting tangential divergence =", translated)
    print("nonlinear two-sided stress minus trace formula =", stress_residual)
    print("receiving amplitude slope with restoring example =", receiving)
    print("receiving slopes at lambda=0,1 =", [receiving.subs(lam, k) for k in (0, 1)])
    if (current_jump != 0 or translated != sp.diff(energy, y) or stress_residual != 0
            or sp.simplify(receiving - (lam-1)*g*(amp-1)) != 0):
        raise ArithmeticError("Interface momentum balance failed")


def signed_work(speed):
    path = archived_packet("mouth_intake_redistribution.json", f"n32_restoring_slope_u{speed:g}.npz")
    with np.load(path, allow_pickle=False) as data:
        amplitude = abs(data["field"][-1])
        current, mouth = data["current"], data["mouth"]
        dx, drive = float(data["dx"]), float(data["speed"])
    x = (np.arange(amplitude.shape[0]) + 0.5 - amplitude.shape[0] / 2) * dx
    axial = x[:, None, None]
    transverse2 = x[None, :, None]**2 + x[None, None, :]**2
    def base(offset):
        return np.maximum(1 - (axial - offset)**2 - transverse2, 0)

    if (drive != speed or not np.array_equal(mouth, base(0) > 0)
            or not np.all(current[mouth] > 0) or np.any(current[~mouth] != 0)):
        raise ValueError("Unexpected mouth geometry, drive or inward current")
    # Odd part across opposite incident speeds on the same symmetric grid.
    # The reflected field is the negative-speed solution checked by the source readout.
    density_work = ((amplitude-1)**2 - (amplitude[::-1]-1)**2) / (4*speed)
    slope = dx**3 * np.sum(density_work * (-6*axial*base(0)**2))
    shift = 1e-4 * dx
    def energy(field, offset):
        return dx**3 * np.sum((field-1)**2 * base(offset)**3) / 2

    def signed_force(field):
        return (energy(field, -shift) - energy(field, shift)) / (2*shift)

    differenced = (signed_force(amplitude) - signed_force(amplitude[::-1])) / (2*speed)
    rate = dx**3 * np.sum(current)
    print(f"abs_U={speed:g} inward_rate={rate:.12g} restoring_drive_per_speed={slope:.12g} "
          f"translated_energy_check={differenced:.12g} fraction_of_continuum_incoming={slope/rate:.12g}",
          flush=True)
    np.testing.assert_allclose(slope, differenced, rtol=2e-6, atol=1e-13)
    np.testing.assert_allclose(rate, 0.25, rtol=0, atol=1e-12)


if __name__ == "__main__":
    trace_balance()
    for speed in (0.01, 0.02):
        signed_work(speed)
