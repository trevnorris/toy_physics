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
    lam, g = sp.symbols("lambda g", real=True)
    restoring = lam * c * g * (sp.sqrt(n) - 1)**2
    amp_jump = sp.simplify(match[right] - left).subs(energy_n, sp.diff(restoring, n)).subs(n, amp**2)
    receiving = sp.factor(-g * (amp-1) + amp_jump)
    print("outgoing minus incoming normal number current, in units hbar/m =", current_jump)
    print("outgoing minus incoming normal log-density slope =", log_jump)
    print("feeding minus receiving x momentum flux =", drive)
    print("after subtracting tangential divergence =", translated)
    print("receiving amplitude slope with restoring example =", receiving)
    print("receiving slopes at lambda=0,1 =", [receiving.subs(lam, k) for k in (0, 1)])
    if current_jump != 0 or translated != sp.diff(energy, y) or receiving.subs(lam, 1) != 0:
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

    if drive != speed or not np.array_equal(mouth, base(0) > 0) or not np.all(current[mouth] > 0):
        raise ValueError("Unexpected mouth geometry, drive or inward current")
    # Odd part across opposite incident speeds on the same symmetric grid.
    # The reflected field is the negative-speed solution checked by the source readout.
    density_work = ((amplitude-1)**2 - (amplitude[::-1]-1)**2) / (4*speed)
    slope = dx**3 * np.sum(density_work * (-6*axial*base(0)**2))
    shift = 1e-4 * dx
    differenced = dx**3 * np.sum(density_work * (base(-shift)**3-base(shift)**3)) / (2*shift)
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
