"""Constituent-fraction viability of the reversible wave/flow/order law in a uniform patch."""
import json
import numpy as np
from scipy.integrate import solve_ivp


def solve(relative_energy, phase, steps):
    # Two transverse standing profiles proportional to cos(kx) and sin(kx)
    # in perpendicular polarizations make their combined energy density uniform.
    # A spatially uniform carrier flow drops out in the material frame.
    width = .04
    inertia = width**2
    well = 3/width
    k = 8.
    stored = well*relative_energy
    if phase == "kinetic":
        q0, p0 = 0., np.sqrt(2*stored)
    else:
        q0, p0 = np.sqrt(2*stored)/k, 0.
    def rates(t, y):
        chi, pi, q, p = y
        return [pi/inertia, -2*well*chi*(1-chi)*(1-2*chi)+.5*((p/chi)**2-(k*q)**2),
                p/chi, -chi*k*k*q]
    state = np.array([1., 0., q0, p0])
    end = 2*np.pi/np.sqrt(2*well/inertia)
    integration = solve_ivp(rates, (0., end), state, method="DOP853", rtol=2e-11, atol=1e-14,
                            t_eval=np.linspace(0, end, 1001), max_step=end/steps)
    if not integration.success:
        raise ArithmeticError(integration.message)
    chi, pi, q, p = integration.y
    energy = pi**2/(2*inertia)+well*(chi*(1-chi))**2 + .5*(p*p/chi + chi*(k*q)**2)
    return dict(phase=phase, wave_energy_over_well=relative_energy, k=k,
                order_period=end, wave_period=2*np.pi/k, initial_chi=1.,
                initial_chi_acceleration=rates(0,state)[1]/inertia,
                min_chi_quartic_continuation=min(chi), max_chi_quartic_continuation=max(chi),
                relative_energy_error=max(abs(energy/energy[0]-1.)), max_step=end/steps)


def main():
    for e in (.001, .002):
        for phase in ("kinetic", "strain"):
            fine = solve(e, phase, 64)
            coarse = solve(e, phase, 32)
            print(json.dumps(fine), flush=True)
            if (fine["relative_energy_error"]>1e-8 or
                    abs(coarse["max_chi_quartic_continuation"]-fine["max_chi_quartic_continuation"])>1e-9):
                raise ArithmeticError("Uniform Hamiltonian control failed")


if __name__ == "__main__":
    main()
