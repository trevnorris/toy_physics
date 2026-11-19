"""Controlled depth relaxation of the previously supported scalar fundamental.

At fixed mouth radius a commanded mass current gives constant axial load in
the archived kernel. Keep the work of that load in the addition energies.
The optional hard cap is tested separately. The response time uses the stored
collective mass at equilibrium and illustrates two stated damping choices.
"""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from mouth_cavity_occupation import S, R0, L0


CFG = S["parameters"]
RADIUS = R0
TRANSVERSE = CFG["mu_in"]**2 + (CFG["c_w"]*CFG["chi_perp"]/RADIUS)**2
AXIAL = (CFG["c_w"]*CFG["chi_DN"])**2
V_PER_DEPTH = 4*np.pi*RADIUS**3/3
WALL_SLOPE = CFG["alpha_volume"]*V_PER_DEPTH + CFG["sigma"]*4*np.pi*RADIUS**2
CAP_ENERGY = CFG["sigma"]*8*np.pi*RADIUS**3/3
REFERENCE_CURRENT = .1*.1*RADIUS**3/L0
ACTION_STEP = CFG["hbar"]
SLAB = CFG["d_slab"]


def frequency(length):
    return np.sqrt(TRANSVERSE + AXIAL/length**2)


def frequency_slope(length):
    return -AXIAL / (frequency(length)*length**3)


def curvature(length):
    w = frequency(length)
    return 3*AXIAL/(w*length**4) - AXIAL**2/(w**3*length**6)


def potential(length, action, axial_coefficient):
    return CAP_ENERGY + axial_coefficient*length + action*frequency(length)


def best_depth(action, slope, cap):
    if slope <= 0:
        return cap if np.isfinite(cap) else None
    root = brentq(lambda z: slope + action*frequency_slope(z), .001, 100.)
    return min(root, cap)


def response(label, slope, cap=np.inf):
    actions = CFG["I_wave"] + np.arange(3)*ACTION_STEP
    lengths = [best_depth(i, slope, cap) for i in actions]
    if None in lengths:
        print(f"shape model={label} axial_slope={slope:.12g} no_interior_minimum=1")
        return
    l = lengths[0]
    spacings = np.diff([potential(z, i, slope) for z, i in zip(lengths,actions)])
    print(f"shape model={label} axial_slope={slope:.12g} L0={l:.12g} L1={lengths[1]:.12g} L2={lengths[2]:.12g} "
          f"first={spacings[0]:.15g} second={spacings[1]:.15g} U={np.diff(spacings)[0]:.15g} "
          f"width_to_resolve={abs(np.diff(spacings)[0])/ACTION_STEP:.12g} "
          f"max_wave_frequency={max(frequency(z) for z in lengths):.12g}")
    assert np.diff(spacings)[0]<2e-13 and all(frequency(z)<CFG["mu_out"] for z in lengths)
    if l == cap:
        return
    stiffness = actions[0]*curvature(l)
    mass = CFG["rho0"]*V_PER_DEPTH*l/4
    omega_mech = np.sqrt(stiffness/mass)
    print(f"static model={label} stiffness={stiffness:.12g} axial_mass={mass:.12g} mechanical_frequency={omega_mech:.12g} "
          f"curvature_in_action={-frequency_slope(l)**2/stiffness:.12g}")
    for damping in (0., CFG["c_damp_sample"]):
        matrix = np.array([[0.,1.], [-stiffness/mass,-damping/mass]])
        # Normalize the small step in radiation force so the static displacement is one.
        def displacement(t):
            return 1-expm(t*matrix)[0,0]
        half_time = brentq(lambda t: displacement(t)-.5, 0., np.pi/omega_mech)
        t_wave = 2*np.pi/frequency(l)
        print(f"step model={label} damping={damping:g} first_half_time={half_time:.12g} "
              f"one_wave_period={t_wave:.12g} fraction_after_period={displacement(t_wave):.12g}")


def main():
    effective_mass_factor = RADIUS**2 / CFG["rho0"]
    print(f"archive reference_inward_current={REFERENCE_CURRENT:.12g} zero_action_open_endpoint=0 "
          f"trial_action_step={ACTION_STEP} held_radius={RADIUS:.12g}")
    response("rest", WALL_SLOPE)
    for sign in (-1, 1):
        slope = WALL_SLOPE - sign * REFERENCE_CURRENT**2 * effective_mass_factor
        response(f"controlled_{sign:+}", slope)
    cap_on = np.sqrt((WALL_SLOPE+CFG["I_wave"]*frequency_slope(SLAB))/effective_mass_factor)
    zero_slope = np.sqrt(WALL_SLOPE/effective_mass_factor)
    print(f"threshold at_action=1 cap_on_current={cap_on:.12g} no_finite_uncapped_above={zero_slope:.12g}")
    # Exercise both sides of the cap using the archive's mass-current units.
    # No identification with the previous unit-density radius-one drain is made.
    illustrative = (cap_on+zero_slope)/2
    slope = WALL_SLOPE - illustrative**2*effective_mass_factor
    response("above_cap_unrestricted", slope)
    response("above_cap_hard_slab", slope, SLAB)
    print(f"held both_first={ACTION_STEP*frequency(L0):.15g} both_second={ACTION_STEP*frequency(L0):.15g} U=0")


if __name__ == "__main__":
    main()
