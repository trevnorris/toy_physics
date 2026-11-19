"""Flat torsional slab: solve a guided branch and integrate positive energy and power.

Both x and w reference slips in the homogeneous receiver obey the same passive
law. The core and its welded faces are a specification, in units of lattice
pitch. One transverse internal angular polarization, no intake or acoustic load.
"""
from math import sqrt

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, root

from mouth_torsion_bulk_power import C, RHO, H, W0, material_modulus

WIDTH = 64.
HALF = WIDTH/2


def surface(k, omega, density, drag):
    modulus = material_modulus(0, omega, drag)
    p = np.sqrt(complex(RHO*omega**2/C-k**2))
    beta = np.sqrt(complex(k**2-density*omega**2/modulus))
    # U(0)=1, U'(0)=0. Match the actual force on the two welded materials.
    residual = -C*p*np.sin(p*HALF) + modulus*beta*np.cos(p*HALF)
    return p, beta, modulus, residual


def lossless_root(omega, density):
    core_k = omega*sqrt(RHO/C)
    modulus = material_modulus(0, omega, 0).real
    if modulus > 0:
        outer_k = omega*sqrt(density/modulus)
        if outer_k >= core_k:
            return None
        # Nodeless oscillatory interior.
        lo = max(outer_k, sqrt(max(0.,core_k**2-(np.pi/WIDTH)**2)))
        hi = core_k
    else:
        # Nodeless interior growing toward either surface. The matching
        # ratio tends to one at large k; subresonant M is negative.
        if abs(modulus) >= C:
            return None
        lo, hi = core_k, max(2*core_k, .05)
        while surface(hi, omega, density, 0)[3].real < 0:
            hi *= 2
    return brentq(lambda k: surface(k, omega, density, 0)[3].real,
                  lo, hi, xtol=1e-14)


def solve_lossy(omega, density, drag, initial):
    def pair(parts):
        f = surface(complex(*parts), omega, density, drag)[3]
        return [f.real, f.imag]
    sol = root(pair, [initial.real, initial.imag], tol=1e-11)
    k = complex(*sol.x)
    assert sol.success and np.linalg.norm(pair(sol.x)) < 1e-10 and k.imag > 0
    return k


def inventory(k, omega, density, drag):
    p, beta, modulus, mismatch = surface(k, omega, density, drag)
    assert beta.real > 0 and abs(mismatch)<1e-10
    core_norm = 2*quad(lambda w: abs(np.cos(p*w))**2, 0, HALF)[0]
    core_grad = 2*quad(lambda w: abs(p*np.sin(p*w))**2, 0, HALF)[0]
    boundary = abs(np.cos(p*HALF))**2
    tails = boundary / beta.real  # both half-spaces
    shear_tail = (abs(k)**2+abs(beta)**2)*tails
    slip = C/(C-H*omega**2-1j*drag*omega)
    energy_core = ((RHO*omega**2+C*abs(k)**2)*core_norm+C*core_grad)/4
    energy_bulk = (density*omega**2*tails+
                   (H*omega**2*abs(slip)**2+C*abs(1-slip)**2)*shear_tail)/4
    power_core = omega*C*k.real*core_norm/2
    power_bulk = omega*(modulus*k).real*tails/2
    power = power_core + power_bulk
    heat = drag*omega**2*abs(slip)**2*shear_tail/2
    normal_both = -omega*(modulus*beta).imag*boundary
    assert abs(heat-2*k.imag*power) < 1e-9 * max(1e-8, abs(power))
    assert abs(normal_both+2*k.imag*power_bulk-heat) < 1e-9 * max(1e-8, abs(power))
    return (energy_bulk/(energy_core+energy_bulk), power_bulk/power, normal_both,
            power/(energy_core+energy_bulk))


def main():
    low_speed = sqrt(C*WIDTH/(RHO*WIDTH+2*sqrt(RHO*H)))
    print(f"two_surface_added_angular_inertia={2*sqrt(RHO*H):.9f} "
          f"predicted_low_speed={low_speed:.12f} bare_speed={sqrt(C/RHO):.12f}")
    for ratio in (.01, .05, .5, .8, 2., 3.):
        omega = ratio*W0
        for density in ((RHO, 2*RHO) if ratio == 2 else (RHO,)):
            k = lossless_root(omega, density)
            if k is None:
                print(f"rho_bulk={density:g} omega/omega0={ratio:g} "
                      "no_lossless_nodeless_bound_solution_in_continuum_model")
                continue
            energy, power, normal, velocity = inventory(k,omega,density,0)
            delta = omega*1e-5
            group = 2*delta/(lossless_root(omega+delta,density)-lossless_root(omega-delta,density))
            print(f"rho_bulk={density:g} omega/omega0={ratio:g} k={k:.12f} "
                  f"phase_speed={omega/k:.9f} group_speed={group:.9f} "
                  f"work_per_energy={velocity:.9f} energy_bulk_fraction={energy:.9f} "
                  f"axial_bulk_power_fraction_signed={power:.9f} normal_power={normal:.3g}")
            assert abs(velocity-group)<1e-7 and 0<energy<1 and abs(normal)<1e-15
            if ratio == .01:
                assert abs(omega/k-low_speed)<1e-4
            if ratio == 2:
                k_complex = complex(k)
                for drag in (2.,8.,32.):
                    k_complex = solve_lossy(omega,density,drag,k_complex)
                    e,p,n,_ = inventory(k_complex,omega,density,drag)
                    print(f"rho_bulk={density:g} omega/omega0=2 drag={drag:g} "
                          f"Re_k={k_complex.real:.9f} Im_k={k_complex.imag:.9f} "
                          f"through_one_wavelength_power={np.exp(-4*np.pi*k_complex.imag/k_complex.real):.9f} "
                          f"two_face_power_per_unit_length={n:.9g} bulk_energy_fraction={e:.9f}")
    # Independent aperture control for a heavier positive-modulus receiver:
    # impose the bare slab tangential phase and measure the OUTGOING normal
    # power on one prescribed moving face; this is not a leaky mode lifetime.
    omega = 2*W0
    k = omega*sqrt(RHO/C)
    modulus = material_modulus(0,omega,0).real
    for density in (RHO,2*RHO):
        q = np.sqrt(complex(density*omega**2/modulus-k**2))
        outward = omega*modulus*q.real/2
        print(f"unit_face_motion omega/omega0=2 conserved_bare_k={k:.9f} rho_bulk={density:g} "
              f"q_real={q.real:.9f} q_imag={q.imag:.9f} normal_power={outward:.12g}")
        assert (outward>0) == (density>RHO)


if __name__ == "__main__":
    main()
