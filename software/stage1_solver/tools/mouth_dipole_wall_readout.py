"""Line response to axial loading, and access to an open odd bulk channel.

The wall-side and eigenmode normalizations correspond to the separable
finite radial well with a reflected axial end. Coordinates are outward
displacements measured in length units. Rate coefficients use a cold bath.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import spherical_jn, spherical_kn

from mouth_dipole_spectroscopy import (RADIUS, HBAR, SUPPORT_ACTION, S_WAVE, DIPOLE,
                                       CFG, AXIAL, waves, surface, optimum)
from mouth_wave_depth_response import WALL_SLOPE


def normalized_side(ell, mode):
    # Normalize a full radial step-well eigenfunction on the interior and its tail.
    x = np.sqrt(mode[0])*RADIUS/CFG["c_w"]
    y = np.sqrt(CFG["mu_out"]**2-CFG["mu_in"]**2-mode[0])*RADIUS/CFG["c_w"]
    inner = quad(lambda u: (u*spherical_jn(ell,x*u)/spherical_jn(ell,x))**2,0,1)[0]
    outer = quad(lambda u: (u*spherical_kn(ell,y*u)/spherical_kn(ell,y))**2,1,np.inf)[0]
    return (RADIUS**3*(inner+outer))**-.5


def loaded_energy(n, force):
    def gradient(length):
        ws,wd = waves(length)
        return WALL_SLOPE - force - AXIAL/length**3*(SUPPORT_ACTION/ws+n*HBAR/wd)
    depth = brentq(gradient,.2,5)
    return surface(depth,n,0)-force*depth


def main():
    length=optimum(0,0)
    ws,wd=waves(length)
    first=loaded_energy(1,0)-loaded_energy(0,0)
    second=loaded_energy(2,0)-loaded_energy(1,0)
    separation=abs(second-first)/HBAR
    shift=optimum(1,0)-length
    for load in (.001,.0005):
        upper=loaded_energy(1,load)-loaded_energy(0,load)
        lower=loaded_energy(1,-load)-loaded_energy(0,-load)
        inferred=(upper-lower)/(2*load)
        print(f"test_force magnitude={load:g} first_plus={upper:.15g} first_minus={lower:.15g} "
              f"odd_energy_slope={inferred:.12g} independent_depth_prediction={-shift:.12g} "
              f"absolute_error={abs(inferred+shift):.3g}")
        assert abs(inferred+shift)<1e-6
    # z integral with axial ground function sqrt(2/L) sin(pi z/(2L)).
    z_overlap=quad(lambda z: 2/length*np.sin(np.pi*z/(2*length))**2*np.sin(np.pi*z/length),0,length)[0]
    angular=quad(lambda u: np.sqrt(3)*u*u/2,-1,1)[0]
    uniform=quad(lambda u: np.sqrt(3)*u/2,-1,1)[0]
    derivative=(CFG["mu_out"]**2-CFG["mu_in"]**2)*RADIUS**2*z_overlap*angular
    derivative*=normalized_side(0,S_WAVE)*normalized_side(1,DIPOLE)
    f=HBAR*derivative/(2*np.sqrt(ws*wd))
    reference_force=HBAR*(2/(np.pi*np.sqrt(3)))*np.sqrt(S_WAVE[1]*DIPOLE[1])/(2*np.sqrt(ws*wd))
    print(f"matrix uniform_overlap={uniform:.12g} axial_overlap={z_overlap:.12g} angular_overlap={angular:.12g} "
          f"mixed_force={f:.12g} independent_error={abs(f-reference_force):.3g}")
    assert abs(f-reference_force)<1e-13 and abs(uniform)<1e-14
    # A supported s mode with mean n_s=I_s/hbar assumes the action is normal ordered.
    # Distinct positive frequencies, weak coupling and a stationary cold bow.
    mean_s=SUPPORT_ACTION/HBAR
    cminus=2/HBAR*f**2*(mean_s+1)
    cplus=2/HBAR*f**2*mean_s
    print(f"cold_rate delta_frequency={wd-ws:.12g} sum_frequency={wd+ws:.12g} "
          f"coeff_down_to_monopole={cminus:.12g} coeff_pair_down={cplus:.12g} "
          f"classical_line_gap_frequency={separation:.12g} "
          f"inverse_line_gap_time={1/separation:.12g} "
          f"Q_if_equal_width={wd/separation:.12g}")
    print(f"isolated_path_bench if_ImChi_equal={separation/(cminus+cplus):.12g} "
          f"if_only_difference={separation/cminus:.12g} if_only_sum={separation/cplus:.12g}")
    # Delta omega from rms length noise in a held cavity is g_L delta L.
    slope=abs(-AXIAL/(wd*length**3))
    print(f"frozen_length line_slope={slope:.12g} length_jitter_for_one_gap={separation/slope:.12g}")


if __name__ == "__main__":
    main()
