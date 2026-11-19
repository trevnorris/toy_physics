"""Conditional addition spectrum of an empty dipole in the finite-side cavity.

Keep the archived radius, one unit of fundamental action and a reflected end.
The axial material coordinate can redistribute the load from the new mode.
An optional action-diagonal quantization of the axial motion tests line strengths.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal, solve_banded
from scipy.optimize import brentq

from mouth_cavity_occupation import radial, MI
from mouth_wave_depth_response import (RADIUS, REFERENCE_CURRENT, WALL_SLOPE, CAP_ENERGY,
                                      CFG, AXIAL, ACTION_STEP, V_PER_DEPTH)

S_WAVE, DIPOLE = [radial(ell, 1, RADIUS, True) for ell in (0,1)]
HBAR = ACTION_STEP
SUPPORT_ACTION = CFG["I_wave"]


def waves(length):
    return tuple(np.sqrt(MI**2 + row[0] + AXIAL/length**2) for row in (S_WAVE,DIPOLE))


def slope_for(sign):
    return WALL_SLOPE - sign*REFERENCE_CURRENT**2*RADIUS**2/CFG["rho0"]


def surface(length, n, sign):
    ws, wd = waves(length)
    return CAP_ENERGY + slope_for(sign)*length + SUPPORT_ACTION*ws + n*HBAR*wd


def optimum(n, sign):
    def derivative(z):
        ws, wd = waves(z)
        return slope_for(sign) - AXIAL/z**3*(SUPPORT_ACTION/ws + n*HBAR/wd)
    return brentq(derivative, .3, 5.)


def classical(sign):
    lengths = np.array([optimum(n,sign) for n in range(3)])
    spacings = np.diff([surface(z,n,sign) for n,z in enumerate(lengths)])
    ws,wd = waves(lengths[0])
    print(f"static sign={sign:+d} L0={lengths[0]:.12g} L1={lengths[1]:.12g} L2={lengths[2]:.12g} "
          f"support_frequency={ws:.12g} dipole_frequency={wd:.12g} rigid_step={HBAR*wd:.14g} "
          f"first={spacings[0]:.14g} second={spacings[1]:.14g} U={np.diff(spacings)[0]:.12g}")
    # Rigid radial motion is rotationally invariant; an odd longitudinal bow is
    # a separately imposed displacement epsilon*cos(theta)*sin(pi*z/L).
    # The step-potential matrix element uses normalized wavefunction values
    # at the radial side. Both angular harmonics have unit L2 norm.
    mixed = (2/(np.pi*np.sqrt(3))) * np.sqrt(S_WAVE[1]*DIPOLE[1])
    force_per_pair = HBAR*mixed/(2*np.sqrt(ws*wd))
    k_sound = []
    for w in (wd-ws,wd+ws,2*wd):
        cs2=5*CFG["K_eos"]*CFG["rho0"]**4
        beta=(CFG["hbar"]/(2*CFG["m"]))**2
        k_sound.append(np.sqrt(2*w**2/(cs2+np.sqrt(cs2**2+4*beta*w**2))))
    print(f"odd_bow sign={sign:+d} derivative_abs={mixed:.12g} "
          f"single_pair_force={force_per_pair:.12g} diff_freq={wd-ws:.12g} sum_freq={wd+ws:.12g} "
          f"twice_dipole={2*wd:.12g} sound_k_diff={k_sound[0]:.12g} sound_k_sum={k_sound[1]:.12g} "
          f"sound_k_twice_dipole={k_sound[2]:.12g}")
    assert all(lengths<CFG["d_slab"]) and wd<CFG["mu_out"]
    return spacings


def quantum(sign, cells, upper, diagnostics=False):
    # Metric ordering in one dimension: choose y so M(L) dL^2=dy^2.
    mass_slope=CFG["rho0"]*V_PER_DEPTH/4
    ymax=(2/3)*np.sqrt(mass_slope)*upper**1.5
    dy=ymax/(cells+1)
    y=dy*np.arange(1,cells+1)
    length=(1.5*y/np.sqrt(mass_slope))**(2/3)
    hop=-HBAR**2/(2*dy**2)
    off=np.full(cells-1,hop)
    diagonals=[surface(length,n,sign)-2*hop for n in range(3)]
    spectra=[eigh_tridiagonal(diag,off,select="i",select_range=(0,10))
             for diag in diagonals]
    energies=np.array([row[0][0] for row in spectra])
    steps=np.diff(energies)
    central=[]
    for n in (0,1):
        psi=spectra[n][1][:,0]
        weights=abs(spectra[n+1][1].T@psi)**2
        central.append(weights[0])
        if diagnostics:
            beyond_cap=np.sum(psi[length>CFG["d_slab"]]**2)
            print(f"sideband sign={sign:+d} start_n={n} zero_wall_excitation_weight={weights[0]:.12g} "
                  f"one_wall_excitation_weight={weights[1]:.12g} missing_weight_above_ten={1-sum(weights):.3g} "
                  f"wall_first_frequency={(spectra[n+1][0][1]-spectra[n+1][0][0])/HBAR:.12g} "
                  f"prob_L_above_optional_cap={beyond_cap:.12g}")
    print(f"dressed sign={sign:+d} cells={cells} L_limit={upper:g} first={steps[0]:.14g} "
          f"second={steps[1]:.14g} U={np.diff(steps)[0]:.12g} width_to_resolve={abs(np.diff(steps)[0])/HBAR:.12g}")
    susceptibility=[]
    if diagnostics:
        # Coupling G(a_d sigma+ + h.c.) independent of L. G here has energy units.
        # The declared instrumental width regularizes the *probe*, not the drain.
        width=HBAR*.001
        for detuning_freq in (.08,-.08):
            ep=steps[0]-HBAR*detuning_freq
            shifts=[]
            for n in (0,1):
                psi=spectra[n][1][:,0]
                channels=[]
                for adjacent, energy, factor in ((n-1,energies[n]-ep,n),(n+1,energies[n]+ep,n+1)):
                    if factor==0:
                        channels.append(0j)
                        continue
                    band=np.zeros((3,cells),complex)
                    band[1]=energy+1j*width-diagonals[adjacent]
                    band[0,1:]=-hop
                    band[2,:-1]=-hop
                    channels.append(factor * (psi @ solve_banded((1,1),band,psi)))
                shifts.append(sum(channels)/2)
            centroid=(shifts[1]-shifts[0]).real
            # The scalar rotating-wave contact is an independent hypothesis.
            # Isolate the transition leaving the geometry in its lowest state.
            d0=steps[0]-ep
            d1=steps[1]-ep
            zpl=central[0]*(1/(d0+1j*width)).real-central[1]*(1/(d1+1j*width)).real
            approx=1/d0-1/d1
            print(f"readout sign={sign:+d} detuning_freq={detuning_freq:g} width_freq=.001 "
                  f"full_shift_per_G2={centroid:.12g} lowest_lines_only={zpl:.12g} "
                  f"unit_overlap_zero_width={approx:.12g}")
            susceptibility.append(centroid)
    return np.r_[steps,central,susceptibility]


if __name__ == "__main__":
    for sign in (0,-1,1):
        classical(sign)
        coarse=quantum(sign,2000,6)
        fine=quantum(sign,4000,6,True)
        domain=quantum(sign,6200,8,True)
        err=np.max(abs(fine[:2]-coarse[:2]))
        ext=np.max(abs(domain-fine))
        print(f"refinement sign={sign:+d} max_step_difference={err:.3g} "
              f"domain_max_step_difference={np.max(abs(domain[:2]-fine[:2])):.3g} "
              f"domain_max_readout_difference={np.max(abs(domain[-2:]-fine[-2:])):.3g}")
        assert err<2e-7 and np.max(abs(domain[:2]-fine[:2]))<5e-8 and ext<.03
