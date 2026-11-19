"""Same-mode pressure and spectral thresholds of the archived standing wave.

The sharp D/N variational frequency and the frozen smooth two-ended gate are
distinct experiments evaluated at the same held (a,L). The latter is solved
on R^3 x R in the s-wave and w-even sector via u=r A. Positive ground-state
symmetry places the actual lowest mode in this sector. All frequencies are
angular. An action quantum for the optional scalar is a trial calibration.
"""
import argparse
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.special import expit


ARCHIVE = Path(__file__).resolve().parents[1] / "reports/pathA_26_derrick_results.yaml"


def mode(geometry, constants, step, half_extent):
    radius, length = geometry
    cells = round(half_extent / step)
    assert np.isclose(cells * step, half_extent)
    r = (1 + np.arange(cells - 1)) * step
    w = (.5 + np.arange(cells)) * step
    # Radial grid includes the u(0)=0 and outer Dirichlet edges implicitly.
    # Axial finite-volume ghost: even at w=0, zero at the outer cell face.
    nr, nw = len(r), len(w)
    radial = sparse.diags((-np.ones(nr-1), 2*np.ones(nr), -np.ones(nr-1)), (-1,0,1))
    center = 2*np.ones(nw)
    center[0], center[-1] = 1, 3
    axial = sparse.diags((-np.ones(nw-1), center, -np.ones(nw-1)), (-1,0,1))
    kinetic = constants["c_w"]**2 / step**2 * sparse.kronsum(axial, radial, format="csc")
    gr = expit(2*(radius-r)/constants["delta_r"])
    gw = expit((length-2*np.sqrt(w*w+constants["eps_w"]**2))/constants["delta_parallel"])
    depth = constants["mu_out"]**2 - constants["mu_in"]**2
    gate = np.outer(gr,gw)
    potential = constants["mu_out"]**2 - depth*gate
    operator = kinetic + sparse.diags(potential.ravel(), format="csc")
    seed = (r[:,None] * np.exp(-r[:,None]**2-w[None,:]**2)).ravel()
    levels, states = eigsh(operator,k=1,sigma=0.,which="LM",ncv=12,tol=1e-11,v0=seed)
    frequency = np.sqrt(levels[0])
    state = states[:,0]
    weight = state.reshape(nr,nw)**2
    dg_radius = np.outer(2*gr*(1-gr)/constants["delta_r"],gw)
    dg_length = np.outer(gr,gw*(1-gw)/constants["delta_parallel"])
    slopes = -depth/(2*frequency)*np.array([np.sum(weight*dg) for dg in (dg_radius,dg_length)])
    residual = np.linalg.norm(operator@state-levels[0]*state)
    return frequency,slopes,np.sum(weight*gate),residual


def ports(label, frequency, slopes, params):
    trial = params["hbar"]
    pressure = -params["I_wave"]*slopes
    pair = -trial/np.sqrt(2)*slopes
    # This is an on-shell kinematic bulk check at homogeneous rho=rho0.
    sound2 = 5*params["K_eos"]*params["rho0"]**4/params["m"]
    dispersion = (params["hbar"]/(2*params["m"]))**2
    print(f"ports mode={label} omega={frequency:.12g} scalar_gap={params['mu_out']:.12g} "
          f"evanescent_kappa={np.sqrt(params['mu_out']**2-frequency**2)/params['c_w']:.12g} "
          f"trial_single_energy={trial*frequency:.12g} trial_double_energy={2*trial*frequency:.12g} "
          f"pair_force_a={pair[0]:.12g} pair_force_L={pair[1]:.12g} "
          f"coherent_Fa_DC_and_at_2omega={pressure[0]:.12g} coherent_FL_DC_and_at_2omega={pressure[1]:.12g}",flush=True)
    for harmonic in (1,2,3):
        omega = harmonic*frequency
        k2 = 2*omega**2/(sound2+np.sqrt(sound2**2+4*dispersion*omega**2))
        relative = abs(sound2*k2+dispersion*k2**2-omega**2)/omega**2
        kg = np.sqrt(max(0.,omega**2-params["mu_out"]**2))/params["c_w"]
        print(f"threshold mode={label} harmonic={harmonic} omega={omega:.12g} "
              f"scalar_propagates={int(omega>params['mu_out'])} "
              f"scalar_k_if_open={kg:.12g} homogeneous_GNLS_k={np.sqrt(k2):.12g} "
              f"GNLS_dispersion_residual={relative:.3g}",flush=True)
        assert relative<5e-14
    assert frequency < params["mu_out"] < 2*frequency


def main(archive):
    data = yaml.safe_load(archive.read_text())
    p = data["parameters"]
    q = np.array(data["phase_A"]["q_star"])
    chi = np.array([np.pi,np.pi/2])
    sharp = np.sqrt(p["mu_in"]**2+sum((p["c_w"]*chi/q)**2))
    slopes = -(p["c_w"]*chi)**2/(sharp*q**3)
    ports("sharp_fixed_mass_surrogate",sharp,slopes,p)
    baseline = None
    for step, extent in ((.12,7.2),(.06,7.2),(.03,7.2),(.015,7.2),(.06,9.0)):
        omega, grad, participation, residual = mode(q,p,step,extent)
        print(f"smooth step={step:g} extent={extent:g} omega={omega:.12g} "
              f"domega_da={grad[0]:.12g} domega_dL={grad[1]:.12g} "
              f"gate_overlap={participation:.12g} eigen_residual={residual:.3g}",flush=True)
        assert residual<3e-7
        if step==.06 and extent==7.2:
            baseline = np.r_[omega,grad]
            delta = .0002
            for index in range(2):
                displacement = delta*np.eye(2)[index]
                left = mode(q-displacement,p,step,extent)[0]
                right = mode(q+displacement,p,step,extent)[0]
                secant = (right-left)/(2*delta)
                print(f"derivative coord={index} HF={grad[index]:.12g} "
                      f"separate_eigenfrequency_slope={secant:.12g} error={abs(secant-grad[index]):.3g}",flush=True)
                assert abs(secant-grad[index]) < 2e-5
        if step==.03:
            coarse = np.r_[omega,grad]
        if step==.015:
            difference = abs(np.r_[omega,grad]-coarse)
            print(f"refinement last_frequency_change={difference[0]:.8g} "
                  f"last_max_slope_change={max(difference[1:]):.8g}",flush=True)
            ports("smooth_fixed_gate_finest",omega,grad,p)
        if extent==9.0:
            change=abs(np.r_[omega,grad]-baseline)
            print(f"domain_change step={step:g} omega_difference={change[0]:.8g} "
                  f"max_slope_difference={max(change[1:]):.8g}",flush=True)
            assert max(change)<1e-5


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--archive",type=Path,default=ARCHIVE)
    main(parser.parse_args().archive)
