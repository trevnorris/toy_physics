"""Local transverse momentum exchange in one material with two structural states.

Read mixed-order probe locations from the actual four-dimensional mouth. Test a
homogeneous patch in its material frame, using independent tag relaxation and
a finite driven chain. No relaxation time is inferred from the carrier.
"""
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse import bmat, diags, eye
from scipy.sparse.linalg import spsolve

from mouth_free_contact import Carrier
from mouth_free_order_diffuse import read


ROOT = Path("_scratch/s11c")
ARCHIVE = Path("research/pde_ledger_v3/_measurements/mouth_free_order_depth.json")


def emit(**values):
    print(json.dumps(values, allow_nan=False), flush=True)


def local_response():
    # Solve the two driven balances directly in complex amplitudes.
    q, rho, mu, b, omega = sp.symbols("q rho mu b omega", positive=True)
    wave_number_sq = sp.Symbol("k_squared")
    vo, vd = sp.symbols("v_o v_d")
    velocity_ratio = sp.solve(-sp.I*omega*rho*(1-q)*vd-b*(vo-vd), vd)[0] / vo
    mass = sp.factor(rho*(q+(1-q)*velocity_ratio))
    dispersion = -sp.I*omega*mass + mu*q*wave_number_sq/(-sp.I*omega)
    k2 = sp.solve(dispersion, wave_number_sq)[0]
    x = sp.symbols("x", positive=True)
    normalized = sp.simplify((mass/rho).subs(b, omega*rho*(1-q)/x))
    circle = (sp.re(normalized)-(1+q)/2)**2 + sp.im(normalized)**2 - ((1-q)/2)**2
    heat = b*abs(1-velocity_ratio)**2
    heat_residual = sp.simplify(heat-omega*sp.im(mass))
    effective_modulus = mu*q*rho/mass
    equivalence = sp.simplify(k2-omega**2*rho/effective_modulus)
    emit(test="symbolic", velocity_ratio=str(velocity_ratio), dynamic_mass=str(mass),
         wave_number_squared=str(k2), normalized_mass=str(normalized),
         circle_residual=str(sp.factor(circle)), heat_residual=str(heat_residual),
         constant_mass_modulus=str(sp.factor(effective_modulus)), dispersion_difference=str(equivalence))
    assert sp.factor(circle) == heat_residual == equivalence == 0


def read_patches():
    hashes = json.loads(ARCHIVE.read_text())
    for profile in ("cubic", "delayed"):
        for mesh in (34, 48):
            carrier_folder = ROOT / f"mouth_flared_{profile}_n{mesh}"
            carrier = Carrier(carrier_folder)
            carrier_file = carrier_folder/"field.npz"
            with carrier_file.open("rb") as stream:
                carrier_hash = file_digest(stream, "sha256").hexdigest()
            assert hashes["input_carriers_sha256"][str(carrier_file)] == carrier_hash
            with np.load(carrier_file) as data:
                carrier_keys = sorted(data.files)
            intake = json.loads((carrier_folder/"stdout").read_text().splitlines()[0])["bulk"]
            for state, seed in (("near", "0.44"), ("deep", "1.55")):
                folder = ROOT / f"mouth_free_order_{profile}_n{mesh}_d04"
                if (profile, mesh, state) == ("cubic", 48, "deep"):
                    folder = folder.with_name(folder.name+"_outer_residual")
                order_file = folder / f"order_seed{seed}.npz"
                with order_file.open("rb") as stream:
                    digest = file_digest(stream, "sha256").hexdigest()
                assert hashes["diffuse_orders_sha256"][str(order_file)] == digest
                with np.load(order_file) as data:
                    chi, order_keys = data["chi"], sorted(data.files)
                r = .27
                def order_at(z):
                    return read(carrier.domain, chi, r, z)[0]
                w = brentq(lambda z: order_at(z)-.25, .50001, 1.98)
                q, qr, qw = read(carrier.domain, chi, r, w)
                psi, pr, pw = read(carrier.domain, carrier.field, r, w)
                ur, uw = (pr/psi).imag, (pw/psi).imag
                material_rate = ur*qr+uw*qw
                emit(test="archive_patch", profile=profile, mesh=mesh, state=state,
                     intake=float(intake), r=r, w=w, chi=q, n=abs(psi)**2,
                     carrier_ur=ur, carrier_uw=uw, order_gradient_length=q/np.hypot(qr, qw),
                     carrier_order_change_time=q/abs(material_rate), carrier_D_chi=material_rate,
                     carrier_file=str(carrier_file), carrier_sha256=carrier_hash, carrier_keys=carrier_keys,
                     order_file=str(order_file), order_sha256=digest, order_keys=order_keys)
                assert abs(intake-.25) < 1e-8 and abs(q-.25) < 1e-9 and material_rate < 0


def calibrate_tags(q, drag):
    # A homogeneous opposite impulse has zero barycentric momentum. Infer
    # the time from the two velocities after removing the drives.
    def rates(t, velocities):
        force = drag*(velocities[0]-velocities[1])
        return -force/q, force/(1-q)
    times = np.array([.25, .75])
    out = solve_ivp(rates, (0., times[-1]), [1-q, -q], t_eval=times, rtol=1e-11, atol=1e-13)
    if not out.success:
        raise ArithmeticError(out.message)
    delta = out.y[0]-out.y[1]
    tau = np.diff(times)[0]/np.log(delta[0]/delta[1])
    emit(test="synthetic_tag_calibration", chi=q, times=times.tolist(),
         ordered_v=out.y[0].tolist(), disordered_v=out.y[1].tolist(),
         inferred_free_slip_time=tau, max_net_momentum=float(max(abs(q*out.y[0]+(1-q)*out.y[1]))))
    return tau


def chain(q, drag, omega, per_length):
    # Moving masses at every site, ordered sites linked by springs. Drive the
    # first ordered site harmonically and clamp the last. Keep the reflected
    # wave: a three-point recurrence identifies the interior wavenumber for
    # either mixture of traveling waves, away from a displacement node.
    cells = 12*per_length
    h = 1/(omega*per_length)
    size = cells-1
    spring = q/h**2
    neighbor = -spring*np.ones(size-1)
    ordered = diags([neighbor, np.full(size, 2*spring-omega**2*q-1j*omega*drag), neighbor], [-1, 0, 1])
    cross = 1j*omega*drag*eye(size)
    disordered = (-omega**2*(1-q)-1j*omega*drag)*eye(size)
    force = np.zeros(2*size, dtype=complex)
    force[0] = spring
    solution = spsolve(bmat([[ordered, cross], [cross, disordered]], format="csc"), force)
    ordered_readings = solution[:size]
    i = 1+np.argmax(abs(ordered_readings[1:size//3]))
    cosine = (ordered_readings[i-1]+ordered_readings[i+1])/(2*ordered_readings[i])
    fitted = np.arccos(cosine)/h
    if fitted.imag < 0:
        fitted = -fitted
    return fitted, solution[size+i]/solution[i]


def test_waves():
    # One illustrative collision time for the material, in local units with
    # rho=mu=1. Fix it at half order in a separate zero-wavenumber preparation.
    q0, gamma = .5, 1.
    tau = calibrate_tags(q0, gamma*q0*(1-q0))
    for q in (.25, .75):
        drag = gamma*q*(1-q)
        for x in (.1, 1., 10.):
            omega = x*q/tau
            tagged_ratio = 1/(1-1j*omega*tau/q)
            mass = q+(1-q)*tagged_ratio
            predicted = omega*np.sqrt(mass/q)
            cp = omega/predicted.real
            amplitude = np.exp(-2*np.pi*predicted.imag/predicted.real)
            errors = []
            for refinement in (24, 48, 96):
                measured, ratio = chain(q, drag, omega, refinement)
                error = abs(measured/predicted-1)
                emit(test="held_out_patch", chi=q, x=x, omega=omega, refinement=refinement,
                     predicted_k=[predicted.real, predicted.imag], fitted_k=[measured.real, measured.imag],
                     tag_velocity_ratio=[ratio.real, ratio.imag], relative_k_error=error,
                     phase_speed=cp, amplitude_after_one_wavelength=amplitude,
                     mass_real=mass.real, mass_imag=mass.imag)
                assert abs(ratio-tagged_ratio) < 1e-8 and measured.real > 0 and measured.imag > 0
                errors.append(error)
            assert errors[0] > errors[1] > errors[2] and errors[2] < 1e-4


if __name__ == "__main__":
    local_response()
    read_patches()
    test_waves()
