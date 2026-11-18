"""Signed, time-dependent neutral-carrier response between the entrance and first cut.

The prescribed straight receiver and both uniform-phase baths follow the same
translation. Its baseline is a previously solved open-sleeve field. A(t) grows
as exp(s*t) from zero in the remote past. All reported changes are first order.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import splu

from open_sleeve_lateral_balance import bath_quantum, cut_currents, side_stress
from stage1_solver.config import TensorGridSpec
from stage1_solver.open_sleeve import OpenSleeve


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "_scratch/s11c/sleeve_response_current_fine01/fields"
HASHES = {
    64: "cc30eecbf1bf148809a229e9825e187e21968baaa05f72f8941e2549d5390bad",
    128: "03531f6219273a88eedbcbd0170345871e90594a6702e63df347ed5b7fa5a9b4",
}


def background(level):
    path = BASE / f"level_{level}.npz"
    if hashlib.sha256(path.read_bytes()).hexdigest() != HASHES[level]:
        raise ValueError(f"Unrecognized steady field at {path}")
    with np.load(path, allow_pickle=False) as data:
        nr, nw = (int(data[key].max()) + 1 for key in ("r_index", "w_index"))
        grid = TensorGridSpec(nr=nr, nw=nw, r_max=nr*data["dr"], w_min=0, w_max=nw*data["dw"])
        domain = OpenSleeve(grid, int(data["mouth_cells"]), int(data["sleeve_cells"]),
                            *(float(data[k]) for k in ("hbar", "constituent_mass", "eos_k", "bath_density")))
        if not all(np.array_equal(getattr(domain, key), data[key]) for key in ("r_index", "w_index")):
            raise ValueError("Cell ordering mismatch")
        field = data["psi"]
        rates = domain.budget(field, float(data["feed_phase"]))[:3]
        if np.max(np.abs(np.array(rates) - 0.125)) > 1e-10:
            raise ValueError("Incorrect baseline current")
    return domain, field, rates[0]


def linear_operator(domain, psi, growth_rate):
    # In the translating frame: i hbar s eta = D(H psi) eta + m A r psi.
    volume = domain.volume
    density = abs(psi)**2
    h = 5*domain.eos_k/4
    angular = domain.kinetic*8*np.pi*domain.dr*domain.dw
    base = domain.operator + diags(angular + volume*h*(density**4 - domain.bath_density**4))
    slope = 8*volume*h*density**3
    mixed = slope*psi.real*psi.imag
    temporal = domain.hbar*growth_rate*volume
    return bmat([[base + diags(slope*psi.real**2), diags(mixed + temporal)],
                 [diags(mixed - temporal), base + diags(slope*psi.imag**2)]], format="csc")


def radial_derivative(values, spacing, *, odd):
    left = -values[:1] if odd else values[:1]
    padded = np.concatenate((left, values, values[-1:]), axis=0)
    return (padded[2:] - padded[:-2])/(2*spacing)


def measure(domain, psi, perturbation, rate, s, acceleration):
    mask = domain.w_index < domain.sleeve_cells
    shape = (domain.mouth_cells, domain.sleeve_cells)
    field = psi[mask].reshape(shape)
    eta = perturbation[mask].reshape(shape)
    volume = domain.volume[mask].reshape(shape)
    end = domain.sleeve_cells//4
    if not end or 4*end != domain.sleeve_cells:
        raise ValueError("Quarter cut must be a grid face")
    m, hb, dr, dw = domain.constituent_mass, domain.hbar, domain.dr, domain.dw
    r = dr*(np.arange(shape[0])[:, None] + 0.5)
    density = abs(field)**2
    logn = np.log(density)
    phi = hb/m*np.unwrap(np.unwrap(np.angle(field), axis=0), axis=1)
    b, g = 2*(eta/field).real, hb/m*(eta/field).imag
    flow, stress = cut_currents(logn, phi, b, g, z=end, spacing=dw, dr=dr, m=m, hbar=hb)
    side = sum(side_stress(logn, b, 0, end, dr=dr, dw=dw, m=m, hbar=hb, eos_k=domain.eos_k))
    feed = bath_quantum(field, eta, domain.bath_density, dr=dr, dw=dw, hbar=hb, m=m, at_exit=False)

    # Angular average of hbar Im(conj(psi) partial_x psi) to first order.
    grad_eta = radial_derivative(eta, dr, odd=True) + 2*eta/r
    grad_field = radial_derivative(field, dr, odd=False)
    momentum_density = hb*(field.conj()*grad_eta + eta.conj()*grad_field).imag/3
    p_relative = np.sum((volume*momentum_density)[:, :end])
    mass = m*np.sum((volume*density)[:, :end])
    velocity = acceleration/s
    p_lab = p_relative + mass*velocity
    export = flow + stress + m*rate*velocity
    incoming = feed + m*rate*velocity
    return dict(cut_flow_comoving=flow, cut_stress=stress, boost_flux=m*rate*velocity,
                cut_lab=export, feed_lab=incoming, side_to_support=side, mass=mass,
                p_relative=p_relative, p_lab=p_lab, storage_rate_lab=s*p_lab,
                balance_residual=s*p_lab + export + side - incoming,
                peak_log_density_perturbation=float(abs(b).max()))


def main():
    for level in HASHES:
        domain, psi, rate = background(level)
        unit_drive = -domain.constituent_mass*domain.dr*(domain.r_index + 0.5)*domain.volume*psi
        drive = np.r_[unit_drive.real, unit_drive.imag]
        for s in (0.02, 0.01):
            operator = linear_operator(domain, psi, s)
            solve = splu(operator)
            signed = []
            for amplitude in (1e-4, -1e-4):
                rhs = amplitude*drive
                vector = solve.solve(rhs)
                eta_real, eta_imag = np.split(vector, 2)
                readings = measure(domain, psi, eta_real + 1j*eta_imag, rate, s, amplitude)
                linear_residual = float(abs(operator @ vector - rhs).max())
                print(json.dumps(dict(level=level, s=s, A=amplitude, rate=rate,
                                      linear_residual=linear_residual, **readings), allow_nan=False), flush=True)
                if linear_residual > 1e-13 or readings["peak_log_density_perturbation"] >= 0.01:
                    raise ArithmeticError("Linear response solve or amplitude check failed")
                signed.append([readings[k] for k in ("cut_lab", "p_lab", "side_to_support", "feed_lab")])
            np.testing.assert_allclose(signed[0], -np.array(signed[1]), atol=1e-15, rtol=1e-12)


if __name__ == "__main__":
    main()
