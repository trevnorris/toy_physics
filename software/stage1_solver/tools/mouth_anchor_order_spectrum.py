"""Does the supplied half-head scalar order spontaneously break up into sites?

Uses the fixed carrier, the driven order Jacobian and regular spherical harmonics
on the brane. No microscopic coordinates, contact parameters or clamp readings
enter the operator. Rates use the model clock of the order equation.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse.linalg import eigsh, LinearOperator

from mouth_carrier_lateral_response import prepare, quadrature_radii
from mouth_conversion_supply import SuppliedOrder
from mouth_free_order_diffuse import NU, positive_factor, read
from mouth_order_wave_pulse import digest, emit


ARCHIVE = Path("research/pde_ledger_v3/_measurements/mouth_offset_brane_readout.json")
SCRATCH = Path("_scratch/s11c")


def mode_rate(hessian, mass, degree):
    factor = positive_factor(hessian)
    if factor is None:
        raise ArithmeticError(f"Could not certify a positive Hessian for degree {degree}")
    inverse = LinearOperator(hessian.shape, matvec=factor.solve, dtype=float)
    values, vectors = eigsh(hessian, k=1, M=mass, sigma=0., which="LM", OPinv=inverse,
                            tol=1e-10, v0=np.ones(hessian.shape[0]))
    rate, vector = float(values[0]), vectors[:, 0]
    mass_vector = mass @ vector
    residual = np.linalg.norm(hessian @ vector - rate*mass_vector)/np.linalg.norm(mass_vector)
    # A second factorization certifies a lower bound on every eigenvalue,
    # including modes that could have little overlap with the eigensolver seed.
    lower = .99*rate
    shifted = positive_factor(hessian - lower*mass)
    emit(test="angular_order", degree=degree, free_dofs=len(vector),
         rate=rate, eigen_residual=residual, lower_test=lower,
         shifted_positive=shifted is not None)
    if shifted is None or residual > 1e-7 or rate <= 0:
        raise ArithmeticError("Lowest rate verification failed")
    return rate


def run(mesh):
    prior = json.loads(ARCHIVE.read_text())["verified_runs"][f"native_offset_brane{mesh}_01"]
    records = [json.loads(line) for line in prior["literal_stdout"]]
    state = next(row for row in records if row["test"] == "state")
    probe = next(row for row in records if row["test"] == "contact_free_bulk")
    supply = SCRATCH / f"contact_supply_n{mesh}_01"
    provenance = [json.loads(line) for line in (supply/"stdout").read_text().splitlines()]
    calibration = next(row for row in provenance if row["test"] == "calibration")
    level = next(row for row in provenance if row["test"] == "supplied_order" and row["multiple"] == .5)
    order_path = Path(level["path"])
    if (state["head"] != level["head"] or state["intake"] != calibration["intake"]
            or state["mobility"] != calibration["mobility"]):
        raise ValueError("Mismatched preparation")
    for name, expected in state["fields"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Changed frozen field {name}")
    if digest(order_path) != state["fields"][str(order_path)]:
        raise ValueError("Half-head order did not match native preparation")
    domain, psi = prepare(SimpleNamespace(source=SCRATCH/f"mouth_flared_delayed_n{mesh}"))
    with np.load(order_path) as saved:
        order = saved["chi"]
    model = SuppliedOrder(domain, psi, calibration["width"], state["mobility"])
    radius = quadrature_radii(domain)
    if not np.all((radius > 0) & (radius <= domain.far)):
        raise ValueError("Unexpected quadrature radii")
    c = (domain.N @ order)/model.quad_sum
    d = (domain.N @ (1-order))/model.quad_sum
    drive = model.wt*model.mobility*level["head"]*model.density
    residual, jac = model.measure(order)
    residual += model.Nf.T @ (drive*6*c*d)
    jac += model.mass_weight(drive*6*(d-c))
    rc, wc = probe["position"]
    observed = float(read(domain, order, rc, wc)[0])
    emit(test="state", mesh=mesh, intake=state["intake"], head=level["head"],
         fields=state["fields"], patch=[rc, wc], local_order=observed,
         max_stationary_weak_residual=float(max(abs(residual))), outer_radius=domain.far,
         nu=NU, minimum_weight=float(min(model.wt)))
    if max(abs(residual)) > 1e-9 or abs(observed-.5) > 1e-10 or min(model.wt) <= 0:
        raise ValueError("Invalid reference order or weights")
    # A regular nonspherical harmonic vanishes on the brane axis r=0.
    regular = domain.points[domain.free, 0] > 0
    angular_mass = model.mass_weight(model.wt*NU/radius**2)
    base = mode_rate(jac, model.mass, 0)
    for degree in (1, 2):
        multiplier = degree*(degree+1)
        hessian = (jac + multiplier*angular_mass)[regular][:, regular]
        rate = mode_rate(hessian, model.mass[regular][:, regular], degree)
        bound = .99*base + NU*multiplier/domain.far**2
        emit(test="positive_angular_cost", degree=degree, rate=rate, lower_bound=bound)
        if rate <= bound:
            raise ArithmeticError("Angular Rayleigh lower bound violated")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", required=True, type=int, choices=(34, 48))
    run(parser.parse_args().mesh)
