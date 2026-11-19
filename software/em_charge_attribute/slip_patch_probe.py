"""Probe the saved collar at the moving zero and the equator of its section."""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Legendre, Polynomial
from scipy.linalg import eigh

from relaxed_neck import Neck


def probe_rows(r, mu, radial_order, angular_order):
    """Value and radial derivative of the inherited outer-element Ritz basis."""
    t = Polynomial([0, 1])
    bubbles = [t*(1-t)*Legendre.basis(i).convert(kind=Polynomial)(2*t-1)
               for i in range(radial_order-1)]
    polynomials = [1-t, t] + [Polynomial([0])] * len(bubbles) + bubbles
    angular = [Legendre.basis(i)(mu) for i in range(angular_order+1)]
    return tuple(np.kron([p.deriv(d)(r-2) for p in polynomials], angular) for d in (0, 1))


def gram(rows, weights):
    return rows.T @ (weights[:, None]*rows)


def probe(label, path, radial, angular, quadrature):
    state = next(s for s in path["states"] if s["beta"] == 2.5)
    model = Neck(*path["shape"], path["free_angle"], radial, angular, quadrature)
    hr, ht = (model.shapes @ np.r_[1., state["beta"], state["coeff"]])[:2]
    volume = np.sqrt(1 + hr**2 + ht**2)
    weights = model.wt * model.chi
    mass = gram(model.val, weights*volume)
    stiffness = (gram(model.dr, weights/volume) + gram(model.dt, weights/volume)
                 + gram(ht[:, None]*model.dr-hr[:, None]*model.dt, weights/volume))
    eigenvalues, modes = eigh(stiffness, mass, subset_by_index=(0, 1))
    omega = np.sqrt(eigenvalues[0])
    mode = modes[:, 0]
    action = path["action"]
    print(f"{label} wave={radial},{angular} quad={quadrature} omega={omega:.10f} "
          f"saved={state['parts'][2]:.10f} next={np.sqrt(eigenvalues[1]):.8f}")
    sensitivities = []
    for name, mu in (("slip", 1.), ("equator", 0.)):
        value, derivative = probe_rows(2.5, mu, radial, angular)
        u, ur = value @ mode, derivative @ mode
        # Positive control. A small inertial load with matrix element u(site)^2
        # is added at fixed stiffness. Only the first variation is interpreted.
        step = 1e-4
        perturbation = step * np.outer(value, value)
        frequencies = [np.sqrt(eigh(stiffness, mass+sign*perturbation,
                                    eigvals_only=True, subset_by_index=(0, 0))[0])
                       for sign in (-1, 1)]
        predicted = -action*omega*u**2/2
        measured = action*(frequencies[1]-frequencies[0])/(2*step)
        sensitivities.append(predicted)
        print(f"  {name}: U^2={u*u:.10f} d_r(U^2)={2*u*ur:.10f} "
              f"dE/d(inertial load)={predicted:.8f} centered={measured:.8f}")
        if not np.isclose(measured, predicted, rtol=2e-5, atol=2e-7):
            raise ValueError("Point-readout and mass-loading derivative disagree")
    # Occupation changes the background energy and the positive control's response.
    for occupied in (action, 0.):
        print(f"  J={occupied:.9f}: wave energy={occupied*omega:.8f}, "
              f"control slip-equator per inertial load={occupied/action*(sensitivities[0]-sensitivities[1]):.8f}")
    if not (0 < eigenvalues[0] < eigenvalues[1]) or abs(mode @ mass @ mode - 1) > 1e-9:
        raise ValueError("Invalid fundamental or normalization")
    if (radial, angular, quadrature) == (*path["wave"], path["quad"]):
        if abs(omega-state["parts"][2]) > 1e-8:
            raise ValueError("Geometry or mode does not reproduce the saved preparation")


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[2] / "docs/_measurements"
    for label, path in json.loads((folder/"relaxed_neck_paths.json").read_text()).items():
        for level in [(*path["wave"], path["quad"]), (11, 10, 36)]:
            probe(label, path, *level)
