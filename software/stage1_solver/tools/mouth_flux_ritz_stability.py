"""Local stability, concentrated starts and omitted modes of the collar trials.

Run after mouth_flux_ritz.py --fields SOURCE. The returned curvature values
are finite-basis Hessian eigenvalues in the producer's scaled coordinates.
"""

import argparse
from math import factorial
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import lpmv

from mouth_flux_ritz import Collar


def concentrated(c):
    # Two non-antipodal spots carrying equal and opposite integrated flux.
    # The analytic spherical sections are peaked at their respective spots,
    # with additional localization and noncollinear complex phase gradients.
    points = np.array([[.71,.05,.70], [-.30,.82,-.49]])
    points /= np.linalg.norm(points, axis=1)[:, None]
    directions = np.vstack((c.sin*np.cos(c.phi), c.sin*np.sin(c.phi), c.x))
    dot = points @ directions
    values = {}
    for h, (point, (labels, B, _, _, _), m) in enumerate(zip(points, c.sections, c.ms)):
        zp = (1+point[2])/2
        phip = np.arctan2(point[1], point[0])
        coherent = (np.sqrt(c.z*zp) + np.sqrt((1-c.z)*(1-zp))*
                    np.exp(1j*np.sign(m)*(c.phi-phip)))**abs(m)
        f = 1.8*coherent*np.exp(2.5*(dot[h]-1) + .4j*directions[1-h])
        projected = B.conj() @ (c.wa*f)
        for (p, j), amp in zip(labels, projected):
            values["re",h,0,p,j] = amp.real
            values["im",h,0,p,j] = amp.imag
    # star d(*dY/sqrt(l(l+1))) = -sqrt(l(l+1))*Y on the unit sphere.
    flux = 4*(np.exp(4*(dot[1]-1))-np.exp(4*(dot[0]-1)))
    for l, q in c.gauge[0][1]:
        k = abs(q)
        angular = np.sin(k*c.phi) if q < 0 else np.cos(k*c.phi)
        norm = np.sqrt((2*l+1)*(1+(k>0))*factorial(l-k)/factorial(l+k))
        Y = norm*lpmv(k, l, c.x)*angular
        values["curl",0,l,q] = -(Y @ (c.wa*flux))/np.sqrt(l*(l+1))
    return c.vector(values)


def local_hessian(c, vector):
    step = 2e-6
    H = np.empty((len(vector), len(vector)))
    for j in range(len(vector)):
        delta = np.zeros_like(vector)
        delta[j] = step
        H[:,j] = (c.energy(vector+delta)[1]-c.energy(vector-delta)[1])/(2*step)
    asymmetry = np.max(abs(H-H.T))
    eig = np.linalg.eigvalsh((H+H.T)/2)
    return asymmetry, eig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for path in sorted(args.source.glob("*_l3.npz")):
        name = path.stem.removesuffix("_l3")
        with np.load(path) as saved:
            n = int(saved["n"])
            gains = tuple(saved["gains"])
            coarse = Collar(n, gains, *saved["level"])
            vector = saved["vector"]
            old_energy = saved["energies"][0]
        previous = coarse.physical(vector)
        asymmetry, eig = local_hessian(coarse, vector)
        print(f"{name} Hessian dof={len(eig)} asymmetry={asymmetry:.2e} "
              f"first7={np.array2string(eig[:7], precision=8, max_line_width=200)} "
              f"negative_below_-1e-5={sum(eig < -1e-5)}", flush=True)
        if asymmetry > 2e-7:
            raise ArithmeticError("Inconsistent symmetry in directional Hessian.")
        for i, level in enumerate(((5,6,3,4), (7,8,3,4))):
            c = Collar(n, gains, *level)
            finer = Collar(n, gains, *level, extra=3)
            transfer = c.vector(previous)
            E0, grad0 = c.energy(transfer)
            print(f"{name} enrich={level} embedding_energy_error={E0-old_energy:.2e} "
                  f"omitted_gradient_inf={np.linalg.norm(grad0, np.inf):.4e}", flush=True)
            if abs(E0-old_energy)>1e-10:
                raise ArithmeticError("Nested prolongation failed.")
            starts = [("transferred", transfer), ("two_spots", concentrated(c))]
            best = None
            for label, x0 in starts:
                res = minimize(c.energy, x0, jac=True, method="L-BFGS-B",
                               options={"ftol": 4e-15, "gtol": 1e-8, "maxcor": 16,
                                        "maxiter": 2000, "maxls": 40})
                grad = np.linalg.norm(res.jac, np.inf)
                values = finer.energy(res.x, True)
                err = abs(values[0]-res.fun)
                print(f"{name} level={level} seed={label} E={res.fun:.12f} "
                      f"shift={res.fun-E0:.4e} dof={len(res.x)} it={res.nit} "
                      f"ok={res.success} grad={grad:.2e} quadrature={err:.2e}", flush=True)
                if not res.success or grad>8e-6 or err>2e-10:
                    raise ArithmeticError(f"Unresolved Ritz check: {res.message}")
                if best is None or res.fun<best[0]:
                    best = (res.fun, res.x, values)
            old_energy, vector, values = best
            previous = c.physical(vector)
            np.savez(args.output/f"{name}_l{4+i}.npz", n=n, gains=gains, level=level,
                     vector=vector, energies=values)


if __name__ == "__main__":
    main()
