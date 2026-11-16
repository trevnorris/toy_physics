"""Hold and probe finite spherical intakes in a shared 3-D potential-flow exterior.

Units: radius, suction speed and displaced mass all equal one. Multipoles
satisfy the normal velocity boundary conditions on every intake. At release
the centers are at rest, so the potential time derivative is the response
to the two center accelerations with the geometry still at its initial value.
"""

from math import pi

import numpy as np
from numpy.polynomial.legendre import leggauss


DENSITY = 3 / (4 * pi)
RATE = 4 * pi
PROBE = 3 / 200


def basis(axial, transverse, centers, order):
    """Potential and axial/transverse gradients for all decaying solid harmonics."""
    phi, ux, ur = [], [], []
    for center in centers:
        offset = axial - center
        radius = np.hypot(offset, transverse)
        mu = offset / radius
        sine = transverse / radius
        previous = np.zeros_like(mu)
        legendre = np.ones_like(mu)
        for ell in range(order + 1):
            # dP/dmu is regular on these open Gauss nodes.
            derivative = ell * (previous - mu * legendre) / (1 - mu**2)
            scale = radius**(-ell - 2)
            phi.append(scale * radius * legendre)
            ux.append(scale * (-(ell + 1) * mu * legendre + (1 - mu**2) * derivative))
            ur.append(-scale * sine * ((ell + 1) * legendre + mu * derivative))
            previous, legendre = legendre, ((2 * ell + 1) * mu * legendre - ell * previous) / (ell + 1)
    return tuple(np.column_stack(columns) for columns in (phi, ux, ur))


def measure(distance, order):
    centers = [0.0] if distance is None else [0.0, distance]
    count = len(centers)
    mu, weights = leggauss(2 * order + 5)
    sine = np.sqrt(1 - mu**2)
    rows, values = [], []
    for i, center in enumerate(centers):
        _, bx, br = basis(center + mu, sine, centers, order)
        rows.append(weights[:, None]**0.5 * (mu[:, None] * bx + sine[:, None] * br))
        rhs = np.zeros((len(mu), count + 1))
        rhs[:, 0], rhs[:, i + 1] = -1, mu
        values.append(weights[:, None]**0.5 * rhs)
    coeff, _, rank, _ = np.linalg.lstsq(np.vstack(rows), np.vstack(values), rcond=None)
    if rank != coeff.shape[0]:
        raise ArithmeticError("singular Neumann fit")

    # Read force on a different, denser grid from the boundary fit.
    mu, weights = leggauss(3 * order + 8)
    sine = np.sqrt(1 - mu**2)
    steady = np.zeros(count)
    added = np.zeros((count, count))
    incident = np.zeros(count)
    error = 0.0
    for i, center in enumerate(centers):
        phi, bx, br = basis(center + mu, sine, centers, order)
        vx, vr = bx @ coeff[:, 0], br @ coeff[:, 0]
        normal = mu * vx + sine * vr
        bernoulli = -(vx**2 + vr**2) / 2
        steady[i] = -2 * pi * DENSITY * np.dot(weights, vx * normal + bernoulli * mu)
        added[i] = -2 * pi * DENSITY * ((weights * mu) @ phi @ coeff[:, 1:])
        expected = np.zeros((len(mu), count))
        expected[:, i] = mu
        normal_rate = (mu[:, None] * bx + sine[:, None] * br) @ coeff[:, 1:]
        error = max(error, np.max(np.abs(normal + 1)), np.max(np.abs(normal_rate - expected)))
        # On a center the regular field from the other sphere is exactly axial.
        for j, other in enumerate(centers):
            if i == j:
                continue
            ell = np.arange(order + 1)
            dx = center - other
            other_coeff = coeff[j * (order + 1):(j + 1) * (order + 1), 0]
            incident[i] += np.dot(other_coeff, -(ell + 1) * np.sign(dx)**(ell + 1) / abs(dx)**(ell + 2))

    inertia = np.eye(count) + added
    mobility = np.linalg.solve(inertia, np.eye(count))
    all_results = []
    for fraction in (0, 1):
        # Controls from the one-intake benchmark: no lateral outlet stress;
        # outward mass flow has mean axial velocity fraction * incident.
        support = DENSITY * RATE * fraction * incident - steady
        direct = np.linalg.solve(inertia, -support)
        measured_cols = []
        for i in range(count):
            probe = PROBE * np.eye(count)[:, i]
            plus = np.linalg.solve(inertia, probe - support)
            minus = np.linalg.solve(inertia, -probe - support)
            measured_cols.append((plus - minus) / (2 * PROBE))
            if i == 0:
                target_probes = (minus[0], plus[0])
        measured_mobility = np.column_stack(measured_cols)
        predicted = measured_mobility @ (-support)
        discrepancy = np.max(np.abs(predicted - direct))
        values = (support[0], *target_probes, direct[0], discrepancy)
        print("k=%d H_target=%+.12g a_minus=%+.12g a_plus=%+.12g a_free=%+.12g prediction_error=%.2g" %
              (fraction, *values), flush=True)
        if discrepancy > 2e-12 or np.max(np.abs(measured_mobility - mobility)) > 2e-12:
            raise ArithmeticError("held-force and probe reconstruction failed")
        if count == 2 and max(abs(support.sum()), abs(direct.sum())) > 2e-12:
            raise ArithmeticError("unforced pair lost exchange symmetry")
        all_results.extend(values[:4])
    print("max_normal_error=%.3g mobility_target_row=%s" %
          (error, " ".join("%.12g" % x for x in mobility[0])), flush=True)
    if error > 2e-9 or np.max(np.abs(inertia - inertia.T)) > 2e-12:
        raise ArithmeticError("unresolved surface or asymmetric added mass")
    if count == 1 and (abs(added[0, 0] - 0.5) > 2e-12 or abs(steady[0]) > 2e-12):
        raise ArithmeticError("isolated spherical control failed")
    return np.asarray(all_results)


if __name__ == "__main__":
    for distance in (None, 10.0, 20.0):
        previous = None
        for order in (8, 14):
            print(f"separation={distance} degree={order}", flush=True)
            result = measure(distance, order)
            if previous is not None:
                change = max(abs(result - previous))
                print("max_refinement_change=%.3g" % change, flush=True)
                if change > 2e-10:
                    raise ArithmeticError("surface-force refinement failed")
            previous = result
