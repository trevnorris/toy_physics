"""Test the isolated density tail and the pressure contrast on a second opening."""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RegularGridInterpolator
import sympy as sp

from lateral_reservoir_readout import depth_sum
from lateral_source_rate import verified_packet


def paper_checks():
    r, q = sp.symbols("r q", positive=True)
    b, c, width_k, grad_k, k = sp.symbols("B C kW D k", real=True, nonzero=True)

    def lap(f):
        return sp.diff(r**2 * sp.diff(f, r), r) / r**2

    # Eliminate the radial phase velocity with r^2 n v_r = -q. Expand
    # the real carrier equation in inverse distance with a free density tail.
    n1, n4, n6 = sp.symbols("n1 n4 n6")
    n = 1 + n1/r + n4/r**4 + n6/r**6
    chem = (n**4 - 1)/4
    real_row = chem + q**2/(2*r**4*n**2) - lap(sp.sqrt(n))/(2*sp.sqrt(n))
    def coeff(expr, power):
        return sp.limit(r**power * expr, r, sp.oo)

    row1 = coeff(real_row, 1)
    row4 = coeff(real_row.subs(n1, 0), 4)
    row6 = coeff(real_row.subs({n1: 0, n4: -q**2/2}), 6)
    print("radial carrier residual at powers 1,4,6 =", row1, row4, row6, flush=True)
    assert row1 == n1 and sp.expand(row4 - n4 - q**2/2) == 0
    assert sp.expand(row6 - n6 - 3*q**2/2) == 0
    # Uniform static density/width slice with freely relaxed densification.
    hessian = sp.Matrix([[b, c], [c, width_k + grad_k*k**2]])
    width_response = sp.factor(hessian.inv()[1, 1])
    print("width response to width load at wave number k =", width_response, flush=True)
    assert sp.simplify(1/width_response - (width_k - c**2/b + grad_k*k**2)) == 0


def isolated_profile(field, x, radial, dz, rate):
    # Compute density before interpolating. Interpolating the complex phase
    # first would itself lower the interpolated density through phase mixing.
    deficit = depth_sum((1 - abs(field)) * (1 + abs(field)), dz)
    interp = RegularGridInterpolator((x, np.r_[-radial[0], radial]),
                                     np.column_stack((deficit[:, 0], deficit)))
    radii = (4, 6, 8, 12, 16)
    rows = []
    for order in (96, 192):
        mu, weights = leggauss(order)
        averages = []
        for radius in radii:
            points = radius * np.column_stack((mu, np.sqrt(1 - mu**2)))
            averages.append(np.dot(weights, interp(points))/2)
        rows.append(np.array(averages))
    analytic_coefficient = (rate/(4*np.pi))**2/2
    for radius, deficit_at_r in zip(radii, rows[-1], strict=True):
        print(f"  r={radius} mean_deficit={deficit_at_r:.12g} "
              f"r4_deficit_over_prediction={radius**4*deficit_at_r/analytic_coefficient:.9g}")
    relative_change = np.max(abs(rows[1] - rows[0])/np.maximum(abs(rows[1]), 1e-16))
    print(f"  max_relative_quadrature_change={relative_change:.6g}")
    if not np.isfinite(rows).all() or relative_change > 0.005:
        raise ValueError("Shell density quadrature did not converge")
    for i, j in ((0, 2), (2, 4)):
        ratio = rows[-1][j]/rows[-1][i]
        exponent = -np.log2(ratio) if ratio > 0 else np.nan
        print(f"  doubling={radii[i]}->{radii[j]} deficit_ratio={ratio:.9g} "
              f"effective_power={exponent:.7g}")


def half_contrast(values, x, radial):
    # The companion lies at negative x. Use the actual 3D volume of
    # each in-brane half-ball on the feed side at w = mouth-face height.
    mask = x[:, None]**2 + radial**2 < 1
    weights = np.broadcast_to(radial, mask.shape)
    def mean(half):
        selection = mask & half[:, None]
        return np.average(values[selection], weights=weights[selection])
    return mean(x < 0) - mean(x > 0)


def read_fields(tag, hashes):
    fields = Path("_scratch/s11c") / f"lateral_bath_{tag}_01" / "fields"
    with verified_packet(fields / "d0.npz", hashes) as packet:
        x, radial, dz = packet["x"], packet["r"], float(packet["dz"])
        f, j = packet["field"], packet["currents"][0]
    dx = x[1] - x[0]
    dv = 2*np.pi*radial*dx**2
    rate = np.sum(j*dv)
    print(f"grid={tag} rate={rate:.12g} width={(f.shape[0]-1)*dz:.12g}", flush=True)
    if not np.isclose(rate, 0.25, atol=1e-12, rtol=0) or (f.shape[0]-1)*dz != 1:
        raise ValueError("Reference preparation changed")
    isolated_profile(f, x, radial, dz, rate)
    n0 = abs(f[-1])**2
    if abs(half_contrast(n0, x, radial)) > 1e-12:
        raise ValueError("The isolated opening lacks reflection symmetry")
    for d in (4, 8):
        with verified_packet(fields / f"d{d}.npz", hashes) as packet:
            if (int(packet["d"]) != d or float(packet["dz"]) != dz
                    or not np.array_equal(x, packet["x"]) or not np.array_equal(radial, packet["r"])):
                raise ValueError("Pair grid differs from reference")
            n = abs(packet["layers"][-1])**2
            rates = (packet["currents"] * dv).sum(axis=(1, 2))
        if rates.shape != (2,) or not np.allclose(rates, rate, atol=1e-12, rtol=0):
            raise ValueError("Pair intake differs from reference")
        dn = half_contrast(n-n0, x, radial)
        dp = half_contrast((n**5-n0**5)/5, x, radial)
        print(f"  d={d} near_minus_far_paired_minus_alone dn={dn:.12g} dP={dp:.12g}")
        if not np.isfinite([dn, dp]).all():
            raise ValueError("Nonfinite opening contrast")


if __name__ == "__main__":
    paper_checks()
    record = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(p): digest for run in record["runs"] for p, digest in run.get("fields", {}).items()}
    for grid in ("e24_c4", "e24_c6", "e48_c4"):
        read_fields(grid, hashes)
