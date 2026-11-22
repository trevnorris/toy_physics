"""Test an affine backing pressure against the eight saved radial-work bands."""

import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import linprog

from mouth_wave_flare_work import radial_profile


ARCHIVE = Path(__file__).resolve().parents[3] / (
    "research/pde_ledger_v3/_measurements/mouth_occupied_radial_balance.json")


def depth_moments(bands, degree):
    nodes, weights = leggauss(degree)
    t = (1 + nodes) / 2
    moments, volumes = [], []
    for lo, hi in bands:
        w = lo + (hi - lo) * t
        r, _ = radial_profile(w)
        dvol = 4 * np.pi * r**2 * 16 * t**2 * (1 - t)**2 * weights * (hi - lo) / 2
        volumes.append(sum(dvol))
        moments.append(w @ dvol / volumes[-1])
    return np.array(moments), np.array(volumes)


def compare(matrix, observed):
    scales = np.linalg.norm(matrix, axis=0)
    design = matrix / scales
    # Least maximum mismatch for all real coefficients strictly enlarges the
    # admissible material class. The dual contrast cancels all four columns.
    upper = np.column_stack([np.vstack([design, -design]), -np.ones(2 * len(observed))])
    fit = linprog([0, 0, 0, 0, 1], A_ub=upper, b_ub=np.r_[observed, -observed],
                  bounds=[(None, None)] * 4 + [(0, None)], method="highs")
    if not fit.success:
        raise ArithmeticError(fit.message)
    coef = fit.x[:-1] / scales
    residual = observed - matrix @ coef
    dual_hi, dual_lo = np.split(fit.ineqlin.marginals, 2)
    contrast = dual_hi - dual_lo
    row = dict(backing_at_origin=coef[0], backing_slope=coef[1], tension=coef[2],
               wave_action=coef[3], max_error=fit.fun, residual=residual.tolist(),
               dual_weights=contrast.tolist(), dual_l1=np.sum(abs(contrast)),
               dual_columns=(contrast @ matrix).tolist(), dual_pressure=contrast @ observed)
    return row


def main():
    archived = json.loads(ARCHIVE.read_text())["runs"]
    for n in (26, 34, 48):
        case = archived[f"occupied_radial_n{n}_01"]
        source = case["literal_stdout"]
        if (case["guarded_exit"] != 0 or case["worker_exit"] != 0
                or case["limits_verified"] is not True or case["swap_max"] != 0
                or hashlib.sha256(source.encode()).hexdigest() != case["stdout_sha256"]):
            raise ValueError("Invalid guarded input")
        records = [json.loads(line) for line in source.splitlines()]
        carrier, = [r for r in records if r["test"] == "carrier_work"]
        wave, = [r for r in records if r["test"] == "wave_state"]
        mid, vol = depth_moments(carrier["bands"], 12)
        mid2, vol2 = depth_moments(carrier["bands"], 18)
        np.testing.assert_allclose(mid, mid2, atol=1e-14, rtol=0)
        np.testing.assert_allclose(vol, vol2, atol=1e-14, rtol=0)
        np.testing.assert_allclose(vol, wave["projected_volumes"], atol=2e-14, rtol=0)
        matrix = np.column_stack((np.ones_like(mid), mid, carrier["curvature"],
                                  -np.array(wave["pressure_per_action"])))
        print(json.dumps(dict(test="preparation", mesh=n, source_sha256=case["stdout_sha256"],
                              bands=carrier["bands"], volume_mean_depths=mid.tolist(),
                              design=matrix.tolist())))
        for method in ("surface", "volume"):
            observed = np.array(carrier[method])
            row = compare(matrix, observed)
            print(json.dumps(dict(test="affine_work", mesh=n, method=method,
                                  carrier_pressure=observed.tolist(), **row)))
            np.testing.assert_allclose(row["dual_columns"], 0, atol=3e-13)
            np.testing.assert_allclose(row["dual_l1"], 1, atol=1e-12)
            np.testing.assert_allclose(row["dual_pressure"], row["max_error"], atol=2e-11)
            np.testing.assert_allclose(np.max(np.abs(row["residual"])), row["max_error"],
                                       rtol=2e-7, atol=2e-11)


if __name__ == "__main__":
    main()
