"""Virtual translation of a diagnostic surface energy on the measured carrier trace.

The energy is zero on the common resting trace. Its strength is left as an
unknown multiplier. All force readings are per unit multiplier in solver units.
"""

import argparse
from pathlib import Path

import numpy as np

from held_crossflow_budget import archived_packet
from mouth_carrier_momentum import solve_tangent
from mouth_feed_weak_check import dx4
from mouth_feed_weak_flow import FeedPatch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--face-directory", type=Path, required=True)
    args = parser.parse_args()
    args.face_directory.mkdir(parents=True, exist_ok=False)
    baseline = archived_packet("mouth_feed_weak_run.json", "n32_restoring_slope_u0.npz")
    with np.load(baseline, allow_pickle=False) as data:
        field = data["field"]
        patch = FeedPatch(field.shape[1], field.shape[0], "restoring_slope")
        dx = patch.dx
        if (float(data["speed"]) != 0 or float(data["dx"]) != dx
                or float(data["dz"]) != patch.dz
                or not np.array_equal(data["mouth"], patch.mouth)
                or not np.array_equal(data["current"], patch.current)):
            raise ValueError("Resting field and tangent solver have different preparations")
    resting = abs(field[-1])**2
    # The saved solver differentiates with respect to mouth speed at fixed
    # feed. The held-mouth feed speed has the opposite sign.
    moving, steps, error = solve_tangent(patch, field)
    density_tangent = -2*(field[-1].conj()*moving[-1]).real
    del moving, field
    coord = (np.arange(resting.shape[0]) + 0.5 - resting.shape[0]/2)*dx
    radius_sq = sum(v**2 for v in np.meshgrid(coord, coord, coord, indexing="ij", sparse=True))
    weight = np.maximum(1 - 4*radius_sq, 0)**3
    cell = dx**3
    rate = cell * patch.current.sum()
    print(f"inward number rate={rate:.12g}; density tangent iterations={steps} last update={error:.4g}")

    # With coordinates attached to the mouth, a virtual shift by dX samples
    # n(y+dX). Only the incident trace varies in this representation of E.
    slope = -cell*np.sum(weight*density_tangent*dx4(resting, dx))
    print(f"virtual-work force per strength per incident speed at U=0: {slope:+.15g}", flush=True)
    np.testing.assert_allclose(rate, 0.25, rtol=0, atol=1e-13)
    np.testing.assert_allclose(resting, resting[::-1], rtol=0, atol=5e-13)
    np.testing.assert_allclose(density_tangent, -density_tangent[::-1], rtol=0, atol=5e-12)
    resting.astype("<f8").tofile(args.face_directory/"rest.f64")
    density_tangent.astype("<f8").tofile(args.face_directory/"tangent.f64")

    def force(density):
        return -cell*np.sum(weight*(density-resting)*dx4(density, dx))

    print(f"unloaded profile force at rest: {force(resting):+.15g}")
    for speed in (0.02, 0.01):
        filename = f"n32_restoring_slope_u{speed:g}.npz"
        path = archived_packet("mouth_intake_redistribution.json", filename)
        with np.load(path, allow_pickle=False) as data:
            density = abs(data["field"][-1])**2
            if (float(data["speed"]) != speed or float(data["dx"]) != dx
                    or float(data["dz"]) != patch.dz
                    or not np.array_equal(data["mouth"], patch.mouth)):
                raise ValueError("Driven preparation changed")
            np.testing.assert_allclose(cell*data["current"].sum(), rate, rtol=0, atol=1e-13)
        plus, minus = force(density), force(density[::-1])
        secant = (plus-minus)/(2*speed)
        print(f"abs_U={speed:g} force(+U)={plus:+.15g} force(-U)={minus:+.15g} "
              f"odd slope={secant:+.15g} slope minus tangent={secant-slope:+.8g}", flush=True)
        density.astype("<f8").tofile(args.face_directory/f"flow{speed:g}.f64")
        np.testing.assert_allclose(plus, -minus, rtol=1e-9, atol=1e-18)


if __name__ == "__main__":
    main()
