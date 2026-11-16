"""Validate and difference the redistributed and fixed-profile saved solutions."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_feed_weak_check import check_file, dx4


def face_fields(path):
    with np.load(path) as packet:
        field = packet["field"]
        speed = float(packet["speed"])
        dx, dz = float(packet["dx"]), float(packet["dz"])
        face = field[-1]
        density = abs(face)**2
        vx = speed + (face.conj()*dx4(face, dx)).imag/density
        normal = (11*face-18*field[-2]+9*field[-3]-2*field[-4])/(6*dz)
        current_fd = (face.conj()*normal).imag
        return density, vx, packet["current"], current_fd, dx, packet["mouth"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--driven", type=Path, required=True)
    parser.add_argument("--level", type=int, required=True)
    parser.add_argument("--reference-manifest", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.reference_manifest.read_text())
    hashes = {Path(path).resolve(): digest for run in manifest["runs"]
              for path, digest in run.get("fields", {}).items()}
    for law in ("flat_density_slope", "restoring_slope"):
        stem = f"n{args.level}_{law}_u"
        ref = [args.reference / (stem + f"{speed:g}.npz") for speed in (0, 0.01, 0.02)]
        for path in ref:
            if hashlib.sha256(path.read_bytes()).hexdigest() != hashes[path.resolve()]:
                raise ValueError(f"Invalid fixed-intake archive {path}")
        base = check_file(ref[0], law)
        n0, v0, j0, _, dx, mouth = face_fields(ref[0])
        x = (np.arange(args.level) - args.level/2 + 0.5)*dx
        downstream_half = mouth & (x[:, None, None] > 0)
        dv = dx**3
        for fixed_path, speed in zip(ref[1:], (0.01, 0.02)):
            free_path = args.driven / (stem + f"{speed:g}.npz")
            fixed, free = (check_file(path, law) for path in (fixed_path, free_path))
            _, vf, jf, jf_fd, _, _ = face_fields(fixed_path)
            n, vr, jr, jr_fd, _, _ = face_fields(free_path)
            # Verify the imposed shape of normal speed against the archive.
            # The one allowed change is a scalar common to all aperture cells.
            speed_ratio = (jr[mouth]/n[mouth]) / (j0[mouth]/n0[mouth])
            if (np.ptp(speed_ratio) > 2e-12 or np.any(jr[mouth] <= 0)
                    or np.any(jr[~mouth] != 0) or not np.array_equal(jf, j0)):
                raise ValueError("Local intake law disagrees with the saved fields")
            def integral(f):
                return dv * f[mouth].sum()
            shift = dv*(jr-jf)[downstream_half].sum()
            shift_fd = dv*(jr_fd-jf_fd)[downstream_half].sum()
            transported = integral((jr-jf)*vf)/speed
            velocity_change = integral(jr*(vr-vf))/speed
            flow_difference = (free[3]-fixed[3])/speed
            if abs(transported+velocity_change-flow_difference) > 2e-13:
                raise ArithmeticError("Carried-momentum difference did not close")
            if abs(free[2]-speed) > 1e-12 or abs(integral(jr-j0)) > 2e-13:
                raise ArithmeticError("Drive or total current changed")
            if max(abs(base[3]), abs(base[4])) > 1e-12:
                raise ArithmeticError("Zero-drive baseline has a net sideways load")
            print(f"n={args.level} law={law} U={speed:g} "
                  f"half_aperture_rate_shift={shift:.12g} shift_fd={shift_fd:.12g} "
                  f"centroid_shift={integral((jr-j0)*x[:,None,None])/integral(jr):.12g} "
                  f"normal_speed_factor={speed_ratio.mean():.12g} "
                  f"fixed_flow={(fixed[3]-base[3])/speed:.12g} "
                  f"responsive_flow={(free[3]-base[3])/speed:.12g} "
                  f"flow_difference={flow_difference:.12g} "
                  f"flow_difference_fd={(free[5]-fixed[5])/speed:.12g} "
                  f"weight_term={transported:.12g} velocity_term={velocity_change:.12g} "
                  f"fixed_stress={(fixed[4]-base[4])/speed:.12g} "
                  f"responsive_stress={(free[4]-base[4])/speed:.12g} "
                  f"stress_difference={(free[4]-fixed[4])/speed:.12g} "
                  f"stress_difference_fd={(free[6]-fixed[6])/speed:.12g} "
                  f"responsive_stress_fd={(free[6]-base[6])/speed:.12g} "
                  f"direct_redistribution_on_baseline_velocity={integral((jr-j0)*v0)/speed:.12g} "
                  f"max_flux_error={free[1]:.3g} max_residual={free[0]:.3g}", flush=True)


if __name__ == "__main__":
    main()
