"""Test whether the brane density deficit supplies an independent mouth trajectory at release."""

import json
from pathlib import Path
import shlex
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import CubicSpline

from lateral_mouth_acceleration import Translation
from lateral_source_rate import verified_packet


def intake_moment(motion, spline, shift):
    slab = motion.slab
    x = slab.x[:, None]
    g = np.maximum(1-(x-shift)**2-slab.r**2, 0)**3
    raw = g*motion.n[-1]/spline(slab.x-shift)
    current = motion.rate*raw/np.sum(slab.dv*raw)
    return np.sum(slab.dv*x*current)


def main():
    archive = json.loads(Path("docs/_measurements/lateral_mouth_acceleration.json").read_bytes())
    for run in archive["runs"]:
        if not run["limits_verified"] or run["worker_exit"] or run["guard_exit"]:
            raise ValueError("Unverified isolated preparation")
        path = Path(shlex.split(run["command"])[-1])
        with verified_packet(path, {path: run["packet_sha256"]}) as packet:
            x, r, dz = packet["x"], packet["r"], float(packet["dz"])
            dx = 2*r[0]
            if (int(packet["d"]) or packet["currents"].shape[0] != 1
                    or not np.allclose(np.diff(x), dx) or not np.allclose(np.diff(r), dx)):
                raise ValueError("Incorrect isolated-grid geometry")
            slab = SimpleNamespace(x=x, r=r, dx=dx, dz=dz, dv=2*np.pi*r*dx**2)
            motion = Translation(slab, packet["field"], packet["currents"][0])
        # Count deficit relative to the uniform density-one reference, over the
        # entire brane cylinder, with a half-weight at each end in depth.
        depth = np.full(len(motion.f), dz)
        depth[[0, -1]] /= 2
        if not np.isclose(depth.sum(), 1):
            raise ValueError("Incorrect slab depth")
        deficit = (1-motion.n)*depth[:, None, None]*slab.dv
        count = np.sum(deficit)
        initial_dipole = np.sum(deficit*x[:, None])
        if count <= 0 or abs(initial_dipole) > 1e-5*count:
            raise ArithmeticError("Cannot normalize a centered density deficit")

        # X=t²/2 per unit acceleration. The first change of the lab field is
        # its third derivative b=i ∂X(normal)/dz, supported on the mouth plane.
        field_third = 1j*motion.displaced_normal/dz
        deficit_third = -2*(motion.face.conj()*field_third).real*depth[-1]*slab.dv
        count_third = np.sum(deficit_third)
        dipole_third = np.sum(x[:, None]*deficit_third)
        centroid_third = dipole_third/count - initial_dipole*count_third/count**2

        # Independent boundary calculation: differentiate the first moment of
        # the normalized incoming current at fixed density, with two shift sizes.
        spline = CubicSpline(x, motion.n[-1], axis=0)
        differences = [(intake_moment(motion, spline, e)-intake_moment(motion, spline, -e))/(2*e)
                       for e in (dx*1e-4, dx*5e-5)]
        if (abs(count_third)>1e-12 or not np.isfinite(centroid_third)
                or not np.allclose(differences, dipole_third, rtol=1e-7, atol=1e-13)):
            raise ArithmeticError("Density-moment and intake-moment derivatives disagree")
        print(json.dumps(dict(field=str(path), L=int(round(r[-1]+r[0])), dx=dx, dz=dz, rate=motion.rate,
                              brane_deficit=count, initial_deficit_dipole=initial_dipole,
                              deficit_count_third_per_A=count_third, deficit_dipole_third_per_A=dipole_third,
                              frozen_intake_moment_slopes=differences, deficit_centroid_third_per_A=centroid_third,
                              deficit_centroid_t3_coefficient_per_A=centroid_third/6),
                         allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
