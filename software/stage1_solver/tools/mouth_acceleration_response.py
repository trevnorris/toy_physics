"""Periodic small-signal brane response at two crossings of the same mouth speed.

In y=x-X(t), with the lab feed at rest, i f_t = H(f) + i V(t) f_y.
Write V=V0+epsilon*cos(omega*t) and f=f0+epsilon*(p*cos+q*sin).
The two real-linear carrier equations are Lp+i*d_x f0=i*omega*q
and Lq=-i*omega*p. Normal speed shape and total intake follow the
archived redistributing mouth. Compute the periodic particular response.
"""

import argparse
from pathlib import Path

import numpy as np

from mouth_carrier_momentum import read_state
from mouth_feed_weak_flow import FeedPatch
from mouth_release_input_check import manifest_fields


def boundary_variation(patch, f, current, h):
    face, change = f[-1], h[-1]
    amplitude = abs(face)
    n = amplitude**2
    dn = 2*(face.conj()*change).real
    mean = np.sum(current*dn/n)/current.sum()
    speed = current/n
    slope = -patch.restoring*(1-1/amplitude)
    dslope = -patch.restoring*dn/(2*amplitude**3)
    return (slope+1j*speed)*change + (dslope-1j*speed*mean)*face


def solve_periodic(patch, base, current, v0, omega):
    density = abs(base)**2
    local = patch.shift - (density**4-1)/4
    lift = patch.elliptic_step(np.zeros_like(base), np.zeros_like(base[-1]))
    drive = 1j*patch.x_derivative(base)

    def advance(h, source):
        dn = 2*(base.conj()*h).real
        rhs = local*h - density**3*dn*base - 1j*v0*patch.x_derivative(h) + source
        return patch.elliptic_step(rhs, boundary_variation(patch, base, current, h)) - lift

    p, q = np.zeros_like(base), np.zeros_like(base)
    for step in range(1, 251):
        new_p = advance(p, 1j*omega*q - drive)
        new_q = advance(q, -1j*omega*p)
        error = max(np.max(abs(new_p-p)), np.max(abs(new_q-q)))
        if not np.isfinite(error):
            raise ArithmeticError("Nonfinite harmonic field update")
        p, q = new_p, new_q
        if error < 2e-13:
            return p, q, step, error
    raise ArithmeticError(f"Unconverged harmonic field update {error}")


def load_preparation(level, law):
    coarse = level == 24
    scratch, docs = Path("_scratch/s11c"), Path("docs/_measurements")
    held_dir = scratch / ("mouth_feed_weak01" if coarse else "mouth_feed_depth01") / "fields"
    driven_dir = scratch / ("mouth_redistribution_coarse01" if coarse else "mouth_redistribution_depth01") / "fields"
    zero, j0, dx, dz = read_state(held_dir, level, law, 0, manifest_fields(docs/"mouth_feed_weak_run.json"))
    base, current, dx2, dz2 = read_state(driven_dir, level, law, 0.01,
                                         manifest_fields(docs/"mouth_intake_redistribution.json"))
    # The stored U=+0.01 envelope represents V=-0.01. Reflecting x
    # gives V0=+0.01 with identical density law, radius and inward rate.
    base = np.flip(base, axis=1).copy()
    current = np.flip(current, axis=0).copy()
    patch = FeedPatch(level, base.shape[0], law)
    raw = j0*abs(base[-1])**2 / abs(np.flip(zero[-1], axis=0))**2
    expected = raw*0.25/(patch.dx**3*raw.sum())
    if ((dx,dz) != (patch.dx,patch.dz) or (dx2,dz2) != (dx,dz)
            or np.max(abs(current-expected)) > 2e-12 or np.max(abs(base[0]-1)) > 1e-12):
        raise ValueError("Base geometry, lab feed or local intake law differs from the archive")
    return patch, base, current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, choices=[24,32], required=True)
    parser.add_argument("--law", choices=["flat_density_slope","restoring_slope"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    patch, base, current = load_preparation(args.level, args.law)
    args.output.mkdir(parents=True, exist_ok=False)
    for omega in (0.1, 0.2):
        p, q, step, err = solve_periodic(patch, base, current, 0.01, omega)
        np.savez_compressed(args.output/f"omega{omega:g}.npz", p=p, q=q, omega=omega,
                            v0=0.01, amplitude=0.002, law=args.law, level=args.level)
        print(f"n={args.level} depth={base.shape[0]-1} law={args.law} V0=0.01 dV=0.002 "
              f"omega={omega:g} crossing_accelerations=+/-{0.002*omega:.7g} "
              f"rate={patch.dx**3*current.sum():.12g} steps={step} update={err:.3g}", flush=True)


if __name__ == "__main__":
    main()
