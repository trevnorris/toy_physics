"""Stored lateral carrier momentum for a steadily translated draining aperture.

Feed at lateral rest in the lab, aperture speed V, coordinate y=x-Vt.
The stationary equation for the periodic envelope has the extra term
+i V d_y f. Thus f_V equals the archived held-aperture envelope at U=-V.
Constituent mass and hbar are one. Normal depth and all brane coordinates
are integrated; the grid ends at the brane side of the aperture.
"""

from hashlib import sha256
from pathlib import Path

import numpy as np

from mouth_feed_weak_check import dx4
from mouth_feed_weak_flow import FeedPatch
from mouth_release_input_check import manifest_fields


def read_state(folder, n, law, drive, hashes):
    path = folder / f"n{n}_{law}_u{drive:g}.npz"
    if sha256(path.read_bytes()).hexdigest() != hashes[path]:
        raise ValueError(f"Unverified field {path}")
    with np.load(path) as packet:
        if float(packet["speed"]) != drive:
            raise ValueError("Incorrect archived drive")
        return packet["field"], packet["current"], float(packet["dx"]), float(packet["dz"])


def integrate(field, dx, dz, mask=None):
    planes = field if mask is None else field[:, mask]
    totals = planes.sum(axis=tuple(range(1, planes.ndim)))
    return float(dx**3 * dz * (totals[1:-1].sum() + (totals[0]+totals[-1])/2))


def solve_tangent(patch, f):
    """Differentiate the discrete moving-frame field and wall law at V=0."""
    n = abs(f)**2
    face = f[-1]
    amp = abs(face)
    j = patch.current
    cell = patch.dx**3
    total = cell*j.sum()
    normal_log_amp = -patch.restoring*(1-1/amp)
    normal_speed = j/n[-1]
    # elliptic_step always sets its feed to 1. Subtract its zero-source
    # lift so its linearized action here has homogeneous feed data.
    lift = patch.elliptic_step(np.zeros_like(f), np.zeros_like(face))
    drive = -1j*patch.x_derivative(f)
    g = np.zeros_like(f)
    for steps in range(1, 251):
        dn = 2*(f.conj()*g).real
        rhs = ((patch.shift + patch.mu0 - patch.mu0*n**4)*g
               - 4*patch.mu0*n**3*dn*f + drive)
        mean_fraction = cell*np.sum(j*dn[-1]/n[-1])/total
        d_log_amp = -patch.restoring*dn[-1]/(2*amp**3)
        derivative = ((normal_log_amp+1j*normal_speed)*g[-1]
                      + (d_log_amp-1j*normal_speed*mean_fraction)*face)
        updated = patch.elliptic_step(rhs, derivative) - lift
        err = np.max(abs(updated-g))
        if not np.isfinite(err):
            raise ArithmeticError("Nonfinite translational tangent")
        g = updated
        if err < 1e-12:
            return g, steps, err
    raise ArithmeticError(f"Unconverged translation tangent: {err}")


def main():
    records = Path("docs/_measurements")
    zero_hashes = manifest_fields(records/"mouth_feed_weak_run.json")
    driven_hashes = manifest_fields(records/"mouth_intake_redistribution.json")
    scratch = Path("_scratch/s11c")
    for n, reference_run, driven_run in (
            (24, "mouth_feed_weak01", "mouth_redistribution_coarse01"),
            (32, "mouth_feed_depth01", "mouth_redistribution_depth01")):
        for law in ("flat_density_slope", "restoring_slope"):
            f, j, dx, dz = read_state(scratch/reference_run/"fields", n, law, 0, zero_hashes)
            patch = FeedPatch(n, f.shape[0], law)
            if dx != patch.dx or dz != patch.dz or not np.array_equal(j, patch.current):
                raise ValueError("Archived and tangent preparations disagree")
            g, steps, err = solve_tangent(patch, f)
            # dx4 differentiates axis 0; swap it with the depth axis.
            def grad(a):
                return dx4(a.swapaxes(0, 1), dx).swapaxes(0, 1)
            dg = (g.conj()*grad(f) + f.conj()*grad(g)).imag
            intrinsic = integrate(dg, dx, dz)
            # Periodic integration by parts gives the equivalent full-cell
            # expression using only the derivative of the stationary field.
            by_parts = 2*integrate((g.conj()*grad(f)).imag, dx, dz)
            coarse_derivative = (g.conj()*patch.x_derivative(f)
                                 + f.conj()*patch.x_derivative(g)).imag
            coord = (np.arange(n)+0.5-n/2)*dx
            r2 = sum(v**2 for v in np.meshgrid(coord, coord, coord, indexing="ij", sparse=True))
            cuts = [integrate(dg, dx, dz, r2 < radius**2) for radius in (1, 1.5)]
            rest = (f.conj()*grad(f)).imag
            print(f"n={n} depth={f.shape[0]-1} law={law} rate={dx**3*j.sum():.12g} "
                  f"rest_P={integrate(rest, dx, dz):.12g} steps={steps} update={err:.3g} "
                  f"tangent_full={intrinsic:.12g} by_parts={by_parts:.12g} "
                  f"second_order_full={integrate(coarse_derivative, dx, dz):.12g} "
                  f"within_r1={cuts[0]:.12g} within_r1.5={cuts[1]:.12g}", flush=True)
            if abs(intrinsic-by_parts) > 1e-12 or np.max(abs(g[0])) > 1e-12:
                raise ArithmeticError("Periodic momentum identity or lab feed check failed")
            for speed in (0.01, 0.02):
                moved, moved_j, field_dx, field_dz = read_state(
                    scratch/driven_run/"fields", n, law, speed, driven_hashes)
                if (dx != field_dx or dz != field_dz or abs(dx**3*moved_j.sum()-0.25) > 1e-12
                        or np.max(abs(moved[0]-1)) > 1e-12):
                    raise ValueError("Frame change altered grid, intake or feed envelope")
                density = abs(moved)**2
                current = (moved.conj()*grad(moved)).imag
                slope = integrate(current-rest, dx, dz)/(-speed)
                # Independently remove the uniform boost from the held-flow
                # momentum. Without it the full background inventory enters.
                held = integrate(current+speed*density, dx, dz)
                boost = speed*integrate(density, dx, dz)
                lab = integrate(current, dx, dz)
                # Odd reflection extracts the finite-speed translation
                # tangent of the full envelope independently of the solve.
                odd = (moved - np.flip(moved, axis=1))/(-2*speed)
                tangent_error = np.max(abs(odd-g))/np.max(abs(g))
                print(f"n={n} law={law} V={-speed:g} held_P={held:.12g} "
                      f"boost={boost:.12g} lab_P={lab:.12g} secant={slope:.12g} "
                      f"relative_field_tangent_difference={tangent_error:.6g}", flush=True)
                if abs((held-boost)-lab) > 2e-12:
                    raise ArithmeticError("Momentum frame transformation failed")


if __name__ == "__main__":
    main()
