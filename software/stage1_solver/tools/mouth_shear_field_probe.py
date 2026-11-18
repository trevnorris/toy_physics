"""Drive independent 0PN shear fields along the archived one-mouth accelerations."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from stage1_solver.mouth_shear import MouthShearField


ARCHIVE = Path("docs/_measurements/mouth_prescribed_acceleration.txt")
ARCHIVE_HASH = "db7ae4e93ffb1ddf105d2afd02194f29009464d9a642721a18e3edd3aeeeae73"
GAPS = (0.000025, 0.0001)
SPACING = 0.000004


def archived_drives():
    raw = ARCHIVE.read_bytes()
    if sha256(raw).hexdigest() != ARCHIVE_HASH:
        raise ValueError("Different prescribed-acceleration benchmark")
    drives = [json.loads(row) for row in raw.splitlines() if row.startswith(b'{"acceleration":')]
    if [d["acceleration"] for d in drives] != [0, .02, -.02, .01, -.01]:
        raise ValueError("Incomplete four-drive family and stationary control")
    for d in drives:
        rows = d["32"]
        np.testing.assert_allclose([r["t"] for r in rows], [.00025, .0005, .001], rtol=1e-14, atol=0)
        for row in rows:
            np.testing.assert_allclose(row["x"], d["acceleration"] * row["t"]**2 / 2,
                                       rtol=1e-14, atol=0)
            np.testing.assert_allclose(row["inward_rate"], .25, rtol=1e-14, atol=0)
    return drives


def check_retarded_field():
    """Derive the two displacement projections and the cut force from a vector potential."""
    a, r, c, rho, t = sp.symbols("a r c rho t", positive=True)
    b = sp.Function("B")(t - (r-a)/c)
    potential = -a**3 * (b + r/c * sp.diff(b, t)) / r**2
    wave = sp.diff(potential, t, 2)/c**2 - (sp.diff(r**2 * sp.diff(potential, r), r)
                                            - 2 * potential)/r**2
    p = 2 * potential/r
    d = sp.diff(r * potential, r)/r
    amplitude = (a/r)**3 * (b + r/c * sp.diff(b, t))
    curl_force = 8 * sp.pi * rho * c**2 * r * (p - sp.diff(r*d, r))/3
    scale = 8 * sp.pi * rho * a**3/3
    expected_force = scale * (sp.diff(b, t, 2) + r/c * sp.diff(b, t, 3))
    residuals = [wave, p + 2*amplitude, d - amplitude - (a/r)**3 * (r/c)**2 * sp.diff(b, t, 2),
                 curl_force - expected_force]
    if any(sp.simplify(expr) != 0 for expr in residuals):
        raise ArithmeticError("Outgoing wave, displacement or signed curl traction failed")
    print("SymPy wave, radial displacement, tangential displacement, projected force residuals: 0 0 0 0")


def run_drive(drive, density, speed, drag):
    acc = drive["acceleration"]
    samples = drive["32"]
    field = MouthShearField(1., density, speed, drag, lambda t: acc*t, samples[-1]["t"])
    rows = []
    for sample in samples:
        t = sample["t"]
        _, rim_d = field.displacement(t, 0.)
        rim_v = field.rim_velocity(t)
        rim_force = field.drag * (acc*t - rim_v)
        # Sample beyond the disturbance as a separate causality control.
        np.testing.assert_array_equal(field.displacement(t, 2*speed*t), (0., 0.))
        shells = []
        for gap in GAPS:
            p, d = field.displacement(t, gap)
            coarse = field.shell_force(t, gap, SPACING)
            fine = field.shell_force(t, gap, SPACING/2)
            # Centered stencils on a smooth part of the outgoing wave.
            force = (4*fine - coarse)/3
            shells.append(dict(gap=gap, p=p, d=d, force_coarse=coarse, force_fine=fine,
                               shear_force=force))
        rows.append(dict(t=t, x=sample["x"], inward_rate=sample["inward_rate"], rim_d=rim_d,
                         rim_v=rim_v, rim_force=rim_force, shells=shells))
    return dict(slip_drag=drag, acceleration=acc, samples=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--density", type=float, required=True)
    parser.add_argument("--wave-speed", type=float, required=True)
    parser.add_argument("--slip-drag", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    drives = archived_drives()
    check_retarded_field()
    if max(GAPS) + SPACING >= min(s["t"] for s in drives[0]["32"]) * args.wave_speed:
        raise ValueError("Selected shell stencils must be inside the wavefront at each sample")
    results = [run_drive(d, args.density, args.wave_speed, drag)
               for drag in args.slip_drag for d in drives]
    payload = dict(archive_sha256=ARCHIVE_HASH, radius=1., density=args.density,
                   wave_speed=args.wave_speed, spacing=SPACING, runs=results)
    with args.output.open("x") as file:
        json.dump(payload, file, indent=2, allow_nan=False)
        file.write("\n")
    print("t=0.001 gap=0.0001; drag acceleration X equatorial_material_displacement force_on_exterior")
    for result in results:
        row = result["samples"][-1]
        shell = row["shells"][-1]
        print(f"{result['slip_drag']:g} {result['acceleration']:+g} {row['x']:+.12g} "
              f"{shell['d']:+.12g} {shell['shear_force']:+.12g}")


if __name__ == "__main__":
    main()
