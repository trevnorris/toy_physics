"""Inspect the signed runs and derive the spherical brane-displacement readout."""

from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

import sympy as sp


def inspect_runs():
    manifest = Path("docs/_measurements/mouth_prescribed_acceleration.txt")
    lines = set(manifest.read_text().splitlines())

    def verified(path):
        if f"{sha256(path.read_bytes()).hexdigest()}  {path}" not in lines:
            raise ValueError(f"File differs from the recorded run: {path}")
        return path

    source = verified(Path("software/stage1_solver/tools/mouth_prescribed_acceleration.py"))
    seed = verified(Path("_scratch/s11c/mouth_feed_depth01/fields/n32_restoring_slope_u0.npz"))
    with ZipFile(seed) as data:
        print("initial packet arrays:", sorted(Path(name).stem for name in data.namelist()))
    fields = set()
    for tag, acceleration in (("zero", 0), ("positive", .02), ("negative", -.02),
                              ("positive_half", .01), ("negative_half", -.01)):
        folder = Path(f"_scratch/s11c/mouth_prescribed_accel_{tag}02")
        outcome = json.loads(verified(folder/"outcome.json").read_bytes())
        if not outcome["limitsVerified"] or outcome["exitCode"] or outcome["childOutcome"]["exitCode"]:
            raise ValueError(f"Unverified trajectory {tag}")
        command = json.loads((folder/"invocation.json").read_bytes())["command"]
        data = json.loads(verified(folder/"stdout").read_bytes())
        if command[1] != str(source) or float(command[2]) != acceleration or data["acceleration"] != acceleration:
            raise ValueError(f"Different driver or acceleration in {tag}")
        for grid in ("16", "32"):
            fields.update(tuple(sorted(row)) for row in data[grid])
        fine = data["32"]
        print(f"A={acceleration:+g} sample_times={[row['t'] for row in fine]} "
              f"commanded_final_X={fine[-1]['x']:+g}")
    print("per-sample keys across all five runs and both time refinements:", sorted(fields))


def reconstruct_shear():
    r, theta, mu = sp.symbols("r theta mu_R", positive=True)
    p, d = (sp.Function(name)(r) for name in ("p", "d"))
    ur, ut = p*sp.cos(theta), -d*sp.sin(theta)
    # Force delivered to the brane exterior. Its inner normal is -e_r.
    curl_phi = (sp.diff(r*ut, r)-sp.diff(ur, theta))/r
    def angular(f):
        return sp.integrate(2*sp.pi*sp.sin(theta)*f, (theta, 0, sp.pi))
    radial = sp.simplify(3*angular(ur*sp.cos(theta))/(4*sp.pi))
    tangential = sp.simplify(-3*angular(ut*sp.sin(theta))/(8*sp.pi))
    force = sp.factor(mu*r**2*angular(sp.sin(theta)*curl_phi))
    print("recovered radial and tangential dipole amplitudes:", radial, tangential)
    print("axial rotational shear force on exterior:", force)
    expected = 8*sp.pi*mu*r*(p-sp.diff(r*d,r))/3
    if radial != p or tangential != d or sp.simplify(force-expected) != 0:
        raise ArithmeticError("Dipole readout or force has incorrect normalization")
    # A genuine rigid translation has zero curl. A potential gradient also
    # has zero curl even when its radial and tangential components differ.
    f = sp.Function("f")(r)
    for normal, tangent in ((1,1), (sp.diff(f,r),f/r)):
        if sp.simplify(force.subs({p:normal,d:tangent}).doit()) != 0:
            raise ArithmeticError("A curl-free displacement produced shear traction")


if __name__ == "__main__":
    inspect_runs()
    reconstruct_shear()
