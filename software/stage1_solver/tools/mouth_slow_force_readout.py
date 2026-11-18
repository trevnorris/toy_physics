"""Project a signed harmonic reaction and test the slow extrapolation on saved carrier fluxes."""

from pathlib import Path
import hashlib
import json

import sympy as sp


def phase_readout():
    phase = sp.symbols("phase", real=True)
    amp, w = sp.symbols("A w", positive=True)
    k, drag, mass, offset = sp.symbols("K G M H", real=True)
    x = amp * sp.cos(phase)
    v = w * x.diff(phase)
    force = offset + k * x + drag * v + mass * w * v.diff(phase)
    cosine, sine = [sp.integrate(force * basis, (phase, 0, 2 * sp.pi)) / (amp * sp.pi)
                    for basis in (sp.cos(phase), sp.sin(phase))]
    print("force cosine and sine per displacement amplitude =", cosine, sine, flush=True)
    print("static-subtracted mass and velocity coefficient =", (k - cosine) / w**2, -sine / w,
          flush=True)
    if sp.simplify((k - cosine) / w**2 - mass) != 0 or sp.simplify(-sine / w - drag) != 0:
        raise ValueError("Phase signs or normalization disagree")


def carrier_readout():
    manifest = json.loads(Path("docs/_measurements/lateral_mouth_acceleration.json").read_text())
    folder = Path("_scratch/s11c/lateral_accel_e24_c6_01")
    record, = [run for run in manifest["runs"] if str(folder) in run["command"].split()]
    raw = (folder / "stdout").read_bytes()
    if hashlib.sha256(raw).hexdigest() != record["stdout_sha256"]:
        raise ValueError("Carrier source does not match the recorded run")
    rows = {}
    for line in raw.decode().splitlines():
        if line.startswith("omega="):
            parts = dict(item.split("=") for item in line.split())
            w = sp.Rational(parts["omega"])
            rows[w] = -sp.Rational(parts["carried_per_A"]) - sp.Rational(parts["stress_per_A"])
    if sorted(rows) != [0, sp.Rational(1, 200), sp.Rational(1, 100)]:
        raise ValueError("Expected the static-limit calculation and both finite rates")
    zero, w1, w2 = sorted(rows)
    for w in (w1, w2):
        print("frequency, carrier drive per acceleration, in-phase force per displacement after static subtraction =",
              float(w), float(rows[w]), float(-w**2 * rows[w]), flush=True)
    predicted = (w2**2 * rows[w1] - w1**2 * rows[w2]) / (w2**2 - w1**2)
    print("carrier slow extrapolation, archived slow tangent, signed difference =",
          float(predicted), float(rows[zero]), float(predicted - rows[zero]), flush=True)


if __name__ == "__main__":
    phase_readout()
    carrier_readout()
