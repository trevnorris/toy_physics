"""Calibrate a 0PN release probe against measured held-mouth loads and derive its readout."""

from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from stage1_solver.mouth_probe import spherical_rim_probe


ARCHIVE = Path("docs/_measurements/held_force_incident_flow.json")


def probe_scales():
    rows = {r["name"]:r for r in json.loads(ARCHIVE.read_text())["runs"]}
    mu, weights = leggauss(24)
    for d in (8, 16):
        run = rows[f"L96C6_d{d}"]
        if run["guard_exit"] or run["worker_exit"] or not run["limits_verified"]:
            raise ValueError("Unverified companion load")
        reading = json.loads(run["stdout"][-1])
        if reading["grid"] != [96, 6, 24] or reading["d"] != d:
            raise ValueError("Wrong release preparation")
        (path, digest), = run["output_fields_sha256"].items()
        if sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError("Paired field checksum mismatch")
        with np.load(path, allow_pickle=False) as p:
            if not np.array_equal(p["centers"], [-d, 0]) or not np.array_equal(p["rates"], [.25, .25]):
                raise ValueError("Mismatched center or intake")
            print(f"d={d} packet_keys={','.join(sorted(p.files))} "
                  f"centers={p['centers'].tolist()} rates={p['rates'].tolist()}")
        load = reading["force"]
        for factor in (0, -1, -.5, .5, 1):
            applied = factor*load
            traction = spherical_rim_probe(1., mu, applied)
            # Polar axis points toward the companion. dA=2*pi*a^2*dmu.
            net = 2*np.pi*np.dot(weights*mu, traction)
            monopole = 2*np.pi*np.dot(weights, traction)
            if abs(net-applied)>2e-14*load or abs(monopole)>2e-14*load:
                raise ArithmeticError("Rim force or zero mean traction check failed")
            print(f"d={d} fraction={factor:g} applied={applied:.12g} "
                  f"integrated={net:.12g} mean_traction_integral={monopole:.4g}")


def release_readout():
    # This is an identifiability calculation. All acceleration readings are symbols.
    f, h = sp.symbols("f h", nonzero=True)
    a0, chi, c2, c3, c4, c5 = sp.symbols("a0 chi c2 c3 c4 c5")
    acceleration = a0+chi*f+c2*f**2+c3*f**3+c4*f**4+c5*f**5
    def slope(step):
        return sp.expand((acceleration.subs(f,step)-acceleration.subs(f,-step))/(2*step))
    full, half = slope(h), slope(h/2)
    combined = sp.expand((4*half-full)/3)
    if combined != chi-c5*h**4/4:
        raise ArithmeticError("Incorrect small-force extrapolation")
    print(f"full_signed_slope={full} half_signed_slope={half} extrapolated_slope={combined}")
    # Define positive equivalent force as the applied force producing the same
    # acceleration increment. The canceling force has the opposite sign.
    drift, fraction = sp.symbols("a_iso eta", real=True)
    load = sp.Symbol("L", nonzero=True)
    drive = chi*f+a0
    equivalent = (a0-drift)/chi
    for expression in (drive.subs(f,-equivalent)-drift,
                       equivalent.subs(a0,drift+chi*fraction*load)/load-fraction):
        if sp.simplify(expression)!=0:
            raise ArithmeticError("Release force equivalence failed")
    print("linear_response: equivalent_force=(a0-a_iso)/chi; "
          "canceling_force=-equivalent_force; eta=equivalent_force/L")
    print("measured_free_center_accelerations=unavailable "
          "measured_force_to_acceleration_response=unavailable")


if __name__ == "__main__":
    probe_scales()
    release_readout()
