"""Predict archived two-drain ideal forces using only one-drain calibrations."""

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re

from mpmath import mp


def archived_stdout(run, record):
    folder = Path("_scratch/s11c") / run
    outcome = json.loads((folder / "outcome.json").read_text())
    if (outcome["exitCode"] or outcome["childOutcome"]["exitCode"]
            or not outcome["limitsVerified"]):
        raise ValueError(f"Unverified archive: {folder}")
    path = folder / "stdout"
    data = path.read_bytes()
    fingerprint = f"Stdout {path} ({len(data)} bytes, SHA256 {sha256(data).hexdigest()})."
    if fingerprint not in record.splitlines():
        raise ValueError(f"Output differs from recorded run: {path}")
    return data.decode().splitlines()


def calibrate(lines):
    rows = {}
    for line in lines:
        if line.startswith("U="):
            speed, values = line.split(" holder [W=0, W=U/2, W=U]: ")
            rows[Fraction(speed[2:])] = tuple(Fraction(x) for x in values.strip("[]").split(", "))
    full = Fraction(1, 50)
    if set(rows) != {full, -full, full/2, -full/2}:
        raise ValueError("Incomplete signed one-drain preparations")
    # H is the external holder force. Brane-only forcing selects A.
    # Common incident speed selects A+B (equivalent to opposite towing
    # speed only when both reservoirs and all other data transform).
    def intrinsic_slope(column):
        return (rows[-full][column]-rows[full][column])/(2*full)

    brane = intrinsic_slope(0)
    common = intrinsic_slope(2)
    bulk = common-brane
    print(f"brane slope={brane}, common slope={common}, inferred bulk slope={bulk}")
    for speed, values in sorted(rows.items()):
        for ratio, held in zip((0, Fraction(1, 2), 1), values, strict=True):
            if abs(speed) == full and ratio in (0, 1):
                continue
            predicted = (brane+bulk*ratio)*speed
            print(f"held check U={speed} W/U={ratio}: predicted pull={predicted}, "
                  f"recorded pull={-held}, residual={predicted+held}")
            if predicted != -held:
                raise ArithmeticError("Signed stream control")
    return brane, bulk


def compare_pairs(lines, brane, bulk):
    mp.dps = 40
    # Independent inputs of the archived two-source quadrature; mass rates
    # in each medium share units, volume densities in 3D and 4D do not.
    source_rate = mp.mpf(3)/8
    rho3, rho4 = 2, 3
    a_br = mp.mpf(brane.numerator)/brane.denominator
    a_bulk = mp.mpf(bulk.numerator)/bulk.denominator
    cuts = []
    for line in lines:
        if not line.startswith('"distance, cut radius, outward brane and bulk flux, total holder: "'):
            continue
        values = line.split("InputForm[{")[1].removesuffix("}]")
        d, r, flux_br, flux_bulk, held = (
            mp.mpf(re.sub(r"`[\d.]+", "", s.strip()).replace("*^", "e"))
            for s in values.split(","))
        cuts.append((d, r))
        incident_br = source_rate/(4*mp.pi*rho3*d**2)
        incident_bulk = -source_rate/(2*mp.pi**2*rho4*d**3)
        predictions = (a_br*incident_br, a_bulk*incident_bulk)
        observed = (-flux_br, -flux_bulk, -held)
        pulls = (*predictions, sum(predictions))
        errors = [abs(calc-obs)/abs(obs) for calc, obs in zip(pulls, observed)]
        print(f"d={int(d)} R={int(r)}: predicted brane,bulk,total="
              f"{','.join(mp.nstr(x,17) for x in pulls)}, "
              f"recorded total={mp.nstr(-held,17)}, "
              f"largest channel/total relative difference={mp.nstr(max(errors),5)}")
        # The original record labels its displayed quadratures at 15 digits.
        if max(errors) > mp.mpf("1e-14"):
            raise ArithmeticError("Independent surface quadrature mismatch")
    if sorted(cuts) != [(8, 1), (8, 3), (16, 1), (16, 3)]:
        raise ValueError("Incomplete two-drain quadratures")


if __name__ == "__main__":
    record = Path("docs/_measurements/mouth_outer_holding_force.txt").read_text()
    calibration = archived_stdout("mouth_outer_hold_sympy01", record)
    brane_slope, bulk_slope = calibrate(calibration)
    compare_pairs(archived_stdout("mouth_outer_hold_pair01", record), brane_slope, bulk_slope)
