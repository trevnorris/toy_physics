"""Check the ideal cross terms and compare the guarded carrier readouts."""

import json
from pathlib import Path

import sympy as sp


def ideal_reference():
    # For a regular point-drain potential the earlier full-sphere projection
    # gives <n_x b_n>=U/3. The sphere mean of its harmonic axial velocity is U.
    # Rebuild both from the normal and tangential pieces of the regular
    # l=1 potential; higher Legendre modes are orthogonal in the linear flux.
    m, u, t, rho, a = sp.symbols("mu U vt rho a")
    radial = u*m
    tangential_x = u*(1-m**2)
    pressure = -rho*t*radial*m
    carried = rho*t*(radial*m + tangential_x + radial*m)
    for name, term in (("pressure", pressure), ("carriage", carried)):
        print(name, "outward linear flux =", sp.factor(2*sp.pi*a**2*sp.integrate(term, (m,-1,1))))
    # Uniform pressure gradient on a small control sphere. Applied traction
    # on its contents is -p*n, so a dipole toward +x above ambient pushes -x.
    grad = sp.symbols("grad")
    print("ambient pressure gradient force =", sp.integrate(-2*sp.pi*a**2*m*(grad*a*m), (m,-1,1)))


def compare():
    by_distance = {}
    for distance in (8, 16):
        path = Path(f"_scratch/s11c/companion_uniform_d{distance}_01")
        outcome = json.loads((path/"outcome.json").read_text())
        if outcome["exitCode"] or outcome["childOutcome"]["exitCode"] or not outcome["limitsVerified"]:
            raise ValueError(f"Invalid comparison at {distance}")
        by_distance[distance] = {}
        for line in (path/"stdout").read_text().splitlines():
            if line.startswith("{"):
                row = json.loads(line)
                by_distance[distance][row["match"]] = row
                print(f"d={distance} match={row['match']} "
                      f"carrier_excess_ppm={1e6*row['combined_difference']/row['scale']:.9g} "
                      f"interior_excess_ppm={1e6*row['difference'][2]/row['scale']:.9g}")
        if set(by_distance[distance]) != {"center", "intake_mean"}:
            raise ValueError("Both calibrations are required")
    for match in ("center", "intake_mean"):
        lo, hi = [by_distance[d][match] for d in (8,16)]
        print(f"match={match} carrier_excess16_over_8={hi['combined_difference']/lo['combined_difference']:.9g} "
              f"stress_excess16_over_8={hi['difference'][1]/lo['difference'][1]:.9g} "
              f"paired_flux16_over_8={sum(hi['pair'][:2])/sum(lo['pair'][:2]):.9g}")


if __name__ == "__main__":
    ideal_reference()
    compare()
