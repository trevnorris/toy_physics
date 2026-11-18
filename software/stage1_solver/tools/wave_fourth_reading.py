"""Fit three real-Laplace readings of the archived reduced rim trajectory."""

import json
from pathlib import Path

import sympy as sp


def main():
    archive = json.loads(Path("docs/_measurements/trapped_wave_push.json").read_text())
    series = archive["runs"]["sympy"]["stdout_claim_lines"]
    # Preserve the source trajectory without executing text from the archive.
    archived = "-alpha*f*cos(z) + alpha*f - f*z**2*(alpha - 1)/2"
    if f"normalized rim trajectory = {archived}" not in series:
        raise ValueError("Archived reduced-model rim trajectory has changed")
    z, p, u = sp.symbols("z p u", positive=True)
    alpha, f = sp.symbols("alpha f", nonzero=True, real=True)
    rim = alpha*f*(1-sp.cos(z)) + f*z**2*(1-alpha)/2
    displacement = sp.laplace_transform(rim, z, p, noconds=True)
    mass = sp.factor((f/p)/(p**2*displacement))
    m0 = sp.limit(mass, p, 0, dir="+")
    m1, m2 = mass.subs(p, u), mass.subs(p, 2*u)
    ratio = sp.factor((m0-m1)/(m0-m2))
    gap_squared = sp.factor(4*u**2*(1-ratio)/(4*ratio-1))
    wave_mass = sp.factor((m0-m1)*(gap_squared+u**2)/u**2)
    rim_mass = sp.factor(m0-wave_mass)
    predicted = sp.factor(rim_mass + wave_mass*gap_squared/(gap_squared+9*u**2))
    print("three readings =", m0, m1, m2)
    print("inferred squared rate in relative-mode units =", gap_squared)
    print("fourth predicted from three =", predicted)
    # The unused rate enters the data reduction only after the prediction.
    measured = sp.factor(mass.subs(p, 3*u))
    residual = sp.factor(measured-predicted)
    print("fourth read from the trajectory =", measured)
    print("fourth residual =", residual)
    flat_deficits = [sp.simplify((m0-m).subs(alpha, 0)) for m in (m1,m2)]
    print("flat-response calibration deficits =", flat_deficits)
    if residual != 0 or flat_deficits != [0,0]:
        raise ArithmeticError("Holdout or unresolved-response control failed")


if __name__ == "__main__":
    main()
