"""Virtual work per unit three-area on a flat patch of the conversion front."""

import sympy as sp


def check(label, measured, reference):
    error = sp.factor(measured - reference)
    print(f"{label}: measured={measured}, reference={reference}, residual={error}")
    if error != 0:
        raise ValueError(label)


# q increases toward bulk. t is the ordered fraction at coincident centers.
# Translating chi or the shear profile through fixed space acts as -d/dq
# on that field alone. The two derivatives below keep the other field live.
t = sp.symbols("t", positive=True)
width, peak = sp.symbols("delta b0", positive=True)
dtdq = -t * (1 - t) / width
measure = -1 / dtdq  # bounds reversed so t runs from 0 to 1
b = sp.Function("b")(t)
front_density = dtdq * b
profile_density = t * sp.diff(b, t) * dtdq
check("common shift is a spatial derivative",
      front_density + profile_density, sp.diff(t * b, t) * dtdq)
check("local force magnitude / local stored energy",
      sp.factor(-front_density / (t * b)), (1 - t) / width)


def integrate_over_patch_depth(density):
    return sp.factor(sp.integrate(sp.cancel(density * measure), (t, 0, 1)))


ratios = []
for exponent in (1, 2):
    # One tangential displacement component with localized normal gradient.
    bare_energy = peak * (4 * t * (1 - t)) ** exponent
    energy = integrate_over_patch_depth(t * bare_energy)
    front = integrate_over_patch_depth(dtdq * bare_energy)
    profile = integrate_over_patch_depth(t * sp.diff(bare_energy, t) * dtdq)
    ratio = sp.factor(-front / energy)
    print(f"p={exponent}: E/A={energy}, F_front/A={front}, "
          f"F_profile/A={profile}, -F_front/E={ratio}")
    check(f"p={exponent} rigid shift", profile, -front)
    check(f"p={exponent} relaxed front load", front.subs(peak, 0), 0)
    ratios.append(ratio)

print(f"equal stored energy, narrow / broad force = {ratios[1] / ratios[0]}")
