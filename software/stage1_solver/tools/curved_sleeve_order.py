"""Bare order force on a curved layer and the M1c sleeve shape.

The diffuse force uses a reversible parcel displacement carrying its order
fraction. The thin interface has constant tension. No current is prescribed.
"""

import json
from pathlib import Path

import sympy as sp


BACKGROUND = Path("software/stage1_solver/frozen/m1c") / (
    "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_10x8.json"
)


def main():
    r, w, z, q, n, slope = sp.symbols("r w z q n g", real=True)
    a, k, eos, width, b = sp.symbols("a k K ell b", positive=True)
    mean = sp.Symbol("H", real=True)
    well = a*q**2*(1-q)**2
    density = eos*n**5/4 + well + k*slope**2/2
    mu_n = sp.diff(density, n)
    print("constituent chemical potential and direct order derivative =", [mu_n, sp.diff(mu_n, q)])

    # Signed distance z increases radially outward; the transverse Laplacian
    # here is d_z^2 + div(N) d_z for profiles depending only on z.
    profile = 1/(1+sp.exp(z/width))
    eps = sp.Symbol("eps")
    perturbed = profile - eps*sp.exp(-(z/b)**2)
    dq = sp.diff(perturbed, z)
    mu_order = sp.diff(well, q).subs(q, perturbed) - k*(sp.diff(dq, z) + mean*dq)
    force_normal = mu_order*dq
    mu_at_surface = sp.simplify(mu_order.subs({z: 0, eps: 0}))
    # Gradient of the profile is unaltered at z=0, so normalize the
    # change of the normal force by that initial signed gradient.
    response = sp.simplify(sp.diff(force_normal, eps).subs({z: 0, eps: 0}) / sp.diff(profile, z).subs(z, 0))
    controls = [sp.simplify(response.subs({b: value, k: 2*a*width**2})/a)
                for value in (width, 3*width)]
    print("order derivative at surface, first order force change per profile slope =", [mu_at_surface, response])
    print("normalized change at b=ell and 3ell =", controls)

    radius = sp.Function("R")(w)
    scale = sp.sqrt(1+sp.diff(radius, w)**2)
    nr, nw = 1/scale, -sp.diff(radius, w)/scale
    # Divergence of the unit level-set normal to r-R(w) in four-space.
    curvature = (sp.diff(r**2*nr, r)/r**2 + sp.diff(nw, w)).subs(r, radius)
    axial = -curvature*nw
    print("curvature sum and axial traction per tension =", [sp.factor(curvature), sp.factor(axial)])
    darea_dw = 4*sp.pi*radius**2*scale
    edge = 4*sp.pi*radius**2/scale
    rim_residual = sp.simplify(axial*darea_dw - sp.diff(edge, w))
    print("axial force minus rim derivative per tension =", rim_residual)

    record = json.loads(BACKGROUND.read_text(), parse_float=sp.Rational)
    cfg = record["config"]["branch"]
    if cfg["geometry_profile"] != "cubic_smoothstep":
        raise ValueError("Unsupported archived shape")
    lo, hi = cfg["w_min"], cfg["w_max"]
    x = (w-lo)/(hi-lo)
    shape = cfg["r_mouth"]+(cfg["r_exit"]-cfg["r_mouth"])*x**2*(3-2*x)
    substitutions = {sp.diff(radius, w, 2): sp.diff(shape, w, 2),
                     sp.diff(radius, w): sp.diff(shape, w), radius: shape}
    axial_at_mid = sp.factor(axial.xreplace(substitutions).subs(w, (lo+hi)/2))
    edge_on_shape = edge.xreplace(substitutions)
    integrated = sp.simplify(edge_on_shape.subs(w, hi)-edge_on_shape.subs(w, lo))
    print("archived [midpoint axial traction, integrated axial side force] per tension =", [axial_at_mid, integrated])
    if rim_residual != 0:
        raise ValueError("Side and edge forms disagree")


if __name__ == "__main__":
    main()
