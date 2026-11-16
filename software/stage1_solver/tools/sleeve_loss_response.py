"""Instantaneous bare relaxation of an ordered-to-disordered curved profile.

Normal distance d increases outwards. The order field initially depends on d
alone. The saved sleeve supplies a reference geometry, not kinetics.
"""

import json

import sympy as sp

from curved_sleeve_order import BACKGROUND


def main():
    d, w, c = sp.symbols("d w c", real=True)
    ell, k, mobility, density, extra = sp.symbols("ell k M n lam", positive=True)
    u, v = sp.symbols("u V", real=True)
    angular, meridian = sp.symbols("k_ang k_mer", real=True)
    jacobian = (1+d*angular)**2*(1+d*meridian)
    curvature = sp.diff(jacobian, d)/jacobian
    # Radial part of the Euler derivative with volume measure J dd.
    q = sp.Function("q")(d)
    well = k*c**2*(1-c)**2/(2*ell**2)
    mu_order = sp.diff(well, c).subs(c, q) - k*sp.diff(jacobian*q.diff(d), d)/jacobian
    gradient = -c*(1-c)/ell
    mu = mu_order.xreplace({q.diff(d, 2): sp.diff(gradient, c)*gradient,
                            q.diff(d): gradient, q: c})
    bare_loss = sp.factor(density*mobility*mu)
    expected = -density*mobility*k*curvature*gradient
    profile_residual = sp.factor(bare_loss-expected)
    force = mu*gradient
    speed = sp.solve((u-v)*gradient + mobility*mu, v)[0]
    speed_residual = sp.factor(speed-u+mobility*k*curvature)
    work_residual = sp.factor(force*(speed-u)-mobility*mu**2)
    print("signed-distance Jacobian =", jacobian)
    print("loss minus n*M*k*H(d)*q*(1-q)/ell =", profile_residual)
    print("contour speed minus material speed plus M*k*H(d) =", speed_residual)
    print("force*(V-u) minus M*mu_order^2 =", work_residual)

    # An illustrative allowed change in the unspecified rate, with the same
    # profile and stress at this instant. It vanishes at the layer midpoint.
    extra_loss = density*extra*c*(1-c)*(1-2*c)**2
    print("extra loss / (n*lam) at q=1/2 and 1/4 =",
          [sp.simplify(extra_loss.subs(c, level)/(density*extra))
           for level in (sp.Rational(1, 2), sp.Rational(1, 4))])

    branch = json.loads(BACKGROUND.read_text(), parse_float=sp.Rational)["config"]["branch"]
    if branch["geometry_profile"] != "cubic_smoothstep":
        raise ValueError("Unsupported sleeve profile")
    s = (w-branch["w_min"])/(branch["w_max"]-branch["w_min"])
    radius = branch["r_mouth"]+(branch["r_exit"]-branch["r_mouth"])*s**2*(3-2*s)
    scale = sp.sqrt(1+radius.diff(w)**2)
    principal = {angular: 1/(radius*scale), meridian: -radius.diff(w, 2)/scale**3}
    midpoint = (branch["w_min"]+branch["w_max"])/2
    h_mid = sp.simplify(curvature.subs(d, 0).subs(principal).subs(w, midpoint))
    coefficient = sp.simplify(bare_loss.subs({d: 0, c: sp.Rational(1, 2)})
                              .subs(principal).subs(w, midpoint)*ell/(density*mobility*k))
    print("archived midpoint H and loss*ell/(n*M*k) =", [h_mid, coefficient])
    if any(value != 0 for value in (profile_residual, speed_residual, work_residual)):
        raise ValueError("Local profile, kinematics or work identity failed")


if __name__ == "__main__":
    main()
