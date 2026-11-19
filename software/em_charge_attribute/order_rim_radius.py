"""Axial boundary work of the order wall, and the finite-width cylindrical trial."""

import json
from pathlib import Path

import sympy as s


def check(label, expr):
    value = s.factor(expr)
    print(f"{label} = {value}", flush=True)
    if value != 0:
        raise ArithmeticError(label)


def sheet_forces():
    # z measures depth from the mouth. Positive axial force is the wall's
    # pull on its near attachment toward increasing z. p,q,j are radius jets.
    r, a, sigma, bend = s.symbols("R a sigma k_b", positive=True)
    p, q, j = s.symbols("p q j", real=True)
    speed = s.sqrt(1+p*p)
    area = 4*s.pi*r*r*speed
    mean = (2/(r*speed)-q/speed**3)/3

    def dz(expr):
        return expr.diff(r)*p + expr.diff(p)*q + expr.diff(q)*j

    def axial(lagrangian):
        moment = lagrangian.diff(q)
        return s.factor(lagrangian - p*(lagrangian.diff(p)-dz(moment)) - q*moment)

    tension = axial(sigma*area)
    curvature = axial(bend*area*mean**2)
    print("tension force =", tension)
    print("bending force at axial tangent =", s.factor(curvature.subs(p, 0)))
    check("integrated tension curvature minus end-flux derivative",
          4*s.pi*sigma*r**2*p*3*mean - dz(tension))
    scale = {r:a*r, q:q/a, j:j/a**2}
    check("similarity tension residual", tension.xreplace(scale)-a*a*tension)
    check("similarity bending residual", curvature.xreplace(scale)-curvature)
    print("straight segment bending force =", curvature.subs({p:0,q:0,j:0}))

    path = Path("software/stage1_solver/frozen/m1c") / (
        "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_10x8.json")
    cfg = json.loads(path.read_text(), parse_float=s.Rational)["config"]["branch"]
    if cfg["geometry_profile"] != "cubic_smoothstep":
        raise ValueError("Unrecognized reference shape")
    b = cfg["r_exit"]/cfg["r_mouth"]
    aspect = (cfg["w_max"]-cfg["w_min"])/cfg["r_mouth"]
    z, length = s.symbols("z L", positive=True)
    t = z/length
    profile = a*(1+(b-1)*(3*t**2-2*t**3))
    rim = {key: profile.diff(z,k).subs(z,0) for k,key in enumerate([r,p,q,j])}
    rim_bend = s.factor(curvature.xreplace(rim))
    print("flare bending force at fixed L =", rim_bend)
    print("similar flare bending / cylinder bending =",
          s.factor(rim_bend.subs(length,aspect*a)/(16*s.pi*bend/9)))
    # Axial virtual displacement eta(z) changes the stretch and both radius
    # derivatives. Integrate its eta'' term by parts once to get the end load.
    density = s.Function("L")(r,p,q)
    stretch_term = density-p*density.diff(p)-2*q*density.diff(q)
    gradient_term = -p*density.diff(q)
    check("axial virtual-work reduction", stretch_term-dz(gradient_term)-axial(density))


def diffuse_cylinder():
    # Normal distance in units of the planar width. q=1/(1+exp(-d/ell)).
    # Relative to its tension the planar energy measure is 6 q(1-q)dq.
    t = s.symbols("t", real=True)
    gen = s.gamma(2+t)*s.gamma(2-t)  # normalized integral of exp(t*d/ell)
    variance = s.simplify(gen.diff(t,2).subs(t,0))
    print("planar normal moments [zero, first, second] =",
          [gen.subs(t,0),gen.diff(t).subs(t,0),variance])
    a, ell, sigma, k = s.symbols("a ell sigma kappa", positive=True)
    force = 4*s.pi*sigma*(a**2 + ell**2*variance)
    print("cylindrical trial axial pull through constant order =", force)
    u, d = s.symbols("u d", real=True)
    slope = u*(1-u)/ell
    second = slope.diff(u)*slope
    well = k*u*u*(1-u)**2/(2*ell**2)
    mu = s.factor(well.diff(u)-k*(second+2*slope/(a+d)))
    print("curved trial order residual / planar slope =", s.factor(mu/slope))
    check("planar kink first integral", k*slope**2/2-well)
    # A symmetric O(ell/a) width distortion changes the tension at second
    # order, which enters the integrated end load at the same constant order.
    eps, c = s.symbols("eps c", real=True)
    width_ratio = 1+c*eps
    rel_tension = (width_ratio+1/width_ratio)/2
    print("trial tension relative correction, coefficient of (ell/a)^2 =",
          s.simplify(rel_tension.diff(eps,2).subs(eps,0)/2))


if __name__ == "__main__":
    sheet_forces()
    diffuse_cylinder()
