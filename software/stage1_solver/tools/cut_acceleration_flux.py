"""Receiving-side momentum current on an accelerating mouth's transverse cut."""

from decimal import Decimal
from pathlib import Path
import re

import sympy as sp


def verify(label, expr, expected):
    value = sp.simplify(expr)
    print(f"{label} = {value}", flush=True)
    if sp.simplify(value - expected) != 0:
        raise ArithmeticError(label)


def local_receiver():
    # Differences between two forced solutions obey the homogeneous equations.
    # Y is a regular l=1 eigenfunction with transverse Laplacian -Y.
    w, r, a = sp.symbols("w r a", positive=True)
    eps, y, yx = sp.symbols("eps Y Y_x", real=True)
    beta, cs2 = sp.Rational(1, 24), sp.Rational(35, 24)
    h = -3*sp.exp(-2*w) + 8*sp.exp(-3*w)
    g = -2*sp.exp(-2*w) + 3*sp.exp(-3*w)
    verify("continuity difference", h.diff(w) + g.diff(w, 2) - g, 0)
    verify("Bernoulli difference", g.diff(w) + cs2*h - beta*(h.diff(w, 2)-h), 0)
    verify("decaying downstream data", sp.limit(h**2+g**2, w, sp.oo), 0)
    density = 1 + eps*h*y
    normal_current = density*(1+eps*g.diff(w)*y)
    adv = normal_current*eps*g*yx
    stress = -beta*density*sp.diff(sp.log(density), w, y)*yx
    # Take w -> 0 from within the receiver. Both signs of eps are allowed.
    def at_cut(expr):
        return sp.diff(expr, eps).subs({eps: 0, w: 0})
    verify("normal-current difference on cut", at_cut(normal_current), 0)
    verify("carried sideways coefficient", at_cut(adv), yx)
    verify("mixed-stress coefficient", at_cut(stress), 3*yx/4)
    # Normalize Y on the spherical rim to cos(theta); divergence theorem.
    mu = sp.symbols("mu", real=True)
    area_integral = 2*sp.pi*a**2*sp.integrate(mu**2, (mu, -1, 1))
    verify("integrated difference per unit perturbation",
           (at_cut(adv)+at_cut(stress)).subs(yx, area_integral), 7*sp.pi*a**2/3)


def archived_cuts():
    record = Path("docs/_measurements/open_sleeve_lateral_balance.txt").read_text()
    # Use the recorded field results at the nearest saved interior cut.
    for drive, heading in [(0, "sleeve_lateral_balance_fine01"), (1, "sleeve_lateral_outlet_control01")]:
        text = record.split("--log-directory _scratch/s11c/"+heading+" --", 1)[1]
        text = text.split("drain_minus_zero level=64", 1)[0]
        pattern = (r"^rate=(\S+) level=(64|128) rate_in/out/far=[^\n]+\n"
                   r"  sleeve_quarter: carried=(\S+) quantum=(\S+) ")
        cases = {(Decimal(rate), int(level)): (Decimal(flow), Decimal(stress))
                 for rate, level, flow, stress in re.findall(pattern,text,re.MULTILINE)}
        if len(cases) != 4:
            raise ValueError("Missing archived cut or control")
        for level in (64, 128):
            flow, stress = cases[Decimal("0.125"), level]
            total = flow+stress
            ref = sum(cases[Decimal(0), level])
            print(f"outlet_tau={drive} level={level}: flow={flow} stress={stress} "
                  f"total={total} drain_minus_zero={total-ref}", flush=True)


if __name__ == "__main__":
    local_receiver()
    archived_cuts()
