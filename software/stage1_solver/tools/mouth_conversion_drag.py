"""Momentum and normal-traction compatibility of the two ideal exteriors."""

import sympy as sp


def check(label, actual, expected):
    residual = sp.factor(actual - expected)
    print(f"{label}: {sp.factor(actual)} | {sp.factor(expected)} | residual {residual}", flush=True)
    if residual != 0:
        raise ArithmeticError(label)


def balances():
    rate = sp.Symbol("D", positive=True)
    brane, bulk, mouth, boost = sp.symbols("U_br U_bulk V boost", real=True)
    carried, interstress, other_feed, other_bulk = sp.symbols("C T K_feed K_bulk", real=True)
    force_feed, force_bulk = sp.symbols("F_feed F_bulk", real=True)
    # Supplied ideal outer integrals, in the translating mouth frame.
    # C is the same carried rate at the two sides of local conversion.
    # T is inter-medium stress transfer, distinct from force on the mouth.
    # K_feed and K_bulk are transfers to independently supported surroundings.
    incoming = rate*(brane-mouth)
    outgoing = rate*(bulk-mouth)
    equations = [incoming-carried-interstress-force_feed-other_feed,
                 carried+interstress-outgoing-force_bulk-other_bulk]
    solution = sp.solve(equations, (force_feed, force_bulk), dict=True)[0]
    print("required feed stress reaction =", sp.factor(solution[force_feed]), flush=True)
    print("required bulk stress reaction =", sp.factor(solution[force_bulk]), flush=True)
    total = sp.factor((force_feed+force_bulk).subs(solution))
    check("total reaction including other receiver", total,
          rate*(brane-bulk)-other_feed-other_bulk)
    check("dependence on local conversion carriage", sp.diff(total, carried), 0)
    check("dependence on internal stress transfer", sp.diff(total, interstress), 0)
    check("common change of reference speed",
          total.subs({brane:brane-boost, bulk:bulk-boost, mouth:mouth-boost})-total, 0)
    check("ideal tow through resting media",
          total.subs({brane:0, bulk:0, other_feed:0, other_bulk:0}), 0)
    check("ideal held brane stream",
          total.subs({bulk:0, mouth:0, other_feed:0, other_bulk:0}), rate*brane)

    # Keep intake redistribution as well as the changed velocity. Premise 3
    # uses these same material quantities on both sides of each transfer.
    eps, j0, j1, vx0, vx1 = sp.symbols("eps j0 j1 vx0 vx1", real=True)
    carry_density = (j0+eps*j1)*(vx0+eps*vx1)
    variation = sp.diff(carry_density, eps).subs(eps, 0)
    check("first variation of transferred carriage", variation, j0*vx1+j1*vx0)
    area0, area1 = sp.symbols("area0 area1", real=True)
    live_transfer = carry_density*(area0+eps*area1)
    check("live carried first variation including measure",
          sp.diff(live_transfer, eps).subs(eps, 0),
          area0*(j0*vx1+j1*vx0)+area1*j0*vx0)
    # Carriage in a uniformly shifted frame changes by local transferred mass.
    shifted = (j0+eps*j1)*(vx0+eps*vx1-boost)
    check("carriage frame change", shifted-carry_density, -boost*(j0+eps*j1))


def pressure_projection():
    r, radius = sp.symbols("r a", positive=True)
    mu = sp.Symbol("mu", real=True)
    p0, p1, h0, h1 = (sp.Function(s)(r) for s in ("p0", "p1", "h0", "h1"))
    # Upper bulk above graph w=h, n dA=(-grad h,1) d^3x.
    # Pressure on the material below is -p n; its x projection is p*h_x.
    angular = p1*h0.diff(r)*mu**2 + p0*(h1.diff(r)*mu**2+h1*(1-mu**2)/r)
    shell = sp.factor(2*sp.pi*r**2*sp.integrate(angular, (mu, -1, 1)))
    print("projected pressure per radial interval =", shell, flush=True)
    boundary = 4*sp.pi*p0*r**2*h1/3
    interior = 4*sp.pi*r**2*(p1*h0.diff(r)-p0.diff(r)*h1)/3
    check("live pressure and shape, integration by parts", shell, boundary.diff(r)+interior)
    check("flat unchanged face", shell.subs({h0:0, h1:0}).doit(), 0)
    # An admissible normal-force moment, with no asserted fluid solution.
    slope, pressure_gradient = sp.symbols("k_h k_p", real=True)
    example = shell.subs({h0:slope*r**2/2, p1:pressure_gradient*r, h1:0}).doit()
    check("curved face pressure moment", sp.integrate(example, (r, 0, radius)),
          4*sp.pi*slope*pressure_gradient*radius**5/15)


def radius_requirement():
    d, rho, dt, ds, mt, ms = sp.symbols("d rho D_t D_s M_t M_s", positive=True)
    common_rate = sp.Symbol("lambda", positive=True)
    # Independent M values denote the measured scalar slow response.
    # This branch additionally assumes all source intake reaches the far monopole.
    incident = ds/(4*sp.pi*rho*d**2)
    acceleration = dt*incident/mt
    coefficient = sp.factor(d**2*acceleration/ms)
    print("mass-normalized inverse-square coefficient =", coefficient, flush=True)
    check("common intake per inertia condition",
          coefficient.subs({dt:common_rate*mt, ds:common_rate*ms}),
          common_rate**2/(4*sp.pi*rho))


if __name__ == "__main__":
    print("SymPy version =", sp.__version__, flush=True)
    balances()
    pressure_projection()
    radius_requirement()
