"""Static far field and bend work of the previously tested lateral order load."""

import sympy as sp


def emit(label, value):
    print(f"{label} = {value}", flush=True)


def exterior():
    # Invert the transverse and longitudinal Fourier operators separately.
    # G is I/(4 pi mu r) + (1/B - 1/mu) times Hessian(r)/(8 pi).
    x, y, z, d, b = sp.symbols("x y z d b", positive=True)
    B, mu = sp.symbols("B mu", positive=True)
    xyz = (x, y, z)
    radius = sp.sqrt(sum(t**2 for t in xyz))
    G = (sp.eye(3)/(4*sp.pi*mu*radius)
         + (1/B-1/mu)*sp.hessian(radius, xyz)/(8*sp.pi))
    axial = sp.simplify(G[0, 0].subs({x:d, y:0, z:0}))
    transverse = sp.simplify(G[1, 1].subs({x:d, y:0, z:0}))
    # Work bilinear for two collinear pairs of unit, opposite forces spaced b.
    # d > b. Differentiate its magnitude only to infer force scaling.
    pairs = sp.factor(2*axial-axial.subs(d,d-b)-axial.subs(d,d+b))
    pair_gradient = sp.diff(pairs, d)
    def leading_power(expr):
        return sp.limit(d*sp.diff(expr, d)/expr, d, sp.oo)
    powers = [leading_power(f) for f in (axial, sp.diff(axial,d), pairs, pair_gradient)]
    emit("unit net force kernel, longitudinal and transverse", [axial, transverse])
    emit("two balanced finite pairs, bilinear work", pairs)
    emit("first nonzero coefficient of d^-3", sp.limit(d**3*pairs,d,sp.oo))
    emit("powers, net work, gradient, balanced work, gradient", powers)
    return powers


def bend():
    h, Y, Z, s, q0 = sp.symbols("h Y Z s q0", real=True)
    C, ke, c = [sp.Function(name)(h) for name in ("C", "ke", "c")]
    # Hold the geometric rim at zero. C is the one-sided layer compliance.
    E = (Y-Z)**2/(2*C) + (ke+c)*Z**2/2
    Zstar = sp.solve(sp.diff(E,Z), Z)[0]
    Ered = sp.factor(E.subs(Z,Zstar))
    P = sp.factor(sp.diff(Ered,Y))
    S = ke*Zstar
    H = -c*Zstar
    Q = -sp.diff(Ered,h)
    Ceff = C+1/(ke+c)
    residuals = [sp.factor(v) for v in (
        P+H-S, Q-P**2*sp.diff(Ceff,h)/2,
        sp.diff(Q,Y)+sp.diff(P,h))]
    emit("on-shell energy", Ered)
    emit("balance, compliance and reciprocity residuals", residuals)
    gauge_example = Ered-s*q0*h
    ambiguity = [sp.simplify(sp.diff(gauge_example-Ered,v)) for v in (Y,h)]
    emit("added baseline term, changes in lateral drive force and bend source",
         [ambiguity[0], -ambiguity[1]])

    # Stretch a layer whose tension profile stays fixed in its material labels.
    # alpha=0 translates both ends and alpha=1 holds the far depth fixed.
    u = sp.symbols("u", nonnegative=True)
    L, ell, T = sp.symbols("L ell T", positive=True)
    alpha = sp.symbols("alpha", real=True)
    length = L-alpha*s*ell*h
    compliance = sp.integrate(length/(T*(1+u)**2),(u,0,1))
    emit("compliance of variable tension example", compliance)
    trial = {C:compliance, ke:sp.Integer(2), c:sp.Integer(2)}
    qtrial = -sp.diff(Ered.subs(trial),h)
    rows = []
    for side, anchor, shift in [
        (1,1,"1/10"), (1,1,"-1/10"), (-1,1,"1/10"), (1,1,"1/20"),
        (1,1,"-1/20"), (1,1,"0"), (1,0,"1/10")
    ]:
        case = {s:side, alpha:anchor, Y:sp.Rational(shift), h:0, L:3, ell:1, T:2}
        rows.append([sp.cancel(f.subs(trial).subs(case)) for f in (P,S,H,qtrial)])
    emit("control rows P S H Q", rows)
    assert residuals == [0,0,0] and ambiguity == [0,-s*q0]
    return rows


if __name__ == "__main__":
    powers = exterior()
    rows = bend()
    print("audit " + " ".join(map(str, powers + [v for row in rows for v in row])), flush=True)
    assert powers == [-1,-2,-3,-4]
    assert rows[0][3] and rows[1][3] == rows[0][3] == -rows[2][3]
    assert rows[3][3] == rows[4][3] == rows[0][3]/4
    assert rows[5][3] == rows[6][3] == 0
