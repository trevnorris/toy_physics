"""Local conservation of a candidate rim charge and steady carrier background."""

import sympy as sp


def check(label, actual, expected):
    result = sp.factor(actual)
    print(f"{label} = {result}", flush=True)
    if sp.simplify(result-expected) != 0:
        raise AssertionError(label)


# Cross sections are followed in a meridian coordinate pointing feed -> bulk.
# B is cumulative signed turning, chi drops from 1 to 0.
u, t = sp.symbols("u t", real=True)
B, chi = (sp.Function(name)(u,t) for name in ("B","chi"))
f = sp.Function("f")(B,chi)
line_density = chi*B.diff(u) + f*chi.diff(u)
current = -chi*B.diff(t) - f*chi.diff(t)
wedge = chi.diff(t)*B.diff(u)-chi.diff(u)*B.diff(t)
check("uncompensated continuity coefficient", sp.factor(
    ((chi*B.diff(u)).diff(t)-(chi*B.diff(t)).diff(u))/wedge), 1)
check("local completion continuity coefficient",
      sp.factor((line_density.diff(t)+current.diff(u))/wedge), 1-f.diff(B))

K, s = sp.symbols("K s", real=True)
potential = chi*B+(1-chi)*K
check("excess over neutral completion, density",
      (potential-potential.subs(K,0)).diff(u), -K*chi.diff(u))
check("excess over neutral completion, current",
      -(potential-potential.subs(K,0)).diff(t), K*chi.diff(t))
check("fixed-sector endpoint charge",
      potential.subs(chi,0)-potential.subs({B:0,chi:1},simultaneous=True), K)
check("charge added by conversion on a vertical segment",
      sp.integrate(s-K,(u,1,0)), K-s)
sector = sp.Function("K")(t)
check("deep current if the bound sector varies",
      -potential.subs(K,sector).diff(t).subs({chi.diff(t):0,chi:0}), -sector.diff(t))
shift = sp.Function("Y")(t)
front = sp.Function("c")(u-shift)
check("layer rides the front under translation",
      front.diff(t)+shift.diff(t)*front.diff(u), 0)

# In a broad test transition choose D=1-chi as coordinate. Shape changes at
# fixed lip sign. The free parameter modulates how turning overlaps order loss.
D, alpha = sp.symbols("D alpha",real=True)
turn = s*(D+alpha*D*(1-D))
ordered = sp.integrate((1-D)*turn.diff(D),(D,0,1))
correction = sp.integrate(K-turn,(D,0,1))
check("ordered count when shape is varied",ordered,s*(sp.Rational(1,2)+alpha/6))
check("total through all shape variations", ordered+correction, K)
g, probe = sp.symbols("g probe",real=True)
check("constant probe charge for proposed bound coupling",
      -sp.diff(-probe*g*(ordered+correction),probe),g*K)

# The signed angular area can be constructed from the direction of the
# horizontal normal alone. Use an axisymmetric stretch of its unit sphere.
# With x=cos(theta), signed area per dx dphi is d(n_z)/dx for this map.
x = sp.symbols("x",real=True)
k, a = sp.symbols("k a",positive=True)
nz = k*x/sp.sqrt(1+(k**2-1)*x**2)
check("signed angular Jacobian magnitude", nz.diff(x), k/(1+(k**2-1)*x**2)**sp.Rational(3,2))
for orientation in (1,-1):
    for stretch in (sp.Rational(1,2),1,2):
        degree = orientation*sp.integrate(nz.diff(x).subs(k,stretch),(x,-1,1))/2
        check(f"director degree orientation {orientation:+d}, stretch {stretch}", degree, orientation)
areal_weight=1/(4*sp.pi*a**2)
check("areal weight at double radius relative to original",areal_weight.subs(a,2*a)/areal_weight,
      sp.Rational(1,4))

# Through-flow check in a distinct 1D tube coordinate. Both preparations
# have identical carrier throughput and endpoint densities.
q, F, epsilon = sp.symbols("q F epsilon",real=True)
L, v0 = sp.symbols("L v0",positive=True)
inverse_speed=(1+epsilon*sp.sin(sp.pi*u/L)**2)/v0
check("steady transit charge excess", sp.integrate(q*F*(inverse_speed-1/v0),(u,0,L)),
      q*F*epsilon*L/(2*v0))
qbulk = sp.symbols("q_bulk",real=True)
check("layer accumulation without another outgoing current",q*F-qbulk*F,(q-qbulk)*F)
