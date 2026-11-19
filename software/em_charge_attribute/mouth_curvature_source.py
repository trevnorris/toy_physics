"""Check the signed curvature density and its boundary and time currents on a graph."""

import itertools

import sympy as sp


def check(label, actual, expected):
    print(f"{label} = {sp.factor(actual)}", flush=True)
    if sp.simplify(actual-expected) != 0:
        raise AssertionError(label)


# Graph area times the determinant of the shape operator, with the normal
# pointing toward the same +w reference on either face.
p = sp.Matrix(sp.symbols("p0:3", real=True))
metric = sp.eye(3)+p*p.T
check("graph metric determinant", metric.det(), 1+p.dot(p))
# det(second form) = det(Hessian)/(1+p.p)^(3/2); divide by sqrt(det(metric)).
weight = (1+p.dot(p))**-2/sp.pi**2

# Spatial Piola identity for a general Hessian, with fully symmetric third
# derivatives. This makes the following identity valid even at zero slope.
x = sp.symbols("x y z", real=True)
h = sp.Function("H")(*x)
cof = sp.hessian(h,x).cofactor_matrix()
for i in range(3):
    check(f"cofactor row {i} divergence",
          sum(sp.diff(cof[i,j], x[j]) for j in range(3)).expand(), 0)
m = sp.Matrix(3,3,sp.symbols("M0:9"))
residual = m*m.cofactor_matrix().T - m.det()*sp.eye(3)
check("cofactor contraction squared residual", sum(v.expand()**2 for v in residual), 0)

# Radial primitive F(p) in slope space: div_p F = (1+|p|^2)^-2 / pi^2.
# The signed slope v gives the flux through a radially oriented spherical
# boundary in real space, increasing radius counted positive.
v = sp.symbols("v", real=True)
boundary = 4/sp.pi * sp.integrate(v**2/(1+v**2)**2,v)
check("signed boundary function", boundary, 2/sp.pi*(sp.atan(v)-v/(1+v**2)))
check("slope-space primitive regular at zero", sp.limit(boundary/(4*sp.pi*v**3),v,0),
      weight.subs(dict.fromkeys(p,0))/3)
check("flat-center to flat-infinity control", boundary.subs(v,0)-boundary.subs(v,0), 0)
for angle in (sp.pi/4, sp.pi/3):
    inward = -sp.tan(angle)
    result = boundary.subs(v,0)-boundary.subs(v,inward)
    check(f"upper mouth with lip angle {angle}", result, (2*angle-sp.sin(2*angle))/sp.pi)

# Each sphere is oriented outward from the exterior domain; excised lip
# boundaries therefore enter with a minus sign. Bounded remote slope and
# angular derivatives do not change the limiting radial direction at the lip.
for signs in [(1,),(-1,),*itertools.product((1,-1),repeat=2)]:
    terms = [-sp.limit(boundary,v,-s*sp.oo) for s in signs]
    check(f"enclosed exterior curvature for {signs}", sum(terms), sum(signs))

# The two boundary terms distinguish the integrated exterior density from
# its local potential evaluated on the outer sphere. C sets the far bend.
r, C = sp.symbols("r C",positive=True)
outer = boundary.subs(v,-C/r**2)
check("large-sphere geometric boundary flux", sp.limit(outer,r,sp.oo), 0)
check("upper tail times r^6", sp.limit(r**6*outer,r,sp.oo), -4*C**3/(3*sp.pi))

# Move the graph in time: rho = w(p) det(Dp) and j = -w(p) cof(Dp)^T dp/dt.
# Jacobi's determinant derivative and the cofactor row identity cancel the
# second-derivative terms. The following is the remaining weight derivative.
velocity = sp.Matrix(sp.symbols("v0:3"))
grad_weight = sp.Matrix([weight.diff(pi) for pi in p])
transport = (m.T*grad_weight).dot(m.cofactor_matrix().T*velocity)
check("time-current weight contraction residual",
      sp.factor(transport-m.det()*grad_weight.dot(velocity)), 0)

# Independent mechanical control. The scalar exterior has compliance 1/a
# in these units; minimize over outward depth and reflector position. The
# radius is controlled in this test; it is not selected by wave energy.
z, tip = sp.symbols("rim tip", real=True)
for a, action, expected_length in [(1,sp.Integer(2),2),(2,sp.Rational(25,6),sp.Rational(5,2))]:
    cavity = tip-z
    energy = a*z**2/2+cavity+action/cavity+(tip-3)**2/2
    load = 1-action/expected_length**2
    state = {z:load/a, tip:3-load}
    residuals = [energy.diff(q).subs(state) for q in (z,tip)]
    check(f"stationary energy residual at radius {a}", sum(v**2 for v in residuals), 0)
    check(f"reflecting length at radius {a}", cavity.subs(state), expected_length)
    print("radius, action, wave energy, scalar flux =",
          (a, action, (action/cavity).subs(state), a*state[z]), flush=True)
