"""Quasistatic work of a broad imposed bend, relative to a slab without a mouth."""

import sympy as sp


def check(label, result, expected):
    print(f"{label} = {sp.factor(result)}", flush=True)
    if sp.simplify(result - expected) != 0:
        raise AssertionError(label)


# Spherically symmetric midplane displacement on a < r < R. The same driver
# moves the sphere at R in both preparations. Everything outside it cancels.
a, R, kappa, r = sp.symbols("a R kappa r", positive=True)
z, b, Y, rest, eta = sp.symbols("z b Y rest eta", real=True)
shape = b + (z-b)*(1/r-1/R)/(1/a-1/R)
shell = sp.integrate(2*sp.pi*kappa*r**2*sp.diff(shape, r)**2, (r, a, R))
K = 4*sp.pi*kappa/(1/a-1/R)
check("shell after blank subtraction", shell, K*(z-b)**2/2)
check("blank constant interior gradient", sp.diff(b, r), 0)

# Wall and fixed-action wave act across the same length. The positive rear
# spring is a test support, whose reference can follow the imposed bend.
ke, kb, T, length, J, nu, gap = sp.symbols("k k_b T L J nu gap", positive=True)
wave = J*sp.sqrt(gap**2+(nu/length)**2)
core = T*length + wave
P, kc = sp.symbols("P k_c", real=True)
core_function = sp.Function("V")
potential = ke*(z-b)**2/2 + core_function(Y-z) + kb*(Y-rest-eta*b)**2/2
print("wave length stiffness =", sp.factor(wave.diff(length, 2)), flush=True)
check("wave curvature times positive factor",
      wave.diff(length, 2)*length**3*(gap**2*length**2+nu**2)**sp.Rational(3, 2),
      J*nu**2*(3*gap**2*length**2+2*nu**2))
check("net tension plus wave compression force", core.diff(length)-wave.diff(length), T)

# Use independent length during differentiation to express the constitutive
# derivatives without binding a specific cavity spectrum to the response test.
coords = (z, Y)
hessian = sp.hessian(potential, coords)
core_second = hessian[0, 0] - ke
hessian = hessian.subs(core_second, kc)
coupling = sp.Matrix([potential.diff(q, b) for q in coords])
rates = -hessian.inv()*coupling
ks = 1/(1/ke+1/kc+1/kb)
check("rim motion per unit bend", rates[0], 1-(1-eta)*ks/ke)
check("reflector motion per unit bend", rates[1], eta+(1-eta)*ks/kb)
check("cavity length change per unit bend", rates[1]-rates[0], -(1-eta)*ks/kc)
on_branch = {z: b+P/ke, Y: rest+eta*b-P/kb}
slope = potential.diff(b).subs(on_branch, simultaneous=True)
check("blank-subtracted work slope at equilibrium", slope, (eta-1)*P)
check("relaxed work curvature", potential.diff(b, 2)+coupling.dot(rates), ks*(1-eta)**2)
check("shared translation work slope", slope.subs(eta, 1), 0)

# Illustrative dimensionless equilibrium: k = kb = J*nu = 1, gap = 0,
# T = 2, rest = 3, eta = 0. Both equilibrium forces are checked at each drive.
example_core = core.subs({T:2, J:1, nu:1, gap:0})
check("example core energy at positive length", example_core, 2*length+1/length)
energy = (z-b)**2/2 + example_core.subs(length, Y-z) + (Y-3)**2/2
forces = [energy.diff(q) for q in coords]
baseline = energy.subs({b:0,z:1,Y:2})
print("b       rim_shift      reflector_shift  wave_change    excess_work", flush=True)
for drive in map(sp.Rational, ("0.02", "-0.02", "0.01", "-0.01")):
    rim, tip = sp.nsolve([f.subs(b, drive) for f in forces], coords, (1, 2), prec=45)
    values = {b:drive,z:rim,Y:tip}
    change = energy.subs(values)-baseline
    readings = [drive, rim-1, tip-2, 1/(tip-rim)-1, change]
    print(" ".join(f"{float(x):+.10f}" for x in readings), flush=True)
    if tip <= rim or any(abs(f.subs(values)) > 1e-35 for f in forces):
        raise ArithmeticError("unbalanced or nonphysical root")
