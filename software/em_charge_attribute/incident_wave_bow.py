"""Project a standing wave's confining-wall force and check two-mode preparations."""

import sympy as sp


def emit(label, value):
    print(f"{label} = {value}", flush=True)


# Instantaneous real scalar displacement on the wall.  nu is the direction
# cosine from the wave axis, z = 2 k r.  The 1/2 in the wave Hamiltonian is
# applied after this angular integral.
nu, phase = sp.symbols("nu phase", real=True)
z = sp.symbols("z", positive=True)
# Split the square into even and odd terms before quadrature.
intensity = (1 + sp.cos(z * nu) * sp.cos(2 * phase)
             - sp.sin(z * nu) * sp.sin(2 * phase)) / 2
angular = sp.factor(2 * sp.pi * sp.integrate(nu * intensity, (nu, -1, 1)))
emit("angular integral of direction cosine times squared wave", angular)
emit("angular integral at z=pi and phases +/-pi/4, 0",
     [angular.subs({z: sp.pi, phase: p}) for p in (sp.pi/4, -sp.pi/4, 0)])
emit("angular integral at zero lateral wave number", sp.limit(angular, z, 0))

# Uniform effective wall, joined smoothly at w=0, zero axial end traction.
# L and I are the half-length and the physical dipole inertia per depth.
y = sp.symbols("y", real=True)  # w/L
L, Iwall = sp.symbols("L Iwall", positive=True)
peak = sp.symbols("S", real=True)
bow = sp.cos(sp.pi * y)
shapes = sp.Matrix([1, bow])
def axial_average(expr):
    return sp.integrate(sp.expand_trig(expr), (y, -1, 1)) / 2

mass = 2 * L * Iwall * (shapes * shapes.T).applyfunc(axial_average)
def generalized_load(profile):
    return 2 * L * peak * (shapes * profile).applyfunc(axial_average)

load = generalized_load((1 + bow) / 2)
accels = sp.simplify(mass.inv() * load)
emit("physical generalized forces for common shift and bow", list(load))
emit("initial modal accelerations", list(accels))
contrast = (bow.subs(y, 0) - (bow.subs(y, -1) + bow.subs(y, 1)) / 2) * accels[1]
emit("midplane minus mean end acceleration", contrast)
emit("constant-depth control accelerations", list(mass.inv() * generalized_load(1)))
emit("bow derivatives at both ends", [sp.diff(bow, y).subs(y, v) for v in (-1, 1)])

# A separate quantized, number-conserving diagnostic.  Build creation/annihilation
# matrices on two number ladders.  Three levels suffice for the N=2 restriction
# of a_plus^dagger a_minus and for lowering to vacuum.
a = sp.zeros(3)
for n in (1, 2):
    a[n - 1, n] = sp.sqrt(n)
right = sp.kronecker_product(a, sp.eye(3))
left = sp.kronecker_product(sp.eye(3), a)
transfer = right.H * left
sector = [6, 4, 2]  # |2,0>, |1,1>, |0,2>
quad = ((transfer - transfer.H) / (2 * sp.I)).extract(sector, sector)
pair = sp.Matrix([0, 1, 0])
def pair_mean(operator):
    return sp.simplify((pair.H * operator * pair)[0])

emit("one each: mean, second moment of odd fringe quadrature",
     [pair_mean(quad), pair_mean(quad**2)])
weights = {}
for value, _, vectors in quad.eigenvects():
    weights[value] = sum(sp.simplify(abs((v.normalized().H * pair)[0])**2)
                         for v in vectors)
emit("quadrature eigenvalues and one-each probabilities", sorted(weights.items()))
theta = sp.symbols("theta", real=True)
positive_frequency = sp.exp(sp.I * theta) * right + sp.exp(-sp.I * theta) * left
annihilation = sp.expand((positive_frequency**2)[0, 4])
emit("vacuum to one-each lowering matrix element (bra vacuum)", annihilation)
emit("angular dipole of the two-photon annihilation amplitude",
     2 * sp.pi * sp.integrate(nu * annihilation, (nu, -1, 1)))

# Structural checks, including symmetry and the local equation at release.
assert sp.trigsimp(intensity - sp.cos(z*nu/2 + phase)**2, method="fu") == 0
assert sp.simplify(angular + angular.subs(phase, -phase)) == 0
assert sp.simplify(Iwall * (shapes.dot(accels)) - peak * (1 + bow) / 2) == 0
assert sum(weights.values()) == 1 and quad == quad.H
