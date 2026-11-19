"""First-order virtual work on translated shear contrast and a scattering control."""

import math

import mpmath as mp
import sympy as sp


def emit(label, value):
    print(f"{label} = {value}", flush=True)


# One linear transverse polarization. Time phase t is in radians, z=2 k R.
# dr and dm are the outward jumps in transverse inertia and stiffness.
nu, t, phi = sp.symbols("nu t phi", real=True)
z, k, omega, radius, amp = sp.symbols("z k omega R A", positive=True)
dr, dm = sp.symbols("dRho dMu", real=True)
theta = z * nu / 2 + phi
wave = amp * sp.cos(theta) * sp.cos(t)
velocity = omega * sp.diff(wave, t)
strain = k * sp.diff(wave, phi)
traction_density = (-dr * velocity**2 + dm * strain**2) / 2
average = sp.integrate(traction_density, (t, 0, 2 * sp.pi)) / (2 * sp.pi)
# Resolve the squared displacement and strain into second harmonics before
# integrating. This keeps symbolic phase terms out of the numerical form factor.
fringe = sp.cos(z*nu)*sp.cos(2*phi) - sp.sin(z*nu)*sp.sin(2*phi)
harmonic_average = amp**2 * (-dr*omega**2*(1+fringe) + dm*k**2*(1-fringe)) / 8
assert sp.simplify(sp.expand((average - harmonic_average).rewrite(sp.exp))) == 0
odd = sp.factor((harmonic_average - harmonic_average.subs(nu, -nu)) / 2)
section = sp.factor(2 * sp.pi * radius**2 * sp.integrate(nu * odd, (nu, -1, 1)))
emit("section force per depth, time averaged", section)
mu, rho, eps, m = sp.symbols("mu rho eps m", real=True, positive=True)
weak_jump = section.subs({dr: eps * m * rho, dm: eps * mu, omega**2: k**2 * mu / rho},
                         simultaneous=True)
emit("section force coefficient of weak order depletion", sp.factor(weak_jump / eps))
emit("signed phase controls (+pi/4,-pi/4,0) divided by positive-phase force",
     [sp.simplify(section.subs(phi, v) / section.subs(phi, sp.pi/4))
      for v in (sp.pi/4, -sp.pi/4, 0)])
emit("matched-impedance first-order load", sp.simplify(weak_jump.subs(m, -1)))
emit("one traveling input, averaged closed-sphere dipole",
     sp.integrate(nu, (nu, -1, 1)) * (dm*k**2-dr*omega**2)*amp**2*sp.pi*radius**2/2)

# Strip only the common material, wave and phase factors. Keep the sphere size.
angular_kernel = sp.cancel(section / (sp.pi * amp**2 * (dr * omega**2 + dm * k**2)
                                      * sp.sin(2*phi) / 2))
radial_kernel = angular_kernel.subs(z, 2 * k * radius)
assert radial_kernel.free_symbols == {radius, k}
ratio = radial_kernel.subs(radius, sp.Rational(3, 2)) / radial_kernel.subs(radius, 1)
emit("limit of end-to-midplane load ratio for long wavelengths", sp.limit(ratio, k, 0))
mp.mp.dps = 36
h = sp.lambdify(radius, radial_kernel.subs(k, sp.Rational(1, 2)), "mpmath")


def flare(s):
    return 1 + (3*s**2 - 2*s**3) / 2


f0, f1 = h(1), h(mp.mpf("1.5"))
common = mp.quad(lambda s: h(flare(s)), [0, 1])
bending = mp.quad(lambda s: h(flare(s))*mp.cos(mp.pi*s), [0, 1])
emit("k=1/2 normalized midplane, end, end/midplane, bow work/common work",
     " ".join(mp.nstr(x, 14) for x in (f0, f1, f1/f0, bending/common)))
emit("constant-radius cosine overlap", sp.integrate(sp.cos(sp.pi*nu), (nu, 0, 1)))


def slab_flux(depletion, exponent, phase):
    """Per-area wave momentum flux difference, solving for both outgoing amplitudes.

    Units A=k=omega=rho_out=mu_out=1. Slab faces x=+/-pi/4.
    Plane-wave coefficients inside the slab use the left face as their origin.
    """
    half_width = mp.pi / 4
    stiffness = 1 - depletion
    density = stiffness**exponent
    inner_k = mp.sqrt(density / stiffness)
    impedance = stiffness * inner_k
    propagation = mp.exp(2j * inner_k * half_width)
    left_in = mp.exp(1j*(phase-half_width)) / 2
    right_in = mp.exp(-1j*(phase+half_width)) / 2
    # Unknowns are left_out, inner_rightgoing, inner_leftgoing, right_out.
    matching = mp.matrix([[1, -1, -1, 0], [-1, -impedance, impedance, 0],
                          [0, propagation, 1/propagation, -1],
                          [0, impedance*propagation, -impedance/propagation, -1]])
    rhs = mp.matrix([-left_in, -left_in, right_in, -right_in])
    left_out, _, _, right_out = mp.lu_solve(matching, rhs)
    flux = (abs(left_in)**2+abs(left_out)**2-abs(right_in)**2-abs(right_out)**2)/2
    power_residual = abs(left_in)**2+abs(right_in)**2-abs(left_out)**2-abs(right_out)**2
    return flux, power_residual


residuals = []
for exponent in (0, 1, -1):
    readings = []
    for strength in map(mp.mpf, ("0.02", "0.01", "0.001")):
        pos, power = slab_flux(strength, exponent, mp.pi/4)
        neg, _ = slab_flux(strength, exponent, -mp.pi/4)
        readings.append(pos/strength)
        residuals.extend([abs(power), abs(pos+neg)])
    emit(f"exact planar flux per depletion m={exponent} at .02,.01,.001",
         " ".join(mp.nstr(x, 14) for x in readings))
emit("largest planar power-balance or reflection residual", mp.nstr(max(residuals), 6))
assert max(residuals) < mp.mpf("1e-28")
assert all(math.isfinite(float(value)) for value in (f0, f1, common, bending))
