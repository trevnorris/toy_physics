#!/usr/bin/env python3
"""
Stage V2-14: outgoing l=2 fingerprint audit.

This script derives the compact outgoing l=2 low-frequency fingerprint from the
3D exterior spherical Hankel solution.  The key convention is that the quoted
fingerprint is the normalized Neumann-to-Dirichlet compliance, not the raw
Dirichlet-to-Neumann admittance.

Frequency convention: exp(-i omega t), outgoing radial wave h_l^(1)(k r).
"""

import sympy as sp

z = sp.symbols("z")
omega, a, c_s = sp.symbols("omega a c_s", positive=True)
ell = sp.Integer(2)
I = sp.I

checks = []

def record(name: str, expr) -> None:
    """Record whether a symbolic expression simplifies to zero."""
    simplified = sp.simplify(expr)
    checks.append((name, simplified == 0, simplified))

# Exact closed form for spherical Hankel h_2^(1)(z).
# h_l^(1)(z)=(-i)^(l+1) exp(i z)/z * sum_{m=0}^l i^m (l+m)!/[m!(l-m)!(2z)^m]
h2 = I * sp.exp(I*z) * (z**2 + 3*I*z - 3) / z**3

# Static exterior l=2 singular behavior is h2 ~ -3 i / z^3.
static_singular = -3*I/z**3
static_normalized_amplitude = sp.simplify(h2/static_singular)

# Raw normalized outgoing Dirichlet-to-Neumann map:
#   f(r)=h_2(k r), z=k a.
#   a f'(a)/f(a) = z h_2'(z)/h_2(z).
# Static l=2 exterior field goes as r^-3, so z h'/h -> -3.
raw_DtN_normalized = sp.simplify(-z/(ell + 1) * sp.diff(h2, z) / h2)

# The fingerprint used in the PN bridge is the inverse normalized compliance:
#   Yhat = 1/raw_DtN_normalized = -(ell+1) h_2/(z h_2').
Yhat = sp.simplify(1/raw_DtN_normalized)
Yhat_alt = sp.simplify(-(ell+1)*h2/(z*sp.diff(h2, z)))

record("Yhat equals inverse normalized raw DtN", raw_DtN_normalized*Yhat - 1)
record("Yhat equals direct compliance formula", Yhat - Yhat_alt)
record("static-normalized outgoing amplitude tends to 1", sp.limit(static_normalized_amplitude, z, 0) - 1)
record("raw normalized DtN tends to 1", sp.limit(raw_DtN_normalized, z, 0) - 1)
record("normalized compliance tends to 1", sp.limit(Yhat, z, 0) - 1)

series_raw = sp.series(raw_DtN_normalized, z, 0, 8).removeO().expand()
series_Y = sp.series(Yhat, z, 0, 8).removeO().expand()

target_through_omega5 = 1 + z**2/sp.Integer(9) + 4*z**4/sp.Integer(81) + I*z**5/sp.Integer(27)
record("Yhat matches target through O(z^5)",
       sp.expand(series_Y - target_through_omega5).series(z, 0, 6).removeO())

record("z^1 coefficient vanishes", series_Y.coeff(z, 1))
record("z^2 coefficient is 1/9", series_Y.coeff(z, 2) - sp.Rational(1, 9))
record("z^3 coefficient vanishes", series_Y.coeff(z, 3))
record("z^4 coefficient is 4/81", series_Y.coeff(z, 4) - sp.Rational(4, 81))
record("z^5 coefficient is i/27", series_Y.coeff(z, 5) - I*sp.Rational(1, 27))

# Convert z = a omega / c_s and extract the leading odd coefficient.
Y_omega_series = sp.expand(series_Y.subs(z, a*omega/c_s))
Gamma5 = sp.simplify(Y_omega_series.coeff(omega, 5) / I)
record("Gamma5 equals a^5/(27 c_s^5)", Gamma5 - a**5/(27*c_s**5))

# Dimensional/channel discrimination.
# For a d-spatial-dimensional outgoing partial wave, the leading absorptive
# small-z power scales as z^(2*l+d-2), since J_nu/Y_nu ~ z^(2 nu)
# with nu = l + (d-2)/2.
d = sp.symbols("d", integer=True, positive=True)
absorptive_power = sp.simplify(2*ell + d - 2)
power_3d = absorptive_power.subs(d, 3)
power_4d = absorptive_power.subs(d, 4)
record("3D spatial l=2 absorptive power is 5", power_3d - 5)
record("4D spatial l=2 absorptive power is 6", power_4d - 6)

# Closed rational form for display.
Yhat_closed = sp.factor(Yhat)
raw_closed = sp.factor(raw_DtN_normalized)

print("STAGE V2-14: OUTGOING l=2 FINGERPRINT AUDIT")
print("=" * 72)
print()
print("Frequency convention: exp(-i omega t)")
print("Dimensionless frequency z = k a = a*omega/c_s")
print()
print("Exact spherical Hankel h_2^(1)(z):")
print("  h2 =", sp.sstr(h2))
print()
print("Static-normalized outgoing amplitude h2 / (-3 i/z^3):")
print("  ", sp.sstr(sp.simplify(static_normalized_amplitude)))
print("  series =", sp.sstr(sp.series(static_normalized_amplitude, z, 0, 8)))
print()
print("Raw normalized outgoing DtN:")
print("  Dhat_raw = - z/3 * h2'(z)/h2(z)")
print("  Dhat_raw =", sp.sstr(raw_closed))
print("  series =", sp.sstr(sp.series(raw_DtN_normalized, z, 0, 8)))
print()
print("Normalized outgoing compliance / inverse DtN:")
print("  Yhat = -3*h2/(z*h2')")
print("  Yhat =", sp.sstr(Yhat_closed))
print("  series =", sp.sstr(sp.series(Yhat, z, 0, 8)))
print()
print("Target through O(z^5):")
print("  1 + z^2/9 + 4*z^4/81 + i*z^5/27 + O(z^6)")
print()
print("Substituting z=a*omega/c_s gives the leading odd coefficient:")
print("  Im coefficient Gamma5 =", sp.sstr(Gamma5))
print()
print("Dimensional/channel discriminator:")
print("  leading absorptive power in d spatial dimensions = 2*l + d - 2")
print("  for d=3, l=2:", power_3d)
print("  for d=4, l=2:", power_4d)
print()
print("CHECKS")
print("-" * 72)
passed = 0
for name, ok, residual in checks:
    status = "PASS" if ok else "FAIL"
    if ok:
        passed += 1
    print(f"{status:4s} | {name}")
    if not ok:
        print("       residual =", sp.sstr(residual))
print("-" * 72)
print(f"checks_total: {len(checks)}")
print(f"checks_passed: {passed}")
print(f"checks_failed: {len(checks)-passed}")
if passed != len(checks):
    raise SystemExit(1)
