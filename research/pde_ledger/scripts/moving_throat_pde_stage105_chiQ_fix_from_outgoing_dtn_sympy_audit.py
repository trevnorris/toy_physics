#!/usr/bin/env python3
"""
Stage 105 SymPy audit.

Provenance notes
----------------
- `Omega_Q = 3 c_s / (2 a)` is the same minimal-isotropic pole scale fixed by
  the Stage 088 conservative quadrupole module.
- `sigma_Q^can = 4 a^5 / (27 c_s^5)` is the canonical outgoing `l=2` DtN odd
  coefficient inherited from the Stage 104 exact fingerprint.
"""

from __future__ import annotations
import sympy as sp

def banner(title: str) -> None:
    line = "=" * 88
    print("\n" + line)
    print(title)
    print(line)

def expect_zero(name: str, expr: sp.Expr) -> None:
    expr = sp.simplify(sp.expand(expr))
    print(f"{name} = {expr}")
    if expr != 0:
        raise AssertionError(f"{name} is not zero")

banner("STAGE 105 — EXACT FIXING OF chi_Q")

a, c_s, omega, chi_Q = sp.symbols("a c_s omega chi_Q", positive=True, real=True)
# Stage 088 minimal isotropic pole scale, carried into the retarded lane.
Omega_Q = sp.simplify(3*c_s/(2*a))
sigma_can = sp.simplify(sp.Rational(9, 8) / Omega_Q**5)

print("Omega_Q =", Omega_Q)
print("sigma_Q^can =", sigma_can)
expect_zero("sigma_Q^can - 4 a^5/(27 c_s^5)", sigma_can - 4*a**5/(27*c_s**5))

Yret = sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2/Omega_Q**2 - sp.I*chi_Q*sigma_can*omega**5)
Yret_series = sp.expand(sp.series(Yret, omega, 0, 6).removeO())

print("\nY_Q^ret(omega) =")
sp.pprint(Yret_series)

expect_zero("omega^2 coefficient", Yret_series.coeff(omega, 2) - a**2/(9*c_s**2))
expect_zero("omega^4 coefficient", Yret_series.coeff(omega, 4) - 4*a**4/(81*c_s**4))
expect_zero("imag omega^5 coefficient", Yret_series.coeff(omega, 5)/sp.I - chi_Q*a**5/(27*c_s**5))

z_out = sp.symbols("z", real=True)
j2_out = (sp.Integer(3)/z_out**3 - 1/z_out)*sp.sin(z_out) - 3*sp.cos(z_out)/z_out**2
y2_out = -((sp.Integer(3)/z_out**3 - 1/z_out)*sp.cos(z_out) + 3*sp.sin(z_out)/z_out**2)
h2_out = sp.simplify(j2_out + sp.I*y2_out)
Lambda_out = sp.simplify(z_out * sp.diff(h2_out, z_out) / h2_out)
Lambda_out_series = sp.expand(sp.series(Lambda_out, z_out, 0, 6).removeO())
Lambda_out_expected = sp.expand(
    -3 + z_out**2/sp.Integer(3) + z_out**4/sp.Integer(9) + sp.I*z_out**5/sp.Integer(9)
)
print("\nLambda_2^out(z) from exact h_2^(1) =")
sp.pprint(Lambda_out_series)
expect_zero("Lambda_2^out series fingerprint", Lambda_out_series - Lambda_out_expected)

Y_out = sp.simplify(sp.Integer(-3) / Lambda_out)
Y_out_series = sp.expand(sp.series(Y_out, z_out, 0, 6).removeO())
outgoing_z5_coeff = sp.simplify(Y_out_series.coeff(z_out, 5) / sp.I)
expect_zero("Y_2^out imag z^5 coefficient - 1/27", outgoing_z5_coeff - sp.Rational(1, 27))
outgoing_omega5_coeff = sp.simplify(outgoing_z5_coeff * (a/c_s)**5)

chi_sol = sp.solve(sp.Eq(Yret_series.coeff(omega, 5)/sp.I, outgoing_omega5_coeff), chi_Q)[0]
print("\nchi_Q from exact outgoing match =", chi_sol)
if sp.simplify(chi_sol - 1) != 0:
    raise AssertionError("chi_Q did not match to 1.")

# General deformed DtN coefficient xi_Q.
z, xi_Q = sp.symbols("z xi_Q", real=True)
Lam_def = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + sp.I*xi_Q*z**5/sp.Integer(9)
Y_def = sp.expand(sp.series(sp.Integer(-3)/Lam_def, z, 0, 6).removeO())
print("\nDeformed DtN normalized branch =")
sp.pprint(Y_def)

expect_zero("deformed z^2 coefficient", Y_def.coeff(z, 2) - sp.Rational(1, 9))
expect_zero("deformed z^4 coefficient", Y_def.coeff(z, 4) - sp.Rational(4, 81))
expect_zero("deformed imag z^5 coefficient", Y_def.coeff(z, 5)/sp.I - xi_Q/sp.Integer(27))

print("\nRESULT:")
print("  On the canonical compact passive/outgoing grouped-P2 branch, chi_Q = 1 exactly.")
print("  More generally, chi_Q is just the leading imaginary DtN deformation coefficient xi_Q.")
