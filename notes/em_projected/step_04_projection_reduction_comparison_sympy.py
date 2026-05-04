#!/usr/bin/env python3
"""Compare projection-first and reduction-first Maxwell couplings in the zero-mode regime.

Model assumptions for this comparison:
- zero-mode field strength F^{mu nu}(x,w) = f^{mu nu}(x)
- mixed sector suppressed: F^{w nu} = 0, J^w = 0
- projected gauge-driver term set to zero by gauge choice / consistent source sector
- source profile factorizes as J^nu(x,w) = j^nu(x) S(w)

Then the projected inhomogeneous law becomes
    I_WZ * partial_mu f^{mu nu} = mu0 * I_WS * j^nu,
with
    I_WZ = ∫ W(w) Z(w) dw,
    I_WS = ∫ W(w) S(w) dw.

So the projection-first effective coupling to the measured zero-mode field is
    mu0_eff^(proj) = mu0 * I_WS / I_WZ,
while the action reduction in the EM paper gives
    mu0_eff^(red) = mu0 / Z_int.
"""
from __future__ import annotations

import sympy as sp


def line(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.simplify(expr)
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    residue = sp.simplify(expr)
    if residue == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


w, lam, mu0 = sp.symbols("w lam mu0", positive=True)
Z = sp.exp(-w**2 / lam**2)
Z_int = sp.simplify(sp.integrate(Z, (w, -sp.oo, sp.oo)))
Z2_int = sp.simplify(sp.integrate(Z**2, (w, -sp.oo, sp.oo)))

line("1) Generic zero-mode projection formula")
print("Assume F^{mu nu}(x,w) = f^{mu nu}(x), J^nu(x,w) = j^nu(x) S(w), and ∫W(w)dw = 1.")
print("Then measured zero-mode fields satisfy Proj_W[F^{mu nu}] = f^{mu nu}.")
print("The projected inhomogeneous law becomes:")
print("  I_WZ * partial_mu f^{mu nu} = mu0 * I_WS * j^nu")
print("with")
print("  I_WZ := ∫ W(w) Z(w) dw")
print("  I_WS := ∫ W(w) S(w) dw")
print()
print("So the projection-first effective coupling to the measured field is:")
print("  mu0_eff^(proj) = mu0 * I_WS / I_WZ")
print()
print("The reduction-first coupling from the localized action is:")
print("  mu0_eff^(red)  = mu0 / Z_int")
print()
print("These coincide only if I_WS / I_WZ = 1 / Z_int.")

x, k = sp.symbols("x k", real=True)
sigma, tau = sp.symbols("sigma tau", positive=True)
eta = sp.symbols("eta", real=True, nonzero=True)
f_test = sp.sin(k * x) + x**2
j_test = sp.cos(k * x) + x
df_test = sp.diff(f_test, x)
W_smooth = sp.exp(-w**2 / sigma**2) / (sp.sqrt(sp.pi) * sigma)
S_smooth = sp.exp(-w**2 / tau**2) / (sp.sqrt(sp.pi) * tau)
I_WZ_smooth = sp.simplify(sp.integrate(W_smooth * Z, (w, -sp.oo, sp.oo)))
I_WS_smooth = sp.simplify(sp.integrate(W_smooth * S_smooth, (w, -sp.oo, sp.oo)))
projected_zero_mode_residual = sp.simplify(
    sp.integrate(W_smooth * (Z * df_test - mu0 * S_smooth * j_test), (w, -sp.oo, sp.oo))
)
assert_zero(
    "smooth I_WZ overlap",
    I_WZ_smooth - lam / sp.sqrt(lam**2 + sigma**2),
)
assert_zero(
    "smooth I_WS overlap",
    I_WS_smooth - 1 / (sp.sqrt(sp.pi) * sp.sqrt(sigma**2 + tau**2)),
)
assert_zero(
    "smooth zero-mode projected residual",
    projected_zero_mode_residual - (I_WZ_smooth * df_test - mu0 * I_WS_smooth * j_test),
)

F_w_dependent = f_test + eta * x * w**2
field_mutation_lhs = sp.simplify(sp.integrate(W_smooth * Z * sp.diff(F_w_dependent, x), (w, -sp.oo, sp.oo)))
field_mutation_delta = sp.simplify(field_mutation_lhs - I_WZ_smooth * df_test)
assert_zero(
    "w-dependent field mutation amplitude",
    field_mutation_delta - eta * lam**3 * sigma**2 / (2 * (lam**2 + sigma**2) ** sp.Rational(3, 2)),
)
assert_nonzero("w-dependent field mutation breaks zero-mode projection", field_mutation_delta)

J_w_dependent = S_smooth * j_test + eta * x * w**2 * S_smooth
source_mutation_rhs = sp.simplify(sp.integrate(W_smooth * mu0 * J_w_dependent, (w, -sp.oo, sp.oo)))
source_mutation_delta = sp.simplify(source_mutation_rhs - mu0 * I_WS_smooth * j_test)
assert_zero(
    "w-dependent source mutation amplitude",
    source_mutation_delta
    - eta * mu0 * x * sigma**2 * tau**2 / (2 * sp.sqrt(sp.pi) * (sigma**2 + tau**2) ** sp.Rational(3, 2)),
)
assert_nonzero("w-dependent source mutation breaks factorized source projection", source_mutation_delta)
print("Concrete smooth-profile projection residual checks pass, including a w-dependent-field mutation guard.")

line("2) Gaussian localization data")
print("Take Z(w) = exp(-w^2/lambda^2).")
print("  Z_int  =", sp.sstr(Z_int))
print("  Z2_int =", sp.sstr(Z2_int))
print("  Z2_int / Z_int =", sp.sstr(sp.simplify(Z2_int / Z_int)))
assert_zero("Gaussian Z_int", Z_int - sp.sqrt(sp.pi) * lam)
assert_zero("Gaussian Z2_int", Z2_int - sp.sqrt(2 * sp.pi) * lam / 2)

line("3) Case A: matched projection kernel W = Z / Z_int and a brane-localized source S = delta(w)")
W_match = sp.simplify(Z / Z_int)
S_delta = sp.DiracDelta(w)
I_WZ_match = sp.simplify(sp.integrate(W_match * Z, (w, -sp.oo, sp.oo)))
I_WS_match = sp.simplify(W_match.subs(w, 0))  # use the defining delta-sampling property explicitly
mu0_proj_match = sp.simplify(mu0 * I_WS_match / I_WZ_match)
mu0_red = sp.simplify(mu0 / Z_int)

print("  I_WZ(match) =", sp.sstr(I_WZ_match))
print("  I_WS(match) =", sp.sstr(I_WS_match))
print("  mu0_eff^(proj, match) =", sp.sstr(mu0_proj_match))
print("  mu0_eff^(red)         =", sp.sstr(mu0_red))
print("  ratio proj/reduction  =", sp.sstr(sp.simplify(mu0_proj_match / mu0_red)))
assert_zero("matched observer overlap equals Z2_int / Z_int", I_WZ_match - Z2_int / Z_int)
assert_zero("matched observer I_WZ", I_WZ_match - sp.sqrt(2) / 2)
assert_zero("delta-source ratio", sp.simplify(mu0_proj_match / mu0_red) - sp.sqrt(2))
assert_nonzero("mutated delta-source ratio should fail", sp.simplify(mu0_proj_match / mu0_red) - 1)
print()
print("So with the natural matched observer kernel and a delta-localized source,")
print("projection gives mu0 / Z2_int, not mu0 / Z_int.")

line("4) Case B: regularized sharp observer and regularized sharp source")
eps = sp.symbols("epsilon", positive=True)
W_eps = sp.exp(-w**2 / eps**2) / (sp.sqrt(sp.pi) * eps)
S_eps = W_eps
I_WZ_eps = sp.simplify(sp.integrate(W_eps * Z, (w, -sp.oo, sp.oo)))
I_WS_eps = sp.simplify(sp.integrate(W_eps * S_eps, (w, -sp.oo, sp.oo)))
mu0_proj_eps = sp.simplify(mu0 * I_WS_eps / I_WZ_eps)

print("Use normalized Gaussians W_eps = S_eps = exp(-w^2/eps^2)/(sqrt(pi) eps).")
print("  I_WZ(eps) =", sp.sstr(I_WZ_eps))
print("  I_WS(eps) =", sp.sstr(I_WS_eps))
print("  mu0_eff^(proj, eps) =", sp.sstr(mu0_proj_eps))
print("  I_WZ(eps -> 0) =", sp.sstr(sp.limit(I_WZ_eps, eps, 0, dir='+')))
print("  I_WS(eps) diverges as eps -> 0, showing why exact delta/delta sampling is not a finite coupling.")
assert_zero("regularized sharp observer samples Z(0)", sp.limit(I_WZ_eps, eps, 0, dir='+') - 1)
assert_zero("regularized source overlap", I_WS_eps - sp.sqrt(2) / (2 * sp.sqrt(sp.pi) * eps))
assert_nonzero(
    "mutated regularized observer sample should fail",
    sp.limit(I_WZ_eps, eps, 0, dir='+') - sp.sqrt(2) / 2,
)
print()
print("This replaces the distribution-undefined delta squared case with an explicit regulator.")

line("5) Compact comparison")
print("Reduction-first:")
print("  mu0_eff = mu0 / Z_int")
print()
print("Projection-first, zero-mode limit:")
print("  mu0_eff = mu0 * (∫ W S) / (∫ W Z)")
print()
print("Matched Gaussian observer + delta source:")
print("  mu0_eff =", sp.sstr(mu0_proj_match), "= sqrt(2) * mu0 / (sqrt(pi)*lambda)")
print()
print("So projection does not simply reproduce the reduced Maxwell coupling.")
print("It yields a family of observer-dependent effective laws unless extra closure conditions are imposed.")
print("STATUS: PASS")
