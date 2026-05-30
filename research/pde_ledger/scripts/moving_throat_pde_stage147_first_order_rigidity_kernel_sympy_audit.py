#!/usr/bin/env python3
from __future__ import annotations
import sympy as sp

def banner(t):
    print("\n" + "="*88)
    print(t)
    print("="*88)

Pi = sp.symbols("Pi", positive=True, real=True)
x, eps = sp.symbols("x eps", real=True)
kap = sp.pi/2
c = sp.cos(sp.pi*x/2)
Kq = sp.cosh(kap*(1-x))/sp.cosh(kap)

gPi = 2*Pi*(2*Pi*sp.exp(Pi)+sp.pi)/((4*Pi**2+sp.pi**2)*(sp.exp(Pi)-1))
Sformula = sp.simplify(
    Pi*(kap*sp.tanh(kap) + Pi*(sp.exp(-Pi)*sp.sech(kap) - 1))
    / ((1-sp.exp(-Pi))*(kap**2-Pi**2))
)

rF1 = sp.sqrt(12*sp.Rational(37,20)**2/sp.pi**2 - 1)
gminus = sp.simplify(rF1 - sp.sqrt(1+rF1**2)/2)
Pi_star = sp.N(sp.nsolve(gPi - gminus, 1.5), 40)
g_star = sp.N(gPi.subs(Pi, Pi_star), 40)
S_star = sp.N(Sformula.subs(Pi, Pi_star), 40)
gp_star = sp.N(sp.diff(gPi, Pi).subs(Pi, Pi_star), 40)
Sp_star = sp.N(sp.diff(Sformula, Pi).subs(Pi, Pi_star), 40)

Sigma_star = sp.N(Pi_star/(1-S_star/4), 40)
T_star = sp.N(sp.sqrt(9*Sigma_star/20), 40)

AT = sp.N(
    -(sp.Rational(9,1)/(40*T_star)) * (
        1/(gp_star*(1-S_star/4)) + Pi_star*Sp_star/(4*gp_star*(1-S_star/4)**2)
    ),
    30
)
BT = sp.N(
    (sp.Rational(9,1)/(40*T_star)) * Pi_star/(4*(1-S_star/4)**2),
    30
)

banner("Stage 147 audit: first-order rigidity coefficients")
print("Pi_* =", Pi_star)
print("Sigma_* =", Sigma_star)
print("T_* =", T_star)
print("A_T =", AT)
print("B_T =", BT)
print("|A_T|/B_T =", sp.N(abs(AT)/BT, 20))

# --- Audit assertions: numerical anchor against paper-quoted literals ---
AT_paper = sp.Float("-4.27263956256927", 30)
BT_paper = sp.Float("0.134875005736706", 30)
ratio_crosscheck = sp.Float("31.6785", 20)
assert abs(sp.N(AT) - AT_paper) < sp.Float("1e-12", 30), \
    f"A_T deviates from paper-quoted value: {sp.N(AT)} vs {AT_paper}"
print("PASS: A_T matches paper-quoted -4.27263956256927 to 1e-12")
assert abs(sp.N(BT) - BT_paper) < sp.Float("1e-12", 30), \
    f"B_T deviates from paper-quoted value: {sp.N(BT)} vs {BT_paper}"
print("PASS: B_T matches paper-quoted 0.134875005736706 to 1e-12")
assert abs(sp.N(sp.Abs(AT)/BT) - ratio_crosscheck) < sp.Float("1e-3", 30), \
    f"|A_T|/B_T computed ratio cross-check deviates from 31.6785: {sp.N(sp.Abs(AT)/BT)}"
print("PASS: |A_T|/B_T computed ratio cross-check (not a paper literal) matches 31.6785 to 1e-3")

# --- Audit assertion: A_T from automatic differentiation of T_m(Pi) (independent) ---
# Independent primitive: build T_m as a SYMBOLIC function of the free symbol Pi from
# Sformula (NOT from the cached numeric S_star/Sp_star/T_star), and let SymPy derive
# dT_m/dPi automatically. The retuning identity (appendix eq:app-part04-deltaPi-firstorder)
# gives A_T = -(dT_m/dPi)/(dg/dPi) at Pi_*. sp.diff regenerates the (1-S/4)^2 factor and
# the S' factor on its own, so a hand-written sign/power error in the closed-form AT
# (py:33-38) cannot be reproduced here.
Tm_of_Pi = sp.sqrt(sp.Rational(9, 20) * (Pi / (1 - Sformula/4)))
dTm_dPi = sp.diff(Tm_of_Pi, Pi)
dg_dPi = sp.diff(gPi, Pi)
AT_autodiff = sp.N(-(dTm_dPi.subs(Pi, Pi_star)) / (dg_dPi.subs(Pi, Pi_star)), 30)
AT_30 = sp.N(AT, 30)
assert abs(AT_autodiff - AT_30) < sp.Float("1e-20", 30), \
    f"A_T autodiff route disagrees with closed form: {AT_autodiff} vs {AT_30}"
print("PASS: A_T closed form agrees with autodiff of T_m(Pi) (residual < 1e-20)")

# Weight kernel representation
gbar, Sbar = sp.symbols("gbar Sbar", real=True)
dT = sp.simplify(eps*(AT*(gbar-gminus) + BT*(Sbar-Sformula.subs(Pi, Pi_star))))
print("delta T_m =")
sp.pprint(dT)

# Centered kernel form: constant pieces cancel against normalization.
Wcenter = sp.simplify(AT*(c-gminus) + BT*(Kq-Sformula.subs(Pi, Pi_star)))
print("Centered rigidity kernel W_*(x) =")
sp.pprint(Wcenter)

# --- Audit assertion: rigidity-kernel projection identity (independent quadrature) ---
# Notes sec.2 / appendix eq:app-part04-deltaT-firstorder: the traction shift equals the
# projection of W_*(x) against (smallsigma - Sigma_*). Test it for a concrete, normalized,
# NON-canonical deformation smallsigma(x) = 2x (integral_0^1 2x dx = 1, positive, and a
# different shape from Sigma_*). LHS by numerical quadrature of the kernel projection;
# RHS by the two-moment coefficient formula. The two sides share NO derivation step:
# LHS integrates W_* numerically, RHS uses the algebraic A_T, B_T and the moment
# integrals. This also exercises the "integrates to zero" centering claim, since the
# -g_*, -S_* constants inside W_* drop out only because integral(smallsigma-Sigma_*)=0.
Sigma_star_x = Pi_star * sp.exp(-Pi_star * x) / (1 - sp.exp(-Pi_star))
smallsigma_x = 2*x
# normalization sanity (must integrate to 1 each):
norm_s = sp.N(sp.integrate(smallsigma_x, (x, 0, 1)), 40)
norm_Sigma = sp.N(sp.integrate(Sigma_star_x, (x, 0, 1)), 40)
assert abs(norm_s - 1) < sp.Float("1e-30", 40) and abs(norm_Sigma - 1) < sp.Float("1e-30", 40), \
    f"deformation/source not normalized: {norm_s}, {norm_Sigma}"
Wstar_x = AT*(c - g_star) + BT*(Kq - S_star)
lhs_proj = sp.N(sp.integrate(Wstar_x * (smallsigma_x - Sigma_star_x), (x, 0, 1)), 30)
gbar_s = sp.N(sp.integrate(smallsigma_x * c, (x, 0, 1)), 40)
Sbar_s = sp.N(sp.integrate(smallsigma_x * Kq, (x, 0, 1)), 40)
rhs_moment = sp.N(AT*(gbar_s - g_star) + BT*(Sbar_s - S_star), 30)
assert abs(lhs_proj - rhs_moment) < sp.Float("1e-22", 30), \
    f"kernel projection != two-moment formula: {lhs_proj} vs {rhs_moment}"
print("PASS: kernel projection of W_* reproduces two-moment traction shift (residual < 1e-22)")
# --- Audit assertion: source-centering of the rigidity kernel (independent) ---
# CONSULT Q5 (batch 6): the projection identity above is BLIND to the centering
# constants -A_T*g_*, -B_T*S_* -- they vanish against (smallsigma - Sigma_*) because
# both integrate to 1. The x-independence check (R3) only proves the offset is CONSTANT,
# not that it equals -A_T*g_*-B_T*S_*. The kernel's defining centering condition is
# orthogonality to the canonical source: integral_0^1 Sigma_*(x) W_*(x) dx == 0. This
# DOES test the constants: dropping them leaves integral Sigma_*(A_T c + B_T Kq) =
# A_T*g_* + B_T*S_* != 0, so a missing-centering bug now fails here.
center_resid = sp.N(sp.integrate(Sigma_star_x * Wstar_x, (x, 0, 1)), 30)
assert abs(center_resid) < sp.Float("1e-22", 30), \
    f"kernel not centered against Sigma_*: integral Sigma_* W_* = {center_resid}"
print("PASS: rigidity kernel W_* is source-centered (integral Sigma_* W_* = 0, residual < 1e-22)")

# --- Audit assertion: g_*, S_* equal their source-moment integrals (independent quadrature) ---
# Appendix eq:app-part04-gbar-Sbar/c-Kq define g_* = integral_0^1 Sigma_*(x) c(x) dx and
# S_* = integral_0^1 Sigma_*(x) Kq(x) dx, with Sigma_*(x) = Pi_* e^(-Pi_* x)/(1-e^(-Pi_*)).
# The closed forms gPi, Sformula are the ANALYTIC evaluations of those integrals. Compute
# the integrals directly by quadrature and compare to gPi(Pi_*), Sformula(Pi_*): a
# transcription error in gPi/Sformula would be caught (the old resub check could not).
Sigma_star_x = Pi_star * sp.exp(-Pi_star * x) / (1 - sp.exp(-Pi_star))
g_star_moment = sp.N(sp.integrate(Sigma_star_x * c, (x, 0, 1)), 40)
S_star_moment = sp.N(sp.integrate(Sigma_star_x * Kq, (x, 0, 1)), 40)
assert abs(g_star_moment - g_star) < sp.Float("1e-25", 40), \
    f"g_* moment integral != gPi(Pi_*): {g_star_moment} vs {g_star}"
assert abs(S_star_moment - S_star) < sp.Float("1e-25", 40), \
    f"S_* moment integral != Sformula(Pi_*): {S_star_moment} vs {S_star}"
print("PASS: g_*, S_* equal their source-moment integrals (residual < 1e-25)")

print("\nStage 147 complete.")
