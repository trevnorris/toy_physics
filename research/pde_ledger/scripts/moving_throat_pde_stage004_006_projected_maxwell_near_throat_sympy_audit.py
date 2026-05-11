
#!/usr/bin/env python3
"""Near-throat extension of the projected-Maxwell derivation.

This script extends the projection-first Maxwell analysis in the direction most
relevant to the moving-throat PDE:

1. keep the exact projected inhomogeneous Maxwell law with mixed-sector terms,
2. rewrite it on a finite/half-line throat domain,
3. derive two local observation limits:
   - a symmetric narrow kernel around an interior slice,
   - a one-sided mouth kernel anchored at the throat opening,
4. derive the corresponding effective-coupling / effective-gauge expansions in
   the zero-mode limit,
5. evaluate a concrete Gaussian localization example.

The goal is to make the near-throat structure auditable, especially the fact
that the projected theory keeps the mixed flux derivative ∂_w(Z F^{wν}) alive,
while the far-field zero-mode reduction suppresses it.
"""
from __future__ import annotations

import sympy as sp

def line(title: str) -> None:
    print("\n" + "=" * 96)
    print(title)
    print("=" * 96)

def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


# symbols
ell, sigma, lam, xi, mu0 = sp.symbols("ell sigma lambda xi mu0", positive=True, nonzero=True)
m2, m4 = sp.symbols("m2 m4", real=True)
w, u = sp.symbols("w u", real=True)

# generic Taylor coefficients at the mouth/interior point
z0, z1, z2, z3, z4 = sp.symbols("z0 z1 z2 z3 z4")
h0, h1, h2 = sp.symbols("h0 h1 h2")
s0, s1, s2 = sp.symbols("s0 s1 s2")
q0, q1, q2, q3, q4, q5 = sp.symbols("q0 q1 q2 q3 q4 q5")

# -----------------------------------------------------------------------------
line("1) Exact projected inhomogeneous Maxwell law on a throat domain")

print("Take the generalized weighted bulk equation")
print("  ∂_M(Z F^{MN}) + (1/ξ) ∂^N(H B) = μ0 J^N")
print("and project a brane component ν over a throat interval D ⊂ w using a normalized kernel W(w):")
print("  <Q>_W := ∫_D W(w) Q(w) dw,     ∫_D W(w) dw = 1.")
print()
print("Integrating the w-derivative by parts gives the exact projected law")
print("  ∂_μ <Z F^{μν}>")
print("  + [ W Z F^{wν} ]_{∂D}")
print("  - < W' Z F^{wν} >")
print("  + (1/ξ) ∂^ν <H B>")
print("  = μ0 <J^ν>.")
print()
print("Equivalently, by recombining the boundary and W'-term,")
print("  ∂_μ <Z F^{μν}> + < ∂_w(Z F^{wν}) > + (1/ξ) ∂^ν <H B> = μ0 <J^ν>.")
print()
print("This form is exact on a finite throat interval or a half-line throat chart.")
print("The mixed-sector derivative term <∂_w(Z F^{wν})> is the part that disappears in the far-field zero-mode limit.")

# -----------------------------------------------------------------------------
line("2) Symmetric narrow-kernel expansion around an interior slice")

# formal Taylor series average with an even kernel
Q_series = q0 + q1*(sigma*u) + q2*(sigma*u)**2/sp.Integer(2) + q3*(sigma*u)**3/sp.Integer(6) + q4*(sigma*u)**4/sp.Integer(24) + q5*(sigma*u)**5/sp.Integer(120)
# even normalized kernel moments: <1>=1, <u>=0, <u^2>=m2, <u^3>=0, <u^4>=m4, <u^5>=0
subs_even = {sp.integrate(1, (u, -sp.oo, sp.oo)): 1}  # not used directly

# perform moment substitution by replacing powers of u
avg_Q_even = q0 + m2*sigma**2*q2/sp.Integer(2) + m4*sigma**4*q4/sp.Integer(24)
avg_dQ_even = q1 + m2*sigma**2*q3/sp.Integer(2) + m4*sigma**4*q5/sp.Integer(24)

W_even_unit = sp.exp(-u**2 / 2) / sp.sqrt(2 * sp.pi)
avg_Q_even_gauss = sp.integrate(W_even_unit * Q_series.subs(sigma, 1), (u, -sp.oo, sp.oo))
avg_dQ_even_gauss = sp.integrate(
    W_even_unit * sp.diff(Q_series.subs(sigma, 1), u),
    (u, -sp.oo, sp.oo),
)
assert_zero("Gaussian even-kernel Q moments", avg_Q_even_gauss - avg_Q_even.subs({sigma: 1, m2: 1, m4: 3}))
assert_zero("Gaussian even-kernel derivative moments", avg_dQ_even_gauss - avg_dQ_even.subs({sigma: 1, m2: 1, m4: 3}))

print("For an even normalized kernel W_σ(w) centered at an interior slice w0, odd moments vanish.")
print("Writing Q(w0+σu) in a Taylor series, the exact moment expansion begins")
print("  <Q>_σ =", sp.sstr(avg_Q_even))
print("  <∂_w Q>_σ =", sp.sstr(avg_dQ_even))
print()
print("So the projected Maxwell equation near an interior slice is")
print("  ∂_μ[ M0^{μν} + (m2 σ^2/2) M2^{μν} + O(σ^4) ]")
print("  + [ Q1^ν + (m2 σ^2/2) Q3^ν + O(σ^4) ]")
print("  + (1/ξ) ∂^ν[ G0 + (m2 σ^2/2) G2 + O(σ^4) ]")
print("  = μ0 [ J0^ν + (m2 σ^2/2) J2^ν + O(σ^4) ],")
print("where M^{μν}=Z F^{μν},  Q^ν=Z F^{wν},  G=H B.")
print()
print("Important consequence: for a symmetric observer kernel, the first projection-width corrections are even, O(σ^2).")

# -----------------------------------------------------------------------------
line("3) One-sided mouth kernel: exact boundary cancellation and local throat expansion")

# define truncated Taylor series around w=0
Qw = q0 + q1*w + q2*w**2/sp.Integer(2) + q3*w**3/sp.Integer(6) + q4*w**4/sp.Integer(24)
Wexp = sp.exp(-w/ell)/ell

avg_Q_half = sp.simplify(sp.integrate(Wexp*Qw, (w, 0, sp.oo)))
avg_dQ_half = sp.simplify(sp.integrate(Wexp*sp.diff(Qw, w), (w, 0, sp.oo)))
boundary_half = sp.simplify(sp.limit(Wexp*Qw, w, sp.oo) - (Wexp*Qw).subs(w, 0))
minus_Wp_half = sp.simplify(-sp.integrate(sp.diff(Wexp, w)*Qw, (w, 0, sp.oo)))
recombined_half = sp.simplify(boundary_half + minus_Wp_half)
assert_zero("half-line boundary recombination", recombined_half - avg_dQ_half)
assert_zero("half-line derivative expansion", avg_dQ_half - (q1 + ell*q2 + ell**2*q3 + ell**3*q4))

print("Choose the normalized mouth-anchored kernel on the half-line w ≥ 0")
print("  W_ell(w) = exp(-w/ell)/ell.")
print()
print("For a smooth Q(w) at the mouth, SymPy gives")
print("  <Q>_ell =", sp.sstr(sp.expand(avg_Q_half)))
print("  <∂_w Q>_ell =", sp.sstr(sp.expand(avg_dQ_half)))
print()
print("If we keep the boundary split, the two pieces are")
print("  [W Q]_{0}^{∞} =", sp.sstr(sp.expand(boundary_half)))
print("  -<W' Q>       =", sp.sstr(sp.expand(minus_Wp_half)))
print("and they recombine exactly to")
print("  [W Q]_{0}^{∞} - <W' Q> =", sp.sstr(sp.expand(recombined_half)))
print()
print("So the apparent 1/ell boundary singularity cancels.")
print("The mouth-local projected derivative is finite and expands as")
print("  <∂_w Q>_ell = q1 + ell q2 + ell^2 q3 + ell^3 q4 + ...")
print()
print("Therefore the near-mouth projected Maxwell equation is")
print("  ∂_μ[ M0^{μν} + ell M1^{μν} + ell^2 M2^{μν} + ... ]")
print("  + [ Q1^ν + ell Q2^ν + ell^2 Q3^ν + ... ]")
print("  + (1/ξ) ∂^ν[ G0 + ell G1 + ell^2 G2 + ... ]")
print("  = μ0 [ J0^ν + ell J1^ν + ell^2 J2^ν + ... ].")
print()
print("Important consequence: at the mouth, asymmetry turns the first width correction into O(ell), not O(σ^2).")

# -----------------------------------------------------------------------------
line("4) Zero-mode effective parameters in the near-throat limits")

# symmetric kernel zero-mode expansions
mu_eff_sym = sp.simplify(mu0 * (s0 + m2*sigma**2*s2/2) / (z0 + m2*sigma**2*z2/2))
xi_eff_sym = sp.simplify(xi * (z0 + m2*sigma**2*z2/2) / (h0 + m2*sigma**2*h2/2))
mu_eff_sym_series = sp.series(mu_eff_sym, sigma, 0, 3).removeO()  # up to sigma^2
xi_eff_sym_series = sp.series(xi_eff_sym, sigma, 0, 3).removeO()

# one-sided mouth expansions
mu_eff_half = sp.simplify(mu0 * (s0 + ell*s1 + ell**2*s2) / (z0 + ell*z1 + ell**2*z2))
xi_eff_half = sp.simplify(xi * (z0 + ell*z1 + ell**2*z2) / (h0 + ell*h1 + ell**2*h2))
mu_eff_half_series = sp.series(mu_eff_half, ell, 0, 2).removeO()  # up to ell
xi_eff_half_series = sp.series(xi_eff_half, ell, 0, 2).removeO()

print("Now impose the zero-mode ansatz")
print("  A_μ(x,w)=a_μ(x),   A_w=0,   F^{wν}=0,   J^ν(x,w)=j^ν(x) S(w).")
print("Then the projected equation becomes")
print("  <Z> ∂_μ f^{μν} + (<H>/ξ) ∂^ν(∂·a) = μ0 <S> j^ν,")
print("so the local effective parameters are")
print("  μ_eff = μ0 <S>/<Z>,    ξ_eff = ξ <Z>/<H>.")
print()
print("Symmetric interior kernel (first nontrivial correction O(σ^2)):")
print("  μ_eff^(sym) =", sp.sstr(mu_eff_sym_series))
print("  ξ_eff^(sym) =", sp.sstr(xi_eff_sym_series))
print()
print("One-sided mouth kernel (first correction O(ell)):")
print("  μ_eff^(mouth) =", sp.sstr(mu_eff_half_series))
print("  ξ_eff^(mouth) =", sp.sstr(xi_eff_half_series))
print()
print("So the first mismatch measures local shape-slippage:")
print("  symmetric: (S2/S0 - Z2/Z0) and (Z2/Z0 - H2/H0),")
print("  mouth:     (S1/S0 - Z1/Z0) and (Z1/Z0 - H1/H0).")
print()
print("Immediate corollaries:")
print("  - If H = Z, then ξ_eff = ξ exactly, not just approximately.")
print("  - If S is locally proportional to Z at the mouth/interior slice, the first μ_eff correction vanishes.")

# exact cancellations for H=Z and S=C Z
C = sp.symbols("C")
mu_eff_proportional = sp.simplify(mu0 * (C*(z0 + ell*z1 + ell**2*z2)) / (z0 + ell*z1 + ell**2*z2))
xi_eff_exact_HZ = sp.simplify(xi * (z0 + ell*z1 + ell**2*z2) / (z0 + ell*z1 + ell**2*z2))
W_conc = sp.exp(-w / ell) / ell
Z_conc = z0 + z1*w + z2*w**2 / 2
IZ_conc = sp.integrate(W_conc * Z_conc, (w, 0, sp.oo))
H_sym = sp.Function("H")(w)
S_sym = sp.Function("S")(w)
IH_conc = sp.Integral(W_conc * H_sym, (w, 0, sp.oo))
IS_conc = sp.Integral(W_conc * S_sym, (w, 0, sp.oo))
assert_zero("H=Z before effective-gauge cancellation", IH_conc.subs(H_sym, Z_conc) - IZ_conc)
assert_zero("S=CZ before effective-coupling cancellation", IS_conc.subs(S_sym, C * Z_conc) - C * IZ_conc)
assert_zero("H=Z effective gauge", sp.simplify((xi * IZ_conc / IH_conc).subs(H_sym, Z_conc)) - xi)
assert_zero("S=CZ effective coupling", sp.simplify((mu0 * IS_conc / IZ_conc).subs(S_sym, C * Z_conc)) - C * mu0)
print()
print("Exact profile-alignment checks:")
print("  if S = C Z globally, then μ_eff =", sp.sstr(mu_eff_proportional))
print("  if H = Z globally, then ξ_eff =", sp.sstr(xi_eff_exact_HZ))

# -----------------------------------------------------------------------------
line("5) Concrete Gaussian localization examples")

# Symmetric Gaussian kernel and Gaussian localizer
Ws = sp.exp(-w**2/(2*sigma**2)) / (sp.sqrt(2*sp.pi) * sigma)
Zg = sp.exp(-w**2/lam**2)
IWZ_sym_gauss = sp.simplify(sp.integrate(Ws * Zg, (w, -sp.oo, sp.oo)))
IWZ_sym_gauss_series = sp.series(IWZ_sym_gauss, sigma, 0, 5).removeO()

# One-sided exponential mouth kernel with Gaussian localizer
Wm = sp.exp(-w/ell)/ell
IWZ_mouth_gauss = sp.simplify(sp.integrate(Wm * Zg, (w, 0, sp.oo)))
r = sp.symbols("r", positive=True)
IWZ_mouth_gauss_r = sp.sqrt(sp.pi) * lam * r * sp.erfc(lam * r / 2) * sp.exp(lam**2 * r**2 / 4) / 2
IWZ_mouth_gauss_series = sp.simplify(
    sp.series(IWZ_mouth_gauss_r, r, sp.oo, 8).removeO().subs(r, 1 / ell)
)
IWZ_mouth_taylor_integral = sp.integrate(
    Wm * sp.series(Zg, w, 0, 8).removeO(),
    (w, 0, sp.oo),
)
assert_zero("mouth Gaussian asymptotic from erfc closed form", IWZ_mouth_gauss_series - (1 - 2*ell**2/lam**2 + 12*ell**4/lam**4 - 120*ell**6/lam**6))
assert_zero("mouth Gaussian asymptotic from Taylor integration", IWZ_mouth_gauss_series - IWZ_mouth_taylor_integral)

print("Take Z(w)=exp(-w^2/λ^2).")
print()
print("Symmetric Gaussian observer kernel:")
print("  W_σ(w) = exp(-w^2/(2σ^2)) / (sqrt(2π) σ)")
print("  <Z> =", sp.sstr(IWZ_sym_gauss))
print("  series =", sp.sstr(IWZ_sym_gauss_series))
print()
print("So the symmetric kernel samples the throat localizer as")
print("  <Z> = 1 - σ^2/λ^2 + 3σ^4/(2λ^4) + O(σ^6).")
print()
print("One-sided exponential mouth kernel:")
print("  W_ell(w) = exp(-w/ell)/ell,   w ≥ 0")
print("  <Z> =", sp.sstr(IWZ_mouth_gauss))
print("  series =", sp.sstr(IWZ_mouth_gauss_series))
print()
print("So the mouth-local kernel samples the localizer as")
print("  <Z> = 1 - 2 ell^2/λ^2 + 12 ell^4/λ^4 - 120 ell^6/λ^6 + O(ell^8).")
print("The odd correction is absent here because the Gaussian localizer has Z'(0)=0,")
print("but a generic asymmetric source or gauge profile would generate O(ell) corrections.")

# -----------------------------------------------------------------------------
line("6) PDE reading of the result")

print("Near the throat, projection-first Maxwell does NOT collapse immediately to the far-field brane law.")
print("The leading extra structure is the mixed-sector derivative term")
print("  <∂_w(Z F^{wν})>,")
print("which survives exactly on the throat domain and is lost only when the far-field zero-mode suppressions")
print("  F_{μw} ≈ 0,   J^w ≈ 0,   ∂_w A_μ ≈ 0")
print("are imposed.")
print()
print("This is why the near-throat limit is a better match to the moving-throat PDE language:")
print("the Maxwell/mixed block contributes local transverse-gradient data before any reduction to ordinary 3+1 Maxwell.")
print()
print("The main structural distinction is:")
print("  - symmetric interior observation  -> first width correction O(σ^2),")
print("  - mouth-anchored observation     -> first width correction O(ell).")
print()
print("So the mouth itself is the place where projection-first and reduction-first differ most sharply.")
print("That makes it the natural regime to compare against the moving-throat Maxwell/mixed bundle.")
print("STATUS: PASS")
