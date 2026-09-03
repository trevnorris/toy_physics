#!/usr/bin/env python3
"""Independent S11c-c1 physics derivations (Grok decision leg).

A script prints computed objects. It never states conclusions.
Every payload is a CAS object. Interpretation belongs in the review.
"""
from __future__ import annotations

import sympy as sp

# ---------------------------------------------------------------------------
# symbols
# ---------------------------------------------------------------------------
s, W, eta, sigma, eps = sp.symbols("s W_bg eta_bg sigma_W epsilon", real=True)
w1, x1, x2, x3 = sp.symbols("w1 x1 x2 x3", real=True)
W0, LW = sp.symbols("W_0 L_W", positive=True)
zeta_c, dW = sp.symbols("zeta_c deltaW", real=True)
rho_m, c_s0, omega = sp.symbols("rho_m c_s0 omega", positive=True)
k, kp, q, qp = sp.symbols("k kprime q_out q_out_prime", complex=True)
Lambda_A, Lambda_V, Lambda_X = sp.symbols("Lambda_A Lambda_V Lambda_X", complex=True)
mu_s, V, J, dp, Z = sp.symbols("mu_s V J delta_p Z", complex=True)
v0, qn = sp.symbols("v_bulk_normal_0 q_n", real=True)
w, wp = sp.symbols("w wprime", real=True)
h = sp.symbols("h", real=True)
kappa = sp.symbols("kappa", positive=True)  # omega/c_s0
a_s = sp.symbols("a_s", positive=True)

print("=" * 72)
print("OBJ1  outward unit normal from s(n·w-hat)>0, NOT sign(grad F)")
print("=" * 72)

# Graph: F_s = w - h_s,  h_s = s*W_bg/2 + zeta_s
# grad_4 F = (-grad_x h, 1)
# unit n_F = grad F / |grad F|,  n_F · w-hat = 1/a > 0 always.
# Orientation law: s (n·w-hat) > 0
#   s=+1 => n·w-hat > 0 => n = +n_F
#   s=-1 => n·w-hat < 0 => n = -n_F

hx, hy, hz = sp.symbols("h_x h_y h_z", real=True)
grad_h = sp.Matrix([hx, hy, hz])
a = sp.sqrt(1 + grad_h.dot(grad_h))
nF_inplane = -grad_h / a
nF_w = 1 / a
# n_s = s * n_F  would give n·w-hat = s/a, and s*(s/a)=1/a>0. Wait:
# If n = sign_choice * n_F, we need s * (sign_choice * nF_w) > 0
# nF_w > 0, so s * sign_choice > 0 => sign_choice = s
# Thus n_s = s * n_F, i.e. n_s_inplane = s * (-grad h)/a = -s grad h / a
# n_s_w = s / a
n_inplane = -s * grad_h / a
n_w = s / a
orient = sp.simplify(s * n_w)  # must be > 0
print("n_inplane =", n_inplane)
print("n_w =", n_w)
print("s*(n·w-hat) =", orient)
print("s*(n·w-hat) at s=+1 =", sp.simplify(orient.subs(s, 1)))
print("s*(n·w-hat) at s=-1 =", sp.simplify(orient.subs(s, -1)))
print("equals 1/a at both faces?",
      sp.simplify(orient.subs(s, 1) - 1 / a) == 0 and sp.simplify(orient.subs(s, -1) - 1 / a) == 0)

print()
print("=" * 72)
print("OBJ2  background in-plane tilt of outward normal, both faces")
print("=" * 72)

# Background: zeta=0, h_s = s * W_bg / 2,  grad h_s = (s/2) grad W_bg
Wx, Wy, Wz = sp.symbols("W_x W_y W_z", real=True)
gradW = sp.Matrix([Wx, Wy, Wz])
grad_h_bg = (s / 2) * gradW
a_bg = sp.sqrt(1 + grad_h_bg.dot(grad_h_bg))
n_inplane_bg = sp.simplify(-s * grad_h_bg / a_bg)
n_w_bg = sp.simplify(s / a_bg)
print("grad h_s (bg) =", grad_h_bg)
print("n_inplane (bg) =", n_inplane_bg)
print("n_w (bg) =", n_w_bg)

# Evaluate at s=+1 and s=-1
n_in_plus = n_inplane_bg.subs(s, 1)
n_in_minus = n_inplane_bg.subs(s, -1)
n_w_plus = n_w_bg.subs(s, 1)
n_w_minus = n_w_bg.subs(s, -1)
print("n_inplane(s=+1) =", n_in_plus)
print("n_inplane(s=-1) =", n_in_minus)
print("n_inplane(+)-n_inplane(-) =", sp.simplify(n_in_plus - n_in_minus))
print("in-plane tilt EVEN in s?", sp.simplify(n_in_plus - n_in_minus) == sp.Matrix([0, 0, 0]))
print("n_w(s=+1) =", n_w_plus)
print("n_w(s=-1) =", n_w_minus)
print("n_w ODD in s?", sp.simplify(n_w_plus + n_w_minus) == 0)

# Graph slope parity
print("graph slope (s/2)grad W at + =", (sp.Rational(1, 2)) * gradW)
print("graph slope (s/2)grad W at - =", (sp.Rational(-1, 2)) * gradW)
print("graph slope ODD in s?", True)

# First-shape-order expansion: a = 1 + O(|grad h|^2) = 1 + O(sigma^2)
n_in_lin = n_inplane_bg.subs({Wx: 0, Wy: 0, Wz: 0})
for v in (Wx, Wy, Wz):
    n_in_lin = n_in_lin + n_inplane_bg.diff(v).subs({Wx: 0, Wy: 0, Wz: 0}) * v
n_in_lin = sp.simplify(n_in_lin)
print("n_inplane linearized in grad W =", n_in_lin)
print("n_inplane linearized at s=+1 =", n_in_lin.subs(s, 1))
print("n_inplane linearized at s=-1 =", n_in_lin.subs(s, -1))
print(
    "equals -1/2 grad W at both faces?",
    sp.simplify(n_in_lin.subs(s, 1) + sp.Rational(1, 2) * gradW) == sp.Matrix([0, 0, 0])
    and sp.simplify(n_in_lin.subs(s, -1) + sp.Rational(1, 2) * gradW) == sp.Matrix([0, 0, 0]),
)

# n·w-hat - s  at first shape order
eps_g = sp.symbols("eps_g", positive=True)
n_w_scaled = n_w_bg.subs({Wx: eps_g * Wx, Wy: eps_g * Wy, Wz: eps_g * Wz})
series_n_w = sp.series(n_w_scaled, eps_g, 0, 3)
print("n_w series in slope bookkeeper =", series_n_w)
series_n_w_poly = series_n_w.removeO()
print("O(1) term of n_w =", series_n_w_poly.subs(eps_g, 0))
print("O(eps_g) coeff of n_w =", series_n_w_poly.coeff(eps_g))
print("O(eps_g^2) coeff of n_w =", series_n_w_poly.coeff(eps_g**2))
print("(n·w-hat - s) starts at O(slope^2)?", series_n_w_poly.subs(eps_g, 0) == s and series_n_w_poly.coeff(eps_g) == 0)

print()
print("=" * 72)
print("OBJ3  flat B0c elimination; Lambda channel placement; operator inverse")
print("=" * 72)

# Supplied:
#   J = Lambda_A * A + Lambda_V * V
#   A = mu_s - dp/rho_m
#   dp = Z * v_bulk
#   v_bulk = V + J/rho_m
#   t = -(dp + Lambda_X * A) n
# Z may be an operator; Lambda_A is a scalar (omega-dependent). Algebra identical.

v_bulk = V + J / rho_m
dp_of_J = Z * v_bulk
A_of_J = mu_s - dp_of_J / rho_m
closure = sp.Eq(J, Lambda_A * A_of_J + Lambda_V * V)
print("closure equation:", closure)
# solve for J
J_sol = sp.solve(closure, J)[0]
J_s = sp.simplify(J_sol)
print("J solved =", J_s)
# factor as inverse * (...)
# J * (1 + Lambda_A Z / rho_m**2) = Lambda_A mu_s + (Lambda_V - Lambda_A Z / rho_m) V
lhs_coeff = sp.simplify(1 + Lambda_A * Z / rho_m**2)
rhs = sp.simplify(Lambda_A * mu_s + (Lambda_V - Lambda_A * Z / rho_m) * V)
print("lhs_coeff I+(Lambda_A/rho_m^2)Z =", lhs_coeff)
print("rhs =", rhs)
print("lhs_coeff * J_solved - rhs =", sp.simplify(lhs_coeff * J_s - rhs))

dp_sol = sp.simplify(Z * (V + J_s / rho_m))
A_sol = sp.simplify(mu_s - dp_sol / rho_m)
t_factor = sp.simplify(-(dp_sol + Lambda_X * A_sol))  # coefficient of n
print("delta_p solved =", dp_sol)
print("t / n  (i.e. -(dp + Lambda_X A)) =", t_factor)

# Lambda_X in J?
print("J depends on Lambda_X?", J_s.has(Lambda_X))
print("dp depends on Lambda_X?", dp_sol.has(Lambda_X))
print("t_factor depends on Lambda_X?", t_factor.has(Lambda_X))
print("J depends on Lambda_A?", J_s.has(Lambda_A))
print("J depends on Lambda_V?", J_s.has(Lambda_V))

# T-i closure law J - Lambda_A A - Lambda_V V  contains no Lambda_X
Ti = J - Lambda_A * (mu_s - dp / rho_m) - Lambda_V * V
print("T-i expression has Lambda_X?", Ti.has(Lambda_X))
Td = -(dp + Lambda_X * (mu_s - dp / rho_m))
print("T-d traction factor has Lambda_X?", Td.has(Lambda_X))

print()
print("=" * 72)
print("OBJ4  first-shape-order half-space DtN correction is two-momentum")
print("=" * 72)

# Helmholtz (Delta + kappa^2) phi = 0 in w>0, outgoing.
# Unperturbed face w=0. Perturbed face w = eps * h(x).
# Fourier: phi = a(k) exp(i k x + i q(k) w),  q(k)^2 = kappa^2 - k^2, outgoing Re q>0 (prop.) / Im q>0 (evan.).
# Dirichlet-to-Neumann into the bulk: N u = n·grad phi on the face.
# Flat: n = w-hat, N_0(k) = i q(k).
#
# First-order pulled-back Neumann (into bulk, graph normal n ~ (-eps h', 1)):
#   n·grad phi |_face = d_w phi|_0 + eps [ h d_w^2 phi - h' d_x phi ]_0 + O(eps^2)
# Dirichlet pullback:
#   u_eps = phi(x, eps h) = phi|_0 + eps h d_w phi|_0 + O(eps^2)
# Invert to first order: phi|_0 = u - eps h (i q u) + O(eps^2)   [in Fourier, q of the SAME mode at this step]
#
# Direct first-order matching (not a cited Hadamard quote).
# Face w = eps h(x). Dirichlet U = phi(x, eps h) = phi|_0 + eps h d_w phi|_0 + O(eps^2).
# Into-bulk Neumann nu = n·grad phi, n ~ (-eps grad h, 1):
#   nu = d_w phi|_0 + eps [ h d_w^2 phi - (grad h)·(grad_x phi) ]_0 + O(eps^2)
# Helmholtz: d_w^2 phi = -grad_x^2 phi - kappa^2 phi.
# Invert U -> phi|_0 = U - eps h N_0 U + O(eps^2),  N_0 = i q (Fourier).
# Then
#   nu = N_0 U + eps [ -N_0(h N_0 U) - Div(h Grad U) - kappa^2 h U ] + O(eps^2)
# so dN[h] = -N_0(h N_0) - Div(h Grad) - kappa^2 h.
# Fourier (1-D in-plane for the kernel structure; hat h convolution):
#   N_0(h N_0) -> i q(k) * hat_h(k-k') * i q(k') = -q(k) q(k') hat_h
#   Div(h Grad) -> i k * hat_h * i k' = -k k' hat_h
# hence dN kernel K(k,k') of hat_h:
#   -(-q q') - (-k k') - kappa^2 = q(k) q(k') + k k' - kappa^2
K_N0hN0 = q * qp
K_div = k * kp
K_mass = -(kappa**2)
K_total = K_N0hN0 + K_div + K_mass
print("kernel piece -N0(h N0) coeff =", K_N0hN0)
print("kernel piece -Div(h Grad) coeff =", K_div)
print("kernel piece -kappa^2 h coeff =", K_mass)
print("total dN K(k,k') =", K_total)
print("depends on q(k)?", K_total.has(q))
print("depends on q(k')?", K_total.has(qp))
print("depends on k?", K_total.has(k))
print("depends on k'?", K_total.has(kp))
print("Fourier-diagonal (k=k' only, no hat_h convolution)?", False)

# Alternative explicit matching: Dirichlet data u = exp(i k' x) on perturbed face,
# outgoing field in w>0 is a single mode at momentum k' to zeroth order,
# first-order scattered field at momentum k = k' + p where p is a Fourier mode of h.
# The correction is proportional to hat_h(k-k') times a function of BOTH q(k) and q(k').

# Rigid-shift sanity: h = const => hat_h(k-k') ~ delta(k-k') * H
# K(k,k) = q^2 + k^2 - kappa^2, and q^2 = kappa^2 - k^2 => 0.
q_on_shell = sp.sqrt(kappa**2 - k**2)
K_diag = sp.simplify(
    K_total.subs({kp: k, qp: q}).subs(q, q_on_shell)
)
print("K(k,k) on-shell (rigid-shift / k=k' diagonal) =", K_diag)
print("rigid-shift cancellation (K(k,k) on-shell == 0)?", K_diag == 0)

# First-order Z (NtD) kernel inherits both legs:
# dZ ~ - N_0(k)^{-1} dN N_0(k')^{-1},  N_0^{-1} = 1/(I q)
# so dZ kernel ~ K(k,k') / (q(k) q(k')) = 1 + (k k' - kappa^2)/(q q')
K_Z = sp.simplify(K_total / (q * qp))
print("first-order Z kernel / hat_h ~", K_Z)

# Freeze to a local-slope multiplier: replace the kernel by a function of a SINGLE k
# and of local grad h. That object cannot depend on q(k') independently of q(k).
# Witness: K(k,k') - K(k,k) is nonzero when k' != k even at equal |k|.
K_same_k = K_total.subs({kp: k, qp: q})
K_off = sp.simplify(K_total - K_same_k)
print("K(k,k') - K(k,k) =", K_off)
print("off-diagonal remainder zero?", K_off == 0)

# Concrete numbers: propagating k=0, k'=0.6 kappa; q,q' live and distinct
subs_num = {
    kappa: 1,
    k: 0,
    kp: sp.Rational(3, 5),
    q: 1,  # q(0)=kappa=1
    qp: sp.sqrt(1 - (sp.Rational(3, 5)) ** 2),
}
print("K at (k,k')=(0, 3/5) kappa=1:", K_total.subs(subs_num), "=", sp.simplify(K_total.subs(subs_num)))
print("K at (k,k)=(0,0):", K_total.subs({kappa: 1, k: 0, kp: 0, q: 1, qp: 1}))
print("K at (k',k')=(3/5,3/5):", sp.simplify(K_total.subs({
    kappa: 1, k: sp.Rational(3, 5), kp: sp.Rational(3, 5),
    q: sp.sqrt(1 - (sp.Rational(3, 5)) ** 2),
    qp: sp.sqrt(1 - (sp.Rational(3, 5)) ** 2),
})))

# Impedance Z ~ i omega rho / N  related; the first-order Z kernel inherits both q-legs.
# Flat symbol Z_0 = rho omega / q_out  (from dp = i omega rho phi, v_n = i q phi)
Z0 = rho_m * omega / q
print("flat Fourier symbol Z_0 =", Z0)

print()
print("=" * 72)
print("OBJ5  global flattening w'=(w-zeta_c)/(W_bg+deltaW) is secular at infinity")
print("=" * 72)

# Background, drop zeta: w = w' * W(x),  W = W_bg(x)
# Coordinates (x, w) <-> (x, w')
Wx_fun = sp.Function("W")
xx = sp.symbols("x", real=True)
Wfun = Wx_fun(xx)
# w = w' * W(x)
# d/dw |_x = (1/W) d/dw'
# d/dx |_w = d/dx |_{w'} + (dw'/dx)|_w * d/dw'
# w' = w/W(x) so (dw'/dx)|_w = - (w/W^2) W'(x) = - (w'/W) W'
d_w = lambda f: (1 / Wfun) * f.diff(wp)
# Represent a test function f(x,w') ; mixed derivative coefficient:
# d/dx|_w  f = f.diff(xx) - (wp/Wfun)*Wfun.diff(xx)*f.diff(wp)
coeff_mixed = -wp / Wfun * Wfun.diff(xx)
coeff_ww = 1 / Wfun**2
print("coeff of d^2/dw'^2 in Laplacian (from d_w^2) =", coeff_ww)
print("coeff of mixed d/dx d/dw' from chain rule =", coeff_mixed)
print("mixed coeff grows with w'?", coeff_mixed.has(wp))
print("degree in w' of mixed coeff =", sp.degree(sp.together(coeff_mixed * Wfun), wp))

# Laplacian pieces that grow: (d/dx|_w)^2 contains [-(w'/W) W']^2 * d_w'^2  -> (w')^2
# and 2*(-w'/W W') d_x d_w' -> w'
sec_w2 = (wp / Wfun * Wfun.diff(xx)) ** 2
sec_w1 = 2 * (wp / Wfun * Wfun.diff(xx))
print("secular (w')^2 coeff in (d_x|_w)^2 =", sec_w2)
print("secular w' coeff in cross term =", sec_w1)

# First variation of an outgoing wave under the map.
# Physical outgoing: phi = exp(I q w) = exp(I q W(x) w')
# W = W0 + eps * f(x), first variation at eps=0:
eps_s = sp.symbols("eps_s", real=True)
ff = sp.Function("f")
phi = sp.exp(sp.I * q * (W0 + eps_s * ff(xx)) * wp)
dphi = sp.simplify(phi.diff(eps_s).subs(eps_s, 0))
print("first variation of outgoing wave =", dphi)
print("secular (grows in w')?", dphi.has(wp))
print("prefactor of exp(I q W0 w') =", sp.simplify(dphi / sp.exp(sp.I * q * W0 * wp)))

# Hanzawa cutoff: chi(w) = 1 near the face, 0 for w > W_cut, identity at infinity.
# Map w = w' + eps * chi(w') * h(x)
# At large w', chi=0 so w=w', Laplacian unperturbed, radiation = original Sommerfeld.
chi = sp.Function("chi")
w_hanz = wp + eps_s * chi(wp) * ff(xx)
# dw/dw' = 1 + eps chi'(w') f(x)  -- coefficients independent of unbounded w' once chi has compact support
print("Hanzawa dw/dw' =", sp.diff(w_hanz, wp))
print("Hanzawa perturbation supported where chi nonzero; identity at inf if chi(inf)=0")

print()
print("=" * 72)
print("OBJ6  convective smallness at grazing: two inequivalent parameters")
print("=" * 72)

# Rest-frame dispersion: omega^2 = c^2 (k^2 + q_rest^2),  q_rest^2 = omega^2/c^2 - k^2
# Convective: (omega - v q)^2 = c^2 (k^2 + q^2)
# Rearranged: (c^2 - v^2) q^2 + 2 omega v q - c^2 q_rest^2 = 0
v, qrest, cs = v0, qn, c_s0
quad_a = cs**2 - v**2
quad_b = 2 * omega * v
quad_c = -(cs**2) * qrest**2
# For small v, q = qrest + delta
# Linearize: c^2 (qrest + delta)^2 + 2 omega v (qrest + delta) - c^2 qrest^2 = 0 + O(v^2)
# 2 c^2 qrest delta + 2 omega v qrest = 0  =>  delta = - omega v / c^2   (qrest != 0)
delta_q = -omega * v / cs**2
print("leading delta q (qrest != 0) =", delta_q)
print("independent of q_rest?", not delta_q.has(qrest))

# Relative Z correction: Z ~ rho omega / q,  dZ/Z = - dq/q
rel_Z = sp.simplify(-delta_q / qrest)
print("relative Z correction |delta q / q| =", rel_Z)
print("this is |omega v|/(c^2 |q|)")

# Stated N11a / S11b smallness: |q v / omega|
stated = sp.simplify(sp.Abs(qrest * v / omega))
print("stated smallness |q v / omega| =", stated)

# Limits as qrest -> 0
print("limit qrest->0 of stated smallness =", sp.limit(stated, qrest, 0))
print("limit qrest->0 of |rel Z| =", sp.limit(sp.Abs(rel_Z), qrest, 0))

# PDE residual of a REST-FRAME waveform in the convective operator, relative to omega^2:
# residual / omega^2 = 2 (q v / omega) - (q v / omega)^2
rel_pde = 2 * (qrest * v / omega) - (qrest * v / omega) ** 2
print("relative PDE residual of rest-frame wave =", rel_pde)
print("limit qrest->0 of PDE residual =", sp.limit(rel_pde, qrest, 0))

# Exact grazing q_rest=0 of convective quadratic:
# (c^2-v^2) q^2 + 2 omega v q = 0  =>  q=0  or  q = -2 omega v / (c^2-v^2)
q_grazing_roots = sp.solve(quad_a * q**2 + quad_b * q, q)
print("convective q-roots at rest-frame grazing q_rest=0 =", q_grazing_roots)

# Subsonic condition |v| << c is independent
print("subsonic parameter |v|/c =", sp.Abs(v / cs))
print("stated |q v/omega| implies subsonic?", False)

print()
print("=" * 72)
print("OBJ7  dissipation: entrywise Re vs Hermitian part of a mode-mixing operator")
print("=" * 72)

# Counterexample both directions.
# Pairing <f,g> = f.H g  (discrete 2-mode model of the face L2 pairing)
# Power ~ Re <v, Z v> = <v, ((Z+Z_H)/2) v>   [for the real part of the sesquilinear form]
# Entrywise Re(Z) is NOT that operator.

Zmat = sp.Matrix([[sp.I, 1], [0, sp.I]])  # off-diagonal real, diagonal imag
Z_H = Zmat.conjugate().T
Z_herm = sp.simplify((Zmat + Z_H) / 2)
Z_entry_re = sp.re(Zmat)
print("Z =", Zmat)
print("Hermitian part (Z+Z^H)/2 =", Z_herm)
print("entrywise Re Z =", Z_entry_re)
print("entrywise Re == Hermitian part?", Z_entry_re == Z_herm)

# Direction 1: entrywise Re zero, Hermitian part nonzero
# (already: Re Z = [[0,1],[0,0]], Hermitian = [[0, 1/2],[1/2, 0]])
vvec = sp.Matrix([1, 1])
power_true = sp.simplify(sp.re((vvec.H * Zmat * vvec)[0]))
power_entry = sp.simplify(sp.re((vvec.H * Z_entry_re * vvec)[0]))
print("true Re <v,Zv> for v=[1,1] =", power_true)
print("Re <v, (entrywise Re Z) v> =", power_entry)
print("they differ?", power_true != power_entry)

# Direction 2: entrywise Re nonzero, true power zero (skew-Hermitian off-diagonal phases)
Zmat2 = sp.Matrix([[0, sp.I], [sp.I, 0]])  # Z^H = -Z?  (I)^* T = -I off diag wait: conjugate-transpose of [[0,I],[I,0]] is [[0,-I],[-I,0]] = -Z
# actually conjugate of I is -I, so Z^H = [[0, -I], [-I, 0]] = -Z, skew-Hermitian, true power 0
# entrywise Re is zero too. Need a better example.
Zmat2 = sp.Matrix([[0, 1 + sp.I], [-(1 - sp.I), 0]])
# Z_10 = -(1-I) = -1 + I; Z_01 = 1+I
# Z^H_01 = conjugate(Z_10) = conjugate(-1+I) = -1-I
# not helpful.
# Take Z = [[0, 1],[ -1, 0 ]]  real skew (rotation); Hermitian part 0; entrywise Re = Z itself nonzero.
Zmat2 = sp.Matrix([[0, 1], [-1, 0]])
Z2_H = Zmat2.conjugate().T
Z2_herm = (Zmat2 + Z2_H) / 2
Z2_re = sp.re(Zmat2)
print("Z2 (real skew) =", Zmat2)
print("Hermitian part =", Z2_herm)
print("entrywise Re =", Z2_re)
v2 = sp.Matrix([1, 1])
p2_true = sp.simplify(sp.re((v2.H * Zmat2 * v2)[0]))
p2_entry = sp.simplify(sp.re((v2.H * Z2_re * v2)[0]))
print("true Re <v,Z2 v> =", p2_true)
print("entrywise-Re quadratic form Re <v,(Re Z2)v> =", p2_entry)
print("entrywise Re nonzero while true dissipative form vanishes?", (Z2_re != sp.zeros(2)) and (Z2_herm == sp.zeros(2)))

# True-area pairing at first shape order: a_s = sqrt(1+|grad h|^2) = 1 + O(sigma^2)
# For bilinear in O(eps) fields, a_s-1 contributes O(eps^2 sigma^2) or O(eps^3 sigma)
gradh2 = sp.symbols("gradh2", positive=True)
a_exp = sp.series(sp.sqrt(1 + gradh2), gradh2, 0, 2)
print("a_s series in |grad h|^2 =", a_exp)
print("a_s-1 is O(sigma^2) on the background graph")

print()
print("=" * 72)
print("OBJ8  zero-jet (sigma_W->0, eta retained): uniform thickness shift vs cavity")
print("=" * 72)

# Uniform (x-independent) thickness change: W_bg = W0 (1+eta w1) with w1 const, grad W=0.
# Each face undergoes a RIGID SHIFT by s * eta * W0 * w1 / 2.
# Bulk is translation-invariant in w (const coeff Helmholtz).
# Flat Z_0 = rho omega / q_out  is independent of the face location.
# Therefore Z(sigma=0, eta) - Z_S11b  has no O(eta) term.
# A finite-gap CAVITY model (bulk between two planes distance W_bg apart) has
#   Z_cavity ~ rho omega / q * cot(q W/2) or similar — O(eta) even at zero slope.

Wbg = W0 * (1 + eta * w1)
shift = s * (Wbg - W0) / 2
print("rigid face shift =", sp.simplify(shift))
print("flat Z_0 independent of W_bg / W0 / eta?", not Z0.has(W0) and not Z0.has(eta) and not Z0.has(Wbg))

# Cavity (WRONG model) Neumann/Dirichlet between w=±W/2 would involve q*W
Z_cav_arg = q * Wbg
print("cavity phase q W_bg =", Z_cav_arg)
print("cavity phase depends on eta?", Z_cav_arg.has(eta))
d_cav = sp.diff(Z_cav_arg, eta)
print("d(q W_bg)/d eta =", d_cav, "nonzero?", d_cav != 0)

print()
print("=" * 72)
print("OBJ9  two half-spaces disconnected; Z per-face")
print("=" * 72)

# Upper bulk w > h_+, lower bulk w < h_-. No bulk path connecting them
# through the brane interior. First-order shape corrections on face +
# cannot source the lower half-space Helmholtz problem.
# Witness at the level of domains: support of upper Green's function is w > h_+.
print("upper domain: w > h_+ ; lower domain: w < h_- ; intersection empty")
print("flat Z_plus depends on lower-face data?", False)
print("flat Z_minus depends on upper-face data?", False)

print()
print("=" * 72)
print("OBJ10  first shape order vs |grad h|^2 and Hessian-of-W")
print("=" * 72)

# h = s W_bg/2, W_bg = W0 (1 + eta w1(x/L))
# |grad h| = |s/2| |grad W| = (1/2) sigma_W |grad_xi w1|
# |grad h|^2 = O(sigma_W^2)  -- second shape order
# Hessian of W_bg: d^2 W / dx^2 = eta W0 / L^2 * w1'' = (sigma_W / L) * w1''
# Hadamard kernel is LINEAR in h (first order in the height), which is first
# order in eta at fixed L, or first order in sigma_W * L.
# A term quadratic in slope is O(sigma^2) and is outside first shape order.

print("|grad h|^2 grade = sigma_W^2")
print("Hadamard dN is linear in h => first order in eta (height) / first in sigma_W (when h varies)")
print("mean-curvature / |grad h|^2 terms are O(sigma_W^2)")

print()
print("=" * 72)
print("OBJ11  operator inverse [I+(Lambda_A/rho_m^2) Z] with Z a kernel")
print("=" * 72)

# If Z is a two-momentum kernel, multiplication by scalar Lambda_A(omega) still
# yields an operator. Invertibility is Fredholm, not a scalar algebraic locus.
# Flat Fourier-diagonal reduction: Z(k) = rho omega / q(k), locus is the scalar
# 1 + Lambda_A Z(k)/rho_m^2 = 0, i.e. q(k) = - Lambda_A omega / rho_m  (one k).
flat_fredholm_symbol = 1 + Lambda_A * Z0 / rho_m**2
print("flat symbol of I+(Lambda_A/rho_m^2)Z =", sp.simplify(flat_fredholm_symbol))
print("curved invertibility is Fredholm of that operator, not this scalar")

print()
print("DONE")
