#!/usr/bin/env python3
"""c1 energy-route teeth: does 1/2 Re(delta_p V*) vs bulk flux see a t_s sign error?

Prints CAS objects only.
"""
import sympy as sp

rho, omega, q, A, w, Wface = sp.symbols("rho_m omega q_out A w W_face", real=True)
# Outgoing upper-half-space mode, harmonic e^{-i omega t} stripped.
# phi = A exp(I q (w - W_face)), q real > 0 (propagating, energy +w).
phi = A * sp.exp(sp.I * q * (w - Wface))
# v = grad phi; v_w = I q phi
v_w = sp.I * q * phi
# delta_p = -rho d_t phi = I omega rho phi   (d_t -> -I omega)
dp = sp.I * omega * rho * phi
# Face at w = W_face: impermeable V = v_w (outward = +w on upper face)
dp_f = sp.simplify(dp.subs(w, Wface))
vw_f = sp.simplify(v_w.subs(w, Wface))
V = vw_f

# Face pairing written in the spec: 1/2 Re(delta_p V*)
# For a complex amplitude, the period-average is 1/2 Re(dp V.conjugate())
face_spec = sp.simplify(sp.Rational(1, 2) * sp.re(dp_f * sp.conjugate(V)))

# Independent bulk flux: acoustic intensity in +w, period average
# I_w = 1/2 Re(dp * v_w.conjugate())  at any w (no dissipation in the bulk).
flux_w = sp.simplify(sp.Rational(1, 2) * sp.re(dp * sp.conjugate(v_w)))
flux_face = flux_w.subs(w, Wface)
flux_infty = flux_w  # w-independent for a single propagating mode
print("dp_face =", dp_f)
print("V_outward =", V)
print("spec face pairing 1/2 Re(dp V*) =", face_spec)
print("bulk flux 1/2 Re(dp v_w*) at face =", flux_face)
print("bulk flux 1/2 Re(dp v_w*) at generic w =", flux_infty)
print("flux independent of w?", sp.simplify(flux_w - flux_face) == 0)
print("spec pairing minus bulk flux =", sp.simplify(face_spec - flux_face))

# Traction: t = -s_t * dp * n,  s_t = +1 is the supplied law t = -(dp) n (Lambda_X=0)
# Face virtual-power density t · v_face.  v_face = V * n, n = +w-hat on upper face.
# t · v_face = -s_t dp V
s_t = sp.symbols("s_t", real=True)
t_dot_v = -s_t * dp_f * V
face_from_t = sp.simplify(sp.Rational(1, 2) * sp.re(t_dot_v))
# Energy leaving the SLAB through the upper face is - t·v_face  if t is force-on-slab
# (S11b: d(T+U)/dt|_pressure = -Σ P_bulk with P_bulk = 1/2 Re(dp V*)).
# Using t = -s_t dp n, force-on-slab pairing:
#   work_on_slab = t · v_face = -s_t dp V
#   energy LEAVING slab = - work_on_slab = s_t dp V
leave_slab = sp.simplify(sp.Rational(1, 2) * sp.re(s_t * dp_f * sp.conjugate(V)))
print("energy leaving slab from traction, 1/2 Re(s_t dp V*) =", leave_slab)
print("leave_slab - bulk_flux at s_t=+1 =", sp.simplify(leave_slab.subs(s_t, 1) - flux_face))
print("leave_slab - bulk_flux at s_t=-1 =", sp.simplify(leave_slab.subs(s_t, -1) - flux_face))

# Spec residual uses dp V*, not t. It does not contain s_t:
print("spec residual contains s_t?", (face_spec - flux_face).has(s_t))
print("traction residual contains s_t?", (leave_slab - flux_face).has(s_t))

# Left-quantized vs two-sided composition, 1-D Fourier witness (repeat of OBJ4 structure)
k, kp, qk, qkp, kappa = sp.symbols("k kprime q qprime kappa", complex=True)
hhat = sp.symbols("hhat")
# N0 M_h N0 kernel: i q(k) * hhat(k-k') * i q(k') = -q q' hhat
K_two_sided = qk * qkp  # after the overall sign of dN; both legs
K_left_quant = qk ** 2  # left-quantized a(x,k)=h(x) N0(k)^2 / missing q(k')
print("two-sided N0 h N0 kernel ~ q(k) q(k') : depends on q(k')?", K_two_sided.has(qkp))
print("left-quantized h(x) * N0(k)^2 kernel ~ q(k)^2 : depends on q(k')?", K_left_quant.has(qkp))
print("difference (two-sided - left) =", sp.simplify(K_two_sided - K_left_quant))
# Concrete: k=0, k'=3/5, kappa=1, q(0)=1, q(3/5)=4/5
num = {qk: 1, qkp: sp.Rational(4, 5)}
print("two-sided at (0,3/5) =", K_two_sided.subs(num))
print("left-quantized at (0,3/5) using q(k) only =", K_left_quant.subs(num))
print("they differ?", K_two_sided.subs(num) != K_left_quant.subs(num))
