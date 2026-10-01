#!/usr/bin/env python3
"""Leg script 06: small symbolic model of the supplied S11c-a/S11c-b face laws (spec
S11c_a §3b: n_s from the graph h_s, v_face = d_t R_s, V_s = n.v_face,
J_s = rho_m (v_bulk - v_face).n, A_s = mu_theta/rho_br - dp/rho_m,
t_s = -(dp + Lambda_X A_s) n, virtual work a_s t_s . delta_v x_s), on R1 (profile depends on
x1 only) and class P (fields independent of x3), first order in the wave amplitude.

Part A -- control coverage.  The directive's FORM controls K1 = a1 e_i d_k u_i d_k theta,
e = (sin b, 0, cos b), and K2 = a2 theta eps_ijk g_i d_j u_k enter only the stored energy,
hence only mu_theta = dU/dtheta (variational).  Print, for baseline / K1 / K2 / a
background-datum control G3 (a fixed direction-3 first jet g3 of W_bg, NOT transformed by
the mirror), the entries of the odd<->even blocks of the face map:
   d(V_s)/d(u3_t), d(J_s)/d(vb3), d(J_s)/d(u3 jets), d(A_s)/d(u3 jets),
   d(face force on delta_v u_3)/d(dp), d(face force on delta_v u_3)/d(vb_w).
Part B -- pullback evaluation height: the engine's affine trace model
   f(h) = f_ref + (h - s W0/2) f_w,ref   (S11c_a audit :578-592)
evaluated with the pullback Phi taken at the FLAT reference face (w = s W0/2) versus at the
BACKGROUND face (w = s W_bg/2).  Print the difference of the composed trace.
Part C -- weak pairing of a class-P test field: integrate v * d2 u over one x2 period with
(i) test field in class P (same exponential), (ii) conjugate exponential.
Prints computed objects only.
"""
import sympy as sp

x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
s = sp.Symbol("s")  # face sign +-1
eps, eta, W0, rho_m, rho_br, LX, LA, LV, B = sp.symbols("epsilon eta W0 rho_m rho_br Lambda_X Lambda_A Lambda_V B_rho")
a1, a2, beta, g3 = sp.symbols("a_K1 a_K2 beta g3", real=True)
k2, om = sp.symbols("k2 omega", positive=True)

w1 = sp.Function("w1")(x1)                 # R1 profile, x1 only
# fields: class P -> functions of (x1, x2, t) only
U = [sp.Function(f"u{i}")(x1, x2, t) for i in (1, 2, 3)]
TH = sp.Function("theta")(x1, x2, t)
ZC = sp.Function("zeta_c")(x1, x2, t)
DP = sp.Symbol("dp_s")
VB = [sp.Symbol(f"vb_s{i}") for i in (1, 2, 3, 4)]       # bulk velocity trace components (4 = w)
DVU = [sp.Symbol(f"dvu{i}") for i in (1, 2, 3)]          # virtual in-plane displacement
DVZ = sp.Symbol("dvzeta")

def model(control):
    # background thickness: R1 profile; control G3 adds a FIXED direction-3 slope g3*x3
    Wbg = W0 * (1 + eta * w1) + (g3 * x3 if control == "G3" else 0)
    X = (x1, x2, x3)
    # stored energy pieces that matter for mu_theta (variational derivative in theta)
    gradU = [[sp.diff(U[i], X[k]) for k in range(3)] for i in range(3)]
    gradTH = [sp.diff(TH, X[k]) for k in range(3)]
    dens = B * TH**2 / 2
    if control == "K1":
        e = (sp.sin(beta), 0, sp.cos(beta))
        dens += a1 * sum(e[i] * gradU[i][k] * gradTH[k] for i in range(3) for k in range(3))
    if control == "K2":
        g = [sp.diff(Wbg, X[i]) for i in range(3)]
        dens += a2 * TH * sum(sp.LeviCivita(i, j, k) * g[i] * gradU[k][j]
                              for i in range(3) for j in range(3) for k in range(3))
    # variational derivative dU/dtheta
    mu = sp.diff(dens, TH) - sum(sp.diff(sp.diff(dens, gradTH[k]), X[k]) for k in range(3) if gradTH[k] != 0)
    # LAB_HELD face graph h_s = s Wbg/2 + zeta_c (+ s dW/2 omitted: thickness DOF not needed here)
    h = s * Wbg / 2 + eps * ZC
    gh = [sp.diff(h, X[k]) for k in range(3)]
    den = sp.sqrt(1 + sum(c**2 for c in gh))
    n = [-s * c / den for c in gh] + [s / den]
    # face velocity at fixed material label (LAB_HELD): in-plane u_t, vertical zeta_t + s u_t.grad Wbg/2
    ut = [sp.diff(eps * U[i], t) for i in range(3)]
    gW = [sp.diff(Wbg, X[k]) for k in range(3)]
    vface = ut + [sp.diff(eps * ZC, t) + s * sum(ut[k] * gW[k] for k in range(3)) / 2]
    vb = [eps * c for c in VB]
    V = sum(n[i] * vface[i] for i in range(4))
    J = rho_m * sum((vb[i] - vface[i]) * n[i] for i in range(4))
    A = eps * mu / rho_br - eps * DP / rho_m
    trac = [-(eps * DP + LX * A) * n[i] for i in range(4)]
    dvx = DVU + [DVZ + s * sum(DVU[k] * gW[k] for k in range(3)) / 2]
    work = den * sum(trac[i] * dvx[i] for i in range(4))
    lin = lambda e_: sp.expand(sp.series(e_, eps, 0, 2).removeO().coeff(eps, 1))
    return {"V": lin(V), "J": lin(J), "A": lin(A), "WORK": lin(work)}

def coeff_of(expr, target):
    return sp.simplify(sp.diff(expr, target))

for control in ("BASELINE", "K1", "K2", "G3"):
    M = model(control)
    u3t = sp.Derivative(U[2], t)
    u3_1 = sp.Derivative(U[2], x1)
    u3_2 = sp.Derivative(U[2], x2)
    u3_11 = sp.Derivative(U[2], (x1, 2))
    u3_22 = sp.Derivative(U[2], (x2, 2))
    def d_by(expr, target):
        # coefficient of a derivative object: replace it by a dummy then differentiate
        dmy = sp.Dummy()
        return sp.simplify(sp.diff(expr.subs(target, dmy), dmy))
    face_u3 = sp.diff(M["WORK"], DVU[2])   # generalized force on delta_v u_3
    entries = {
        "dV/du3_t": d_by(M["V"], u3t),
        "dJ/dvb3": sp.simplify(sp.diff(M["J"], VB[2])),
        "dJ/du3_t": d_by(M["J"], u3t),
        "dJ/d(u3_x2)": d_by(M["J"], u3_2),
        "dA/d(u3_x2)": d_by(M["A"], u3_2),
        "dA/d(u3_x1x1)": d_by(M["A"], u3_11),
        "dA/d(u3_x2x2)": d_by(M["A"], u3_22),
        "dFaceForce_u3/d(dp)": sp.simplify(sp.diff(face_u3, DP)),
        "dFaceForce_u3/d(vb_w)": sp.simplify(sp.diff(face_u3, VB[3])),
        "dFaceForce_u3/d(theta)": d_by(face_u3, TH),
    }
    for key, val in entries.items():
        print("PART_A", control, key, val)

# ---- Part B: pullback evaluation height --------------------------------------------
Phi = sp.Function("Phi")
zc = sp.Symbol("zeta")
h_phys = s * W0 * (1 + eta * w1) / 2 + eps * zc        # physical face height
w_ref = s * W0 / 2                                     # engine flat reference face
def affine_trace(f_ref, f_w_ref):
    return f_ref + (h_phys - w_ref) * f_w_ref
exact = Phi(x1, h_phys)
ref_eval = affine_trace(Phi(x1, w_ref), sp.diff(Phi(x1, w), w).subs(w, w_ref))
w_bgface = s * W0 * (1 + eta * w1) / 2
bg_eval = affine_trace(Phi(x1, w_bgface), sp.diff(Phi(x1, w), w).subs(w, w_bgface))
def first_order(e_):
    e_ = sp.series(e_.subs(eps, 0), eta, 0, 2).removeO()
    return sp.simplify(e_.doit())
print("PART_B EXACT_MINUS_REFERENCE_FACE_PULLBACK_O(eta)", first_order(exact - ref_eval))
print("PART_B EXACT_MINUS_BACKGROUND_FACE_PULLBACK_O(eta)", first_order(exact - bg_eval))
# in-plane derivative of a trace taken at the background face vs trace of the derivative
tr_bg = Phi(x1, s * W0 * (1 + eta * w1) / 2)
d_of_trace = sp.diff(tr_bg, x1)
trace_of_d = sp.Subs(sp.Derivative(Phi(x1, w), x1), w, s * W0 * (1 + eta * w1) / 2)
print("PART_B D1_OF_BGFACE_TRACE_MINUS_TRACE_OF_D1", sp.simplify((d_of_trace - trace_of_d.doit())))

# ---- Part C: weak pairing with class-P test functions over one x2 period ----------------
Uamp, Vamp = sp.symbols("U_amp V_amp")
Lp = 2 * sp.pi / k2
u_trial = Uamp * sp.exp(sp.I * k2 * x2)
v_same = Vamp * sp.exp(sp.I * k2 * x2)
v_conj = Vamp * sp.exp(-sp.I * k2 * x2)
print("PART_C PAIRING_SAME_CLASS_P", sp.simplify(sp.integrate(v_same * sp.diff(u_trial, x2), (x2, 0, Lp))))
print("PART_C PAIRING_CONJUGATE", sp.simplify(sp.integrate(v_conj * sp.diff(u_trial, x2), (x2, 0, Lp))))
