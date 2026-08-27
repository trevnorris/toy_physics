"""Independent derivation of §3c density background normal jet from §2b members."""
import sympy as sp

# §2a: W_bg depends on in-plane anchor y, not on normal w
y1, y2, y3, w = sp.symbols("y1 y2 y3 w", real=True)
W0, rho_br = sp.symbols("W0 rho_br", positive=True)
eta, w1 = sp.symbols("eta_bg w1_profile", real=True)
# W_bg(y) = W0*(1 + eta*w1(y)); treat w1 as a symbol standing for the profile value at y
W_bg = W0 * (1 + eta * w1)

# §2b representatives — both functions of the in-plane anchor only
rho4_ref = rho_br / W0
rho4_bg_RHO4 = rho4_ref                          # constant in y and w
rho4_bg_RHOBR = rho_br / W_bg                    # depends on y through W_bg, not on w

# §3c law: delta[f(x,h_s)] = delta f(x,h_s0) + delta h_s * d_w f0(x,h_s0)
# Background normal derivative must come from differentiating the §2b member
dw_RHO4 = sp.diff(rho4_bg_RHO4, w)
dw_RHOBR = sp.diff(rho4_bg_RHOBR, w)

# Face values (background) — also w-independent, so evaluation at h_s0 is identity
print("RHO4_CONSTANT: rho4_bg =", rho4_bg_RHO4)
print("RHO4_CONSTANT: d_w rho4_bg =", dw_RHO4)
print("RHOBR_CONSTANT: rho4_bg =", rho4_bg_RHOBR)
print("RHOBR_CONSTANT: d_w rho4_bg =", dw_RHOBR)
print("both_dw_are_zero =", (dw_RHO4 == 0) and (dw_RHOBR == 0))

# Contrast: an ungrounded free premise rhoBulkBackground(x,w,t) keeps a live d/dw
rhoBulkBackground = sp.Function("rhoBulkBackground")
x1, x2, x3, t = sp.symbols("x1 x2 x3 t")
free = rhoBulkBackground(x1, x2, x3, w, t)
print("ungrounded free premise d_w =", sp.diff(free, w))
