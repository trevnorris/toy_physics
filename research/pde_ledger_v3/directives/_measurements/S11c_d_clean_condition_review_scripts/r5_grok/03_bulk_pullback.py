#!/usr/bin/env python3
"""Independent derivation of the class-P rest-frame potential-flow pullback.

Harmonic convention of S11b: exp(i(k·x - omega t)).
Acoustics: v = grad_4 phi, dp = -rho_m d_t phi, d_t^2 phi = c^2 lap_4 phi.
Mass: d_t rho + rho_m div v = 0 in the rest-frame bulk (rho_4D^0 = rho_m).
"""
import sympy as sp

x1, x2, x3, w, t = sp.symbols("x1 x2 x3 w t", real=True)
k2, omega, rho_m, c_s0 = sp.symbols("k2 omega rho_m c_s0", nonzero=True, real=True)
Phi = sp.Function("Phi")(x1, w)
phase = sp.exp(sp.I * (k2 * x2 - omega * t))
phi = Phi * phase
coords = (x1, x2, x3, w)


def amp(expr):
    return sp.simplify(sp.expand(expr / phase))


v = [sp.diff(phi, c) for c in coords]
dp = -rho_m * sp.diff(phi, t)
wave = sp.diff(phi, t, 2) - c_s0**2 * sum(sp.diff(phi, c, 2) for c in coords)
phi_ww = sp.solve(sp.Eq(sp.simplify(wave / phase), 0), sp.diff(Phi, w, 2))[0]

print("DELTA_P", amp(dp))
print("V_BULK", [amp(c) for c in v])
print("V3_IS_ZERO", amp(v[2]))
print("D_W_DELTA_P", amp(sp.diff(dp, w)))
print(
    "D_W_V",
    [
        sp.simplify(amp(sp.diff(c, w)).subs(sp.diff(Phi, w, 2), phi_ww))
        for c in v
    ],
)

# Continuity for delta_rho
dr = sp.Function("dr")(x1, w) * phase
cont = sp.diff(dr, t) + rho_m * sum(sp.diff(v[i], coords[i]) for i in range(4))
dr_amp = sp.solve(sp.Eq(sp.simplify(cont / phase), 0), sp.Function("dr")(x1, w))[0]
dr_amp = sp.simplify(dr_amp.subs(sp.diff(Phi, w, 2), phi_ww))
print("DELTA_RHO", dr_amp)
print("DELTA_RHO_MINUS_DP_OVER_C2", sp.simplify(dr_amp - amp(dp) / c_s0**2))
print("D_T_DELTA_RHO_ENVELOPE", sp.simplify(-sp.I * omega * dr_amp))
print("DELTA_J", [sp.simplify(rho_m * amp(c)) for c in v])

# Pairing over one x2 period
A, B = sp.symbols("trial_amplitude test_amplitude")
period = 2 * sp.pi / k2
trial = A * sp.exp(sp.I * k2 * x2)
test_same = B * sp.exp(sp.I * k2 * x2)
test_conj = B * sp.exp(-sp.I * k2 * x2)
print("L2_SAME_PHASE", sp.simplify(sp.integrate(test_same * trial, (x2, 0, period))))
print("L2_CONJUGATE", sp.simplify(sp.integrate(test_conj * trial, (x2, 0, period))))

# Affine trace: engine coordinates at w = s W0/2
s, W0, eta = sp.symbols("s W0 eta", real=True)
w1 = sp.Function("w1")(x1)
PhiF = sp.Function("Phi")
w_ref = s * W0 / 2
w_bg = s * W0 * (1 + eta * w1) / 2
exact = PhiF(x1, w_bg)
aff_ref = PhiF(x1, w_ref) + (w_bg - w_ref) * sp.diff(PhiF(x1, w), w).subs(w, w_ref)
aff_bg = PhiF(x1, w_bg) + (w_bg - w_ref) * sp.diff(PhiF(x1, w), w).subs(w, w_bg)


def first_eta(expr):
    return sp.simplify(sp.series(expr, eta, 0, 2).removeO().doit())


print("EXACT_MINUS_REF_AFFINE", first_eta(exact - aff_ref))
print("EXACT_MINUS_BG_AFFINE", first_eta(exact - aff_bg))
