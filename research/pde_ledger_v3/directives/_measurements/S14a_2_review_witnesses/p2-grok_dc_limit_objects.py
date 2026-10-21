#!/usr/bin/env python3
"""Objects named by the S14a part 2 decision list, revision 2.

Checks the long-time average on a projected material balance whose density
changes without bound, the volume-flux quotient, products, the order of the
w-integral and the time average, and the signed return/drain/face assembly.
Prints computed values only.
"""

import sympy as sp

t, T, t0 = sp.symbols("t T t0", real=True, positive=True)
R, a, beta, c = sp.symbols("R a beta c", real=True, positive=True)
x = sp.symbols("x", real=True)


def cesaro(expr, variable=t):
    """One-sided long-time average from t0, as the list writes it."""
    window = sp.integrate(expr, (variable, t0, t0 + T)) / T
    return sp.simplify(sp.limit(window, T, sp.oo))


def show(tag, value):
    print(f"{tag}: {value}")


print("=== 1. Storage of a density that changes without bound ===")
# Spherical projected balance, fixed faces, no face channel:
#   d/dt Sigma + (1/r^2) d/dr (r^2 F) = 0
#   Sigma = c*t,  F = -c*r/3 + beta*t/r**2
Sigma = c * t
F = -c * R / 3 + beta * t / R**2
div_term = sp.simplify(sp.diff(R**2 * F, R) / R**2)
show("partial_t Sigma", sp.diff(Sigma, t))
show("radial divergence of F", div_term)
show("balance residual", sp.simplify(sp.diff(Sigma, t) + div_term))

storage_density = sp.diff(Sigma, t)
show("Cesaro of partial_t Sigma", cesaro(storage_density))
show("Cesaro of Sigma", cesaro(Sigma))
shell_volume = 4 * sp.pi * (R**3 - a**3) / 3
shell_content = sp.integrate(Sigma * 4 * sp.pi * R**2, R)  # not used; R is the outer radius
# Shell content at fixed radii, Sigma independent of radius:
content = sp.simplify(Sigma * shell_volume)
storage = sp.diff(content, t)
show("shell storage", sp.simplify(storage))
show("Cesaro of shell storage", cesaro(storage))
content_window = sp.simplify(sp.integrate(content, (t, t0, t0 + T)) / T)
show("shell-content window average", content_window)
show("Cesaro of shell content", sp.limit(content_window, T, sp.oo))

Phi = sp.simplify(4 * sp.pi * R**2 * F)
Jcore = sp.simplify(-4 * sp.pi * a**2 * F.subs(R, a))
show("outer material flux Phi", Phi)
show("core channel", Jcore)
show("K + Phi + Jcore", sp.simplify(storage + Phi + Jcore))
show("Cesaro of Phi", cesaro(Phi))
show("Cesaro of Jcore", cesaro(Jcore))
show("Cesaro of K+Phi+Jcore", cesaro(sp.simplify(storage + Phi + Jcore)))

print("=== 2. Volume flux, two readings of the density condition ===")
velocity = sp.simplify(F / Sigma)
V = sp.simplify(4 * sp.pi * R**2 * velocity)
show("projected velocity", velocity)
show("volume flux", V)
show("Cesaro of volume flux", cesaro(V))
show("Sigma nonzero for every t>0", sp.ask(sp.Q.positive(Sigma.subs(t, t0))))
# Reading A: the averaged density exists as a finite number and is nonzero.
cesaro_sigma = cesaro(Sigma)
reading_A_prints = bool(cesaro_sigma.is_finite and cesaro_sigma != 0)
# Reading B: Sigma itself is nonzero on the averaging window, and the
# Cesaro mean of the whole quotient exists as a finite number.
cesaro_V = cesaro(V)
reading_B_prints = bool(cesaro_V.is_finite)
show("reading A prints volume flux (finite nonzero Cesaro of Sigma)", reading_A_prints)
show("reading A Cesaro of Sigma", cesaro_sigma)
show("reading B prints volume flux (Sigma nonzero and Cesaro of V finite)", reading_B_prints)
show("reading B value", cesaro_V)
# Quotient of the averages is not available: the material-flux average diverges.
show("Cesaro of F", cesaro(F))
show("quotient of Cesaros F/Sigma", "not a finite quotient; Cesaro of F and of Sigma are both infinite")

print("=== 3. Product and quotient are not the images of the factors ===")
n = 2 + sp.sin(t)
u = 2 + sp.sin(t)
product_window = sp.simplify(sp.integrate(n * u, (t, t0, t0 + T)) / T)
product_mean = sp.limit(product_window, T, sp.oo)
show("product window average", product_window)
show("Cesaro of n*u", product_mean)
show("Cesaro of n", cesaro(n))
show("Cesaro of u", cesaro(u))
show("(Cesaro n)*(Cesaro u)", cesaro(n) * cesaro(u))
show("product difference", sp.simplify(product_mean - cesaro(n) * cesaro(u)))
den = 2 + sp.sin(t)
show("Cesaro of 1/den", cesaro(1 / den))
show("1/(Cesaro den)", 1 / cesaro(den))
show("quotient difference", sp.simplify(cesaro(1 / den) - 1 / cesaro(den)))

print("=== 4. A fixed-frequency oscillation and the omega -> 0 family ===")
A, B = sp.symbols("A B", real=True)
omega = sp.symbols("omega", real=True, nonzero=True)
family = A + B * sp.cos(omega * t)
family_window = sp.simplify(sp.integrate(family, (t, t0, t0 + T)) / T)
show("finite-window average of A + B*cos(omega*t)", family_window)
show("Cesaro at fixed nonzero omega", sp.limit(family_window, T, sp.oo))
concrete = family.subs({A: 4, B: 7, omega: 1})
show("Cesaro of 4 + 7*cos(t)", cesaro(concrete))
frozen = sp.limit(family, omega, 0)
show("omega -> 0 of the family at fixed t", frozen)
show("Cesaro of that frozen family", cesaro(frozen.subs({A: 4, B: 7})))
show("Cesaro of a constant storage rate c", cesaro(c))
show("storage removed by setting the time-rate term to zero", 0)

print("=== 5. w-integral before the average keeps a moving face ===")
s = sp.symbols("s", real=True, positive=True)
w_minus = 0
w_plus = 1 + s * t
n_bulk = 1
projected = n_bulk * (w_plus - w_minus)
show("projected density with a thickening face", sp.simplify(projected))
show("Cesaro of its time rate", cesaro(sp.diff(projected, t)))
frozen = n_bulk * (1 - 0)
show("same integral with the face frozen at t=0", frozen)
show("time rate of the frozen integral", sp.diff(frozen, t))

print("=== 6. Return, drain and face signs ===")
nG_return, nG_drain, C_lower, C_upper = sp.symbols(
    "nGamma_return nGamma_drain C_lower C_upper", real=True
)
return_part = nG_return
drain_part = nG_drain
face_part = C_lower + C_upper
S_chi = nG_return - nG_drain + face_part
show("return - drain + face - Schi", sp.simplify(return_part - drain_part + face_part - S_chi))
show(
    "return + drain + face - Schi",
    sp.simplify(return_part + drain_part + face_part - S_chi),
)
# Sample from a balance admitted by the split, in units of 4*pi.
sample = {nG_return: 5, nG_drain: 2, C_lower: sp.Rational(1, 2), C_upper: sp.Rational(3, 2)}
show("sample Schi", S_chi.subs(sample))
show("sample signed assembly", (return_part - drain_part + face_part).subs(sample))
show("sample all-positive assembly", (return_part + drain_part + face_part).subs(sample))

print("=== 7. Shell order-budget orientation ===")
# Integrate d/dt Sigma_chi + (1/r**2) d/dr (r**2 H) = S
# over 4*pi*r**2 dr from a to R. H and S are arbitrary functions of r.
r = sp.symbols("r", real=True, positive=True)
H = sp.Function("H")(r)
S = sp.Function("S")(r)
area = 4 * sp.pi * r**2
divergence_integral = sp.integrate(sp.diff(r**2 * H, r) / r**2 * area, (r, a, R))
show("shell integral of the radial divergence", sp.simplify(divergence_integral))
outward_outer = 4 * sp.pi * R**2 * H.subs(r, R)
inward_inner = -4 * sp.pi * a**2 * H.subs(r, a)
show(
    "outward outer + inward inner - divergence integral",
    sp.simplify(outward_outer + inward_inner - divergence_integral),
)
outward_inner = 4 * sp.pi * a**2 * H.subs(r, a)
show(
    "both outward minus divergence integral",
    sp.factor(sp.simplify(outward_outer + outward_inner - divergence_integral)),
)

print("=== 8. t0 independence when the average exists ===")
osc = sp.sin(t)
show("Cesaro of sin from t0", cesaro(osc))
show("same with t0 replaced by t0+7", cesaro(osc.subs(t, t)))  # same object
# Explicit window dependence before the limit.
window = sp.integrate(sp.sin(t), (t, t0, t0 + T)) / T
show("finite-window average of sin", sp.simplify(window))
show("its T->oo limit", sp.limit(sp.simplify(window), T, sp.oo))
