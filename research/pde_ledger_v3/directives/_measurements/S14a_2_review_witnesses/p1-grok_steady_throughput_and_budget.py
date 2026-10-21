#!/usr/bin/env python3
"""Objects for the S14a part 2 decision-list review.

Prints residuals and values. No verdicts.
"""

import sympy as sp

r, w, t = sp.symbols("r w t", real=True)
h0, s, n0, C0 = sp.symbols("h0 s n0 C0", positive=True)
b, gamma0 = sp.symbols("b gamma0", real=True)
omega, B = sp.symbols("omega B", real=True)


def material_residual(n, ur, uw):
    return sp.simplify(
        sp.diff(n, t)
        + (1 / r**2) * sp.diff(r**2 * n * ur, r)
        + sp.diff(n * uw, w)
    )


def order_residual(n, ur, uw, chi, gamma_b):
    return sp.simplify(
        sp.diff(chi * n, t)
        + (1 / r**2) * sp.diff(r**2 * (chi * n * ur), r)
        + sp.diff(chi * n * uw, w)
        - n * gamma_b
    )


def projected_storage(n, chi, h):
    sigma = sp.integrate(n, (w, -h / 2, h / 2))
    sigma_chi = sp.integrate(chi * n, (w, -h / 2, h / 2))
    return sp.simplify(sigma), sp.simplify(sp.diff(sigma, t)), sp.simplify(sigma_chi), sp.simplify(
        sp.diff(sigma_chi, t)
    )


h = h0 + s * t
xi = w / h
n = n0
uw = (s / h) * w
ur = -(s / h) * r / 3 + C0 / r**2

print("SECTION 1 homologous spreading flow")
print("material_residual", material_residual(n, ur, uw))

# Family S: rates are functions of the spreading coordinate only.
chi_s = sp.Integer(1) / 2
gamma_return_s = xi**2
gamma_drain_s = xi**2
gamma_b_s = sp.simplify(gamma_return_s - gamma_drain_s)
print("family_S_gamma_B", gamma_b_s)
print("family_S_order_residual", order_residual(n, ur, uw, chi_s, gamma_b_s))
dt_n_gamma_return_s = sp.simplify(sp.diff(n * gamma_return_s, t))
print("family_S_dt_n_Gamma_return", dt_n_gamma_return_s)
print(
    "family_S_dt_n_Gamma_return_at_w_h0_t0",
    sp.simplify(dt_n_gamma_return_s.subs({w: h0, t: 0})),
)
sigma_s, dt_sigma_s, sigma_chi_s, dt_sigma_chi_s = projected_storage(n, chi_s, h)
print("family_S_projected_material", sigma_s)
print("family_S_dt_projected_material", dt_sigma_s)
print("family_S_dt_projected_order", dt_sigma_chi_s)
print("family_S_face_speed_minus_interior_uw_lower", sp.simplify(uw.subs(w, -h / 2) - sp.diff(-h / 2, t)))
print("family_S_ur", ur)
print("family_S_uw", uw)

# Family E: Eulerian-constant source densities, secular order, same spreading flow.
chi_e = sp.Integer(1) / 2 + b * t
gamma_return_e = b + gamma0
gamma_drain_e = gamma0
gamma_b_e = sp.simplify(gamma_return_e - gamma_drain_e)
print("SECTION 2 Eulerian-constant source densities")
print("family_E_order_residual", order_residual(n, ur, uw, chi_e, gamma_b_e))
print("family_E_dt_n_Gamma_return", sp.simplify(sp.diff(n * gamma_return_e, t)))
print("family_E_dt_n_Gamma_drain", sp.simplify(sp.diff(n * gamma_drain_e, t)))
print("family_E_dt_chi", sp.simplify(sp.diff(chi_e, t)))
_, dt_sigma_e, _, dt_sigma_chi_e = projected_storage(n, chi_e, h)
print("family_E_dt_projected_material", dt_sigma_e)
print("family_E_dt_projected_order", dt_sigma_chi_e)
# Range on a finite interval, concrete numbers.
chi_e_num = sp.Lambda(t, chi_e.subs(b, sp.Rational(1, 10)))
print("family_E_chi_at_0", chi_e_num(0))
print("family_E_chi_at_1", chi_e_num(1))

# Family A: faces fixed, Eulerian-stationary, flow only the core piece.
h_a = h0
ur_a = C0 / r**2
uw_a = 0
chi_a = sp.Integer(1) / 2
gamma_return_a = gamma0
gamma_drain_a = gamma0
print("family_A_material_residual", material_residual(n, ur_a, uw_a))
print(
    "family_A_order_residual",
    order_residual(n, ur_a, uw_a, chi_a, sp.simplify(gamma_return_a - gamma_drain_a)),
)
print("family_A_dt_n_Gamma_return", sp.diff(n * gamma_return_a, t))
sigma_a, dt_sigma_a, _, dt_sigma_chi_a = projected_storage(n, chi_a, h_a)
print("family_A_dt_projected_material", dt_sigma_a)
print("family_A_dt_projected_order", dt_sigma_chi_a)
print("family_A_minus_family_S_dt_projected_material", sp.simplify(dt_sigma_a - dt_sigma_s))

# Product versus rate.
n_var, gamma_var = sp.symbols("n_var gamma_var", cls=sp.Function)
print("SECTION 3 product versus rate")
n1 = 2 + t
g1 = 3
g2 = 1 / (2 + t)
print("dt_Gamma_when_product_changes", sp.diff(g1, t))
print("dt_n_Gamma_when_rate_constant", sp.diff(n1 * g1, t))
print("dt_n_Gamma_when_product_constant", sp.simplify(sp.diff(n1 * g2, t)))
print("dt_Gamma_when_product_constant", sp.simplify(sp.diff(g2, t)))

# omega -> 0 of a modulation on the spreading profile.
print("SECTION 4 omega to 0 on the spreading profile")
gamma_family = xi**2 * (1 + B * sp.sin(omega * t))
gamma_limit = sp.limit(gamma_family, omega, 0)
print("omega_limit_of_Gamma_return", sp.simplify(gamma_limit))
print("dt_of_that_limit", sp.simplify(sp.diff(n * gamma_limit, t)))
print(
    "dt_of_that_limit_at_w_h0_t0",
    sp.simplify(sp.diff(n * gamma_limit, t).subs({w: h0, t: 0})),
)
# A storage term that a frequency limit drops, and the piece it keeps.
sigma_mode = sp.symbols("Sigma0")
sigma_rate, amplitude = sp.symbols("sigma_rate A", real=True)
storage_mode = sp.diff(sigma_mode + sigma_rate * t + amplitude * sp.sin(omega * t), t)
print("storage_before_limit", sp.simplify(storage_mode))
print("storage_omega_limit", sp.simplify(sp.limit(storage_mode, omega, 0)))
print("storage_dropped_term", sp.simplify(storage_mode - sp.limit(storage_mode, omega, 0)))

# Material balance has no source densities. Expression difference under the constraint is empty.
print("SECTION 5 material balance contains neither rate")
n_sym, ur_sym, uw_sym = sp.symbols("n u_r u_w", cls=sp.Function)
material_expr = sp.diff(n_sym(r, w, t), t) + (1 / r**2) * sp.diff(
    r**2 * n_sym(r, w, t) * ur_sym(r, w, t), r
) + sp.diff(n_sym(r, w, t) * uw_sym(r, w, t), w)
print("Gamma_return_in_material", material_expr.has(sp.Symbol("Gamma_return")))
print("material_free_symbols", sorted(str(s) for s in material_expr.free_symbols))
print("material_expression_minus_itself", sp.simplify(material_expr - material_expr))

# Order-budget assemblies.
print("SECTION 6 order budget assemblies")
fna_m, b2_m, b1, fna_p, b2_p, b3 = sp.symbols(
    "fna_minus b2_minus b1_minus fna_plus b2_plus b3_plus"
)
int_n_return, int_n_drain = sp.symbols("int_n_Gamma_return int_n_Gamma_drain")
lower_before = -(fna_m + b2_m)
upper_before = fna_p + b2_p
lower_after = -b1
upper_after = b3
c_minus = fna_m + b2_m - b1
c_plus = -fna_p + b3 - b2_p
int_n_gamma = int_n_return - int_n_drain
# Interior integration of the order equation, per unit area:
# dt + div + before_sum - int_n_gamma = 0.
# Re-derived full-source form uses after_sum and face coefficients.
before_sum = lower_before + upper_before
after_sum = lower_after + upper_after
full_source = int_n_gamma + c_minus + c_plus
# If the interior identity holds, this residual is the full balance's failure
# when crossing is taken to be the before terms and the source includes faces.
residual_before_plus_full_source = sp.simplify(before_sum - after_sum + (full_source - int_n_gamma) - (c_minus + c_plus) + (after_sum - before_sum - (c_minus + c_plus)))
# Directly: full identity is dt+div+after_sum - full_source = dt+div+before_sum - int_n_gamma
#   + (after_sum - before_sum - c_minus - c_plus)
gap_after_versus_before = sp.simplify(after_sum - before_sum - c_minus - c_plus)
print("C_minus", sp.simplify(c_minus))
print("C_plus", sp.simplify(c_plus))
print("after_minus_before_minus_C", gap_after_versus_before)
print("before_minus_after", sp.simplify(before_sum - after_sum))
# Concrete nonzero face jump.
concrete = {fna_m: 5, b2_m: 7, b1: 3, fna_p: 11, b2_p: 13, b3: 17}
print("concrete_C_minus", c_minus.subs(concrete))
print("concrete_C_plus", c_plus.subs(concrete))
print("concrete_before_minus_after", (before_sum - after_sum).subs(concrete))
print("concrete_gap_after_before_C", gap_after_versus_before.subs(concrete))

# Sphere orientation.
a, R = sp.symbols("a R", positive=True)
H = sp.Function("H_r")
divergence_integral = sp.integrate(sp.diff(r**2 * H(r), r), (r, a, R))
outward_from_origin = lambda radius: radius**2 * H(radius)
split_that_matches = sp.simplify(
    outward_from_origin(R) - outward_from_origin(a) - divergence_integral
)
split_both_plus = sp.simplify(
    outward_from_origin(R) + outward_from_origin(a) - divergence_integral
)
print("orientation_residual_outer_minus_inner", split_that_matches)
print("orientation_residual_both_plus", split_both_plus)
print(
    "orientation_residual_both_plus_value",
    split_both_plus.subs({H(a): 4, H(R): 9, a: 2, R: 5}),
)

# Totals against part 1's order-lost sign.
order_lost = -(int_n_gamma + c_minus + c_plus)
total_plain_difference = int_n_return - int_n_drain
total_order_lost_sign = -int_n_return - (-int_n_drain)
print("order_lost_minus_drain_plus_return_plus_faces", sp.simplify(order_lost - (-int_n_return + int_n_drain - c_minus - c_plus)))
print(
    "plain_totals_minus_order_lost_without_faces",
    sp.simplify(total_plain_difference + order_lost),
)
print(
    "plain_totals_relation_residual_without_faces_concrete",
    sp.simplify(total_plain_difference + order_lost).subs(
        {int_n_return: 8, int_n_drain: 3, **concrete}
    ),
)
print("order_lost_sign_totals_difference_from_plain", sp.simplify(total_order_lost_sign - total_plain_difference))

# Dropping J_chi from the in-plane order flux.
print("SECTION 7 order flux and far radius")
Jr, Jw = sp.symbols("J_r J_w", cls=sp.Function)
chi_f = sp.Function("chi")
n_f = sp.Function("n")
ur_f = sp.Function("u_r")
uw_f = sp.Function("u_w")
full_radial = chi_f(r, w, t) * n_f(r, w, t) * ur_f(r, w, t) + Jr(r, w, t)
advective_radial = chi_f(r, w, t) * n_f(r, w, t) * ur_f(r, w, t)
print(
    "radial_order_flux_minus_advective",
    sp.simplify(full_radial - advective_radial),
)
# Shell face channel depends on the outer radius when the face density is live.
face_density = sp.Function("face_density")
area = 4 * sp.pi * r**2
face_channel = sp.integrate(area * face_density(r, t), (r, a, R))
print("d_face_channel_d_outer", sp.simplify(sp.diff(face_channel, R)))
print(
    "d_face_channel_d_outer_at_sample",
    sp.simplify(sp.diff(face_channel, R).subs({face_density(R, t): 1, R: 5, sp.pi: sp.pi})),
)
# Lower limit of the shell integral is the inner radius, not the origin.
shell_integral = sp.Integral(area * face_density(r, t), (r, a, R))
print("shell_integral", shell_integral)
print("shell_limits", shell_integral.limits)
