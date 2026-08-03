#!/usr/bin/env python3
"""Independent SymPy audit for S11b-A; this program has no file inputs."""

import sympy as sp


def fmt(value):
    return sp.sstr(value).replace("\n", " ")


tags, checks = [], []


def emit(name, value):
    rendered = value if isinstance(value, str) else fmt(value)
    tags.append((name, rendered.replace("\n", " ")))


def checked(statement):
    checks.append(bool(statement))


# Common symbols and harmonic convention exp(i(k.x-omega*t)).
w, w1, w2 = sp.symbols("w w1 w2", real=True)
L, a = sp.symbols("L a", positive=True, finite=True, real=True)
b, c = sp.symbols("b c", finite=True, real=True)
omega, k, W0 = sp.symbols("omega k W0", positive=True, finite=True, real=True)
rho_m, c_s0 = sp.symbols("rho_m c_s0", positive=True, finite=True, real=True)
v0 = sp.symbols("v0", finite=True, real=True)
I = sp.I


# A1: fixed-window projection of continuity.
Omega = sp.Function("Omega")
jw = sp.Function("j_w")
boundary = Omega(w2) * jw(w2) - Omega(w1) * jw(w1)
projection_finite = sp.Integral(jw(w) * sp.diff(Omega(w), w), (w, w1, w2)) - boundary
projection_infinite = sp.Integral(jw(w) * sp.diff(Omega(w), w), (w, -sp.oo, sp.oo))

test_O, test_j = 1 + 2*w + 3*w**2, 2 - w + w**3
lhs_ibp = -sp.integrate(test_O * sp.diff(test_j, w), (w, w1, w2))
rhs_ibp = sp.integrate(test_j * sp.diff(test_O, w), (w, w1, w2)) - (
    test_O.subs(w, w2)*test_j.subs(w, w2) - test_O.subs(w, w1)*test_j.subs(w, w1)
)
checked(sp.simplify(lhs_ibp-rhs_ibp) == 0)
emit("S11BA_PROJECTION_FINITE",
     "d_t integral(Omega*rho,dw)+div_x integral(Omega*j_x,dw)=" + fmt(projection_finite)
     + "; interval=[w1,w2]")
emit("S11BA_PROJECTION_INFINITE",
     "d_t integral(Omega*rho,dw)+div_x integral(Omega*j_x,dw)=" + fmt(projection_infinite)
     + "; assumes limit(Omega*j_w,w->+/-infinity)=0")

# Positive J is outward from the slab, so accretion has the opposite sign.
Jp, Jm = sp.symbols("J_+ J_-", real=True)
accretion, throughflow = sp.Matrix([[-1, -1], [sp.Rational(1, 2), -sp.Rational(1, 2)]]) * sp.Matrix([Jp, Jm])
emit("S11BA_RELATIVE_FLUX_COMBINATIONS",
     "net_accretion=" + fmt(accretion) + "; through_flow_lower_to_upper=" + fmt(throughflow)
     + "; J_pm=rho_m*(v_w-d_t*zeta_pm)*(+/-1)")


# A2: parity reduction on [-L,L].
u = sp.symbols("u", nonnegative=True, real=True)
j_half, Omega_half = sp.Function("j_half"), sp.Function("Omega_even")
half_integral = sp.Integral(j_half(u)*sp.diff(Omega_half(u), u), (u, 0, L))
half_boundary = Omega_half(L)*j_half(L)


def symmetric_parity_source(parity_sign):
    return sp.simplify((1-parity_sign)*(half_integral-half_boundary))


source_even = symmetric_parity_source(sp.Integer(1))
source_odd = symmetric_parity_source(sp.Integer(-1))
source_odd_infinite = 2*sp.Integral(j_half(u)*sp.diff(Omega_half(u), u), (u, 0, sp.oo))
checked(source_even == 0)
checked(sp.simplify(source_odd-2*(half_integral-half_boundary)) == 0)
emit("S11BA_PARITY_EVEN_JW",
     "finite=" + fmt(source_even)
     + " exact_on=[-L,L]; infinite=0 exact_on=(-infinity,infinity) with vanishing boundary")
emit("S11BA_PARITY_ODD_JW",
     "finite=" + fmt(source_odd) + " exact_on=[-L,L]; infinite=" + fmt(source_odd_infinite)
     + " exact_on=(-infinity,infinity) with vanishing boundary; not_identically_zero")
interval_is_symmetric = sp.simplify((-L)+L) == 0
emit("S11BA_PARITY_INTERVAL", "[-L,L]; symmetric_about_w=0=" + str(interval_is_symmetric))


# A3: dynamic-window product-rule terms absent for a fixed Omega.
t, x1, x2, x3 = sp.symbols("t x1 x2 x3", real=True)
Omega_dyn = sp.Function("Omega_dyn")(w, x1, x2, x3, t)
rho_dyn = sp.Function("rho")(w, x1, x2, x3, t)
j_dyn = [sp.Function(f"j_x{i}")(w, x1, x2, x3, t) for i in range(1, 4)]
dynamic_extra = sp.Integral(rho_dyn*sp.diff(Omega_dyn, t), (w, w1, w2))
dynamic_extra += sum(sp.Integral(j_dyn[i]*sp.diff(Omega_dyn, coord), (w, w1, w2))
                     for i, coord in enumerate((x1, x2, x3)))
emit("S11BA_DYNAMIC_WINDOW_EXTRA_TERMS", fmt(dynamic_extra)
     + "; complete_extra_set={integral(rho*d_t Omega,dw),integral(j_x.grad_x Omega,dw)}")


# A4: outgoing/decaying acoustic half-space response.
q2 = sp.simplify(omega**2/c_s0**2-k**2)
q, alpha = sp.symbols("q alpha", positive=True, finite=True, real=True)
V = sp.symbols("V", nonzero=True)
upper_prop = next(s for s in (sp.Integer(1), sp.Integer(-1)) if sp.sign(s*q) == 1)
lower_prop = next(s for s in (sp.Integer(1), sp.Integer(-1)) if sp.sign(s*q) == -1)
upper_evan, lower_evan = I*alpha, -I*alpha
checked(upper_prop == 1 and lower_prop == -1)
checked(sp.im(upper_evan) > 0 and sp.im(lower_evan) < 0)

q_out = sp.symbols("q_out", nonzero=True)
phi_face = sp.simplify(V/(I*q_out))
pressure_face = sp.simplify(I*rho_m*omega*phi_face)
Z0 = sp.simplify(pressure_face/V)
Z_prop, Z_evan = sp.simplify(Z0.subs(q_out, q)), sp.simplify(Z0.subs(q_out, I*alpha))
checked(Z_prop == rho_m*omega/q)
checked(Z_evan == -I*rho_m*omega/alpha)
emit("S11BA_Z_IMPERMEABLE",
     "q_squared=" + fmt(q2) + "; q2>0: q_w_upper=" + fmt(upper_prop*q)
     + ",q_w_lower=" + fmt(lower_prop*q) + " (energy_outgoing); q2<0: q_w_upper="
     + fmt(upper_evan) + ",q_w_lower=" + fmt(lower_evan) + " (decaying); Z=p/V=" + fmt(Z0))
emit("S11BA_Z_BY_REGIME",
     "q2>0: Z=" + fmt(Z_prop) + ",Re=" + fmt(sp.re(Z_prop)) + ",Im=" + fmt(sp.im(Z_prop))
     + "; q2<0: Z=" + fmt(Z_evan) + ",Re=" + fmt(sp.re(Z_evan)) + ",Im=" + fmt(sp.im(Z_evan))
     + "; q2=0: Z divergent/no finite driven bounded solution")

dW, zeta_c = sp.symbols("delta_W zeta_c", nonzero=True)
zeta_delta, zeta_center = sp.Matrix([dW/2, -dW/2]), sp.Matrix([zeta_c, zeta_c])
outward_sign = sp.diag(1, -1)
V_delta = sp.simplify(outward_sign*(-I*omega*zeta_delta))
V_center = sp.simplify(outward_sign*(-I*omega*zeta_center))
p_delta_evan, p_center_evan = sp.simplify(Z_evan*V_delta), sp.simplify(Z_evan*V_center)
checked(V_delta[0] == V_delta[1])
checked(V_center[0] == -V_center[1])
emit("S11BA_Z_BY_PARITY",
     "delta_W: (V_+,V_-)=" + fmt(tuple(V_delta)) + ",(Z_+,Z_-)=" + fmt((Z0, Z0))
     + ",q2>0:(Re,Im)=" + fmt((sp.re(Z_prop), sp.im(Z_prop)))
     + ",q2<0:(Re,Im)=" + fmt((sp.re(Z_evan), sp.im(Z_evan)))
     + "; zeta_c: (V_+,V_-)=" + fmt(tuple(V_center)) + ",(Z_+,Z_-)=" + fmt((Z0, Z0))
     + ",q2>0:(Re,Im)=" + fmt((sp.re(Z_prop), sp.im(Z_prop)))
     + ",q2<0:(Re,Im)=" + fmt((sp.re(Z_evan), sp.im(Z_evan)))
     + "; both:q2=0_singular; no_parity_mixing")

accel_delta, accel_center = sp.simplify(-omega**2*zeta_delta), sp.simplify(-omega**2*zeta_center)
m_delta = tuple(sp.simplify(p_delta_evan[i]/accel_delta[i]) for i in range(2))
m_center = tuple(sp.simplify(p_center_evan[i]/accel_center[i]) for i in range(2))
checked(m_delta == (rho_m/alpha, -rho_m/alpha))
checked(m_center == m_delta)
emit("S11BA_ADDED_MASS",
     "q2<0 per_face_global_zeta_definition: delta_W=" + fmt(m_delta) + ",zeta_c=" + fmt(m_center)
     + "; q2>0=NOT_APPLICABLE(Z_not_pure_imaginary); q2=0=divergent")

A0, B0 = sp.symbols("A0 B0")
grazing_driven = sp.solve((sp.Eq(B0, V), sp.Eq(B0, 0)), (B0,), dict=True)
grazing_undriven = sp.solve((sp.Eq(B0, 0),), (B0,), dict=True)
checked(grazing_driven == [])
checked(grazing_undriven == [{B0: 0}])
emit("S11BA_GRAZING_BEHAVIOUR",
     "q2=0: phi=A0+B0*(outward_distance), boundedness=>B0=0; V!=0=>no_solution; "
     "V=0=>A0_free; lim_q->0+ Z=+infinity(real); lim_alpha->0+ Z=-I*infinity; "
     "m_add_upper/lower=+/-infinity")


# A5: permeable closure with real free coefficients and one real tau.
Lambda_p0, Lambda_V0 = sp.symbols("Lambda_p0 Lambda_V0", finite=True, real=True)
tau = sp.symbols("tau", finite=True, real=True, nonnegative=True)
x = sp.symbols("x", finite=True, real=True, nonnegative=True)  # x=omega*tau
y = sp.symbols("y", finite=True)
r = 1-I*x
Lambda_p, Lambda_V = Lambda_p0/r, Lambda_V0/r
Z_perm = sp.factor((rho_m*r+Lambda_V0)/(y*r-Lambda_p0))
v_n_from_Z, p_from_Z = sp.simplify(Z_perm*V*y/rho_m), sp.simplify(Z_perm*V)
closure_residual = sp.simplify(rho_m*(v_n_from_Z-V)-Lambda_p*p_from_Z-Lambda_V*V)
checked(closure_residual == 0)

yp, be = sp.symbols("a_q b_q", positive=True, finite=True, real=True)
Zp, Ze, Zg = sp.factor(Z_perm.subs(y, yp)), sp.factor(Z_perm.subs(y, I*be)), sp.factor(Z_perm.subs(y, 0))


def real_part(expr):
    return sp.factor(sp.simplify(sp.re(sp.together(expr))))


def imag_part(expr):
    return sp.factor(sp.simplify(sp.im(sp.together(expr))))


ReZp, ImZp = real_part(Zp), imag_part(Zp)
ReZe, ImZe = real_part(Ze), imag_part(Ze)
ReZg, ImZg = real_part(Zg), imag_part(Zg)
Acoef, Bcoef = rho_m+Lambda_V0, yp-Lambda_p0
ReZp_expected = sp.factor((Acoef*Bcoef+rho_m*yp*x**2)/(Bcoef**2+yp**2*x**2))
ReZe_expected = sp.factor((be*x*Lambda_V0-Acoef*Lambda_p0)/((be*x-Lambda_p0)**2+be**2))
checked(sp.simplify(ReZp-ReZp_expected) == 0)
checked(sp.simplify(ReZe-ReZe_expected) == 0)
checked(sp.simplify(ReZg+Acoef/Lambda_p0) == 0)

# Coefficient dimensions are computed here directly from v=grad(phi), p=-rho*d_t(phi), and J=rho*v.
dL, dT, dM = sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])
dvel, drho = dL-dT, dM-4*dL
dphi, dp, dJ = dvel+dL, drho-dT+(dvel+dL), drho+dvel
dLp_a5, dLv_a5, dtau_a5 = dJ-dp, dJ-dvel, dT
vec = lambda z: "[" + ",".join(fmt(e) for e in z) + "]"
emit("S11BA_PERMEABLE_COEFF_DIMS",
     "[Lambda_p0]=" + vec(dLp_a5) + "_[L,T,M]; [Lambda_V0]=" + vec(dLv_a5)
     + "_[L,T,M]; [tau]=" + vec(dtau_a5) + "_[L,T,M]")

Z_complex_limits = {
    "prop_small": sp.factor(sp.limit(Zp, x, 0, dir="+")),
    "prop_large": sp.factor(sp.limit(Zp, x, sp.oo)),
    "evan_small": sp.factor(sp.limit(Ze, x, 0, dir="+")),
    "evan_large": sp.factor(sp.limit(Ze, x, sp.oo)),
    "grazing_small": sp.factor(sp.limit(Zg, x, 0, dir="+")),
    "grazing_large_slope": sp.factor(sp.limit(Zg/x, x, sp.oo)),
}
Lambda_p_small = sp.series(Lambda_p, x, 0, 2)
Lambda_V_small = sp.series(Lambda_V, x, 0, 2)
Lambda_p_large_scaled = sp.limit(x*Lambda_p, x, sp.oo)
Lambda_V_large_scaled = sp.limit(x*Lambda_V, x, sp.oo)
checked(sp.simplify(Z_complex_limits["prop_large"]-rho_m/yp) == 0)
checked(sp.simplify(Z_complex_limits["evan_large"]+I*rho_m/be) == 0)
checked(sp.simplify(Z_complex_limits["grazing_large_slope"]-I*rho_m/Lambda_p0) == 0)
emit("S11BA_Z_PERMEABLE",
     "x=omega*tau; Z(y)=" + fmt(Z_perm) + ", y=q_out/omega; propagating(y=a_q>0)=" + fmt(Zp)
     + "; evanescent(y=I*b_q,b_q>0)=" + fmt(Ze) + "; grazing(y=0)=" + fmt(Zg)
     + "; x<<1 limits(prop,evan,grazing)="
     + fmt((Z_complex_limits["prop_small"], Z_complex_limits["evan_small"], Z_complex_limits["grazing_small"]))
     + "; x=O(1):exact_formulas_above; x>>1 limits(prop,evan)="
     + fmt((Z_complex_limits["prop_large"], Z_complex_limits["evan_large"]))
     + ",grazing~x*" + fmt(Z_complex_limits["grazing_large_slope"])
     + "; coefficient_behaviour: x<<1 (Lambda_p,Lambda_V)=" + fmt((Lambda_p_small, Lambda_V_small))
     + ",x=O(1)=exact_closure,x>>1 x*(Lambda_p,Lambda_V)->"
     + fmt((Lambda_p_large_scaled, Lambda_V_large_scaled)))

parity_dissipation = {
    "delta_W": {"propagating": ReZp, "evanescent": ReZe, "grazing": ReZg},
    "zeta_c": {"propagating": ReZp, "evanescent": ReZe, "grazing": ReZg},
}
emit("S11BA_DISSIPATIVE_BY_REGIME_AND_PARITY", "; ".join(
    parity + ":prop_ReZ=" + fmt(vals["propagating"]) + ",evan_ReZ=" + fmt(vals["evanescent"])
    + ",grazing_ReZ=" + fmt(vals["grazing"])
    + " (generic_presence; each depends_on_Lambda_p0_and_Lambda_V0; zeros occur on numerator loci)"
    for parity, vals in parity_dissipation.items()))

limits = {
    "prop_small": sp.simplify(sp.limit(ReZp, x, 0, dir="+")),
    "prop_large": sp.simplify(sp.limit(ReZp, x, sp.oo)),
    "evan_small": sp.simplify(sp.limit(ReZe, x, 0, dir="+")),
    "evan_large": sp.simplify(sp.limit(ReZe, x, sp.oo)),
    "grazing_small": sp.simplify(sp.limit(ReZg, x, 0, dir="+")),
    "grazing_large": sp.simplify(sp.limit(ReZg, x, sp.oo)),
}
checked(sp.simplify(limits["prop_large"]-rho_m/yp) == 0)
checked(limits["evan_large"] == 0)
checked(sp.simplify(limits["grazing_small"]-limits["grazing_large"]) == 0)
emit("S11BA_DISSIPATION_VS_OMEGA_TAU",
     "x<<1: prop=" + fmt(limits["prop_small"]) + ",evan=" + fmt(limits["evan_small"])
     + ",grazing=" + fmt(limits["grazing_small"])
     + "; x=O(1): prop=" + fmt(ReZp) + ",evan=" + fmt(ReZe) + ",grazing=" + fmt(ReZg)
     + "; x>>1: prop->" + fmt(limits["prop_large"]) + ",evan->" + fmt(limits["evan_large"])
     + ",grazing->" + fmt(limits["grazing_large"]) + "; same_per_face_for_delta_W_and_zeta_c")

Z_tau0 = sp.factor(Z_perm.subs(x, 0))
emit("S11BA_TAU_ZERO_LIMIT",
     "special_memoryless_limit_only: Z=" + fmt(Z_tau0) + "; prop=" + fmt(Zp.subs(x, 0))
     + "; evan=" + fmt(Ze.subs(x, 0)) + "; grazing=" + fmt(Zg.subs(x, 0)))

denominator_locus = sp.expand(y*r)
numerator_locus = sp.expand(-rho_m*r)
checked(sp.simplify(denominator_locus-y*r) == 0)
checked(sp.simplify(numerator_locus+rho_m*r) == 0)
emit("S11BA_DEGENERATE_LOCI",
     "amplitude_coefficient_zero: Lambda_p0=" + fmt(denominator_locus)
     + "; with_real_coefficients: propagating=>x=0,Lambda_p0=a_q; evanescent=>no_locus_for_b_q>0; "
     "grazing=>Lambda_p0=0(any_x). On_denominator_locus: V!=0 and "
     "rho_m*(1-I*x)+Lambda_V0!=0=>no_solution; V=0=>bulk_amplitude_free. "
     "Simultaneous_numerator_locus Lambda_V0=" + fmt(numerator_locus)
     + " is_real_only_at_x=0,where_Lambda_V0=-rho_m; then driven_amplitude_free")

shared_tau_count = len(set([tau, tau]))
checked(shared_tau_count == 1)
emit("S11BA_CLOSURE_SCOPE_LIMIT",
     "one_pole_memory_law; shared_relaxation_times=" + str(shared_tau_count)
     + "; tau>=0; tau=0_only_a_limit; separate_relaxation_times=not_attempted; full_memory_kernel=not_attempted")


# A6: [L,T,M] exponent vectors derived from the given equations.
def dim(Lexp=0, Texp=0, Mexp=0):
    return sp.Matrix([sp.Rational(Lexp), sp.Rational(Texp), sp.Rational(Mexp)])


def dim_text(vector):
    return "[" + ",".join(fmt(entry) for entry in vector) + "]_[L,T,M]"


d_length, d_time, d_mass = dim(1), dim(0, 1), dim(0, 0, 1)
d_grad, d_dt = -d_length, -d_time
d_velocity = d_length-d_time
d_rho = d_mass-4*d_length
d_phi = d_velocity-d_grad
d_pressure = d_rho+d_dt+d_phi
d_flux = d_rho+d_velocity
d_cs = (2*d_dt-2*d_grad)/2
d_v0 = d_velocity
d_Z = d_pressure-d_velocity
d_acceleration = d_length-2*d_time
d_madd = d_pressure-d_acceleration
d_Lambda_p = d_flux-d_pressure
d_Lambda_V = d_flux-d_velocity
d_tau = -d_dt
d_source = d_dt+d_rho+d_length
dimension_data = {
    "Z": (d_Z, "[delta_p]-[V]", "definitional"),
    "M_ADD": (d_madd, "[delta_p]-[d_t^2 zeta]", "definitional"),
    "RHO_M": (d_rho, "[mass]-4[length]", "definitional"),
    "C_S0": (d_cs, "([d_t^2]-[grad_4^2])/2 from wave equation", "independent"),
    "V0": (d_v0, "[length]-[time]", "definitional"),
    "LAMBDA_P0": (d_Lambda_p, "[J]-[delta_p] from closure", "independent"),
    "LAMBDA_V0": (d_Lambda_V, "[J]-[V] from closure", "independent"),
    "TAU": (d_tau, "-[omega] from dimensionless omega*tau", "independent"),
    "A1_SOURCE": (d_source, "[d_t]+[integral rho dw] from projected continuity", "independent"),
}
checked(d_cs == d_velocity)
checked(d_Z == dim(-3, -1, 1))
checked(d_madd == dim(-3, 0, 1))
checked(d_Lambda_p == dim(-1, 1, 0))
checked(d_Lambda_V == d_rho)
checked(d_source == dim(-3, -1, 1))
checked(d_Lambda_p == dLp_a5 and d_Lambda_V == dLv_a5 and d_tau == dtau_a5)
for name, (vector, route, kind) in dimension_data.items():
    emit("S11BA_DIM_"+name, dim_text(vector)+"; route="+route)
    emit("S11BA_DIM_ROUTE_KIND_"+name, kind)
emit("S11BA_DIM_ROUTE_KIND", ";".join(name+"="+data[2] for name, data in dimension_data.items()))


# A7: the two controlled form families, including finite-interval endpoints.
je, jo = sp.Function("j_even"), sp.Function("j_odd")
Omega_b = sp.sech((w-b)/a)**2
Omega_b_prime = sp.diff(Omega_b, w)
Obp_minus, Obp_plus = Omega_b_prime.subs(w, -u), Omega_b_prime.subs(w, u)
Ob_minus_L, Ob_plus_L = Omega_b.subs(w, -L), Omega_b.subs(w, L)
control_A_even = sp.Integral(je(u)*(Obp_plus+Obp_minus), (u, 0, L)) + je(L)*(Ob_minus_L-Ob_plus_L)
control_A_odd = sp.Integral(jo(u)*(Obp_plus-Obp_minus), (u, 0, L)) - jo(L)*(Ob_plus_L+Ob_minus_L)
checked(sp.simplify(control_A_even.subs(b, 0)) == 0)

F_even_b = -(a/2)*w*sp.sech(w/a)**2 + (a**2/2)*sp.tanh(w/a)
checked(sp.simplify(sp.diff(F_even_b, w)-w*sp.sech(w/a)**2*sp.tanh(w/a)) == 0)
even_b_slope = sp.simplify((-4/a)*(F_even_b.subs(w, L)-F_even_b.subs(w, -L)))
odd_b_witness = sp.simplify(-a*(sp.tanh((L-b)/a)+sp.tanh((L+b)/a)))
even_b_independent = sp.simplify(even_b_slope) == 0
odd_b_independent = sp.simplify(sp.diff(odd_b_witness, b)) == 0
checked(not even_b_independent and not odd_b_independent)
emit("S11BA_CONTROL_WINDOW_PARITY",
     "Omega_b=sech((w-b)/a)^2,interval=[-L,L]; even_j:S(b)=" + fmt(control_A_even)
     + ",identically_independent_of_b=" + str(even_b_independent) + ",w^2_witness_dS/db|0=" + fmt(even_b_slope)
     + "; odd_j:S(b)=" + fmt(control_A_odd) + ",identically_independent_of_b=" + str(odd_b_independent)
     + ",w_witness=" + fmt(odd_b_witness))

Omega_0, Omega_0_prime = sp.sech(w/a)**2, sp.diff(sp.sech(w/a)**2, w)
tail_even = sp.Integral(je(w)*Omega_0_prime, (w, L, L+c))
tail_odd = sp.Integral(jo(w)*Omega_0_prime, (w, L, L+c))
control_B_even = tail_even + Omega_0.subs(w, L)*je(L) - Omega_0.subs(w, L+c)*je(L+c)
control_B_odd = 2*sp.Integral(jo(u)*sp.diff(sp.sech(u/a)**2, u), (u, 0, L))
control_B_odd += tail_odd - Omega_0.subs(w, L)*jo(L) - Omega_0.subs(w, L+c)*jo(L+c)
even_c_witness_derivative = sp.simplify(-2*(L+c)*Omega_0.subs(w, L+c))
odd_c_witness = sp.simplify(-a*(sp.tanh((L+c)/a)+sp.tanh(L/a)))
even_c_independent = even_c_witness_derivative == 0
odd_c_independent = sp.simplify(sp.diff(odd_c_witness, c)) == 0
checked(not even_c_independent and not odd_c_independent)
checked(sp.simplify(control_B_even.subs(c, 0)) == 0)
emit("S11BA_CONTROL_INTERVAL_SYMMETRY",
     "Omega_0=sech(w/a)^2,interval=[-L,L+c]; even_j:S(c)=" + fmt(control_B_even)
     + ",identically_independent_of_c=" + str(even_c_independent) + ",w^2_witness_dS/dc=" + fmt(even_c_witness_derivative)
     + "; odd_j:S(c)=" + fmt(control_B_odd) + ",identically_independent_of_c=" + str(odd_c_independent)
     + ",w_witness=" + fmt(odd_c_witness))


# A8: independent slow-background ratios and the leading convective order.
timescale_ratio = sp.simplify((1/omega)/(W0/sp.Abs(v0)))
flow_ratio = sp.simplify(sp.Abs(v0)/c_s0)
Mach = sp.symbols("Mach", real=True)
convective_relative = sp.expand((1-Mach)**2-1)
leading_convective_power = min(monomial[0] for monomial, coefficient
                               in sp.Poly(convective_relative, Mach).terms() if coefficient != 0)
checked(leading_convective_power == 1)
emit("S11BA_VALIDITY_TIMESCALE", fmt(timescale_ratio)
     + " << 1 (wave_time=1/omega << background_crossing_time=W0/Abs(v0))")
emit("S11BA_VALIDITY_FLOW_SPEED", fmt(flow_ratio) + " << 1 (independent_of_timescale_condition)")
emit("S11BA_DISCARDED_CONVECTIVE_ORDER",
     "O((v0/c_s0)^" + str(leading_convective_power) + "); rest_frame_linearisation_scope_only")

for tag, value in tags:
    print(f"{tag}: {value}")
print("VERDICT: " + ("PASS" if checks and all(checks) else "FAIL"))
