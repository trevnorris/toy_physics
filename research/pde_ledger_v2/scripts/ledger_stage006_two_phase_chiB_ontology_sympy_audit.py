#!/usr/bin/env python3
"""Ledger stage006 SymPy audit: two-phase chi_B ontology.

Print-only, standalone, no arguments, no file output.  This is a fresh-authored
ledger stage: the chi_B action is postulated, while the closure identities,
wall admission, recovery reduction, and theta-as-phi no-go predicates are
checked exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0

EXPECTED_POST_D16_DRIFT_TOKEN = "DRIFT(5)"
PRE_D16_DRIFT_MEMBERS = (
    "chi_B",
    "a_B",
    "kappa_B",
    "alpha_aniso",
    "Gamma_B",
    "gating structure",
)
RETIRED_D16_DRIFT_MEMBERS = frozenset({"alpha_aniso"})


class AuditFailure(AssertionError):
    pass


def banner(title: str) -> None:
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def subbanner(title: str) -> None:
    print("")
    print(title)
    print("-" * len(title))


def compact(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.cancel(sp.simplify(expr))))


def assert_no_float(name: str, expr: Any) -> None:
    clean = sp.sympify(expr)
    floats = clean.atoms(sp.Float)
    if floats:
        raise AuditFailure(f"{name}: Float atom(s) found in exact audit expression: {floats}")


def _record_pass(message: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(message)


def _record_fail(message: str) -> None:
    global FAIL_COUNT
    FAIL_COUNT += 1
    print(message)


def expect_zero(name: str, residual: sp.Expr) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean == 0:
        _record_pass(f"PASS  {name}")
        return
    _record_fail(f"FAIL  {name}: residual = {compact(clean)}")
    raise AuditFailure(f"{name} residual was not zero")


def expect_bool(name: str, condition: bool) -> None:
    expect_zero(name, sp.Integer(0) if condition else sp.Integer(1))


def expect_nonzero(name: str, residual: sp.Expr) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {compact(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} produced required FAIL (residual = {compact(clean)})")
        return
    _record_fail(f"FAIL  {name}: required mutation/ablation did not fire")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def q(value: int | str | sp.Rational) -> sp.Rational:
    return sp.Rational(value)


@dataclass(frozen=True)
class Dim:
    """Exact exponent vector for base dimensions in {L, T, M} order."""

    l: sp.Rational | int = 0
    t: sp.Rational | int = 0
    m: sp.Rational | int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "l", q(self.l))
        object.__setattr__(self, "t", q(self.t))
        object.__setattr__(self, "m", q(self.m))

    def __mul__(self, other: "Dim") -> "Dim":
        return Dim(self.l + other.l, self.t + other.t, self.m + other.m)

    def __truediv__(self, other: "Dim") -> "Dim":
        return Dim(self.l - other.l, self.t - other.t, self.m - other.m)

    def __pow__(self, power: int | sp.Rational) -> "Dim":
        p = q(power)
        return Dim(self.l * p, self.t * p, self.m * p)

    def components(self) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
        return (self.l, self.t, self.m)

    def __str__(self) -> str:
        pieces: list[str] = []
        for label, exponent in (("L", self.l), ("T", self.t), ("M", self.m)):
            if exponent == 0:
                continue
            pieces.append(label if exponent == 1 else f"{label}^{exponent}")
        return "1" if not pieces else " ".join(pieces)


DIMENSIONLESS = Dim()
LENGTH = Dim(1, 0, 0)
TIME = Dim(0, 1, 0)
MASS = Dim(0, 0, 1)
VELOCITY = LENGTH / TIME
ACTION = MASS * (LENGTH**2) / TIME
ENERGY = ACTION / TIME
F_DENSITY_4D = ENERGY / (LENGTH**4)
BRANE_LAGRANGIAN_DENSITY = ENERGY / (LENGTH**3)
RHO4 = LENGTH**-4
M_GNLS = MASS
K_EOS = F_DENSITY_4D / (RHO4**5)
HBAR = ACTION


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(sum((have - want) ** 2 for have, want in zip(actual.components(), expected.components())))


def expect_dim(name: str, actual: Dim, expected: Dim) -> None:
    expect_zero(name, dim_residual(actual, expected))


def homogeneity_residual(terms: dict[str, Dim]) -> sp.Expr:
    if not terms:
        raise AuditFailure("homogeneity check requires at least one term")
    dims = list(terms.values())
    reference = dims[0]
    return sp.simplify(sum(dim_residual(actual, reference) for actual in dims[1:]))


@dataclass(frozen=True)
class A1Term:
    dim: Dim
    symbols: frozenset[str]


def operative_a1_absence_residual(terms: dict[str, A1Term]) -> sp.Integer:
    live_symbols = set().union(*(term.symbols for term in terms.values()))
    retired_is_present = "alpha_aniso" in live_symbols or any("alpha_aniso" in label for label in terms)
    return sp.Integer(1 if retired_is_present else 0)


@dataclass(frozen=True)
class Symbols:
    chi: sp.Symbol
    c: sp.Symbol
    n: sp.Symbol
    gamma_b: sp.Symbol
    div_j_chi: sp.Symbol
    div_nu: sp.Symbol
    n_t: sp.Symbol
    chi_t: sp.Symbol
    u_grad_chi: sp.Symbol
    boundary_j: sp.Symbol
    integral_wp_j: sp.Symbol
    boundary_jchi: sp.Symbol
    integral_wp_jchi: sp.Symbol
    convert: sp.Symbol
    k: sp.Symbol
    rho_br: sp.Symbol
    mu_R: sp.Symbol
    J: sp.Symbol
    rho_B0: sp.Symbol
    chi_c: sp.Symbol
    kappa_phase: sp.Symbol
    K_theta: sp.Symbol
    B_eff: sp.Symbol
    m_theta2: sp.Symbol
    C_J: sp.Symbol
    u_L: sp.Symbol
    dot_u_L: sp.Symbol
    theta: sp.Symbol
    dot_theta: sp.Symbol
    p_u: sp.Symbol
    pi_theta: sp.Symbol
    kappa4: sp.Symbol


def make_symbols() -> Symbols:
    return Symbols(
        chi=sp.Symbol("chi_B"),
        c=sp.Symbol("c", nonzero=True),
        n=sp.Symbol("n", nonzero=True),
        gamma_b=sp.Symbol("Gamma_B", nonzero=True),
        div_j_chi=sp.Symbol("div_J_chi"),
        div_nu=sp.Symbol("div_nu"),
        n_t=sp.Symbol("partial_t_n"),
        chi_t=sp.Symbol("partial_t_chi_B"),
        u_grad_chi=sp.Symbol("u_dot_grad_chi_B"),
        boundary_j=sp.Symbol("B_j"),
        integral_wp_j=sp.Symbol("I_Wp_j"),
        boundary_jchi=sp.Symbol("B_Jchi"),
        integral_wp_jchi=sp.Symbol("I_Wp_Jchi"),
        convert=sp.Symbol("I_W_n_Gamma_B", nonzero=True),
        k=sp.Symbol("k", positive=True),
        rho_br=sp.Symbol("rho_br", positive=True),
        mu_R=sp.Symbol("mu_R", positive=True),
        J=sp.Symbol("J", positive=True),
        rho_B0=sp.Symbol("rho_B0", positive=True),
        chi_c=sp.Symbol("chi_c", positive=True),
        kappa_phase=sp.Symbol("kappa_phase", positive=True),
        K_theta=sp.Symbol("K_theta"),
        B_eff=sp.Symbol("B_eff"),
        m_theta2=sp.Symbol("m_theta2", positive=True),
        C_J=sp.Symbol("C_J"),
        u_L=sp.Symbol("u_L"),
        dot_u_L=sp.Symbol("dot_u_L"),
        theta=sp.Symbol("theta"),
        dot_theta=sp.Symbol("dot_theta"),
        p_u=sp.Symbol("p_u"),
        pi_theta=sp.Symbol("pi_theta"),
        kappa4=sp.Symbol("kappa4", positive=True),
    )


def print_pin_block() -> None:
    subbanner("Pinned postulated modeling choices P1-P13")
    pins = [
        ("P1", "POSTULATED", "n(X,t)=total conserved 4D constituent number density [n]=L^-4; KINETIC_MASS_FACTOR_PINNED: kinetic term is 1/2*m_GNLS*n*|u|^2."),
        ("P2", "POSTULATED", "chi_B dimensionless in [0,1]; chi_B=1 brane-ordered, chi_B=0 bulk-like; n_B=chi_B*n."),
        ("P3", "POSTULATED", "f_B(chi_B)=a_B*chi_B^2*(1-chi_B)^2, a_B>0, minima {0,1}; new double-well input because cited U(rho) is single-well."),
        ("P4", "POSTULATED", "SAME_DENSITY_DEGENERACY_POSTULATED: f_B is n-independent and f_B(n,0)=f_B(n,1)=0."),
        ("P5", "POSTULATED", "gradient/interface term (kappa_B/2)*|grad_4 chi_B|^2 with [kappa_B]=M T^-2."),
        ("P6", "POSTULATED", "shear gate chi_B*f_shear with displacement u_d distinct from velocity u; brane mu_R projection is dim-consistent only here."),
        ("P7", "POSTULATED", "chi_B is THE wall order parameter: a postulated independent scalar field, NOT currently built as |P_parallel|^2. Decision 16 retires alpha_aniso and the carried P field. FUTURE_GATE_CHI_B_EQ_ABS_P_PARALLEL_SQ remains named high-risk/Part-VII-adjacent and requires a NEW T0 freeze; obsolete as a carried route, not foreclosed."),
        ("P8", "POSTULATED ADJUNCT", "D_t chi_B=-M_chi*mu_chi+Gamma_B; HANDOFF_P_ORDER_N_PLACEMENT_CORRECTED: P_order=int mu_chi*D_t chi_B d4X."),
        ("P9", "POSTULATED DEFAULT", "J_chi=0 default; J_chi!=0 deferred with dim row only."),
        ("P10", "CONVENTION", "recovery target is frozen old-ledger S_leak with j^w=n*u^w and unit-normalized W, [W]=L^-1."),
        ("P11", "POSTULATED GLOBAL", "global returns R_0=-M_0 and R_1=-D_1 printed here, NOT locally asserted."),
        ("P12", "POSTULATED ONTOLOGY", "throat=phase-conversion admittance/outlet driven by stress/mu gradients, NOT pairwise suction; Gate-L delta w=u_w DEFERRED."),
        ("P13", "CONVENTION", "free-energy +1/2*kappa_phase*(grad theta)^2 maps to Lagrangian K_theta=-kappa_phase; Maxwell needs K_theta=C_J^2/rho_br>0."),
    ]
    for label, provenance, text in pins:
        print(f"  {label}  {provenance}: {text}")


def run_consumed_citation_checks(s: Symbols) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    subbanner("Consumed citations and exact-value integrity checks")
    print("  CITED ledger_stage004 (I-1): {L,T,M}; [n]=[rho0]=L^-4; [m_GNLS]=M; [K]=M L^18 T^-2; [hbar]=M L^2 T^-1; U(rho)=(K/4)*rho^5 single-well.")
    print("  CITED ledger_stage005 (I-2): c_s0^2=5*K*rho0^4/m_GNLS.")
    print("  CITED ledger_stage003 (III): c_gamma^2=mu_R/rho_br; C_J=-J*rho_B0; B_eff=rho_B0^2/chi_c; second-class classification rule.")

    K, rho0, m_gnls, c_s0_sq = sp.symbols("K rho0 m_GNLS c_s0_sq", positive=True)
    consumed_c_s0_sq = 5 * K * rho0**4 / m_gnls
    consumed_c_j = -s.J * s.rho_B0
    consumed_b_eff = s.rho_B0**2 / s.chi_c
    expect_zero("CITED exact-value check C_J+J*rho_B0=0", consumed_c_j + s.J * s.rho_B0)
    expect_zero("CITED exact-value check B_eff-rho_B0^2/chi_c=0", consumed_b_eff - s.rho_B0**2 / s.chi_c)
    expect_zero("CITED exact-value check c_s0^2-5*K*rho0^4/m_GNLS=0", consumed_c_s0_sq - 5 * K * rho0**4 / m_gnls)

    expect_dim("CITED stage004 [n]=L^-4", RHO4, LENGTH**-4)
    expect_dim("CITED stage004 [m_GNLS]=M", M_GNLS, MASS)
    expect_dim("CITED stage004 [K]=M L^18 T^-2", K_EOS, Dim(18, -2, 1))
    expect_dim("CITED stage004 [hbar]=M L^2 T^-1", HBAR, Dim(2, -1, 1))
    expect_dim("CITED stage005 [c_s0^2]=L^2 T^-2", K_EOS * (RHO4**4) / M_GNLS, VELOCITY**2)

    d_mu_R = BRANE_LAGRANGIAN_DENSITY
    d_rho_br = MASS * (LENGTH**-3)
    d_rho_B0 = d_rho_br
    d_J = (LENGTH**2) / TIME
    d_CJ = MASS * (LENGTH**-1) / TIME
    d_Ktheta = MASS * LENGTH / (TIME**2)
    d_chi_c = MASS * (LENGTH**-5) * (TIME**2)
    d_B = BRANE_LAGRANGIAN_DENSITY
    expect_dim("CITED stage003 [mu_R]=M L^-1 T^-2", d_mu_R, Dim(-1, -2, 1))
    expect_dim("CITED stage003 [rho_br]=M L^-3", d_rho_br, Dim(-3, 0, 1))
    expect_dim("CITED stage003 [rho_B0]=M L^-3", d_rho_B0, Dim(-3, 0, 1))
    expect_dim("CITED stage003 [J]=L^2 T^-1", d_J, Dim(2, -1, 0))
    expect_dim("CITED stage003 [C_J]=M L^-1 T^-1", d_CJ, Dim(-1, -1, 1))
    expect_dim("CITED stage003 [K_theta]=M L T^-2", d_Ktheta, Dim(1, -2, 1))
    expect_dim("CITED stage003 [chi_c]=M L^-5 T^2", d_chi_c, Dim(-5, 2, 1))
    expect_dim("CITED stage003 [B]=M L^-1 T^-2", d_B, Dim(-1, -2, 1))
    expect_dim("CITED stage003 [c_gamma^2]=L^2 T^-2", d_mu_R / d_rho_br, VELOCITY**2)
    return consumed_c_j, consumed_b_eff, consumed_c_s0_sq


def run_leg_a_dimensions() -> None:
    subbanner("Leg A1 dimensional audit: OPERATIVE post-Decision-16 action and adjunct rows")
    chi = DIMENSIONLESS
    n = RHO4
    u_velocity = VELOCITY
    u_displacement = LENGTH
    grad4 = LENGTH**-1
    kappa_b = MASS / (TIME**2)
    a_b = F_DENSITY_4D
    mu_R_4d = F_DENSITY_4D
    f_throat = F_DENSITY_4D
    f_mix = F_DENSITY_4D

    terms = {
        "POSTULATED P1 kinetic 1/2*m_GNLS*n*|u|^2": A1Term(M_GNLS * n * (u_velocity**2), frozenset({"m_GNLS", "n", "u"})),
        "CITED I-1 U(n)=(K/4)*n^5": A1Term(K_EOS * (n**5), frozenset({"K", "n"})),
        "POSTULATED P3 f_B=a_B*chi_B^2*(1-chi_B)^2": A1Term(a_b, frozenset({"a_B", "chi_B"})),
        "POSTULATED P5 (kappa_B/2)*|grad_4 chi_B|^2": A1Term(kappa_b * ((grad4 * chi) ** 2), frozenset({"kappa_B", "chi_B"})),
        "POSTULATED P6 chi_B*f_shear": A1Term(chi * mu_R_4d * ((grad4 * u_displacement) ** 2), frozenset({"chi_B", "mu_R^(4)", "u"})),
        "DEFERRED_PLACEHOLDER f_throat": A1Term(f_throat, frozenset({"f_throat"})),
        "DEFERRED_PLACEHOLDER f_mix": A1Term(f_mix, frozenset({"f_mix"})),
    }
    for name, term in terms.items():
        expect_dim(f"{name} has 4D free-energy-density dim M L^-2 T^-2", term.dim, F_DENSITY_4D)
    expect_zero("POSTULATED operative F integrand homogeneity", homogeneity_residual({name: term.dim for name, term in terms.items()}))
    expect_zero("Decision-16 operative A1 surface excludes retired symbol alpha_aniso", operative_a1_absence_residual(terms))

    print("  RETIRED-HISTORICAL (not in operative A1): alpha_aniso*chi_B*(P.w_hat)^2 had [alpha_aniso]=M L^-2 T^-2; retired with P by Decision 16, not by a dimensional defect.")
    expect_dim(
        "RETIRED-HISTORICAL P7 alpha_aniso*chi_B*(P.w_hat)^2 was dimensionally homogeneous",
        F_DENSITY_4D * chi * (DIMENSIONLESS**2),
        F_DENSITY_4D,
    )
    reinjected_terms = dict(terms)
    reinjected_terms["REINJECTED alpha_aniso*chi_B*(P.w_hat)^2"] = A1Term(
        F_DENSITY_4D * chi * (DIMENSIONLESS**2),
        frozenset({"alpha_aniso", "chi_B", "P", "w_hat"}),
    )
    expect_fail(
        "Decision-16 A1 tooth: re-inject alpha_aniso term trips operative retired-symbol absence",
        operative_a1_absence_residual(reinjected_terms),
    )
    expect_zero("Decision-16 baseline operative A1 surface remains alpha_aniso-free after copy mutation", operative_a1_absence_residual(terms))

    pde_terms = {
        "partial_t(chi_B*n)": chi * n / TIME,
        "div_4(chi_B*n*u)": chi * n * u_velocity / LENGTH,
        "div_4 J_chi": (LENGTH**-3 / TIME) / LENGTH,
        "n*Gamma_B": n / TIME,
    }
    for name, dim in pde_terms.items():
        expect_dim(f"POSTULATED balance row [{name}]=L^-4 T^-1", dim, LENGTH**-4 / TIME)
    expect_zero("POSTULATED balance PDE homogeneity implies [Gamma_B]=T^-1", homogeneity_residual(pde_terms))

    mu_chi = F_DENSITY_4D
    M_chi = (LENGTH**2) * TIME / MASS
    P_order = mu_chi * (TIME**-1) * (LENGTH**4)
    M_n = (LENGTH**-4) * TIME / MASS
    number_mu = ENERGY
    expect_dim("POSTULATED ADJUNCT [mu_chi]=M L^-2 T^-2", mu_chi, F_DENSITY_4D)
    expect_dim("POSTULATED ADJUNCT [M_chi]=L^2 T M^-1", M_chi, Dim(2, 1, -1))
    expect_dim("POSTULATED ADJUNCT P_order=int mu_chi*D_t chi_B d4X = M L^2 T^-3", P_order, ENERGY / TIME)
    expect_dim("POSTULATED P12 [M_n]=L^-4 T M^-1", M_n, Dim(-4, 1, -1))
    expect_dim("POSTULATED P12 J_repair=-M_n*grad(mu) gives L^-3 T^-1", M_n * number_mu / LENGTH, LENGTH**-3 / TIME)

    W = LENGTH**-1
    rho_B_projected = W * RHO4 * LENGTH
    S_leak = rho_B_projected / TIME
    mu_R_projected = F_DENSITY_4D * LENGTH
    expect_dim("CONVENTION P10 [W]=L^-1", W, LENGTH**-1)
    expect_dim("CONVENTION P10 projected [rho_B]=int W*chi_B*n dw = L^-4", rho_B_projected, LENGTH**-4)
    expect_dim("CONVENTION P10 [S_leak]=L^-4 T^-1", S_leak, LENGTH**-4 / TIME)
    expect_dim("POSTULATED/PENDING P6 int chi_B*mu_R^(4) dw has brane [mu_R]=M L^-1 T^-2", mu_R_projected, BRANE_LAGRANGIAN_DENSITY)


def run_leg_a_structural_closure(s: Symbols) -> None:
    subbanner("Leg A2 structural-closure identities via SymPy Functions (EARNED)")
    x, w, t = sp.symbols("x w t", real=True)
    n = sp.Function("n")(x, w, t)
    chi = sp.Function("chi_B")(x, w, t)
    ux = sp.Function("u_x")(x, w, t)
    uw = sp.Function("u_w")(x, w, t)
    jx = sp.Function("J_chi_x")(x, w, t)
    jw = sp.Function("J_chi_w")(x, w, t)
    gamma = sp.Function("Gamma_B")(x, w, t)

    total_balance = sp.diff(n, t) + sp.diff(n * ux, x) + sp.diff(n * uw, w)
    order_left = (
        sp.diff(chi * n, t)
        + sp.diff(chi * n * ux + jx, x)
        + sp.diff(chi * n * uw + jw, w)
    )
    order_residual = order_left - n * gamma
    div_j_chi = sp.diff(jx, x) + sp.diff(jw, w)
    advective_form = sp.diff(chi, t) + ux * sp.diff(chi, x) + uw * sp.diff(chi, w)
    advective_residual = sp.expand((order_residual - chi * total_balance).doit() / n) - (
        advective_form + div_j_chi / n - gamma
    )
    expect_zero(
        "EARNED advective form D_t chi_B=Gamma_B-(1/n)div J_chi",
        sp.factor(sp.expand(advective_residual.doit())),
    )

    disorder_left = (
        sp.diff((1 - chi) * n, t)
        + sp.diff((1 - chi) * n * ux - jx, x)
        + sp.diff((1 - chi) * n * uw - jw, w)
    )
    disorder_residual = disorder_left + n * gamma
    expect_zero(
        "EARNED disorder complement is total minus order residual",
        sp.expand((disorder_residual - (total_balance - order_residual)).doit()),
    )
    expect_zero(
        "EARNED order+disorder sum reproduces total conservation with Gamma_B cancelling",
        sp.expand((order_residual + disorder_residual - total_balance).doit()),
    )
    source_coeff = sp.diff(order_residual, gamma)
    expect_nonzero("EARNED order balance is genuinely sourced by n*Gamma_B", source_coeff)


def run_leg_a_wall_admission() -> None:
    subbanner("Leg A3 wall admission (EARNED relative to P3/P5)")
    w = sp.Symbol("w", real=True)
    chi_sym = sp.Symbol("chi_B")
    a_b, kappa_b, delta = sp.symbols("a_B kappa_B delta", positive=True)
    f_b = a_b * chi_sym**2 * (1 - chi_sym) ** 2
    profile = sp.Rational(1, 2) * (1 + sp.tanh(w / (2 * delta)))

    coeff = sp.factor(kappa_b * sp.diff(profile, w, 2) / (profile * (1 - profile) * (1 - 2 * profile)))
    fprime_coeff = sp.factor(sp.diff(f_b, chi_sym) / (chi_sym * (1 - chi_sym) * (1 - 2 * chi_sym)))
    delta_sq = sp.solve(sp.Eq(coeff, fprime_coeff), delta**2)[0]
    delta_derived = sp.sqrt(delta_sq)
    expect_zero("EARNED kink width solved from kappa_B*chi''=f_B'", delta_derived - sp.sqrt(kappa_b / (2 * a_b)))

    profile_derived = profile.subs(delta, delta_derived)
    residual = sp.simplify(kappa_b * sp.diff(profile_derived, w, 2) - sp.diff(f_b, chi_sym).subs(chi_sym, profile_derived))
    expect_zero("EARNED kink EL residual vanishes by substitution", residual)

    sigma = sp.integrate(kappa_b * sp.diff(profile_derived, w) ** 2, (w, -sp.oo, sp.oo))
    sigma = sp.simplify(sigma)
    expect_zero("EARNED sigma_wall exact closed form", sigma - sp.sqrt(2 * a_b * kappa_b) / 6)
    expect_dim("EARNED [sigma_wall]=M L^-1 T^-2", (F_DENSITY_4D * (MASS / TIME**2)) ** sp.Rational(1, 2) * LENGTH**0, BRANE_LAGRANGIAN_DENSITY)


def run_leg_b_projection_profiles() -> None:
    subbanner("Leg B1 projected two-source law via sign-sensitive IBP (EARNED)")
    w = sp.Symbol("w", real=True)
    lam, A, B = sp.symbols("lambda A B", positive=True)

    gaussian_W = sp.exp(-w**2 / lam**2) / (lam * sp.sqrt(sp.pi))
    gaussian_Q = A * w * sp.exp(-w**2 / lam**2)
    rational_W = lam / (sp.pi * (w**2 + lam**2))
    rational_Q = B * w / (w**2 + lam**2)

    for label, W_expr, Q_expr in (
        ("Gaussian profile family", gaussian_W, gaussian_Q),
        ("Rational Cauchy profile family", rational_W, rational_Q),
    ):
        boundary = sp.Integer(0)
        int_w_qprime = sp.integrate(W_expr * sp.diff(Q_expr, w), (w, -sp.oo, sp.oo))
        int_wprime_q = sp.integrate(sp.diff(W_expr, w) * Q_expr, (w, -sp.oo, sp.oo))
        s_flux = -boundary + int_wprime_q
        expect_zero(f"EARNED IBP identity {label}: -int W*Q' = S_flux", sp.simplify(-int_w_qprime - s_flux))
        expect_fail(f"Leg-B ablation {label}: flip W' sign breaks IBP identity", sp.simplify(-int_w_qprime + int_wprime_q))
        expect_nonzero(f"EARNED {label}: W' contribution is not accidentally zero", int_wprime_q)


def boundary_jump(expr: sp.Expr, w: sp.Symbol) -> sp.Expr:
    return sp.simplify(sp.limit(expr, w, sp.oo) - sp.limit(expr, w, -sp.oo))


def assemble_projected_two_source(
    *,
    W_expr: sp.Expr,
    n_expr: sp.Expr,
    u_w_expr: sp.Expr,
    chi_expr: sp.Expr,
    gamma_expr: sp.Expr,
    j_chi_w_expr: sp.Expr,
    w: sp.Symbol,
) -> dict[str, sp.Expr]:
    q_order = sp.expand(chi_expr * n_expr * u_w_expr + j_chi_w_expr)
    raw_projected_flux = -sp.integrate(W_expr * sp.diff(q_order, w), (w, -sp.oo, sp.oo))
    boundary = boundary_jump(W_expr * q_order, w)
    s_flux = -boundary + sp.integrate(sp.diff(W_expr, w) * q_order, (w, -sp.oo, sp.oo))
    s_convert = sp.integrate(W_expr * n_expr * gamma_expr, (w, -sp.oo, sp.oo))
    return {
        "q_order": q_order,
        "raw_projected_flux": sp.simplify(raw_projected_flux),
        "boundary": sp.simplify(boundary),
        "s_flux": sp.simplify(s_flux),
        "s_convert": sp.simplify(s_convert),
        "s_total": sp.simplify(s_flux + s_convert),
    }


def run_leg_b_recovery_and_anchor(s: Symbols) -> None:
    subbanner("Leg B2/B3 recovery reduction and frozen Gaussian anchor (EARNED)")
    w = sp.Symbol("w", real=True)
    lam, n0, u0 = sp.symbols("lambda n0 u0", positive=True)
    chi0, chi_amp, gamma_amp, jchi_amp = sp.symbols("chi0 chi_amp gamma_amp jchi_amp", real=True)

    W_general = sp.exp(-w**2 / lam**2) / (lam * sp.sqrt(sp.pi))
    n_general = n0 * sp.exp(-w**2 / lam**2)
    u_w_general = u0 * w
    chi_general = chi0 + chi_amp * sp.exp(-w**2 / lam**2)
    gamma_general = gamma_amp * sp.exp(-w**2 / lam**2)
    j_chi_general = jchi_amp * w * sp.exp(-w**2 / lam**2)

    assembled = assemble_projected_two_source(
        W_expr=W_general,
        n_expr=n_general,
        u_w_expr=u_w_general,
        chi_expr=chi_general,
        gamma_expr=gamma_general,
        j_chi_w_expr=j_chi_general,
        w=w,
    )
    expect_zero(
        "EARNED B2 general chi_B/Gamma_B/J_chi projected order-balance IBP identity",
        assembled["raw_projected_flux"] - assembled["s_flux"],
    )
    expect_nonzero("EARNED B2 general chi_B profile is live in S_flux", sp.diff(assembled["s_flux"], chi_amp))
    expect_nonzero("EARNED B2 general Gamma_B profile is live in S_convert", sp.diff(assembled["s_convert"], gamma_amp))
    expect_nonzero("EARNED B2 general J_chi^w profile is live in S_flux", sp.diff(assembled["s_flux"], jchi_amp))

    limit_rules = {chi0: 1, chi_amp: 0, gamma_amp: 0, jchi_amp: 0}
    s_two_source_limit = sp.simplify(assembled["s_total"].subs(limit_rules))

    W_frozen = sp.exp(-w**2 / lam**2) / (lam * sp.sqrt(sp.pi))
    j_frozen = n0 * u0 * w * sp.exp(-w**2 / lam**2)
    s_leak_frozen = sp.simplify(
        -boundary_jump(W_frozen * j_frozen, w)
        + sp.integrate(sp.diff(W_frozen, w) * j_frozen, (w, -sp.oo, sp.oo))
    )
    expect_zero("FROZEN stage_243 S_leak target equals independent closed form on B2 profile", s_leak_frozen + sp.sqrt(2) * n0 * u0 / 4)
    expect_zero("EARNED recovery reduction at chi_B=1,Gamma_B=0,J_chi=0 against frozen stage_243 S_leak", s_two_source_limit - s_leak_frozen)

    residual_chi_c = sp.simplify(assembled["s_total"].subs({chi0: s.c, chi_amp: 0, gamma_amp: 0, jchi_amp: 0}) - s_leak_frozen)
    expect_zero("Leg-B conditionality chi_B=c residual computes to (c-1)*S_leak", residual_chi_c - (s.c - 1) * s_leak_frozen)
    expect_nonzero("Leg-B conditionality chi_B=c!=1 leaves (c-1)*S_leak residual", residual_chi_c)
    residual_gamma = sp.simplify(assembled["s_total"].subs({chi0: 1, chi_amp: 0, jchi_amp: 0}) - s_leak_frozen)
    expect_zero("Leg-B conditionality Gamma_B residual computes to S_convert integral", residual_gamma - assembled["s_convert"].subs({chi0: 1, chi_amp: 0, jchi_amp: 0}))
    expect_nonzero("Leg-B conditionality Gamma_B!=0 leaves computed S_convert residual", residual_gamma)
    residual_jchi = sp.simplify(assembled["s_total"].subs({chi0: 1, chi_amp: 0, gamma_amp: 0}) - s_leak_frozen)
    expected_jchi_flux = sp.simplify(assembled["s_flux"].subs({chi0: 1, chi_amp: 0, gamma_amp: 0}) - s_two_source_limit)
    expect_zero("Leg-B conditionality J_chi^w residual computes to its flux terms", residual_jchi - expected_jchi_flux)
    expect_nonzero("Leg-B conditionality J_chi^w!=0 leaves computed flux residual", residual_jchi)
    expect_fail(
        "Leg-B ablation corrupt frozen S_leak transcription alone breaks B2 reduction",
        s_two_source_limit - (s_leak_frozen + sp.sqrt(2) * n0 * u0 / 8),
    )
    expect_fail(
        "Leg-B ablation corrupt general assembly alone breaks B2 reduction",
        (s_two_source_limit + sp.sqrt(2) * n0 * u0 / 8) - s_leak_frozen,
    )

    mu_w, rho0, E0 = sp.symbols("mu_w rho0 E0", positive=True)
    W_lambda = sp.exp(-w**2 / lam**2) / (lam * sp.sqrt(sp.pi))
    phi_lambda = (2 * w / (sp.sqrt(sp.pi) * lam**3)) * sp.exp(-w**2 / lam**2)
    E_w = -E0 * phi_lambda
    j_w = mu_w * rho0 * E_w
    derived = sp.simplify(sp.integrate(sp.diff(W_lambda, w) * j_w, (w, -sp.oo, sp.oo)))
    frozen = sp.sqrt(2) * mu_w * rho0 * E0 / (2 * sp.sqrt(sp.pi) * lam**3)
    expect_zero("EARNED stage_244 Gaussian one-mode S_leak direct integration", derived - frozen)
    corrupt_phi = (2 * w / (sp.sqrt(sp.pi) * lam**2)) * sp.exp(-w**2 / lam**2)
    corrupt_j = mu_w * rho0 * (-E0 * corrupt_phi)
    corrupt = sp.simplify(sp.integrate(sp.diff(W_lambda, w) * corrupt_j, (w, -sp.oo, sp.oo)))
    expect_fail("Leg-B ablation corrupt Gaussian phi lambda^2 vs lambda^3 mismatches stage_244", corrupt - frozen)


def poisson(expr_a: sp.Expr, expr_b: sp.Expr, pairs: tuple[tuple[sp.Symbol, sp.Symbol], ...]) -> sp.Expr:
    out = sp.Integer(0)
    for coord, mom in pairs:
        out += sp.diff(expr_a, coord) * sp.diff(expr_b, mom) - sp.diff(expr_a, mom) * sp.diff(expr_b, coord)
    return sp.simplify(out)


def token_from_discriminators(
    *,
    bracket_zero: bool,
    B_eff_zero: bool,
    m_theta_zero: bool,
    provenance_flags: dict[str, bool],
) -> tuple[str, str]:
    if not m_theta_zero:
        return "FAIL_SECOND_CLASS_NOT_MAXWELL", "NOT_MAXWELL"
    if bracket_zero and B_eff_zero and m_theta_zero:
        flag = "WITH_PROVENANCE" if all(provenance_flags.values()) else "BY_TUNING"
        return "C5_RESOLVED_MAXWELL_BY_TUNING", flag
    if bracket_zero and (not B_eff_zero or not m_theta_zero):
        return "FAIL_SECOND_CLASS_NOT_MAXWELL", "NOT_MAXWELL"
    if B_eff_zero:
        return "FAIL_C5_LONGITUDINAL_ZERO_MODE", "PROVENANCE_NO_GO"
    return "FAIL_CAUCHY_STRAY_LONGITUDINAL", "PROVENANCE_NO_GO"


def algebraic_zero(expr: sp.Expr) -> bool:
    return sp.simplify(expr) == 0


def run_leg_c_no_go(s: Symbols, consumed_c_j: sp.Expr, consumed_b_eff: sp.Expr) -> tuple[str, str, str, str]:
    subbanner("Leg C theta-as-Maxwell-phi no-go (OWNED predicates, CONSUMED classification rule)")
    K_provenance = -s.kappa_phase
    K_maxwell = consumed_c_j**2 / s.rho_br
    expect_bool("OWNED C1 stable kappa_phase>0 gives K_theta=-kappa_phase<0", K_provenance.is_negative is True)
    expect_bool("OWNED C1 Maxwell K_theta=C_J^2/rho_br>0", K_maxwell.is_positive is True)
    expect_bool("OWNED C1 sign-lock predicates are OPPOSITE", (K_provenance.is_negative is True) and (K_maxwell.is_positive is True))
    expect_dim("OWNED C1 [kappa_phase]=[K_theta]=M L T^-2", MASS * LENGTH / (TIME**2), Dim(1, -2, 1))

    c_j_general = s.C_J
    lagrangian = (
        sp.Rational(1, 2) * s.rho_br * s.dot_u_L**2
        - c_j_general * s.k * s.u_L * s.dot_theta
        + sp.Rational(1, 2) * s.K_theta * s.k**2 * s.theta**2
        - sp.Rational(1, 2) * s.m_theta2 * s.theta**2
        - sp.Rational(1, 2) * s.B_eff * s.k**2 * s.u_L**2
    )
    p_u = sp.diff(lagrangian, s.dot_u_L)
    pi_theta = sp.diff(lagrangian, s.dot_theta)
    expect_zero("OWNED C2 p_u=rho_br*dot_u_L from specified probe Lagrangian", p_u - s.rho_br * s.dot_u_L)
    expect_zero("OWNED C2 primary constraint Phi_1=pi_theta-J*k*rho_B0*u_L", pi_theta.subs(s.C_J, consumed_c_j) - s.J * s.k * s.rho_B0 * s.u_L)

    H = (
        s.p_u**2 / (2 * s.rho_br)
        - sp.Rational(1, 2) * s.K_theta * s.k**2 * s.theta**2
        + sp.Rational(1, 2) * s.m_theta2 * s.theta**2
        + sp.Rational(1, 2) * s.B_eff * s.k**2 * s.u_L**2
    )
    phi1 = s.pi_theta + s.C_J * s.k * s.u_L
    phi2 = poisson(phi1, H, ((s.u_L, s.p_u), (s.theta, s.pi_theta)))
    phi2_prov = sp.simplify(phi2.subs({s.C_J: consumed_c_j, s.K_theta: -s.kappa_phase, s.m_theta2: 0}))
    expect_zero(
        "OWNED C2 secondary constraint Phi_2=-k*(J*p_u*rho_B0+k*kappa_phase*rho_br*theta)/rho_br",
        phi2_prov + s.k * (s.J * s.p_u * s.rho_B0 + s.k * s.kappa_phase * s.rho_br * s.theta) / s.rho_br,
    )
    bracket_general = poisson(phi1, phi2, ((s.u_L, s.p_u), (s.theta, s.pi_theta)))
    expect_zero(
        "OWNED C2 bracket route includes m_theta^2 algebraically",
        bracket_general - (s.C_J**2 * s.k**2 / s.rho_br - s.K_theta * s.k**2 + s.m_theta2),
    )
    bracket_prov = sp.simplify(bracket_general.subs({s.C_J: consumed_c_j, s.K_theta: -s.kappa_phase, s.m_theta2: 0}))
    expected_bracket = s.k**2 * (s.J**2 * s.rho_B0**2 + s.kappa_phase * s.rho_br) / s.rho_br
    expect_zero("OWNED C2 bracket {Phi_1,Phi_2}=k^2*(J^2*rho_B0^2+kappa_phase*rho_br)/rho_br", bracket_prov - expected_bracket)
    expect_nonzero("OWNED C2 provenance branch bracket is nonzero", bracket_prov)

    kappa_free = sp.Symbol("kappa_phase_free", real=True)
    bracket_free = bracket_prov.subs(s.kappa_phase, kappa_free)
    kappa_locus = sp.solve(sp.Eq(bracket_free, 0), kappa_free)[0]
    Ktheta_locus = sp.simplify(-kappa_locus)
    expect_zero("OWNED C2 zero locus solved: K_theta=+J^2*rho_B0^2/rho_br", Ktheta_locus - s.J**2 * s.rho_B0**2 / s.rho_br)

    bracket_tuned = sp.simplify(bracket_general.subs({s.C_J: consumed_c_j, s.K_theta: Ktheta_locus, s.m_theta2: 0}))
    expect_zero("OWNED C2 tuned Maxwell fixture bracket computes to zero", bracket_tuned)

    provenance_flags = {
        "K_theta_forced_by_frozen_defs": False,
        "B_eff_removed_by_frozen_defs": False,
        "m_theta_zero_forced": False,
    }
    provenance_branch = {s.C_J: consumed_c_j, s.K_theta: -s.kappa_phase, s.B_eff: consumed_b_eff, s.m_theta2: 0}
    zero_mode_branch = {s.C_J: consumed_c_j, s.K_theta: -s.kappa_phase, s.B_eff: 0, s.m_theta2: 0}
    tuned_branch = {s.C_J: consumed_c_j, s.K_theta: Ktheta_locus, s.B_eff: 0, s.m_theta2: 0}
    provenance_token, provenance_flag = token_from_discriminators(
        bracket_zero=algebraic_zero(bracket_general.subs(provenance_branch)),
        B_eff_zero=algebraic_zero(provenance_branch[s.B_eff]),
        m_theta_zero=algebraic_zero(provenance_branch[s.m_theta2]),
        provenance_flags=provenance_flags,
    )
    zero_mode_token, _ = token_from_discriminators(
        bracket_zero=algebraic_zero(bracket_general.subs(zero_mode_branch)),
        B_eff_zero=algebraic_zero(zero_mode_branch[s.B_eff]),
        m_theta_zero=algebraic_zero(zero_mode_branch[s.m_theta2]),
        provenance_flags=provenance_flags,
    )
    tuned_token, tuned_flag = token_from_discriminators(
        bracket_zero=algebraic_zero(bracket_general.subs(tuned_branch)),
        B_eff_zero=algebraic_zero(tuned_branch[s.B_eff]),
        m_theta_zero=algebraic_zero(tuned_branch[s.m_theta2]),
        provenance_flags=provenance_flags,
    )
    expect_zero("CONSUMED classification emits FAIL_CAUCHY_STRAY_LONGITUDINAL on provenance branch", 0 if provenance_token == "FAIL_CAUCHY_STRAY_LONGITUDINAL" else 1)
    expect_zero("CONSUMED classification emits FAIL_C5_LONGITUDINAL_ZERO_MODE when B_eff=0", 0 if zero_mode_token == "FAIL_C5_LONGITUDINAL_ZERO_MODE" else 1)
    expect_zero("CONSUMED classification emits C5_RESOLVED_MAXWELL_BY_TUNING on tuned fixture", 0 if tuned_token == "C5_RESOLVED_MAXWELL_BY_TUNING" else 1)
    expect_zero("OWNED provenance flag is BY_TUNING NOT WITH_PROVENANCE", 0 if tuned_flag == "BY_TUNING" else 1)

    detune_k_branch = {s.C_J: consumed_c_j, s.K_theta: -s.kappa_phase, s.B_eff: 0, s.m_theta2: 0}
    detune_b_branch = {s.C_J: consumed_c_j, s.K_theta: Ktheta_locus, s.B_eff: consumed_b_eff, s.m_theta2: 0}
    detune_m_branch = {s.C_J: consumed_c_j, s.K_theta: Ktheta_locus, s.B_eff: 0, s.m_theta2: s.m_theta2}
    detune_K, _ = token_from_discriminators(
        bracket_zero=algebraic_zero(bracket_general.subs(detune_k_branch)),
        B_eff_zero=algebraic_zero(detune_k_branch[s.B_eff]),
        m_theta_zero=algebraic_zero(detune_k_branch[s.m_theta2]),
        provenance_flags=provenance_flags,
    )
    detune_B, _ = token_from_discriminators(
        bracket_zero=algebraic_zero(bracket_general.subs(detune_b_branch)),
        B_eff_zero=algebraic_zero(detune_b_branch[s.B_eff]),
        m_theta_zero=algebraic_zero(detune_b_branch[s.m_theta2]),
        provenance_flags=provenance_flags,
    )
    detune_m, _ = token_from_discriminators(
        bracket_zero=algebraic_zero(bracket_general.subs(detune_m_branch)),
        B_eff_zero=algebraic_zero(detune_m_branch[s.B_eff]),
        m_theta_zero=algebraic_zero(detune_m_branch[s.m_theta2]),
        provenance_flags=provenance_flags,
    )
    expect_zero("Leg-C detuning K_theta off locus re-fires FAIL_C5_LONGITUDINAL_ZERO_MODE", 0 if detune_K == "FAIL_C5_LONGITUDINAL_ZERO_MODE" else 1)
    expect_zero("Leg-C detuning B_eff!=0 on locus re-fires FAIL_SECOND_CLASS_NOT_MAXWELL", 0 if detune_B == "FAIL_SECOND_CLASS_NOT_MAXWELL" else 1)
    expect_zero("Leg-C detuning m_theta^2!=0 on locus re-fires FAIL_SECOND_CLASS_NOT_MAXWELL", 0 if detune_m == "FAIL_SECOND_CLASS_NOT_MAXWELL" else 1)
    counterfactual = {
        "K_theta_forced_by_frozen_defs": True,
        "B_eff_removed_by_frozen_defs": True,
        "m_theta_zero_forced": True,
    }
    _, counterfactual_flag = token_from_discriminators(
        bracket_zero=algebraic_zero(bracket_general.subs(tuned_branch)),
        B_eff_zero=algebraic_zero(tuned_branch[s.B_eff]),
        m_theta_zero=algebraic_zero(tuned_branch[s.m_theta2]),
        provenance_flags=counterfactual,
    )
    expect_zero("Leg-C provenance-flag counterfactual flips BY_TUNING to WITH_PROVENANCE", 0 if counterfactual_flag == "WITH_PROVENANCE" else 1)

    kappa_real = sp.Symbol("kappa_phase_real", real=True)
    qstar = sp.solve(sp.Eq(sp.diff(sp.Rational(1, 2) * (kappa_real * sp.Symbol("q") + s.kappa4 * sp.Symbol("q") ** 2), sp.Symbol("q")), 0), sp.Symbol("q"))[0]
    expect_zero("OWNED C3 Lifshitz finite-k minimum k_star^2=-kappa_phase/(2*kappa4)", qstar + kappa_real / (2 * s.kappa4))
    lifshitz_region = sp.solve_univariate_inequality(qstar > 0, kappa_real)
    expect_bool("OWNED C3 k_star^2>0 iff kappa_phase<0", lifshitz_region == (kappa_real < 0))
    expect_bool("OWNED C3 kappa_phase>0 gives k_star^2<0", sp.simplify(qstar.subs(kappa_real, s.kappa_phase)).is_negative is True)
    expect_dim("CITED pathA_25 k_Rstar^2=40*K*m*rho0^4/hbar^2 has L^-2", K_EOS * M_GNLS * (RHO4**4) / (HBAR**2), LENGTH**-2)

    epsilon = sp.Symbol("epsilon", positive=True)
    u_T, dot_u_T, omega2_sym = sp.symbols("u_T dot_u_T omega2")
    transverse_lagrangian = sp.Rational(1, 2) * epsilon * dot_u_T**2 - sp.Rational(1, 2) * s.mu_R * s.k**2 * u_T**2
    transverse_frequency_probe = transverse_lagrangian.subs(dot_u_T**2, omega2_sym * u_T**2)
    transverse_dispersion_poly = sp.simplify(sp.diff(transverse_frequency_probe, u_T) / u_T)
    omega2 = sp.solve(sp.Eq(transverse_dispersion_poly, 0), omega2_sym)[0]
    c_gamma_sq = s.mu_R / s.rho_br
    baseline_speed = sp.simplify((omega2 / s.k**2).subs(epsilon, s.rho_br))
    disturbed_speed = sp.simplify((omega2 / s.k**2).subs(epsilon, 2 * s.rho_br))
    expect_zero("CONSUMED C4 transverse dispersion from L_T gives omega^2=(mu_R/epsilon)*k^2", omega2 - (s.mu_R / epsilon) * s.k**2)
    expect_zero("CONSUMED C4 transverse baseline speed equals consumed c_gamma^2=mu_R/rho_br", baseline_speed - c_gamma_sq)
    expect_zero("CONSUMED C4 epsilon!=rho_br fixture shifts transverse speed to mu_R/(2*rho_br)", disturbed_speed - s.mu_R / (2 * s.rho_br))
    transverse_baseline = "PASS_TRANSVERSE_UNDISTURBED" if algebraic_zero(baseline_speed - c_gamma_sq) else "FAIL_TRANSVERSE_DISTURBED"
    transverse_disturbed = "FAIL_TRANSVERSE_DISTURBED" if not algebraic_zero(disturbed_speed - c_gamma_sq) else "PASS_TRANSVERSE_UNDISTURBED"
    expect_zero("CONSUMED C4 baseline emits PASS_TRANSVERSE_UNDISTURBED", 0 if transverse_baseline == "PASS_TRANSVERSE_UNDISTURBED" else 1)
    expect_zero("CONSUMED C4 epsilon!=rho_br emits FAIL_TRANSVERSE_DISTURBED", 0 if transverse_disturbed == "FAIL_TRANSVERSE_DISTURBED" else 1)
    return provenance_token, zero_mode_token, tuned_token, tuned_flag


def run_ablations(consumed_c_j: sp.Expr, consumed_b_eff: sp.Expr, consumed_c_s0_sq: sp.Expr, s: Symbols) -> None:
    subbanner("Able-to-fail firewall ablations")
    n = RHO4
    expect_fail("Leg-A ablation drop m_GNLS from kinetic term breaks F-density homogeneity", dim_residual(n * (VELOCITY**2), F_DENSITY_4D))
    expect_fail("Leg-A ablation kappa_B*chi_B^2 no-gradient breaks F-density homogeneity", dim_residual(MASS / (TIME**2), F_DENSITY_4D))
    expect_fail("Leg-A ablation Gamma_B in place of n*Gamma_B breaks balance source dim", dim_residual(TIME**-1, LENGTH**-4 / TIME))
    expect_fail("Leg-A ablation handoff P_order=int mu_chi*n*Gamma_B d4X is inhomogeneous", dim_residual(F_DENSITY_4D * RHO4 * (TIME**-1) * (LENGTH**4), ENERGY / TIME))

    expect_fail("Consuming ablation B_eff multiply-vs-divide breaks citation integrity", s.rho_B0**2 * s.chi_c - s.rho_B0**2 / s.chi_c)
    expect_fail("Consuming ablation C_J sign +J*rho_B0 breaks citation integrity", s.J * s.rho_B0 + s.J * s.rho_B0)
    K, rho0, m_gnls = sp.symbols("K rho0 m_GNLS", positive=True)
    expect_fail("Consuming ablation c_s0^2 coefficient 4 breaks citation integrity", 4 * K * rho0**4 / m_gnls - 5 * K * rho0**4 / m_gnls)

    expect_zero("baseline consumed C_J still live after ablations", consumed_c_j + s.J * s.rho_B0)
    expect_zero("baseline consumed B_eff still live after ablations", consumed_b_eff - s.rho_B0**2 / s.chi_c)
    expect_zero("baseline consumed c_s0^2 still live after ablations", consumed_c_s0_sq - 5 * K * rho0**4 / m_gnls)


def drift_partition_residual(
    pre_d16: tuple[str, ...],
    operative: tuple[str, ...],
    retired: tuple[str, ...],
    *,
    verdict_n_delta: int = 0,
) -> sp.Expr:
    pre_set = set(pre_d16)
    operative_set = set(operative)
    retired_set = set(retired)
    residual = sp.Integer(0)
    residual += len(pre_d16) - len(pre_set)
    residual += len(operative) - len(operative_set)
    residual += len(retired) - len(retired_set)
    residual += sp.Integer(0 if retired_set == set(RETIRED_D16_DRIFT_MEMBERS) else 1)
    residual += sp.Integer(0 if operative_set.isdisjoint(retired_set) else 1)
    residual += sp.Integer(0 if operative_set | retired_set == pre_set else 1)
    residual += sp.Integer(0 if operative_set == pre_set - set(RETIRED_D16_DRIFT_MEMBERS) else 1)
    n = len(operative)
    residual += (n - (len(pre_d16) - len(RETIRED_D16_DRIFT_MEMBERS))) ** 2
    token = f"DRIFT({n + verdict_n_delta})"
    residual += sp.Integer(0 if token == EXPECTED_POST_D16_DRIFT_TOKEN else 1)
    return sp.simplify(residual)


def print_carried_tokens_and_drift() -> None:
    subbanner("Carried tokens, deferred items, and computed drift")
    print("  carried no-go tokens verbatim: FAIL_CAUCHY_STRAY_LONGITUDINAL; FAIL_C5_LONGITUDINAL_ZERO_MODE; C5_RESOLVED_MAXWELL_BY_TUNING (BY_TUNING, not WITH_PROVENANCE).")
    print("  pathA_25 lineage carried: FAIL_LIGHT_STARVED finite-k wall/smectic sign-flip route.")
    print("  SECOND_MEDIUM_DRIFT lineage note: chi_B package is Part-I drift; rho_B0 and chi_c also appear in pathA_41 Part-VI drift trio, cross-reference only.")
    print("  THETA_BRANCH_DEAD_NOT_ADMITTED: theta,J,rho_B0,K_theta/kappa_phase,chi_c,B are not live knobs of this stage.")
    print("  DEFERRED: Gate-L delta w=u_w translation-Goldstone hazard; J_chi!=0; f_throat/f_mix; dynamics adjunct.")
    print("  POSTULATED GLOBAL RETURNS: R_0=-M_0, R_1=-D_1 (not locally asserted).")
    pre_d16 = PRE_D16_DRIFT_MEMBERS
    retired = tuple(member for member in pre_d16 if member in RETIRED_D16_DRIFT_MEMBERS)
    operative = tuple(member for member in pre_d16 if member not in RETIRED_D16_DRIFT_MEMBERS)
    n = len(operative)
    token = f"DRIFT({n})"
    expect_zero("pre-Decision-16 drift enumeration computes six members", len(pre_d16) - 6)
    expect_zero("Decision-16 retired drift complement is exactly {alpha_aniso}", 0 if set(retired) == set(RETIRED_D16_DRIFT_MEMBERS) else 1)
    expect_zero("operative drift is the set partition pre_D16_DRIFT6 minus {alpha_aniso}", drift_partition_residual(pre_d16, operative, retired))
    expect_zero("COMPUTED operative DRIFT tally has five live chi_B inputs", n - 5)
    expect_zero("operative DRIFT token is built from computed n", 0 if token == EXPECTED_POST_D16_DRIFT_TOKEN else 1)
    print("  pre-D16: DRIFT(6) incl. alpha_aniso -> operative: DRIFT(5).")
    print(f"  {token} computed {{{'; '.join(operative)}}}.")
    print("  rung_W:140 pre-D16 reconciliation recorded six; Decision 16 removes only alpha_aniso, leaving the operative five-member chi_B package.")

    reinjected = operative + ("alpha_aniso",)
    expect_zero("Decision-16 drift reinjection fixture computes n=6", len(reinjected) - 6)
    expect_fail(
        "Decision-16 drift tooth: re-inject alpha_aniso trips operative DRIFT(5) partition",
        drift_partition_residual(pre_d16, reinjected, retired),
    )
    expect_fail(
        "Decision-16 drift tooth: corrupt computed n before token assembly trips DRIFT(5) equality",
        drift_partition_residual(pre_d16, operative, retired, verdict_n_delta=1),
    )
    expect_zero("Decision-16 baseline drift partition remains valid after copy mutations", drift_partition_residual(pre_d16, operative, retired))


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): CHI_B_ACTION_CLASSIFIED_RECOVERY_VERIFIED  (postulated two-phase chi_B action dimensionally classified; wall admitted; recovery reduction to the frozen S_leak assert-zero; theta-as-phi no-go carried)")
    print("  headline: ACTION_SPECIFIED_CLASSIFIED   (structure; POSTULATED microstructure, all terms labeled)")
    print("  recovery sub-verdict (EARNED rel. to the imposed chi_B split + declared W): RECOVERY_REDUCTION_VERIFIED   (target = frozen stage_243/244 S_leak, incl. the Gaussian one-mode anchor)")
    print("  carried no-go: FAIL_CAUCHY_STRAY_LONGITUDINAL (finite B_eff) / FAIL_C5_LONGITUDINAL_ZERO_MODE (B_eff=0); positive control C5_RESOLVED_MAXWELL_BY_TUNING flagged BY_TUNING NOT WITH_PROVENANCE; the only provenance sign-flip = Lifshitz (pathA_25 wall, killed)")
    print("  drift: pre-D16: DRIFT(6) incl. alpha_aniso -> operative: DRIFT(5) computed {chi_B; a_B; kappa_B; Gamma_B; gating structure}; THETA_BRANCH_DEAD_NOT_ADMITTED; cross-ref rho_B0, chi_c in pathA_41 Part-VI drift")
    print("  consumed: ledger_stage004 {L,T,M}+U(rho) single-well; ledger_stage005 c_s0^2=5*K*rho0^4/m_GNLS; ledger_stage003 c_gamma^2=mu_R/rho_br, C_J=-J*rho_B0, B_eff=rho_B0^2/chi_c, second-class classification rule")
    print("  labeled postulates: P1..P13 (incl. KINETIC_MASS_FACTOR_PINNED, SAME_DENSITY_DEGENERACY_POSTULATED, Decision-16-retired historical alpha_aniso, HANDOFF_P_ORDER_N_PLACEMENT_CORRECTED, global-return R_0=-M_0,R_1=-D_1 NOT locally asserted, throat=phase-conversion ontology)")
    print("  P7 live reframe: chi_B is THE postulated wall OP, not currently |P_parallel|^2; FUTURE_GATE_CHI_B_EQ_ABS_P_PARALLEL_SQ is high-risk/Part-VII-adjacent and requires a NEW T0 freeze (not foreclosed)")
    print("  honest scope: does NOT earn light; dynamics/energy ledger = labeled adjunct; wall-translation Goldstone + J_chi + f_throat/f_mix DEFERRED")


def main() -> None:
    banner("ledger_stage006_two_phase_chiB_ontology SymPy audit")
    s = make_symbols()
    print_pin_block()
    consumed_c_j, consumed_b_eff, consumed_c_s0_sq = run_consumed_citation_checks(s)
    run_leg_a_dimensions()
    run_leg_a_structural_closure(s)
    run_leg_a_wall_admission()
    run_leg_b_projection_profiles()
    run_leg_b_recovery_and_anchor(s)
    run_leg_c_no_go(s, consumed_c_j, consumed_b_eff)
    run_ablations(consumed_c_j, consumed_b_eff, consumed_c_s0_sq, s)
    print_carried_tokens_and_drift()
    print_verdict_labels()
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy classified ledger_stage006 chi_B ontology and recovery/no-go checks exactly")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print("OVERALL FAIL: SymPy stage006 audit did not close")
        raise SystemExit(1)
