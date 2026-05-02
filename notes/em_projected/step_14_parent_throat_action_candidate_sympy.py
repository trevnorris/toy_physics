from __future__ import annotations

import sympy as sp
from sympy.calculus.euler import euler_equations
from sympy.functions.special.spherical_harmonics import Ynm


def assert_zero(label: str, expr: sp.Expr) -> None:
    residue = sp.factor(sp.together(sp.simplify(expr)))
    if residue != 0:
        raise AssertionError(f"{label} failed: {sp.sstr(residue)}")


def assert_nonzero(label: str, expr: sp.Expr) -> None:
    value = sp.factor(sp.together(sp.simplify(expr)))
    if value == 0:
        raise AssertionError(f"{label} unexpectedly vanished")


def boundary_value(expr: sp.Expr, var: sp.Symbol) -> sp.Expr:
    return sp.simplify(sp.limit(expr, var, sp.oo) - sp.limit(expr, var, -sp.oo))


def main() -> None:
    # Exact nonlinear Euler-Lagrange derivation in a local orthonormal angular chart.
    t_exact, w_exact, u_exact, v_exact = sp.symbols('t_exact w_exact u_exact v_exact', real=True)
    Rfield = sp.Function('R')(t_exact, w_exact, u_exact, v_exact)
    muSigma = sp.Function('muSigma')(Rfield, w_exact)
    TwSigma = sp.Function('TwSigma')(Rfield, w_exact)
    TOmegaSigma = sp.Function('TOmegaSigma')(Rfield, w_exact)
    USigma = sp.Function('USigma')(Rfield, w_exact)
    Rt_exact = sp.diff(Rfield, t_exact)
    Rw_exact = sp.diff(Rfield, w_exact)
    Ru_exact = sp.diff(Rfield, u_exact)
    Rv_exact = sp.diff(Rfield, v_exact)
    L_exact = (
        sp.Rational(1, 2) * muSigma * Rt_exact**2
        - sp.Rational(1, 2) * TwSigma * Rw_exact**2
        - sp.Rational(1, 2) * TOmegaSigma * (Ru_exact**2 + Rv_exact**2)
        - USigma
    )
    exact_el = euler_equations(L_exact, Rfield, [t_exact, w_exact, u_exact, v_exact])[0].lhs
    exact_el_lhs = (
        sp.diff(muSigma * Rt_exact, t_exact)
        - sp.diff(TwSigma * Rw_exact, w_exact)
        - sp.diff(TOmegaSigma * Ru_exact, u_exact)
        - sp.diff(TOmegaSigma * Rv_exact, v_exact)
        - sp.diff(muSigma, Rfield) * Rt_exact**2 / 2
        + sp.diff(TwSigma, Rfield) * Rw_exact**2 / 2
        + sp.diff(TOmegaSigma, Rfield) * (Ru_exact**2 + Rv_exact**2) / 2
        + sp.diff(USigma, Rfield)
    )
    assert_zero('exact nonlinear Euler-Lagrange equation', exact_el + exact_el_lhs)
    mutated_exact_el_lhs = exact_el_lhs - 2 * sp.diff(USigma, Rfield)
    assert_nonzero('mutated nonlinear Euler-Lagrange sign should fail', exact_el + mutated_exact_el_lhs)

    w_ibp = sp.symbols('w_ibp', real=True)
    Bfun = sp.Function('B')(w_ibp)
    Afun = sp.Function('A')(w_ibp)
    etaf = sp.Function('eta')(w_ibp)
    linear_boundary_density = -Bfun * etaf
    linear_bulk_density = sp.diff(Bfun, w_ibp) * etaf
    quad_boundary_density = -Afun * etaf**2 / 2
    quad_bulk_density = sp.diff(Afun, w_ibp) * etaf**2 / 2
    assert_zero(
        'generic linear IBP identity',
        -Bfun * sp.diff(etaf, w_ibp) - (sp.diff(linear_boundary_density, w_ibp) + linear_bulk_density),
    )
    assert_zero(
        'generic quadratic IBP identity',
        -Afun * etaf * sp.diff(etaf, w_ibp) - (sp.diff(quad_boundary_density, w_ibp) + quad_bulk_density),
    )
    assert_nonzero(
        'mutated quadratic IBP sign should fail',
        -Afun * etaf * sp.diff(etaf, w_ibp) - (sp.diff(-quad_boundary_density, w_ibp) + quad_bulk_density),
    )
    Bfun_concrete = (1 + w_ibp**2) * sp.exp(-w_ibp**2)
    Afun_concrete = (1 + w_ibp**2 / 2) * sp.exp(-w_ibp**2)
    eta_concrete = sp.exp(-w_ibp**2 / 2)
    linear_boundary_concrete = boundary_value(-Bfun_concrete * eta_concrete, w_ibp)
    linear_cross_concrete = sp.integrate(-Bfun_concrete * sp.diff(eta_concrete, w_ibp), (w_ibp, -sp.oo, sp.oo))
    linear_bulk_concrete = sp.integrate(sp.diff(Bfun_concrete, w_ibp) * eta_concrete, (w_ibp, -sp.oo, sp.oo))
    quad_boundary_concrete = boundary_value(-Afun_concrete * eta_concrete**2 / 2, w_ibp)
    quad_cross_concrete = sp.integrate(
        -Afun_concrete * eta_concrete * sp.diff(eta_concrete, w_ibp),
        (w_ibp, -sp.oo, sp.oo),
    )
    quad_bulk_concrete = sp.integrate(
        sp.diff(Afun_concrete, w_ibp) * eta_concrete**2 / 2,
        (w_ibp, -sp.oo, sp.oo),
    )
    assert_zero('concrete linear IBP boundary discharge', linear_boundary_concrete)
    assert_zero('concrete linear IBP with boundary', linear_cross_concrete - (linear_boundary_concrete + linear_bulk_concrete))
    assert_zero('concrete quadratic IBP boundary discharge', quad_boundary_concrete)
    assert_zero('concrete quadratic IBP with boundary', quad_cross_concrete - (quad_boundary_concrete + quad_bulk_concrete))

    # Small symbolic model for the minimal gauge-fixed nonlinear throat action
    # expanded around a static isotropic background R0(w).
    eps = sp.symbols('eps', real=True)
    eta = sp.symbols('eta', real=True)
    eta_t, eta_w = sp.symbols('eta_t eta_w', real=True)
    grad2 = sp.symbols('grad2', real=True, nonnegative=True)  # |∇_Ω eta|^2
    R0, R0p = sp.symbols('R0 R0p', real=True)
    mu0, Tw0, TO0, U0 = sp.symbols('mu0 Tw0 TO0 U0', real=True)
    TwR0, TwRR0, UR0, URR0 = sp.symbols('TwR0 TwRR0 UR0 URR0', real=True)
    d_TwR_R0p = sp.symbols('d_TwR_R0p', real=True)  # d/dw [Tw_R(R0,w) R0']
    d_Tw_R0p = sp.symbols('d_Tw_R0p', real=True)    # d/dw [Tw(R0,w) R0']

    # Local Taylor expansions around the static isotropic background.
    mu = mu0  # Higher derivatives do not contribute at quadratic order because R0_t = 0.
    Tw = Tw0 + eps * TwR0 * eta + sp.Rational(1, 2) * eps**2 * TwRR0 * eta**2
    TO = TO0  # Higher derivatives do not contribute at quadratic order because ∇_Ω R0 = 0.
    U = U0 + eps * UR0 * eta + sp.Rational(1, 2) * eps**2 * URR0 * eta**2

    Rt = eps * eta_t
    Rw = R0p + eps * eta_w

    L = sp.Rational(1, 2) * mu * Rt**2 - sp.Rational(1, 2) * Tw * Rw**2 - sp.Rational(1, 2) * TO * (eps**2 * grad2) - U
    Lexp = sp.expand(L)
    L1 = sp.expand(sp.diff(Lexp, eps).subs(eps, 0))
    L2_raw = sp.expand(sp.diff(Lexp, eps, 2).subs(eps, 0) / 2)
    assert_zero(
        'linear density before IBP',
        L1 - (-Tw0 * R0p * eta_w - sp.Rational(1, 2) * TwR0 * R0p**2 * eta - UR0 * eta),
    )

    # The linear term becomes a background equation plus a total derivative after integrating by parts.
    linear_bulk = sp.expand(d_Tw_R0p - sp.Rational(1, 2) * TwR0 * R0p**2 - UR0)
    linear_after_ibp = sp.expand((L1 + Tw0 * R0p * eta_w) + d_Tw_R0p * eta)
    assert_zero('linear background coefficient after IBP', linear_after_ibp - linear_bulk * eta)

    # The quadratic cross term -TwR0*R0p*eta*eta_w is integrated by parts.
    A = TwR0 * R0p
    # Replacement after integration by parts:
    cross_after_ibp = sp.Rational(1, 2) * d_TwR_R0p * eta**2
    cross_term = -A * eta * eta_w
    assert_zero('raw quadratic cross term', sp.diff(sp.diff(L2_raw, eta), eta_w) + A)

    K_eta = sp.expand(URR0 - d_TwR_R0p + sp.Rational(1, 2) * TwRR0 * R0p**2)
    L2_ibp = sp.expand(L2_raw - cross_term + cross_after_ibp)
    canonical_L2 = sp.expand(
        sp.Rational(1, 2) * mu0 * eta_t**2
        - sp.Rational(1, 2) * Tw0 * eta_w**2
        - sp.Rational(1, 2) * TO0 * grad2
        - sp.Rational(1, 2) * K_eta * eta**2
    )
    assert_zero('quadratic density after IBP', L2_ibp - canonical_L2)
    canonical_L2_mutated = sp.expand(
        sp.Rational(1, 2) * mu0 * eta_t**2
        - sp.Rational(1, 2) * Tw0 * eta_w**2
        - sp.Rational(1, 2) * TO0 * grad2
        - sp.Rational(1, 2) * (URR0 + d_TwR_R0p + sp.Rational(1, 2) * TwRR0 * R0p**2) * eta**2
    )
    assert_nonzero('mutated K_eta sign should fail', L2_ibp - canonical_L2_mutated)

    th, ph = sp.symbols('th ph', real=True)
    Y20 = sp.expand_func(Ynm(2, 0, th, ph))
    spherical_laplacian_Y20 = sp.simplify(
        sp.diff(sp.sin(th) * sp.diff(Y20, th), th) / sp.sin(th)
        + sp.diff(Y20, ph, 2) / sp.sin(th)**2
    )
    assert_zero('Y20 spherical-Laplacian eigenvalue', spherical_laplacian_Y20 + 6 * Y20)
    assert_nonzero('mutated Y20 spherical-Laplacian eigenvalue should fail', spherical_laplacian_Y20 + 5 * Y20)

    t_mode, w_mode = sp.symbols('t_mode w_mode', real=True)
    q = sp.Function('q')(t_mode, w_mode)
    mu_mode = sp.Function('mu_mode')(w_mode)
    Tw_mode = sp.Function('Tw_mode')(w_mode)
    TO_mode = sp.Function('TO_mode')(w_mode)
    K_mode = sp.Function('K_mode')(w_mode)
    q_t = sp.diff(q, t_mode)
    q_w = sp.diff(q, w_mode)
    angular_norm = sp.simplify(sp.integrate(sp.sin(th) * Y20**2, (ph, 0, 2 * sp.pi), (th, 0, sp.pi)))
    angular_grad = sp.simplify(
        sp.integrate(
            sp.sin(th) * (sp.diff(Y20, th)**2 + sp.diff(Y20, ph)**2 / sp.sin(th)**2),
            (ph, 0, 2 * sp.pi),
            (th, 0, sp.pi),
        )
    )
    projected_modal_density = sp.simplify(
        angular_norm * (mu_mode * q_t**2 / 2 - Tw_mode * q_w**2 / 2 - K_mode * q**2 / 2)
        - angular_grad * TO_mode * q**2 / 2
    )
    closed_modal_density = sp.simplify(
        angular_norm
        * (mu_mode * q_t**2 / 2 - Tw_mode * q_w**2 / 2 - (K_mode + 6 * TO_mode) * q**2 / 2)
    )
    assert_zero('Y20 modal norm', angular_norm - 1)
    assert_zero('Y20 modal angular stiffness', angular_grad - 6)
    assert_zero('Y20 projected modal density', projected_modal_density - closed_modal_density)
    modal_el = euler_equations(projected_modal_density, q, [t_mode, w_mode])[0].lhs
    modal_el_lhs = (
        sp.diff(mu_mode * q_t, t_mode)
        - sp.diff(Tw_mode * q_w, w_mode)
        + (K_mode + 6 * TO_mode) * q
    )
    assert_zero('Y20 fused modal Euler-Lagrange equation', modal_el + modal_el_lhs)

    lines: list[str] = []
    lines.append("Minimal nonlinear parent throat action candidate\n")
    lines.append("Gauge-fixed action density:\n")
    lines.append("  L_Sigma = 1/2 mu(R,w) R_t^2 - 1/2 T_w(R,w) R_w^2 - 1/2 T_Omega(R,w) |∇_Ω R|^2 - U_Sigma(R,w)\n")
    lines.append("with total action S_total = S_psi[psi,A;Sigma] + S_EM[A] + S_Sigma[R].\n")
    lines.append("SymPy also derives the exact nonlinear Euler-Lagrange equation in a local orthonormal angular chart,\n")
    lines.append("which is the coordinate-level proof of the boxed parent equation in the note.\n")
    lines.append("\n")
    lines.append("Background expansion: R = R0(w) + eps*eta(Ω,w,t).\n")
    lines.append("\n")
    lines.append(f"Linear-order density before integration by parts:\n  L1 = {sp.sstr(L1)}\n")
    lines.append("\n")
    lines.append("After integrating the linear derivative term by parts, the bulk linear coefficient is:\n")
    lines.append(f"  E_bg = {sp.sstr(linear_bulk)}\n")
    lines.append("The carried linear boundary term is\n")
    lines.append("  [-T_w(R0,w) R0' eta]_{boundary}.\n")
    lines.append("So the static isotropic background equation is:\n")
    lines.append("  d/dw[T_w(R0,w) R0'] - 1/2 T_{w,R}(R0,w) (R0')^2 - U_{Sigma,R}(R0,w) = 0\n")
    lines.append("\n")
    lines.append("Raw quadratic density before integrating the cross term by parts:\n")
    lines.append(f"  L2_raw = {sp.sstr(L2_raw)}\n")
    lines.append("\n")
    lines.append("The carried quadratic boundary term is\n")
    lines.append("  [-1/2 T_{w,R}(R0,w) R0' eta^2]_{boundary}.\n")
    lines.append("Quadratic density after integrating -T_{w,R}(R0,w) R0' eta eta_w by parts:\n")
    lines.append(f"  L2 = {sp.sstr(L2_ibp)}\n")
    lines.append("\n")
    lines.append("Therefore the audited linear wall action is recovered with\n")
    lines.append("  mu_eta(w) = mu(R0,w)\n")
    lines.append("  T_w(w)    = T_w(R0,w)\n")
    lines.append("  T_Omega(w)= T_Omega(R0,w)\n")
    lines.append(f"  K_eta(w)  = {sp.sstr(K_eta)}\n")
    lines.append("i.e.\n")
    lines.append("  K_eta = U_{Sigma,RR}(R0,w) - d/dw[T_{w,R}(R0,w) R0'] + 1/2 T_{w,RR}(R0,w) (R0')^2\n")
    lines.append("\n")
    lines.append("The quadratic fluctuation density is therefore\n")
    lines.append("  L2 = 1/2 mu_eta eta_t^2 - 1/2 T_w eta_w^2 - 1/2 T_Omega |∇_Ω eta|^2 - 1/2 K_eta eta^2.\n")
    lines.append("\n")
    lines.append("Projecting onto Y_lm then gives the modal operator\n")
    lines.append("  mu_eta q_lm,tt - ∂_w(T_w q_lm,w) + [K_eta + l(l+1) T_Omega] q_lm = S_lm.\n")
    lines.append("For l=0 this reproduces the scalar lane, and the script now also integrates the genuine Y20 harmonic\n")
    lines.append("into the reduced modal Lagrangian itself, where the angular sector contributes the exact +6 T_Omega factor.\n")
    lines.append("\n")
    lines.append("Axisymmetric reduction back to (a,L) uses\n")
    lines.append("  M_AB = 4π ∫ dw mu_eta alpha_A alpha_B\n")
    lines.append("  K_AB = 4π ∫ dw [T_w alpha_A' alpha_B' + K_eta alpha_A alpha_B].\n")
    lines.append("\n")
    lines.append("This verifies: the minimal nonlinear gauge-fixed S_Sigma promotes the throat to an autonomous parent field,\n")
    lines.append("its quadratic limit is the audited S_eta^(2), and the Y20 modal EL is recovered from the genuine S^2 reduction.\n")
    lines.append("\nSTATUS: PASS\n")

    out = "".join(lines)
    print(out)


if __name__ == "__main__":
    main()
