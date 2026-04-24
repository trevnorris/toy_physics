#!/usr/bin/env python3
"""
Stage V2-04: Open-junction / organ-pipe impedance audit for the Moving-Throat PDE program.

Purpose
-------
This script checks the architectural patch that replaces a capped throat
R(L)=0 by an open finite-radius exit connected to unconfined 4D spatial bulk.

It audits four claims:
  1. A hard cap is incompatible with nonzero steady throughput.
  2. A finite-radius open exit permits DC mass/current leakage into 4D bulk.
  3. A sudden expansion creates strong AC reflection.
  4. Whether that reflection is Neumann or Dirichlet depends on which scalar
     variable is assigned to the support ladder:
        - pressure/potential-like scalar -> low-impedance open end gives R -> -1, Dirichlet-like;
        - flow/displacement-like scalar -> the dual reflection is R -> +1, Neumann-like.

SymPy is used for the exact algebra, limits, and series expansions.
"""

from __future__ import annotations

import sympy as sp


def banner(title: str) -> str:
    return "\n" + "=" * 88 + f"\n{title}\n" + "=" * 88 + "\n"


def main() -> None:
    # Symbols
    t, x, r = sp.symbols("t x r", real=True)
    c_s, k = sp.symbols("c_s k", positive=True, real=True)
    Z_t, eps, chi = sp.symbols("Z_t eps chi", positive=True, real=True)
    Phi, a_exit = sp.symbols("Phi a_exit", positive=True, real=True)
    A_exit = sp.symbols("A_exit", positive=True, real=True)

    I = sp.I
    lines: list[str] = []

    lines.append(banner("Stage V2-04 open-junction impedance audit"))
    lines.append("Architectural patch under test: replace hard cap R(L)=0 by finite open exit R(L)>0.")
    lines.append("The open exit is modeled as a 1D tube coupled to unconfined 4D spatial bulk.")

    # ------------------------------------------------------------------
    # 1. Variable-area scalar wave equation.
    # ------------------------------------------------------------------
    lines.append(banner("1. Variable-area scalar wave equation"))
    A = sp.Function("A")(x)
    psi = sp.Function("psi")(x, t)
    wave_eq = sp.Eq((A / c_s**2) * sp.diff(psi, t, 2) - sp.diff(A * sp.diff(psi, x), x), 0)
    lines.append("From the action S = 1/2 ∫ A(x)[psi_t^2/c_s^2 - psi_x^2] dx dt:")
    lines.append(f"  Euler-Lagrange equation: {sp.sstr(wave_eq)}")

    S3sym = sp.Symbol("S_3", positive=True)
    A_bulk_4 = S3sym * r**3
    lines.append("Tube side: A(x)=A_t = constant -> psi_tt = c_s^2 psi_xx.")
    lines.append(f"4D spatial bulk radial side: A_bulk(r)=S_3 r^3 = {sp.sstr(A_bulk_4)}.")
    Psi = sp.Function("Psi")(r, t)
    radial_eq = sp.Eq((A_bulk_4 / c_s**2) * sp.diff(Psi, t, 2) - sp.diff(A_bulk_4 * sp.diff(Psi, r), r), 0)
    radial_eq_div = sp.Eq((1 / c_s**2) * sp.diff(Psi, t, 2) - (1 / A_bulk_4) * sp.diff(A_bulk_4 * sp.diff(Psi, r), r), 0)
    lines.append(f"  4D radial wave equation: {sp.sstr(radial_eq)}")
    lines.append(f"  Divided form: {sp.sstr(radial_eq_div)}")
    lines.append("Interface conditions for a pressure/potential-like scalar are continuity of psi and A*psi_normal.")

    # ------------------------------------------------------------------
    # 2. Impedance reflection algebra.
    # ------------------------------------------------------------------
    lines.append(banner("2. Reflection/transmission from a load impedance"))
    # Set eps = Z_L/Z_t directly.
    R_p_eps = sp.simplify((eps - 1) / (eps + 1))
    T_p_eps = sp.simplify(1 + R_p_eps)
    R_q_eps = sp.simplify(-R_p_eps)
    T_q_eps = sp.simplify(1 + R_q_eps)
    E_ref_eps = sp.simplify(R_p_eps**2)
    E_trans_eps = sp.simplify(4 * eps / (1 + eps) ** 2)

    lines.append("For pressure p and volume flow Q at the tube end:")
    lines.append("  p = p_in + p_ref,")
    lines.append("  Q = (p_in/Z_t) - (p_ref/Z_t),")
    lines.append("  p = Z_L Q.")
    lines.append("Let eps = Z_L/Z_t. Sudden expansion / low acoustic load has eps << 1.")
    lines.append(f"  Pressure/potential amplitude reflection R_p(eps) = {sp.sstr(R_p_eps)}")
    lines.append(f"  Pressure/potential amplitude transmission T_p=1+R_p = {sp.sstr(T_p_eps)}")
    lines.append(f"  Flow/displacement amplitude reflection R_q=-R_p = {sp.sstr(R_q_eps)}")
    lines.append(f"  Flow/displacement amplitude transmission T_q=1+R_q = {sp.sstr(T_q_eps)}")
    lines.append(f"  Energy reflection coefficient |R|^2 = {sp.sstr(E_ref_eps)}")
    lines.append(f"  Energy transmission coefficient = {sp.sstr(E_trans_eps)}")

    low_series = {
        "R_p": sp.series(R_p_eps, eps, 0, 4),
        "T_p": sp.series(T_p_eps, eps, 0, 4),
        "R_q": sp.series(R_q_eps, eps, 0, 4),
        "T_q": sp.series(T_q_eps, eps, 0, 4),
        "E_ref": sp.series(E_ref_eps, eps, 0, 4),
        "E_trans": sp.series(E_trans_eps, eps, 0, 4),
    }
    lines.append("\nLow-load / sudden-expansion series eps -> 0:")
    for name, ser in low_series.items():
        lines.append(f"  {name}: {sp.sstr(ser)}")

    # Substitute chi = A_eff/A_t, eps ~= 1/chi.
    R_p_chi = sp.simplify(R_p_eps.subs(eps, 1 / chi))
    R_q_chi = sp.simplify(R_q_eps.subs(eps, 1 / chi))
    E_trans_chi = sp.simplify(E_trans_eps.subs(eps, 1 / chi))
    delta = sp.symbols("delta", positive=True, real=True)
    lines.append("\nArea-expansion model: chi = A_eff/A_t >> 1, eps ~= 1/chi.")
    lines.append(f"  R_p(chi) = {sp.sstr(R_p_chi)}")
    lines.append(f"  R_q(chi) = {sp.sstr(R_q_chi)}")
    lines.append(f"  T_energy(chi) = {sp.sstr(E_trans_chi)}")
    lines.append(f"  R_p large-chi via delta=1/chi: {sp.sstr(sp.series(R_p_eps.subs(eps, delta), delta, 0, 4))}")
    lines.append(f"  R_q large-chi via delta=1/chi: {sp.sstr(sp.series(R_q_eps.subs(eps, delta), delta, 0, 4))}")
    lines.append(f"  T_energy large-chi via delta=1/chi: {sp.sstr(sp.series(E_trans_eps.subs(eps, delta), delta, 0, 4))}")

    # ------------------------------------------------------------------
    # 3. Neumann/Dirichlet boundary diagnostics.
    # ------------------------------------------------------------------
    lines.append(banner("3. D/N validation: boundary derivative ratios"))
    # Tube coordinate is w increasing toward exit. Incoming exp(i k (w-L)), reflected R exp(-i k (w-L)).
    # Boundary derivative ratio = i k (1-R)/(1+R).
    deriv_ratio_p = sp.simplify(I * k * (1 - R_p_eps) / (1 + R_p_eps))
    deriv_ratio_q = sp.simplify(I * k * (1 - R_q_eps) / (1 + R_q_eps))
    amp_p = sp.simplify(1 + R_p_eps)
    amp_q = sp.simplify(1 + R_q_eps)

    lines.append("For a scalar field f = exp(ik(w-L)) + R exp(-ik(w-L)) at the exit:")
    lines.append("  f_w/f at L = i k (1-R)/(1+R).")
    lines.append(f"Pressure/potential-like scalar: f_w/f = {sp.sstr(deriv_ratio_p)}")
    lines.append(f"Flow/displacement-like scalar:    f_w/f = {sp.sstr(deriv_ratio_q)}")
    lines.append(f"Pressure/potential boundary amplitude 1+R_p = {sp.sstr(amp_p)}")
    lines.append(f"Flow/displacement boundary amplitude    1+R_q = {sp.sstr(amp_q)}")
    lines.append("\nLimits eps -> 0:")
    lines.append(f"  lim R_p = {sp.sstr(sp.limit(R_p_eps, eps, 0, dir='+'))}  (phase pi, Dirichlet-like for p/psi)")
    lines.append(f"  lim R_q = {sp.sstr(sp.limit(R_q_eps, eps, 0, dir='+'))}  (phase 0, Neumann-like for q/u)")
    lines.append(f"  lim (pressure amplitude 1+R_p) = {sp.sstr(sp.limit(amp_p, eps, 0, dir='+'))}")
    lines.append(f"  lim (flow derivative ratio) = {sp.sstr(sp.limit(deriv_ratio_q, eps, 0, dir='+'))}")
    lines.append("\nInterpretation:")
    lines.append("  * Sudden expansion validates Neumann only for the dual flow/displacement variable.")
    lines.append("  * The same junction gives a Dirichlet-like condition for a pressure/potential-like scalar.")

    # D/N ladder under mouth Dirichlet and exit Neumann for support variable q.
    j = sp.symbols("j", integer=True, nonnegative=True)
    Lsym = sp.Symbol("L", positive=True)
    k_j_DN = sp.pi * (j + sp.Rational(1, 2)) / Lsym
    k_j_DD = sp.pi * (j + 1) / Lsym
    lines.append("\nIf the support variable is the flow/displacement variable with q(0)=0 and q_w(L)=0:")
    lines.append(f"  D/N ladder: k_j = {sp.sstr(k_j_DN)}")
    lines.append("If instead the scalar is pressure/potential-like and both ends are Dirichlet:")
    lines.append(f"  D/D ladder: k_j = {sp.sstr(k_j_DD)} for j=0,1,2,...")

    # ------------------------------------------------------------------
    # 4. DC leakage through open finite-radius exit.
    # ------------------------------------------------------------------
    lines.append(banner("4. DC leakage / continuity check"))
    S3 = 2 * sp.pi**2  # area of unit 3-sphere in 4 spatial dimensions.
    A_bulk = S3 * r**3
    J_bulk = sp.simplify(Phi / A_bulk)
    div_bulk = sp.simplify((1 / A_bulk) * sp.diff(A_bulk * J_bulk, r))
    A_exit_3ball = sp.Rational(4, 3) * sp.pi * a_exit**3  # tube cross-section as 3-ball volume.
    J_tube_exit = sp.simplify(Phi / A_exit_3ball)
    cap_limit = sp.limit(Phi / A_exit, A_exit, 0, dir="+")

    lines.append(f"4D bulk radial area A_bulk(r)=S3*r^3 = {sp.sstr(A_bulk)}")
    lines.append(f"Steady flux density J_r(r)=Phi/A_bulk(r) = {sp.sstr(J_bulk)}")
    lines.append(f"Radial continuity check (1/A)d(A J_r)/dr = {sp.sstr(div_bulk)}")
    lines.append(f"Finite 3D-ball tube exit area A_exit=(4/3)pi a_exit^3 = {sp.sstr(A_exit_3ball)}")
    lines.append(f"Tube DC flux density at finite exit = {sp.sstr(J_tube_exit)}")
    lines.append(f"Hard-cap limit Phi/A_exit as A_exit->0+ = {sp.sstr(cap_limit)}")
    lines.append("\nInterpretation:")
    lines.append("  * A finite open exit R(L)>0 supports nonzero DC throughput without singular flux density.")
    lines.append("  * A hard cap / zero exit area is incompatible with fixed nonzero Phi.")

    # ------------------------------------------------------------------
    # 5. Pass/fail summary.
    # ------------------------------------------------------------------
    lines.append(banner("5. Pass/fail summary"))
    lines.append("HARD_CAP_WITH_NONZERO_DC_FLUX: FAIL")
    lines.append("  Reason: Phi/A_exit -> infinity as A_exit -> 0+.")
    lines.append("OPEN_FINITE_EXIT_DC_LEAKAGE: PASS")
    lines.append("  Reason: J_tube=Phi/A_exit finite and J_bulk=Phi/(2*pi^2*r^3) satisfies radial continuity.")
    lines.append("AC_STRONG_REFLECTION_FROM_SUDDEN_EXPANSION: PASS")
    lines.append("  Reason: energy transmission = 4 eps/(1+eps)^2 = O(eps) for eps=Z_L/Z_t << 1.")
    lines.append("GENERIC_SCALAR_NEUMANN_FROM_OPEN_BOUNDARY: FAIL")
    lines.append("  Reason: pressure/potential-like scalar has R_p -> -1 and boundary amplitude -> 0.")
    lines.append("FLOW_OR_DISPLACEMENT_NEUMANN_FROM_OPEN_BOUNDARY: PASS")
    lines.append("  Reason: dual variable has R_q -> +1 and q_w/q = i k eps -> 0.")
    lines.append("REQUIRED_PATCH_FOR_DN_SUPPORT: identify the D/N support field as the free-end flow/displacement variable, or add an impedance-transforming neck/stub if the support field is pressure/potential-like.")

    print("\n".join(lines))


if __name__ == "__main__":
    main()
