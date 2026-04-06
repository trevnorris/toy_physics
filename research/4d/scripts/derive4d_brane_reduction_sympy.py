#!/usr/bin/env python3
"""derive4d_brane_reduction_sympy_v13.py

SymPy harness for Paper-7 hard-mode 4D model derivations.

What changed vs v12
-------------------
This version adds a **paper-style report mode** that prints the *structural* equations
needed to understand where Poisson-like and EM-like sectors can arise **without forcing**.

Key idea: SymPy can easily explode when asked to fully simplify component-expanded
expressions (especially the 4D Euler equation in θ,ρ form). We therefore separate:

  * **Summary / abstract forms** (fast, paper-readable)
  * **Component-expanded forms** (available behind explicit flags)

Nothing here asserts that Poisson or Maxwell must emerge. The script prints:

  * exact identities (projection, continuity, flux/source terms),
  * and clearly labeled *conditional* bridges (quasi-static / longitudinal / linearized).

Usage highlights
----------------

1) Paper-style report (recommended)

  # Gauge ON (default): shows where Lorentz/vorticity structure can enter
  python derive4d_brane_reduction_sympy_v13.py --report

  # Gauge OFF: shows the "pure potential-flow" baseline (Ω=0 in the bulk)
  python derive4d_brane_reduction_sympy_v13.py --report --no_gauge

2) Full brane wrapup (prints more intermediate identities)

  python derive4d_brane_reduction_sympy_v13.py --wrapup --simplify
  python derive4d_brane_reduction_sympy_v13.py --no_gauge --wrapup --simplify

3) Deep bulk Madelung/Euler expansions (can be huge)

  python derive4d_brane_reduction_sympy_v13.py --show_madelung --show_euler --simplify
  python derive4d_brane_reduction_sympy_v13.py --show_madelung --show_euler --full_simplify

If you want the *component form* of the 4D Euler equation in report mode, add:

  --report_with_components

"""

from __future__ import annotations

import argparse
import sympy as sp


# ----------------------------
# Printing / simplification helpers
# ----------------------------


def _light_simplify(expr, *, max_ops: int) -> sp.Expr:
    """A deliberately conservative simplifier.

    Strategy:
      - Skip objects containing Derivative/Integral (common SymPy slow-path triggers).
      - If count_ops(expr) <= max_ops, run sp.simplify.
      - Otherwise return expr unchanged.

    This is a tooling choice (output readability), not a mathematical statement.
    """

    if expr.has(sp.Derivative) or expr.has(sp.Integral):
        return expr

    try:
        ops = int(expr.count_ops()) if hasattr(expr, "count_ops") else int(sp.count_ops(expr))
        if ops > max_ops:
            return expr
        return sp.simplify(expr)
    except Exception:
        return expr


def _maybe_simplify(expr, *, mode: str, max_ops: int) -> sp.Expr:
    if mode == "none":
        return expr
    if mode == "full":
        try:
            return sp.simplify(expr)
        except Exception:
            return expr
    return _light_simplify(expr, max_ops=max_ops)


def latex_print(title: str, obj, *, simplify_mode: str = "none", max_ops: int = 400) -> None:
    """Pretty-print SymPy objects as LaTeX."""

    print("\n" + "=" * 92)
    print(title)
    print("-" * 92)

    def fmt(o):
        if isinstance(o, sp.Equality):
            lhs = _maybe_simplify(o.lhs, mode=simplify_mode, max_ops=max_ops)
            rhs = _maybe_simplify(o.rhs, mode=simplify_mode, max_ops=max_ops)
            return sp.latex(sp.Eq(lhs, rhs))
        return sp.latex(_maybe_simplify(o, mode=simplify_mode, max_ops=max_ops))

    if isinstance(obj, (list, tuple)):
        for k, o in enumerate(obj, start=1):
            print(f"[{k}] {fmt(o)}\n")
    else:
        print(fmt(obj))

    print("=" * 92 + "\n")


def text_note(title: str, lines: list[str]) -> None:
    """Plain-text note blocks (keeps non-LaTeX commentary readable)."""
    print("\n" + "=" * 92)
    print(title)
    print("-" * 92)
    for ln in lines:
        print(ln)
    print("=" * 92 + "\n")


# ----------------------------
# Main
# ----------------------------


def main() -> None:
    ap = argparse.ArgumentParser()

    # Simplification controls
    ap.add_argument("--simplify", action="store_true", help="Light simplification (safe, size-limited).")
    ap.add_argument("--full_simplify", action="store_true", help="Force sympy.simplify() everywhere (can hang).")
    ap.add_argument(
        "--simplify_max_ops",
        type=int,
        default=400,
        help="Operation-count cutoff for --simplify (skip simplify if larger).",
    )

    # Output modes
    ap.add_argument(
        "--report",
        action="store_true",
        help="Paper-style report mode: prints compact structural equations and sector splits (recommended).",
    )
    ap.add_argument(
        "--report_with_integrals",
        action="store_true",
        help="In --report mode, include the integral definitions (ρ_brane, J_brane, Π_ij, source terms).",
    )
    ap.add_argument(
        "--report_with_components",
        action="store_true",
        help="In --report mode, also print the component-expanded 4D Euler equation (can be huge).",
    )

    # Feature flags (bulk)
    ap.add_argument("--show_madelung", action="store_true", help="Show Madelung dictionary + PDEs (expanded).")
    ap.add_argument("--show_euler", action="store_true", help="Show Euler-like equation from Madelung real part (expanded).")
    ap.add_argument(
        "--show_vorticity",
        action="store_true",
        help="Show bulk vorticity identity Ω_ij = ∂i v_j - ∂j v_i (and its gauge reduction).",
    )
    ap.add_argument(
        "--check_madelung_full",
        action="store_true",
        help="Explicitly simplify J_i - ρ v_i residuals (may be slow with gauge).",
    )
    ap.add_argument(
        "--check_madelung_current",
        action="store_true",
        help="Alias for --check_madelung_full (legacy option name).",
    )

    # Brane-observable diagnostics (projection kinematics)
    ap.add_argument(
        "--normalize_weight",
        action="store_true",
        help="Normalize projection weight W̃ = W / ∫_{-Wproj}^{Wproj} W.",
    )
    ap.add_argument(
        "--show_brane_velocity",
        action="store_true",
        help="Print v_brane = J_brane / rho_brane and related brane kinematics.",
    )
    ap.add_argument(
        "--show_brane_vorticity",
        action="store_true",
        help="Print brane-observed 3D vorticity ω = ∇×v_brane (compact current/density form).",
    )
    ap.add_argument(
        "--show_brane_vorticity_expanded",
        action="store_true",
        help="Also print ω components as raw derivatives of J_brane/rho_brane (can be huge).",
    )
    ap.add_argument(
        "--show_poisson_bridge",
        action="store_true",
        help="Print exact div(v_brane) identity from projected continuity and the quasi-static Poisson candidate.",
    )
    ap.add_argument(
        "--check_separable_mode",
        action="store_true",
        help="Check separable ansatz ψ=Ψ3(x,t)χ(w) reduces v_brane to 3D velocity (no gauge).",
    )
    ap.add_argument(
        "--show_brane_momentum",
        action="store_true",
        help="Print brane-projected 3D momentum balance and the extra stress/source sectors from w-projection.",
    )
    ap.add_argument(
        "--show_brane_abstract",
        action="store_true",
        help="Print human-readable brane kinematic identities using abstract rho_brane(x,t) and J_brane(x,t) symbols.",
    )
    ap.add_argument(
        "--show_brane_helmholtz",
        action="store_true",
        help="Print Helmholtz (longitudinal/transverse) reconstruction identities for v_brane (pure math, not physics).",
    )
    ap.add_argument(
        "--show_brane_vorticity_dynamics",
        action="store_true",
        help="Print an abstract curl(momentum) diagnostic showing how transverse/vorticity sector can be sourced (no closure assumed).",
    )
    ap.add_argument(
        "--show_brane_euler_form",
        action="store_true",
        help="Derive brane acceleration (Euler) form from abstract continuity+momentum and show curl/div structure (extra sectors explicit).",
    )
    ap.add_argument(
        "--show_em_analog",
        action="store_true",
        help="Print an emergent-EM dictionary built from brane v and enthalpy and show homogeneous Maxwell identities (no dynamics assumed).",
    )

    ap.add_argument(
        "--wrapup",
        action="store_true",
        help="Enable the consolidated Poisson/EM wrap-up bundle (like v12; does not force physics).",
    )

    # Gauge toggles
    ap.add_argument(
        "--no_gauge",
        action="store_true",
        help="Drop gauge potentials (keeps expressions smaller; velocity becomes ∝∇θ).",
    )
    ap.add_argument(
        "--include_gauge",
        action="store_true",
        help="Alias (gauge is ON by default unless --no_gauge).",
    )

    ap.add_argument("--show_maxwell", action="store_true", help="Print compact 4+1D Maxwell+localization form.")

    args = ap.parse_args()

    # Resolve simplify mode
    if args.full_simplify:
        simplify_mode = "full"
    elif args.simplify:
        simplify_mode = "light"
    else:
        simplify_mode = "none"

    max_ops = int(args.simplify_max_ops)

    # Alias resolution
    if args.check_madelung_current:
        args.check_madelung_full = True

    # Gauge default is ON unless --no_gauge
    use_gauge = (not args.no_gauge) or args.include_gauge

    # Convenience bundle (kept compatible with v12)
    if args.wrapup:
        args.show_brane_velocity = True
        args.show_brane_vorticity = True
        args.show_poisson_bridge = True
        args.show_brane_helmholtz = True
        args.show_brane_momentum = True
        args.show_brane_vorticity_dynamics = True
        args.show_brane_euler_form = True
        args.show_em_analog = True
        args.show_brane_abstract = True
        args.check_separable_mode = True

    # Report mode: force a clean, paper-readable subset
    report_mode = bool(args.report)
    if report_mode:
        # Always enable the core brane outputs that connect to Poisson/EM pathways.
        args.show_brane_velocity = True
        args.show_brane_vorticity = True
        args.show_poisson_bridge = True
        args.show_brane_helmholtz = True
        args.show_brane_momentum = True
        args.show_brane_euler_form = True
        args.show_brane_vorticity_dynamics = True
        args.show_em_analog = True
        args.show_brane_abstract = True
        args.check_separable_mode = True

        # If user didn't explicitly ask, default to including integral defs in report
        if not args.report_with_integrals:
            # keep False unless user asked; we will still print compact abstract defs
            pass

    # ----------------------------
    # Coordinates / parameters
    # ----------------------------

    t = sp.Symbol("t", real=True)
    x, y, z, w = sp.symbols("x y z w", real=True)

    # Use \bar{\h} to match prior project output (purely notational).
    hbar = sp.Symbol(r"\bar{\h}", positive=True, real=True)
    m = sp.Symbol("m", positive=True, real=True)
    q = sp.Symbol("q", real=True)

    # Geometry knobs (enter through V_conf)
    a = sp.Symbol("a", positive=True, real=True)
    L = sp.Symbol("L", positive=True, real=True)

    # ----------------------------
    # Frozen n=5 EOS identities (anchor)
    # ----------------------------

    rho_sym = sp.Symbol("rho", positive=True, real=True)
    K = sp.Symbol("K", real=True)

    P_rho = K * rho_sym**5
    cs2_rho = sp.diff(P_rho, rho_sym)  # = 5K rho^4
    h_rho = sp.Rational(5, 4) * K * rho_sym**4
    U_rho = sp.Rational(1, 4) * K * rho_sym**5

    # In report mode we keep this anchor (short and important); in non-report mode we print anyway.
    latex_print(
        "0) Frozen n=5 EOS anchor",
        [
            sp.Eq(sp.Function("P")(rho_sym), P_rho),
            sp.Eq(sp.Function("c_s^2")(rho_sym), cs2_rho),
            sp.Eq(sp.Function("h")(rho_sym), h_rho),
            sp.Eq(sp.Function("U")(rho_sym), U_rho),
            sp.Eq(sp.Symbol(r"\rho\,h'(\rho)"), sp.simplify(rho_sym * sp.diff(h_rho, rho_sym))),
        ],
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )

    # ----------------------------
    # Fields
    # ----------------------------

    psi = sp.Function(r"\psi")(t, x, y, z, w)
    psib = sp.Function(r"\bar\psi")(t, x, y, z, w)  # interpret as conjugate at the end

    # Gauge potentials
    A0 = sp.Function("A0")(t, x, y, z, w)
    Ax = sp.Function("Ax")(t, x, y, z, w)
    Ay = sp.Function("Ay")(t, x, y, z, w)
    Az = sp.Function("Az")(t, x, y, z, w)
    Aw = sp.Function("Aw")(t, x, y, z, w)

    Vconf = sp.Function(r"V_{\rm conf}")(x, y, z, w, a, L)

    rho = psi * psib

    # ----------------------------
    # Covariant derivatives (minimal coupling)
    # ----------------------------

    # Convention:
    #   iħ ∂_t ψ = [ (1/2m)(-iħ∇ - qA)^2 + qA0 + V + h(ρ) ] ψ

    def Dt_psi(f):
        return sp.diff(f, t) + (sp.I * q / hbar) * A0 * f if use_gauge else sp.diff(f, t)

    def D_psi(coord, Acomp, f):
        return sp.diff(f, coord) - (sp.I * q / hbar) * Acomp * f if use_gauge else sp.diff(f, coord)

    # Conjugate field carries opposite gauge charge
    def D_psib(coord, Acomp, f):
        return sp.diff(f, coord) + (sp.I * q / hbar) * Acomp * f if use_gauge else sp.diff(f, coord)

    def D2_psi(f):
        """Covariant Laplacian on ψ: Σ_i D_i D_i ψ, i∈{x,y,z,w}."""
        out = 0
        for coord, Acomp in [(x, Ax), (y, Ay), (z, Az), (w, Aw)]:
            Di_f = D_psi(coord, Acomp, f)
            out += sp.diff(Di_f, coord)
            if use_gauge:
                out += -(sp.I * q / hbar) * Acomp * Di_f
        return out

    # ----------------------------
    # GNLS equation (compact form)
    # ----------------------------

    h_of_rho = sp.Rational(5, 4) * K * rho**4
    gnls_compact = sp.Eq(
        sp.I * hbar * Dt_psi(psi),
        -(hbar**2 / (2 * m)) * D2_psi(psi) + (Vconf + h_of_rho) * psi,
    )

    latex_print(
        f"1) Canonical {'gauged ' if use_gauge else ''}4D GNLS form (compact)",
        gnls_compact,
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )

    # ----------------------------
    # Current + continuity
    # ----------------------------

    def J_component(coord, Acomp):
        return (hbar / (2 * sp.I * m)) * (psib * D_psi(coord, Acomp, psi) - psi * D_psib(coord, Acomp, psib))

    Jx = J_component(x, Ax)
    Jy = J_component(y, Ay)
    Jz = J_component(z, Az)
    Jw = J_component(w, Aw)

    # In report mode we avoid clutter by printing the *definition* of J (still short).
    # In non-report mode we also print explicit component expressions.
    if (not report_mode) or args.report_with_integrals:
        latex_print(
            "2) Canonical 4D current components (Jx,Jy,Jz,Jw)",
            [Jx, Jy, Jz, Jw],
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )
    else:
        Ji = sp.Symbol("J_i")
        latex_print(
            "2) Canonical 4D current definition (component-wise)",
            sp.Eq(
                Ji,
                (hbar / (2 * sp.I * m)) * (sp.Symbol("\\bar\\psi") * sp.Symbol("D_i\\psi") - sp.Symbol("\\psi") * sp.Symbol("D_i\\bar\\psi")),
            ),
            simplify_mode="none",
            max_ops=max_ops,
        )

    cont4 = sp.Eq(sp.diff(rho, t) + sp.diff(Jx, x) + sp.diff(Jy, y) + sp.diff(Jz, z) + sp.diff(Jw, w), 0)
    latex_print(
        "3) 4D continuity identity: rho_t + div_4 J = 0",
        cont4,
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )

    # ----------------------------
    # Brane projection of continuity
    # ----------------------------

    W = sp.Function("W")(w)
    Wproj = sp.Symbol(r"W_{\rm proj}", positive=True, real=True)

    if args.normalize_weight:
        Nw = sp.Integral(W, (w, -Wproj, Wproj))
        Wt = W / Nw
        W_label = r"\widetilde W(w)=W(w)/\mathcal N_W"
    else:
        Nw = sp.Integer(1)
        Wt = W
        W_label = r"W(w)\ (unnormalized)"

    rho_brane = sp.Integral(Wt * rho, (w, -Wproj, Wproj))
    Jx_brane = sp.Integral(Wt * Jx, (w, -Wproj, Wproj))
    Jy_brane = sp.Integral(Wt * Jy, (w, -Wproj, Wproj))
    Jz_brane = sp.Integral(Wt * Jz, (w, -Wproj, Wproj))

    # Exact source term from J_w (finite window + weight gradient)
    boundary_flux = -(Wt.subs(w, Wproj) * Jw.subs(w, Wproj) - Wt.subs(w, -Wproj) * Jw.subs(w, -Wproj))
    weight_grad_term = sp.Integral(sp.diff(Wt, w) * Jw, (w, -Wproj, Wproj))
    S_rho_brane = boundary_flux + weight_grad_term

    cont_brane_exact = sp.Eq(
        sp.diff(rho_brane, t) + sp.diff(Jx_brane, x) + sp.diff(Jy_brane, y) + sp.diff(Jz_brane, z),
        S_rho_brane,
    )

    # Report-style: always show the structure of projection; optionally show explicit integrals
    latex_print(
        f"4) Brane projection weight: {W_label}",
        sp.Eq(sp.Symbol("W_t"), Wt),
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )
    if args.normalize_weight:
        latex_print(
            "4b) Projection normalization constant 𝒩_W = ∫ W(w) dw over |w|≤Wproj",
            sp.Eq(sp.Symbol(r"\mathcal N_W"), Nw),
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )

    if (not report_mode) or args.report_with_integrals:
        latex_print(
            "5) Brane-projected density rho_brane = ∫ Wt rho dw",
            rho_brane,
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )
        latex_print(
            "6) Brane-projected current components (integrated in w)",
            [Jx_brane, Jy_brane, Jz_brane],
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )
    else:
        latex_print(
            "5–6) Brane observables (definitions)",
            [
                sp.Eq(sp.Symbol(r"\rho_{\rm brane}"), sp.Symbol(r"\int W\,\rho\,dw")),
                sp.Eq(sp.Symbol(r"\mathbf J_{\rm brane}"), sp.Symbol(r"\int W\,\mathbf J\,dw")),
            ],
            simplify_mode="none",
            max_ops=max_ops,
        )

    latex_print(
        "7) Exact projected continuity: (rho_brane)_t + div_3(J_brane) = S_rho_brane",
        cont_brane_exact,
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )
    latex_print(
        "8) Exact brane source term S_rho_brane (purely from J_w)",
        S_rho_brane,
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )

    if report_mode:
        text_note(
            "8-note) Interpretation of S_rho_brane (no claims; just the identity)",
            [
                "S_rho_brane is the *only* way the 3D brane density fails to be conserved.",
                "It arises from bulk flow along w (J_w), through:",
                "  (i) boundary flux at w=±Wproj, and",
                "  (ii) the weight-gradient coupling ∫ (∂w Wt) J_w dw.",
                "If the brane is effectively closed (J_w≈0 in the support of Wt), then S_rho_brane≈0.",
                "If not, the brane is an open system (" "leakage" ") and mass can be injected/removed.",
            ],
        )

    # ----------------------------
    # Brane-observed velocity / vorticity diagnostics
    # ----------------------------

    if args.show_brane_velocity or args.show_brane_vorticity or args.show_poisson_bridge or args.check_separable_mode or args.show_brane_momentum:
        vbx = Jx_brane / rho_brane
        vby = Jy_brane / rho_brane
        vbz = Jz_brane / rho_brane

    if args.show_brane_velocity:
        latex_print(
            "8b) Brane-observed velocity v_brane = J_brane / rho_brane",
            [vbx, vby, vbz],
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )

    if args.show_brane_vorticity:
        # Compact identity: curl(J/ρ) = (curl J)/ρ - (∇ρ × J)/ρ^2
        omega_x = (sp.diff(Jz_brane, y) - sp.diff(Jy_brane, z)) / rho_brane - (
            sp.diff(rho_brane, y) * Jz_brane - sp.diff(rho_brane, z) * Jy_brane
        ) / (rho_brane**2)
        omega_y = (sp.diff(Jx_brane, z) - sp.diff(Jz_brane, x)) / rho_brane - (
            sp.diff(rho_brane, z) * Jx_brane - sp.diff(rho_brane, x) * Jz_brane
        ) / (rho_brane**2)
        omega_z = (sp.diff(Jy_brane, x) - sp.diff(Jx_brane, y)) / rho_brane - (
            sp.diff(rho_brane, x) * Jy_brane - sp.diff(rho_brane, y) * Jx_brane
        ) / (rho_brane**2)

        latex_print(
            "8c) Brane-observed vorticity ω_brane = ∇×v_brane (compact form)",
            [omega_x, omega_y, omega_z],
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )

        if args.show_brane_vorticity_expanded:
            omega_x2 = sp.diff(vbz, y) - sp.diff(vby, z)
            omega_y2 = sp.diff(vbx, z) - sp.diff(vbz, x)
            omega_z2 = sp.diff(vby, x) - sp.diff(vbx, y)

            latex_print(
                "8c-expanded) Same ω components as raw derivatives of J/ρ (often huge)",
                [omega_x2, omega_y2, omega_z2],
                simplify_mode=simplify_mode,
                max_ops=max_ops,
            )

    # ----------------------------
    # Human-readable brane identities (abstract fields)
    # ----------------------------

    if args.show_brane_abstract or args.show_brane_helmholtz:
        rhoB = sp.Function(r"\rho_{\rm brane}")(t, x, y, z)
        JxB = sp.Function(r"J^{\rm brane}_x")(t, x, y, z)
        JyB = sp.Function(r"J^{\rm brane}_y")(t, x, y, z)
        JzB = sp.Function(r"J^{\rm brane}_z")(t, x, y, z)
        SB = sp.Function(r"S_{\rm brane}")(t, x, y, z)

        JB = sp.Matrix([JxB, JyB, JzB])
        grad_rhoB = sp.Matrix([sp.diff(rhoB, x), sp.diff(rhoB, y), sp.diff(rhoB, z)])

        div_JB = sp.diff(JxB, x) + sp.diff(JyB, y) + sp.diff(JzB, z)
        curl_JB = sp.Matrix([
            sp.diff(JzB, y) - sp.diff(JyB, z),
            sp.diff(JxB, z) - sp.diff(JzB, x),
            sp.diff(JyB, x) - sp.diff(JxB, y),
        ])

        vBx = JxB / rhoB
        vBy = JyB / rhoB
        vBz = JzB / rhoB
        vB = sp.Matrix([vBx, vBy, vBz])

        div_vB_kin = div_JB / rhoB - (JB.dot(grad_rhoB)) / (rhoB**2)
        omegaB_kin = curl_JB / rhoB - grad_rhoB.cross(JB) / (rhoB**2)

        if args.show_brane_abstract:
            latex_print(
                "8c-abstract) Definitions: v_brane = J_brane / rho_brane (abstract fields)",
                [
                    sp.Eq(sp.Symbol(r"v_x"), vBx),
                    sp.Eq(sp.Symbol(r"v_y"), vBy),
                    sp.Eq(sp.Symbol(r"v_z"), vBz),
                ],
                simplify_mode="none",
                max_ops=max_ops,
            )

            latex_print(
                "8d-abstract) Kinematic identity: ∇·v = (∇·J)/ρ - (J·∇ρ)/ρ^2",
                sp.Eq(sp.Symbol(r"\nabla\cdot v"), div_vB_kin),
                simplify_mode="none",
                max_ops=max_ops,
            )

            latex_print(
                "8e-abstract) Kinematic identity: ω = ∇×v = (∇×J)/ρ - (∇ρ×J)/ρ^2",
                [
                    sp.Eq(sp.Symbol(r"\omega_x"), omegaB_kin[0]),
                    sp.Eq(sp.Symbol(r"\omega_y"), omegaB_kin[1]),
                    sp.Eq(sp.Symbol(r"\omega_z"), omegaB_kin[2]),
                ],
                simplify_mode="none",
                max_ops=max_ops,
            )

            latex_print(
                "8f-abstract) Open-system brane continuity: ρ_t + ∇·J = S_brane",
                sp.Eq(sp.diff(rhoB, t) + div_JB, SB),
                simplify_mode="none",
                max_ops=max_ops,
            )

            div_vB_cont = (SB - sp.diff(rhoB, t)) / rhoB - (JB.dot(grad_rhoB)) / (rhoB**2)
            latex_print(
                "8g-abstract) Combine: ∇·v = (S_brane - ρ_t)/ρ - (J·∇ρ)/ρ^2",
                sp.Eq(sp.Symbol(r"\nabla\cdot v"), div_vB_cont),
                simplify_mode="none",
                max_ops=max_ops,
            )

        if args.show_brane_helmholtz:
            phi3 = sp.Function(r"\phi_3")(t, x, y, z)
            Ax3 = sp.Function(r"A_x")(t, x, y, z)
            Ay3 = sp.Function(r"A_y")(t, x, y, z)
            Az3 = sp.Function(r"A_z")(t, x, y, z)

            lap_phi3 = sp.diff(phi3, x, 2) + sp.diff(phi3, y, 2) + sp.diff(phi3, z, 2)
            lap_Ax3 = sp.diff(Ax3, x, 2) + sp.diff(Ax3, y, 2) + sp.diff(Ax3, z, 2)
            lap_Ay3 = sp.diff(Ay3, x, 2) + sp.diff(Ay3, y, 2) + sp.diff(Ay3, z, 2)
            lap_Az3 = sp.diff(Az3, x, 2) + sp.diff(Az3, y, 2) + sp.diff(Az3, z, 2)

            latex_print(
                "8h) Helmholtz (brane): choose φ3 so v_L=∇φ3 and ∇²φ3 = ∇·v",
                sp.Eq(lap_phi3, div_vB_kin),
                simplify_mode="none",
                max_ops=max_ops,
            )

            latex_print(
                "8i) Helmholtz (brane): in Coulomb gauge ∇·A=0, curl(v)=ω implies ∇²A = -ω",
                [
                    sp.Eq(lap_Ax3, -omegaB_kin[0]),
                    sp.Eq(lap_Ay3, -omegaB_kin[1]),
                    sp.Eq(lap_Az3, -omegaB_kin[2]),
                    sp.Eq(sp.diff(Ax3, x) + sp.diff(Ay3, y) + sp.diff(Az3, z), 0),
                ],
                simplify_mode="none",
                max_ops=max_ops,
            )

    # ----------------------------
    # Abstract transverse/vorticity dynamics diagnostic: curl of brane momentum equation
    # ----------------------------

    if args.show_brane_vorticity_dynamics:
        rhoB = sp.Function(r"\rho_{\rm brane}")(t, x, y, z)
        JxB = sp.Function(r"J^{\rm brane}_x")(t, x, y, z)
        JyB = sp.Function(r"J^{\rm brane}_y")(t, x, y, z)
        JzB = sp.Function(r"J^{\rm brane}_z")(t, x, y, z)

        Pixx = sp.Function(r"\Pi_{xx}")(t, x, y, z)
        Pixy = sp.Function(r"\Pi_{xy}")(t, x, y, z)
        Pixz = sp.Function(r"\Pi_{xz}")(t, x, y, z)
        Piyx = sp.Function(r"\Pi_{yx}")(t, x, y, z)
        Piyy = sp.Function(r"\Pi_{yy}")(t, x, y, z)
        Piyz = sp.Function(r"\Pi_{yz}")(t, x, y, z)
        Pizx = sp.Function(r"\Pi_{zx}")(t, x, y, z)
        Pizy = sp.Function(r"\Pi_{zy}")(t, x, y, z)
        Pizz = sp.Function(r"\Pi_{zz}")(t, x, y, z)

        SJx = sp.Function(r"S^{J}_x")(t, x, y, z)
        SJy = sp.Function(r"S^{J}_y")(t, x, y, z)
        SJz = sp.Function(r"S^{J}_z")(t, x, y, z)

        JB = sp.Matrix([JxB, JyB, JzB])
        Pi = sp.Matrix([[Pixx, Pixy, Pixz], [Piyx, Piyy, Piyz], [Pizx, Pizy, Pizz]])
        SJ = sp.Matrix([SJx, SJy, SJz])

        def div_tensor_3(T):
            return sp.Matrix([
                sp.diff(T[0, 0], x) + sp.diff(T[0, 1], y) + sp.diff(T[0, 2], z),
                sp.diff(T[1, 0], x) + sp.diff(T[1, 1], y) + sp.diff(T[1, 2], z),
                sp.diff(T[2, 0], x) + sp.diff(T[2, 1], y) + sp.diff(T[2, 2], z),
            ])

        def curl_vec_3(V):
            return sp.Matrix([
                sp.diff(V[2], y) - sp.diff(V[1], z),
                sp.diff(V[0], z) - sp.diff(V[2], x),
                sp.diff(V[1], x) - sp.diff(V[0], y),
            ])

        curlJ = curl_vec_3(JB)
        divPi = div_tensor_3(Pi)

        lhs = sp.diff(curlJ, t) + curl_vec_3(divPi)
        rhs = curl_vec_3(SJ)

        latex_print(
            "8j-dyn) Abstract curl(momentum balance): ∂t(∇×J) + ∇×(∇·Π) = ∇×S_J",
            [sp.Eq(lhs[0], rhs[0]), sp.Eq(lhs[1], rhs[1]), sp.Eq(lhs[2], rhs[2])],
            simplify_mode="none",
            max_ops=max_ops,
        )

        vB = JB / rhoB
        omegaB = curl_vec_3(vB)
        identity = sp.simplify(
            curlJ
            - (
                rhoB * omegaB
                + sp.Matrix([sp.diff(rhoB, x), sp.diff(rhoB, y), sp.diff(rhoB, z)]).cross(vB)
            )
        )

        latex_print(
            "8k-dyn) Identity: ∇×J = ρ ω + (∇ρ×v)  (checked algebraically)",
            [sp.Eq(identity[0], 0), sp.Eq(identity[1], 0), sp.Eq(identity[2], 0)],
            simplify_mode="none",
            max_ops=max_ops,
        )

    # ----------------------------
    # Poisson bridge (exact identity + conditional candidate)
    # ----------------------------

    if args.show_poisson_bridge:
        div_vb = (S_rho_brane - sp.diff(rho_brane, t)) / rho_brane - (
            Jx_brane * sp.diff(rho_brane, x)
            + Jy_brane * sp.diff(rho_brane, y)
            + Jz_brane * sp.diff(rho_brane, z)
        ) / (rho_brane**2)

        latex_print(
            "8d) Exact identity for ∇·v_brane from projected continuity",
            sp.Eq(sp.Symbol(r"\nabla\cdot v_{\rm brane}"), div_vb),
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )

        rho0_sym = sp.Symbol(r"\rho_0", positive=True, real=True)
        Sb_sym = sp.Symbol(r"S_{\rm brane}")
        latex_print(
            "8e) Quasi-static candidate (not forced): ∇·v_brane ≈ S_brane/ρ0",
            sp.Eq(sp.Symbol(r"\nabla\cdot v_{\rm brane}"), Sb_sym / rho0_sym),
            simplify_mode="none",
            max_ops=max_ops,
        )

        phi3_bridge = sp.Function(r"\phi_3")(t, x, y, z)
        lap3_phi3_bridge = sp.diff(phi3_bridge, x, 2) + sp.diff(phi3_bridge, y, 2) + sp.diff(phi3_bridge, z, 2)
        latex_print(
            "8f) If v_L=∇φ3 is the longitudinal part: ∇²φ3 = ∇·v_brane ⇒ ρ0 ∇²φ3 ≈ S_brane",
            sp.Eq(rho0_sym * lap3_phi3_bridge, Sb_sym),
            simplify_mode="none",
            max_ops=max_ops,
        )

    # ----------------------------
    # Separable-mode check (diagnostic)
    # ----------------------------

    if args.check_separable_mode:
        Psi3 = sp.Function(r"\Psi_3")(t, x, y, z)
        Psi3b = sp.Function(r"\overline{\Psi_3}")(t, x, y, z)
        chi = sp.Function(r"\chi")(w)
        chib = sp.Function(r"\bar\chi")(w)

        Jx3 = (hbar / (2 * sp.I * m)) * (Psi3b * sp.diff(Psi3, x) - Psi3 * sp.diff(Psi3b, x))
        Jy3 = (hbar / (2 * sp.I * m)) * (Psi3b * sp.diff(Psi3, y) - Psi3 * sp.diff(Psi3b, y))
        Jz3 = (hbar / (2 * sp.I * m)) * (Psi3b * sp.diff(Psi3, z) - Psi3 * sp.diff(Psi3b, z))

        Cw = sp.Integral(Wt * chib * chi, (w, -Wproj, Wproj))

        rho_b_sep = (Psi3b * Psi3) * Cw
        Jx_b_sep = Jx3 * Cw
        Jy_b_sep = Jy3 * Cw
        Jz_b_sep = Jz3 * Cw

        vbx_sep = sp.simplify(Jx_b_sep / rho_b_sep)
        vby_sep = sp.simplify(Jy_b_sep / rho_b_sep)
        vbz_sep = sp.simplify(Jz_b_sep / rho_b_sep)

        latex_print(
            "8g) Separable-mode check: v_brane reduces to 3D velocity (Cw cancels)",
            [vbx_sep, vby_sep, vbz_sep],
            simplify_mode="none",
            max_ops=max_ops,
        )
        text_note(
            "8g-note)",
            [
                "If --no_gauge and θ3 is single-valued, then in this separable regime v_brane is a 3D gradient field,",
                "so ω_brane=0. Any ω_brane in the full model must then come from:",
                "  • gauge field strength B_ij (bulk), and/or",
                "  • multi-mode w-structure + projection (R_ij, S terms), and/or",
                "  • singular/defect structure in θ (vortices), if allowed.",
            ],
        )

    # ----------------------------
    # Brane-projected momentum balance (flux form) and extra stress sector
    # ----------------------------

    if args.show_brane_momentum:
        rhoB = sp.Function(r"\rho_{\rm brane}")(t, x, y, z)
        JxB = sp.Function(r"J_{x,\rm brane}")(t, x, y, z)
        JyB = sp.Function(r"J_{y,\rm brane}")(t, x, y, z)
        JzB = sp.Function(r"J_{z,\rm brane}")(t, x, y, z)

        if (not report_mode) or args.report_with_integrals:
            latex_print(
                "8h0) Brane shorthand fields (definitions)",
                [
                    sp.Eq(rhoB, rho_brane),
                    sp.Eq(JxB, Jx_brane),
                    sp.Eq(JyB, Jy_brane),
                    sp.Eq(JzB, Jz_brane),
                ],
                simplify_mode="none",
                max_ops=max_ops,
            )

        vbx_s = JxB / rhoB
        vby_s = JyB / rhoB
        vbz_s = JzB / rhoB

        Pi_xx = sp.Function(r"\Pi_{xx}")(t, x, y, z)
        Pi_xy = sp.Function(r"\Pi_{xy}")(t, x, y, z)
        Pi_xz = sp.Function(r"\Pi_{xz}")(t, x, y, z)
        Pi_yy = sp.Function(r"\Pi_{yy}")(t, x, y, z)
        Pi_yz = sp.Function(r"\Pi_{yz}")(t, x, y, z)
        Pi_zz = sp.Function(r"\Pi_{zz}")(t, x, y, z)

        if (not report_mode) or args.report_with_integrals:
            Pi_defs = [
                sp.Eq(Pi_xx, sp.Integral(Wt * (Jx * Jx / rho), (w, -Wproj, Wproj))),
                sp.Eq(Pi_xy, sp.Integral(Wt * (Jx * Jy / rho), (w, -Wproj, Wproj))),
                sp.Eq(Pi_xz, sp.Integral(Wt * (Jx * Jz / rho), (w, -Wproj, Wproj))),
                sp.Eq(Pi_yy, sp.Integral(Wt * (Jy * Jy / rho), (w, -Wproj, Wproj))),
                sp.Eq(Pi_yz, sp.Integral(Wt * (Jy * Jz / rho), (w, -Wproj, Wproj))),
                sp.Eq(Pi_zz, sp.Integral(Wt * (Jz * Jz / rho), (w, -Wproj, Wproj))),
            ]
            latex_print(
                "8h1) Brane momentum-flux tensor Π_{ij} definitions (kept compact)",
                Pi_defs,
                simplify_mode="none",
                max_ops=max_ops,
            )

        # Momentum exchange with bulk through w-direction (open-system term)
        F_xw = Jx * Jw / rho
        F_yw = Jy * Jw / rho
        F_zw = Jz * Jw / rho

        S_Jx = -(Wt.subs(w, Wproj) * F_xw.subs(w, Wproj) - Wt.subs(w, -Wproj) * F_xw.subs(w, -Wproj)) + sp.Integral(
            sp.diff(Wt, w) * F_xw, (w, -Wproj, Wproj)
        )
        S_Jy = -(Wt.subs(w, Wproj) * F_yw.subs(w, Wproj) - Wt.subs(w, -Wproj) * F_yw.subs(w, -Wproj)) + sp.Integral(
            sp.diff(Wt, w) * F_yw, (w, -Wproj, Wproj)
        )
        S_Jz = -(Wt.subs(w, Wproj) * F_zw.subs(w, Wproj) - Wt.subs(w, -Wproj) * F_zw.subs(w, -Wproj)) + sp.Integral(
            sp.diff(Wt, w) * F_zw, (w, -Wproj, Wproj)
        )

        latex_print(
            "8h2) Brane momentum source terms S_{J_i} from w-flux (bulk exchange)",
            [S_Jx, S_Jy, S_Jz],
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )

        # Force-density term from confinement + enthalpy + quantum potential
        Qsym = sp.Function("Q")(t, x, y, z, w)
        h_bulk = sp.Rational(5, 4) * K * rho**4
        mu_eff = Vconf + h_bulk + Qsym

        Fpot_x = -sp.Integral(Wt * (rho / m) * sp.diff(mu_eff, x), (w, -Wproj, Wproj))
        Fpot_y = -sp.Integral(Wt * (rho / m) * sp.diff(mu_eff, y), (w, -Wproj, Wproj))
        Fpot_z = -sp.Integral(Wt * (rho / m) * sp.diff(mu_eff, z), (w, -Wproj, Wproj))

        # Gauge/Lorentz force density (optional)
        if use_gauge:
            v4J = [Jx / rho, Jy / rho, Jz / rho, Jw / rho]
            coords4 = [x, y, z, w]
            A4 = [Ax, Ay, Az, Aw]
            E4 = [-(sp.diff(Acomp, t) + sp.diff(A0, c)) for c, Acomp in zip(coords4, A4)]

            def B(i: int, j: int) -> sp.Expr:
                return sp.diff(A4[j], coords4[i]) - sp.diff(A4[i], coords4[j])

            lor_x = E4[0] + sum(v4J[j] * B(0, j) for j in range(4))
            lor_y = E4[1] + sum(v4J[j] * B(1, j) for j in range(4))
            lor_z = E4[2] + sum(v4J[j] * B(2, j) for j in range(4))

            Fem_x = sp.Integral(Wt * (q * rho / m) * lor_x, (w, -Wproj, Wproj))
            Fem_y = sp.Integral(Wt * (q * rho / m) * lor_y, (w, -Wproj, Wproj))
            Fem_z = sp.Integral(Wt * (q * rho / m) * lor_z, (w, -Wproj, Wproj))
        else:
            Fem_x = sp.Integer(0)
            Fem_y = sp.Integer(0)
            Fem_z = sp.Integer(0)

        eq_mom_x = sp.Eq(
            sp.diff(JxB, t) + sp.diff(Pi_xx, x) + sp.diff(Pi_xy, y) + sp.diff(Pi_xz, z),
            Fpot_x + Fem_x + S_Jx,
        )
        eq_mom_y = sp.Eq(
            sp.diff(JyB, t) + sp.diff(Pi_xy, x) + sp.diff(Pi_yy, y) + sp.diff(Pi_yz, z),
            Fpot_y + Fem_y + S_Jy,
        )
        eq_mom_z = sp.Eq(
            sp.diff(JzB, t) + sp.diff(Pi_xz, x) + sp.diff(Pi_yz, y) + sp.diff(Pi_zz, z),
            Fpot_z + Fem_z + S_Jz,
        )

        latex_print(
            "8i) Brane-projected momentum balance (3D components; compact flux form)",
            [eq_mom_x, eq_mom_y, eq_mom_z],
            simplify_mode="none",
            max_ops=max_ops,
        )

        # Extra brane stress sector (Reynolds-type) from unresolved w-structure
        R_xx = sp.Function(r"R_{xx}")(t, x, y, z)
        R_xy = sp.Function(r"R_{xy}")(t, x, y, z)
        R_xz = sp.Function(r"R_{xz}")(t, x, y, z)
        R_yy = sp.Function(r"R_{yy}")(t, x, y, z)
        R_yz = sp.Function(r"R_{yz}")(t, x, y, z)
        R_zz = sp.Function(r"R_{zz}")(t, x, y, z)

        R_defs = [
            sp.Eq(R_xx, Pi_xx - rhoB * vbx_s * vbx_s),
            sp.Eq(R_xy, Pi_xy - rhoB * vbx_s * vby_s),
            sp.Eq(R_xz, Pi_xz - rhoB * vbx_s * vbz_s),
            sp.Eq(R_yy, Pi_yy - rhoB * vby_s * vby_s),
            sp.Eq(R_yz, Pi_yz - rhoB * vby_s * vbz_s),
            sp.Eq(R_zz, Pi_zz - rhoB * vbz_s * vbz_s),
        ]

        latex_print(
            "8j) Extra stress sector from w-projection: R_{ij} ≡ Π_{ij} - ρ_brane v_i v_j",
            R_defs,
            simplify_mode="none",
            max_ops=max_ops,
        )

        text_note(
            "8j-note)",
            [
                "If v_i(t,x,y,z,w) is effectively w-independent (or a single separable mode), then Π_{ij}=ρ_brane v_i v_j and R_{ij}→0.",
                "Otherwise R_{ij} encodes extra brane physics (multi-mode w-structure / unresolved bulk degrees of freedom).",
            ],
        )

    # ----------------------------
    # Brane Euler (acceleration) form and vorticity dynamics (abstract; no closure forced)
    # ----------------------------

    if args.show_brane_euler_form:
        rhoB = sp.Function(r"\rho_{\rm brane}")(t, x, y, z)
        vxB = sp.Function(r"v^{\rm brane}_x")(t, x, y, z)
        vyB = sp.Function(r"v^{\rm brane}_y")(t, x, y, z)
        vzB = sp.Function(r"v^{\rm brane}_z")(t, x, y, z)
        vB = sp.Matrix([vxB, vyB, vzB])

        S_rhoB = sp.Function(r"S_{\rho,\rm brane}")(t, x, y, z)

        Fx_tot = sp.Function(r"F^{\rm tot}_x")(t, x, y, z)
        Fy_tot = sp.Function(r"F^{\rm tot}_y")(t, x, y, z)
        Fz_tot = sp.Function(r"F^{\rm tot}_z")(t, x, y, z)
        Ftot = sp.Matrix([Fx_tot, Fy_tot, Fz_tot])

        Rxx = sp.Function(r"R_{xx}")(t, x, y, z)
        Rxy = sp.Function(r"R_{xy}")(t, x, y, z)
        Rxz = sp.Function(r"R_{xz}")(t, x, y, z)
        Ryx = sp.Function(r"R_{yx}")(t, x, y, z)
        Ryy = sp.Function(r"R_{yy}")(t, x, y, z)
        Ryz = sp.Function(r"R_{yz}")(t, x, y, z)
        Rzx = sp.Function(r"R_{zx}")(t, x, y, z)
        Rzy = sp.Function(r"R_{zy}")(t, x, y, z)
        Rzz = sp.Function(r"R_{zz}")(t, x, y, z)
        R = sp.Matrix([[Rxx, Rxy, Rxz], [Ryx, Ryy, Ryz], [Rzx, Rzy, Rzz]])

        def div_vec_3(V):
            return sp.diff(V[0], x) + sp.diff(V[1], y) + sp.diff(V[2], z)

        def div_tensor_3(T):
            return sp.Matrix([
                sp.diff(T[0, 0], x) + sp.diff(T[0, 1], y) + sp.diff(T[0, 2], z),
                sp.diff(T[1, 0], x) + sp.diff(T[1, 1], y) + sp.diff(T[1, 2], z),
                sp.diff(T[2, 0], x) + sp.diff(T[2, 1], y) + sp.diff(T[2, 2], z),
            ])

        def curl_vec_3(V):
            return sp.Matrix([
                sp.diff(V[2], y) - sp.diff(V[1], z),
                sp.diff(V[0], z) - sp.diff(V[2], x),
                sp.diff(V[1], x) - sp.diff(V[0], y),
            ])

        cont_eq = sp.Eq(sp.diff(rhoB, t) + div_vec_3(rhoB * vB), S_rhoB)

        conv = sp.Matrix([
            vxB * sp.diff(vxB, x) + vyB * sp.diff(vxB, y) + vzB * sp.diff(vxB, z),
            vxB * sp.diff(vyB, x) + vyB * sp.diff(vyB, y) + vzB * sp.diff(vyB, z),
            vxB * sp.diff(vzB, x) + vyB * sp.diff(vzB, y) + vzB * sp.diff(vzB, z),
        ])

        Dt_v = sp.Matrix([sp.diff(vxB, t), sp.diff(vyB, t), sp.diff(vzB, t)]) + conv
        divR = div_tensor_3(R)

        euler_eqs = [
            sp.Eq(rhoB * Dt_v[0], Ftot[0] - divR[0] - vB[0] * S_rhoB),
            sp.Eq(rhoB * Dt_v[1], Ftot[1] - divR[1] - vB[1] * S_rhoB),
            sp.Eq(rhoB * Dt_v[2], Ftot[2] - divR[2] - vB[2] * S_rhoB),
        ]

        latex_print(
            "8l) Brane continuity (velocity form): ρ_t + ∇·(ρ v) = S_ρ",
            cont_eq,
            simplify_mode="none",
            max_ops=max_ops,
        )

        latex_print(
            "8m) Brane Euler (acceleration form): ρ(∂t+v·∇)v = F_tot - ∇·R - v S_ρ  (abstract)",
            euler_eqs,
            simplify_mode="none",
            max_ops=max_ops,
        )

        omega = curl_vec_3(vB)
        lhs_omega = sp.diff(omega, t) - curl_vec_3(vB.cross(omega))
        rhs_omega = curl_vec_3((Ftot - divR - S_rhoB * vB) / rhoB)
        omega_eqs = [sp.Eq(lhs_omega[i], rhs_omega[i]) for i in range(3)]

        latex_print(
            "8n) Brane vorticity equation: ∂tω - ∇×(v×ω) = ∇×[(F_tot-∇·R-v S_ρ)/ρ]",
            omega_eqs,
            simplify_mode="none",
            max_ops=max_ops,
        )

        text_note(
            "8n-note)",
            [
                "If F_tot is a pure gradient and extra sectors vanish (R=0, S_ρ=0), then RHS=0 and ω is advected.",
                "If not, transverse/vorticity dynamics can be sourced by bulk exchange (S_ρ, R) or non-gradient forces (e.g. gauge/Lorentz).",
            ],
        )

    # ----------------------------
    # Emergent EM-analog dictionary from brane velocity (kinematic identities only)
    # ----------------------------

    if args.show_em_analog:
        rhoB = sp.Function(r"\rho_{\rm brane}")(t, x, y, z)
        vxB = sp.Function(r"v^{\rm brane}_x")(t, x, y, z)
        vyB = sp.Function(r"v^{\rm brane}_y")(t, x, y, z)
        vzB = sp.Function(r"v^{\rm brane}_z")(t, x, y, z)
        vB = sp.Matrix([vxB, vyB, vzB])

        def div_vec_3(V):
            return sp.diff(V[0], x) + sp.diff(V[1], y) + sp.diff(V[2], z)

        def grad_vec_3(f):
            return sp.Matrix([sp.diff(f, x), sp.diff(f, y), sp.diff(f, z)])

        def curl_vec_3(V):
            return sp.Matrix([
                sp.diff(V[2], y) - sp.diff(V[1], z),
                sp.diff(V[0], z) - sp.diff(V[2], x),
                sp.diff(V[1], x) - sp.diff(V[0], y),
            ])

        hB = sp.Rational(5, 4) * K * rhoB**4

        lam_em = sp.Symbol(r"\lambda_{\rm em}", real=True)

        Phi_em = lam_em * (hB + sp.Rational(1, 2) * (vB.dot(vB)))
        A_em = lam_em * vB

        E_em = -grad_vec_3(Phi_em) - sp.diff(A_em, t)
        B_em = curl_vec_3(A_em)
        omega = curl_vec_3(vB)

        # Use MatrixSymbol for vector equalities (SymPy sometimes tries to reduce Eq(Matrix,Matrix) to False)
        AemSym = sp.MatrixSymbol(r"A_{\rm em}", 3, 1)
        EemSym = sp.MatrixSymbol(r"E_{\rm em}", 3, 1)
        BemSym = sp.MatrixSymbol(r"B_{\rm em}", 3, 1)

        latex_print(
            "EM1) Emergent-EM kinematics: Φ_em=λ(h+½|v|^2),  A_em=λ v  ⇒  E=-∇Φ-∂tA,  B=∇×A",
            [
                sp.Eq(sp.Symbol(r"\Phi_{\rm em}"), Phi_em),
                sp.Eq(AemSym, A_em),
                sp.Eq(EemSym, E_em),
                sp.Eq(BemSym, B_em),
            ],
            simplify_mode="none",
            max_ops=max_ops,
        )

        divB = div_vec_3(B_em)
        curlE_plus_dBdt = curl_vec_3(E_em) + sp.diff(B_em, t)

        latex_print(
            "EM2) Homogeneous identities: ∇·B=0,  ∇×E + ∂tB = 0 (true by construction)",
            [
                sp.Eq(divB, 0),
                sp.Eq(curlE_plus_dBdt[0], 0),
                sp.Eq(curlE_plus_dBdt[1], 0),
                sp.Eq(curlE_plus_dBdt[2], 0),
            ],
            simplify_mode="none",
            max_ops=max_ops,
        )

        latex_print(
            "EM3) In this mapping, B is proportional to brane vorticity:  B = λ (∇×v) = λ ω",
            [
                sp.Eq(B_em[0] - lam_em * omega[0], 0),
                sp.Eq(B_em[1] - lam_em * omega[1], 0),
                sp.Eq(B_em[2] - lam_em * omega[2], 0),
            ],
            simplify_mode="none",
            max_ops=max_ops,
        )

        text_note(
            "EM-note)",
            [
                "This block does NOT assume inhomogeneous Maxwell equations. It only shows a kinematic route:",
                "if ω_brane≠0, you can define an EM-like potential A_em∝v whose B-field is vorticity.",
                "Any Maxwell-like *dynamics* would require additional closure assumptions from the brane momentum equation.",
            ],
        )

    # ----------------------------
    # Linearized brane acoustics → forced wave equation → static Poisson candidate (conditional layer)
    # ----------------------------

    rho0 = sp.Symbol(r"\rho_0", positive=True, real=True)
    drho3 = sp.Function(r"\delta\rho_{\rm brane}")(t, x, y, z)
    phi3 = sp.Function(r"\phi_3")(t, x, y, z)

    lap3_phi = sp.diff(phi3, x, 2) + sp.diff(phi3, y, 2) + sp.diff(phi3, z, 2)

    u = sp.Function("u")(t, x, y, z)  # effort variable u ≡ δh
    hprime0 = sp.diff(h_rho, rho_sym).subs(rho_sym, rho0)  # = 5K rho0^3
    cs2_0 = sp.diff(P_rho, rho_sym).subs(rho_sym, rho0)  # = 5K rho0^4

    S3 = sp.Function(r"S_{\rm brane}")(t, x, y, z)

    eq_cont_lin = sp.Eq(sp.diff(drho3, t) + rho0 * lap3_phi, S3)
    eq_effort = sp.Eq(u, hprime0 * drho3)
    eq_mom_lin = sp.Eq(sp.diff(phi3, t) + u / m, 0)

    # In report mode we keep this (it's exactly the controlled bridge to Poisson); otherwise we keep as before.
    latex_print(
        "9) Linearized brane acoustic relations (assumptions stated; do not force Poisson)",
        [
            sp.Eq(sp.Derivative(sp.Function("h")(rho_sym), rho_sym).subs(rho_sym, rho0), hprime0),
            sp.Eq(sp.Symbol(r"c_{s,0}^2"), cs2_0),
            sp.Eq(rho0 * sp.Derivative(sp.Function("h")(rho_sym), rho_sym).subs(rho_sym, rho0), sp.simplify(rho0 * hprime0)),
            eq_cont_lin,
            eq_effort,
            eq_mom_lin,
        ],
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )

    wave_eq = sp.Eq(sp.diff(phi3, t, 2) - (cs2_0 / m) * lap3_phi, -(cs2_0 / (m * rho0)) * S3)
    poisson_static = sp.Eq(rho0 * lap3_phi, S3)

    latex_print("10) Forced wave equation for brane velocity potential φ3", wave_eq, simplify_mode=simplify_mode, max_ops=max_ops)
    latex_print(
        "11) Static-limit candidate (if φ_tt is negligible): ρ0 ∇^2 φ3 = S_brane",
        poisson_static,
        simplify_mode=simplify_mode,
        max_ops=max_ops,
    )

    # ----------------------------
    # Bulk Madelung summary (paper-friendly) and optional full expansions
    # ----------------------------

    # Summary variables
    rho4 = sp.Function(r"\rho")(t, x, y, z, w)
    theta4 = sp.Function(r"\theta")(t, x, y, z, w)
    coords4 = [x, y, z, w]
    A4 = [Ax, Ay, Az, Aw]

    if use_gauge:
        v4 = [(hbar * sp.diff(theta4, c) - q * Acomp) / m for c, Acomp in zip(coords4, A4)]
    else:
        v4 = [(hbar / m) * sp.diff(theta4, c) for c in coords4]

    v4x, v4y, v4z, v4w = v4
    v2 = v4x**2 + v4y**2 + v4z**2 + v4w**2
    div4_rhov = sum(sp.diff(rho4 * vi, ci) for vi, ci in zip(v4, coords4))
    lap4_sqrt_rho = sum(sp.diff(sp.sqrt(rho4), ci, 2) for ci in coords4)
    Q = -(hbar**2 / (2 * m)) * (lap4_sqrt_rho / sp.sqrt(rho4))
    h_rho4 = sp.Rational(5, 4) * K * rho4**4

    madelung_cont = sp.Eq(sp.diff(rho4, t) + div4_rhov, 0)
    madelung_phase = sp.Eq(
        hbar * sp.diff(theta4, t) + (q * A0 if use_gauge else 0) + (m / 2) * v2 + Vconf + h_rho4 + Q,
        0,
    )

    # Always print a compact Madelung summary in report mode.
    if report_mode:
        latex_print(
            "12-summary) Madelung dictionary (bulk) and PDEs (compact)",
            [
                sp.Eq(sp.Symbol(r"\psi"), sp.Symbol(r"\sqrt{\rho}\,e^{i\theta}")),
                sp.Eq(sp.Symbol(r"v_i"), (hbar * sp.Symbol(r"\partial_i\theta") - (q * sp.Symbol("A_i") if use_gauge else 0)) / m),
                sp.Eq(sp.Symbol(r"Q"), Q),
                madelung_cont,
                madelung_phase,
            ],
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )

        # Bulk Euler summary (no component expansion)
        if use_gauge:
            E4 = [-(sp.diff(Acomp, t) + sp.diff(A0, c)) for c, Acomp in zip(coords4, A4)]
        else:
            E4 = [sp.Integer(0) for _ in coords4]

        def B(i: int, j: int) -> sp.Expr:
            if not use_gauge:
                return sp.Integer(0)
            return sp.diff(A4[j], coords4[i]) - sp.diff(A4[i], coords4[j])

        # Abstract Euler component equation
        i_sym = sp.Symbol("i")
        j_sym = sp.Symbol("j")
        latex_print(
            "13-summary) Bulk Euler-like equation (structural form)",
            sp.Eq(
                sp.Symbol(r"m(\partial_t+\mathbf v\cdot\nabla_4)v_i"),
                -sp.Symbol(r"\partial_i(V_{conf}+h(\rho)+Q)")
                + q * sp.Symbol(r"(E_i + v_j B_{ij})"),
            ),
            simplify_mode="none",
            max_ops=max_ops,
        )

        if use_gauge:
            latex_print(
                "13b-summary) Bulk vorticity identity from minimal coupling: Ω_ij = -(q/m) B_ij",
                sp.Eq(sp.Symbol(r"\Omega_{ij}"), -(q / m) * sp.Symbol(r"B_{ij}")),
                simplify_mode="none",
                max_ops=max_ops,
            )
        else:
            latex_print(
                "13b-summary) Bulk vorticity identity for --no_gauge: Ω_ij = 0 (potential flow in bulk)",
                sp.Eq(sp.Symbol(r"\Omega_{ij}"), 0),
                simplify_mode="none",
                max_ops=max_ops,
            )

    # Optional: expanded Madelung/euler blocks (same as v12), only if requested.
    want_madelung = bool(args.show_madelung or args.show_euler or args.show_vorticity or args.check_madelung_full)
    if want_madelung:
        psi_M = sp.sqrt(rho4) * sp.exp(sp.I * theta4)
        psib_M = sp.sqrt(rho4) * sp.exp(-sp.I * theta4)

        if args.show_madelung:
            latex_print(
                "12) Madelung dictionary (definition-level)",
                [
                    sp.Eq(sp.Symbol(r"\psi"), psi_M),
                    sp.Eq(sp.Symbol(r"\rho"), rho4),
                    sp.Eq(sp.Symbol(r"\theta"), theta4),
                    sp.Eq(sp.Symbol(r"v_x"), v4x),
                    sp.Eq(sp.Symbol(r"v_w"), v4w),
                ],
                simplify_mode="none",
                max_ops=max_ops,
            )

            subs_m = {psi: psi_M, psib: psib_M}
            Jx_m = _maybe_simplify(Jx.subs(subs_m), mode="light", max_ops=max_ops)
            target_Jx = _maybe_simplify(rho4 * v4x, mode="light", max_ops=max_ops)
            resJx = _maybe_simplify(Jx_m - target_Jx, mode="light", max_ops=max_ops)

            latex_print(
                "13) Madelung current identity (representative): J_x = ρ v_x",
                [
                    sp.Eq(sp.Symbol(r"J_x"), Jx_m),
                    sp.Eq(sp.Symbol(r"\rho v_x"), target_Jx),
                    sp.Eq(sp.Symbol("residual"), resJx),
                ],
                simplify_mode="none",
                max_ops=max_ops,
            )

            latex_print(
                "14) Madelung PDEs from the 4D GNLS (imag + real parts)",
                [
                    sp.Eq(sp.Symbol(r"\text{(imag)}"), madelung_cont.lhs),
                    madelung_cont,
                    sp.Eq(sp.Symbol(r"\text{(real)}"), madelung_phase.lhs),
                    madelung_phase,
                    sp.Eq(sp.Symbol(r"Q"), Q),
                ],
                simplify_mode=simplify_mode,
                max_ops=max_ops,
            )

        if args.check_madelung_full:
            subs_m = {psi: psi_M, psib: psib_M}
            residuals = []
            for Jcomp, vcomp in [(Jx, v4x), (Jy, v4y), (Jz, v4z), (Jw, v4w)]:
                residuals.append(sp.simplify(Jcomp.subs(subs_m) - rho4 * vcomp))
            latex_print(
                "14b) Full Madelung current check: J_i - ρ v_i (should simplify to 0)",
                [sp.Eq(r, 0) for r in residuals],
                simplify_mode="none",
                max_ops=max_ops,
            )

        if args.show_vorticity:
            pairs = []
            for i in range(4):
                for j in range(i + 1, 4):
                    ci, cj = coords4[i], coords4[j]
                    vi, vj = v4[i], v4[j]
                    Ai, Aj = A4[i], A4[j]
                    Omega_ij = sp.diff(vj, ci) - sp.diff(vi, cj)
                    if use_gauge:
                        Bij = sp.diff(Aj, ci) - sp.diff(Ai, cj)
                        pairs.append(sp.Eq(Omega_ij, -(q / m) * Bij))
                    else:
                        pairs.append(sp.Eq(Omega_ij, 0))
            latex_print(
                "15) Vorticity identity (minimal coupling): Ω_ij = ∂i v_j - ∂j v_i",
                pairs,
                simplify_mode=simplify_mode,
                max_ops=max_ops,
            )

        if args.show_euler or (report_mode and args.report_with_components):
            # Expanded Euler component form (can be huge)
            if use_gauge:
                E4 = [-(sp.diff(Acomp, t) + sp.diff(A0, c)) for c, Acomp in zip(coords4, A4)]
            else:
                E4 = [sp.Integer(0) for _ in coords4]

            def B(i: int, j: int) -> sp.Expr:
                if not use_gauge:
                    return sp.Integer(0)
                return sp.diff(A4[j], coords4[i]) - sp.diff(A4[i], coords4[j])

            def convective_vi(i: int) -> sp.Expr:
                vi = v4[i]
                out = sp.diff(vi, t)
                for vj, cj in zip(v4, coords4):
                    out += vj * sp.diff(vi, cj)
                return out

            eqs = []
            eqs_divm = []
            for i, ci in enumerate(coords4):
                grad_i = sp.diff(Vconf + h_rho4 + Q, ci)
                lor_i = E4[i]
                for j in range(4):
                    lor_i += v4[j] * B(i, j)
                eqs.append(sp.Eq(m * convective_vi(i), -grad_i + q * lor_i))
                eqs_divm.append(sp.Eq(convective_vi(i), -(grad_i / m) + (q / m) * lor_i))

            title = "16) Euler-like equation (component form, 4D space)"
            if report_mode and args.report_with_components and not args.show_euler:
                title = "16) Euler-like equation (component form, 4D space) [requested via --report_with_components]"

            latex_print(title, eqs, simplify_mode=simplify_mode, max_ops=max_ops)
            latex_print("16b) Same equation divided by m (acceleration form)", eqs_divm, simplify_mode=simplify_mode, max_ops=max_ops)

    # ----------------------------
    # (Optional) Maxwell sector (compact form only)
    # ----------------------------

    if args.show_maxwell:
        mu0 = sp.Symbol(r"\mu_0", positive=True, real=True)
        xi_gf = sp.Symbol(r"\xi", positive=True, real=True)
        Z = sp.Function("Z")(t, x, y, z, w, a, L)
        mG = sp.Function(r"m_\gamma")(t, x, y, z, w, a, L)

        coords5 = [t, x, y, z, w]
        A5 = [A0, Ax, Ay, Az, Aw]
        eta = [-1, 1, 1, 1, 1]

        def raise_A(M: int):
            return eta[M] * A5[M]

        def F(M: int, N: int):
            return sp.diff(A5[N], coords5[M]) - sp.diff(A5[M], coords5[N])

        def raise_F(M: int, N: int):
            return eta[M] * eta[N] * F(M, N)

        divA = sum(sp.diff(raise_A(M), coords5[M]) for M in range(5))

        J = [
            sp.Function(r"J0")(t, x, y, z, w),
            sp.Function(r"Jx")(t, x, y, z, w),
            sp.Function(r"Jy")(t, x, y, z, w),
            sp.Function(r"Jz")(t, x, y, z, w),
            sp.Function(r"Jw")(t, x, y, z, w),
        ]

        eqs = []
        for N in range(5):
            term1 = sum(sp.diff(Z * raise_F(M, N), coords5[M]) for M in range(5))
            term2 = (mG**2) * raise_A(N)
            term3 = (1 / xi_gf) * eta[N] * sp.diff(divA, coords5[N])
            eqs.append(sp.Eq(term1 + term2 + term3, mu0 * J[N]))

        latex_print(
            "17) 4+1D Maxwell + localization + Lorenz gauge-fixing (compact component form)",
            eqs,
            simplify_mode="none",
            max_ops=max_ops,
        )

    # ----------------------------
    # Geometry-energy derivatives (optional ledger)
    # ----------------------------

    if not report_mode:
        Pvac = sp.Symbol(r"P_{\rm vac}", real=True)
        sigma = sp.Symbol(r"\sigma", real=True)

        V_geom = (4 * sp.pi / 3) * a**3 * L
        A_geom = (4 * sp.pi * a**2) * L + 2 * (4 * sp.pi / 3) * a**3
        E_geom = Pvac * V_geom + sigma * A_geom

        latex_print(
            "18) Frozen geometry energy E_geom(a,L)=Pvac V + σ A (cylindrical approx)",
            E_geom,
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )
        latex_print(
            "19) Geometry derivatives (for force ledger regression checks): ∂aE, ∂LE",
            [sp.diff(E_geom, a), sp.diff(E_geom, L)],
            simplify_mode=simplify_mode,
            max_ops=max_ops,
        )


if __name__ == "__main__":
    main()
