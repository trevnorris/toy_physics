# g0_closure_card_v0_checks.py — verification for g0_closure_card_v0.md
# Recovered from the Codex card-build (session 019f80ac-765b-7213-b455-dce77e79e5a2) heredocs;
# re-run reproduces the card's admissibility/instantiability checks. This is the HEADLINE block
# (stability margin, generalized wave speeds, transverse factorization + zero-mode norm, dimensions).
# TODO before full trust: (a) tautology/faithfulness review vs the card's ACTUAL operators;
#   (b) Mathematica dual-engine (.wl); (c) fold in remaining build-log blocks (wall factorization, bulk Hessian).

"""Static algebraic prerequisite checks for the G0 closure card.

The original headline calculations are retained and extended below.  This file
does not claim to solve the one-body BVP/eigenproblem or the two-body/far-field
problem.  Run normally for the disposition table and checks, or use
``--self-test`` for the primitive-mutation ablation matrix.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import OrderedDict

import sympy as sp


# Part A: every atomic subclaim in the seven gates and every section 10.4 PASS
# claim/Fredholm-nullspace assertion receives exactly one disposition class.
# Class-(2)/(3) rows name the missing solved quantity rather than fabricating a
# card-level check.
DISPOSITION_ROWS = [
    # Section 10.2, gate 1.
    ("10.2 G1", "localized H zero mode is normalizable and N0>0", 1,
     "PASS_LOCALIZED_H_NORM"),
    ("10.2 G1", "localized H physical kinetic norm M4*N0>0", 1,
     "PASS_LOCALIZED_H_PHYSICAL_NORM"),
    # Gate 2.
    ("10.2 G2", "reduced inertia Mh=N0*M4 equals the declared Mh=1 and is positive", 1,
     "PASS_REDUCED_H_INERTIA"),
    ("10.2 G2", "reduced conservative scalar stiffness is positive", 1,
     "PASS_STABILITY (Sylvester/eigenvalue cross-method predicate)"),
    ("10.2 G2", "parent transverse H operator is nonnegative by factorization", 1,
     "PASS_TRANSVERSE_FACTORIZATION"),
    ("10.2 G2", "planar unconstrained wall Hessian is nonnegative", 1,
     "PASS_PLANAR_WALL_FACTORIZATION"),
    ("10.2 G2", "curved/held sleeve-slab has no negative constrained wall mode", 2,
     "missing lowest eigenvalue of the assembled one-body r_B Hessian with S_hold, BCs, and translation quotient"),
    ("10.2 G2", "uniform bulk density Hessian and its Fourier coefficients are positive", 1,
     "PASS_BULK_U_SECOND_DERIVATIVE + PASS_BULK_FOURIER_DENSITY_HESSIAN"),
    ("10.2 G2", "bulk phase constant mode is explicitly quotiented", 1,
     "PASS_PHASE_CONSTANT_MODE_QUOTIENT"),
    ("10.2 G2", "inhomogeneous one-body quantum/bulk operator has no negative mode after BCs/quotient", 2,
     "missing lowest eigenvalue of the assembled one-body (delta rho, delta phi) Hessian about the sourced solution"),
    # Gate 3.
    ("10.2 G3", "bare projected mouth source is nonzero", 1,
     "PASS_BARE_MOUTH_SOURCE"),
    ("10.2 G3", "dressed response monopole is nonzero", 3,
     "missing 1/R monopole coefficient of the solved dressed two-body/far-field h response"),
    # Gate 4.  No-pinning means no bulk gap/Yukawa mass, not an exact globally
    # unpinned constant in the presence of the localized mouth Robin term.
    ("10.2 G4", "reduced h has no k^0 h^2 term and Qs(0,0)=0", 1,
     "PASS_REDUCED_H_MASSLESSNESS"),
    ("10.2 G4", "parent transverse operator has a zero eigenvalue", 1,
     "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE"),
    # Gate 5.
    ("10.2 G5", "Hadamard force cell (1,0) vanishes on the R_w background", 3,
     "missing solved two-body Hadamard derivative F^(1,0) on the R_w background"),
    ("10.2 G5", "Hadamard force cell (0,1) vanishes on the R_w background", 3,
     "missing solved two-body Hadamard derivative F^(0,1) on the R_w background"),
    # Gate 6.
    ("10.2 G6", "global drain/return mass controller identity", 1,
     "PASS_DRAIN_MASS_CONTROLLER"),
    ("10.2 G6", "D_i and R_0 are normalized and Gamma0>0", 1,
     "PASS_DRAIN_KERNEL_NORMALIZATIONS"),
    ("10.2 G6", "global drain/return momentum controller identity", 1,
     "PASS_DRAIN_MOMENTUM_CONTROLLER"),
    ("10.2 G6", "global drain/return energy controller identity", 1,
     "PASS_DRAIN_ENERGY_CONTROLLER"),
    ("10.2 G6", "one-body finite-volume mass residual closes", 2,
     "missing assembled one-body solved rho_i, phi_i and mass flux through C_i"),
    ("10.2 G6", "one-body finite-volume momentum residual closes", 2,
     "missing assembled one-body solved fields and integrated stress/source momentum residual on C_i"),
    ("10.2 G6", "one-body finite-volume energy residual closes", 2,
     "missing assembled one-body epsilon_i, j_epsilon_i and integrated energy residual on C_i"),
    ("10.2 G6", "drain ablation nulls its attributed local flux piece", 2,
     "missing paired Gamma0 and Gamma0=0 assembled one-body solutions and F_flux,i readout"),
    ("10.2 G6", "pair/two-body finite-volume residual closes", 3,
     "missing solved two-center fields with shared I_ret and P_ret constraints"),
    ("10.2 G6", "outer-return ablation nulls its attributed pair/global flux piece", 3,
     "missing nominal and return-ablated two-center solutions and remote-return momentum/energy ledger"),
    # Gate 7.
    ("10.2 G7", "conservative operator has the action-derived symmetric Hessian prerequisite", 1,
     "PASS_CONSERVATIVE_HESSIAN_SYMMETRY"),
    ("10.2 G7", "solved pair conservative force is integrable before U_11 is written", 3,
     "missing pair-force Jacobian/curl from the solved F_var(X1,X2)"),
    ("10.2 G7", "physical pair total retains F_var+F_flux", 3,
     "missing solved pair readouts F_var and F_flux"),
    # Section 10.4 named PASS claims, split into their atomic content.
    ("10.4 PASS_STABILITY", "Mh=1>0", 1, "PASS_REDUCED_H_INERTIA"),
    ("10.4 PASS_STABILITY", "Beff*Kh-Chu^2=7/4>0", 1, "PASS_STABILITY"),
    ("10.4 PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS",
     "det(K-zM) and the two positive generalized roots have the card values", 1,
     "PASS_STABILITY (same predicate; eigenvalues are cross-method evidence)"),
    ("10.4 PASS_TRANSVERSE_FACTORIZATION", "O_perp=A^dagger A", 1,
     "PASS_TRANSVERSE_FACTORIZATION"),
    ("10.4 PASS_TRANSVERSE_FACTORIZATION", "O_perp f0=0", 1,
     "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE"),
    ("10.4 PASS_WALL_FACTORIZATION", "Lchi^(2)/kappa=Achi^dagger Achi and Lchi^(2)r0'=0", 1,
     "PASS_PLANAR_WALL_FACTORIZATION (planar unconstrained kink only)"),
    ("10.4 PASS_POSITIVE_BULK_HESSIAN", "K*, hbar*, mu* are derived from the section 7 witness", 1,
     "PASS_BULK_U_SECOND_DERIVATIVE"),
    ("10.4 PASS_POSITIVE_BULK_HESSIAN", "U''(rho0)=4>0 is differentiated from U", 1,
     "PASS_BULK_U_SECOND_DERIVATIVE"),
    ("10.4 PASS_POSITIVE_BULK_HESSIAN", "U''+hbar^2 k^2/(4m rho0)>0 for every real k and coefficients positive", 1,
     "PASS_BULK_FOURIER_DENSITY_HESSIAN"),
    ("10.4 PASS_POSITIVE_BULK_HESSIAN", "phase constant mode is the k=0 gradient null and is quotiented", 1,
     "PASS_PHASE_CONSTANT_MODE_QUOTIENT"),
    ("10.4 PASS_DIMENSIONAL_HOMOGENEITY", "wall, H, reduced scalar, and mouth action terms have stated dimensions", 1,
     "PASS_DIMENSIONAL_HOMOGENEITY"),
    ("10.4 PASS_DIMENSIONAL_HOMOGENEITY", "mass source uses [dt rho]=[Gamma0]+[Di]", 1,
     "PASS_DIMENSIONAL_HOMOGENEITY"),
    ("10.4 PASS_DIMENSIONAL_HOMOGENEITY", "momentum source has [dt(m rho v)]", 1,
     "PASS_DIMENSIONAL_HOMOGENEITY"),
    ("10.4 PASS_DIMENSIONAL_HOMOGENEITY", "energy source matches both [ec]+[S] and [R0]+[P_ret]", 1,
     "PASS_DIMENSIONAL_HOMOGENEITY"),
    # Section 10.4 Fredholm/nullspace prose.  These are deliberately not
    # promoted by the static algebra above.
    ("10.4 nullspace prose", "the assembled kernel contains only the bulk-phase constant mode in that sector", 2,
     "missing nullspace basis of the assembled one-body bulk phase operator"),
    ("10.4 nullspace prose", "the assembled wall kernel contains only body translations/wall-shape modes", 2,
     "missing nullspace basis of the constrained assembled one-body wall Hessian"),
    ("10.4 nullspace prose", "the assembled scalar kernel contains the physical massless h mode and no extra modes", 2,
     "missing nullspace basis of the assembled one-body (H,u_L,h) operator"),
    ("10.4 nullspace prose", "the assembled u_L kernel contains only rigid longitudinal translation", 2,
     "missing nullspace basis of the assembled one-body u_L operator"),
    ("10.4 nullspace prose", "phase quotient removes the phase constant from the assembled inverse", 2,
     "missing inverse/domain of the quotiented assembled one-body phase operator"),
    ("10.4 nullspace prose", "S_hold and moduli quotient remove body translation/wall modes from the inverse", 2,
     "missing constrained assembled wall inverse with S_hold and moduli conditions"),
    ("10.4 nullspace prose", "IR condition removes rigid u_L while the physical h Green mode is preserved", 2,
     "missing assembled one-body scalar inverse and its h-pole residue"),
    ("10.4 Fredholm prose", "I_ret,P_ret form finite-dimensional constraints compatible with the assembled operator", 2,
     "missing assembled one-body constraint Jacobian/Schur block"),
    ("10.4 Fredholm prose", "controller constraints do not spoil the positive elliptic principal symbol", 2,
     "missing principal symbol of the fully assembled one-body PDE/constraint system"),
    ("10.4 Fredholm prose", "finite-R_out system is fixed-domain elliptic/Fredholm with fixed nullspace", 2,
     "missing assembled operator domain, boundary map, Fredholm index, and solved nullspace"),
    ("10.4 IR prose", "the required pair/far-field IR extrapolation converges", 3,
     "missing two-center solution family and R_out,R_ret far-field limit of the observables"),
]


# Part D: each runnable predicate has one card-primitive mutation.  Values in
# this table are inputs to the derivation, never expected answers, residuals,
# or booleans.
ABLATIONS = OrderedDict([
    ("PASS_STABILITY", ("C_hu: 1/2 -> 2", {"C_hu": sp.Rational(2)})),
    ("PASS_TRANSVERSE_FACTORIZATION", ("O_perp constant potential coefficient: 4 -> 5", {"O_const": sp.Rational(5)})),
    ("PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE", ("f0 sech exponent: 2 -> 1", {"f0_power": sp.Rational(1)})),
    ("PASS_LOCALIZED_H_NORM", ("parent-H inner-product weight: 2 -> 3", {"H_inner_weight": sp.Rational(3)})),
    ("PASS_LOCALIZED_H_PHYSICAL_NORM", ("M4*: 3/160 -> -3/160", {"M4_star": -sp.Rational(3, 160)})),
    ("PASS_REDUCED_H_INERTIA", ("M4*: 3/160 -> 3/320", {"M4_star": sp.Rational(3, 320)})),
    ("PASS_PLANAR_WALL_FACTORIZATION", ("lambda_chi ell^2/kappa: 2 -> 3", {"wall_lambda_ratio": sp.Rational(3)})),
    ("PASS_BULK_U_SECOND_DERIVATIVE", ("U prefactor K/4 -> K/5", {"U_denominator": sp.Rational(5)})),
    ("PASS_BULK_FOURIER_DENSITY_HESSIAN", ("quantum-gradient energy sign: +1 -> -1", {"quantum_gradient_sign": sp.Rational(-1)})),
    ("PASS_PHASE_CONSTANT_MODE_QUOTIENT", ("declared bulk-phase pinning coefficient: 0 -> 1", {"bulk_phase_mass": sp.Rational(1)})),
    ("PASS_REDUCED_H_MASSLESSNESS", ("declared reduced h^2 coefficient: 0 -> 1", {"reduced_h_mass": sp.Rational(1)})),
    ("PASS_BARE_MOUTH_SOURCE", ("parent J_m*: 1/20 -> 0 (reduced g_chih held at 1)", {"J_m_star": sp.Rational(0)})),
    ("PASS_DRAIN_KERNEL_NORMALIZATIONS", ("B_n polynomial power: 2 -> 1", {"bump_power": sp.Rational(1)})),
    ("PASS_DRAIN_MASS_CONTROLLER", ("return mass gain: 1 -> 2", {"leakage_mass_gain": sp.Rational(2)})),
    ("PASS_DRAIN_MOMENTUM_CONTROLLER", ("I_ret controller sign: -1 -> +1", {"I_ret_sign": sp.Rational(1)})),
    ("PASS_DRAIN_ENERGY_CONTROLLER", ("P_ret controller sign: -1 -> +1", {"P_ret_sign": sp.Rational(1)})),
    ("PASS_CONSERVATIVE_HESSIAN_SYMMETRY", ("u_L Euler cross coefficient: C_hu -> 3/4", {"u_equation_cross": sp.Rational(3, 4)})),
    ("PASS_DIMENSIONAL_HOMOGENEITY", ("[Gamma0] time exponent: -1 -> -2", {"Gamma0_time_power": sp.Rational(-2)})),
])


def print_disposition_table() -> None:
    print("DISPOSITION_TABLE (static checks are class 1; missing solved quantities are classes 2/3)")
    print("| section | atomic subclaim | class | check ID or named missing solved quantity |")
    print("|---|---|---:|---|")
    for section, claim, klass, disposition in DISPOSITION_ROWS:
        print(f"| {section} | {claim} | {klass} | {disposition} |")


def _zero(expr) -> bool:
    return sp.simplify(expr) == 0


def _matrix_zero(matrix: sp.MatrixBase) -> bool:
    return all(_zero(entry) for entry in matrix)


def _positive(expr) -> bool:
    result = sp.simplify(expr > 0)
    return result is sp.true or result is True


def runtime_guard(condition: bool, check_id: str, value: object) -> None:
    """Report a named predicate and force a nonzero exit on failure."""
    if not condition:
        print(f"FIRST_FAILURE={check_id}")
        print(f"FAIL {check_id} value={value}")
        raise SystemExit(1)
    print(f"PASS {check_id} value={value}")


def run_suite(mutations: dict[str, sp.Expr] | None = None) -> None:
    """Derive and evaluate every class-(1) predicate from card primitives."""
    mutations = mutations or {}

    def primitive(name: str, card_value):
        return mutations.get(name, card_value)

    # ------------------------------------------------------------------
    # Retained headline block: scalar stability and generalized speeds.
    Aeff = primitive("A_eff", sp.Rational(1))
    Mh = primitive("M_h_declared", sp.Rational(1))
    Beff = primitive("B_eff", sp.Rational(2))
    Kh = primitive("K_h", sp.Rational(1))
    Chu = primitive("C_hu", sp.Rational(1, 2))
    M = sp.diag(Aeff, Mh)
    K = sp.Matrix([[Beff, Chu], [Chu, Kh]])
    z = sp.symbols("z", real=True)
    char = sp.factor((K - z*M).det())
    roots = sp.solve(char, z)
    margin = sp.factor(Beff*Kh - Chu**2)
    expected_char = z**2 - 3*z + sp.Rational(7, 4)
    expected_roots = {(sp.Rational(3) - sp.sqrt(2))/2,
                      (sp.Rational(3) + sp.sqrt(2))/2}
    # Per Part C, Sylvester margin and positive generalized roots are
    # cross-method evidence for one atomic positivity predicate.
    scalar_stability = (
        Beff > 0 and Kh > 0 and margin == sp.Rational(7, 4)
        and _zero(char - expected_char)
        and set(roots) == expected_roots
        and all(_positive(root) for root in roots)
    )
    runtime_guard(bool(scalar_stability), "PASS_STABILITY",
                  {"margin": margin, "det(K-zM)": char,
                   "c_squared": roots,
                   "numeric": [sp.N(root, 12) for root in roots],
                   "cross_method": "Sylvester margin and generalized eigenvalues"})

    # ------------------------------------------------------------------
    # Retained parent transverse factorization, now split into distinct
    # factorization and zero-eigenvalue predicates.
    w, ell = sp.symbols("w ell", positive=True, real=True)
    f0_power = primitive("f0_power", sp.Rational(2))
    f0 = 1/(ell*sp.cosh(w/ell)**f0_power)
    O_const = primitive("O_const", sp.Rational(4))
    V_H = (O_const - 6/sp.cosh(w/ell)**2)/ell**2
    O_f0 = sp.simplify(-sp.diff(f0, w, 2) + V_H*f0)
    q_H = 2*sp.tanh(w/ell)/ell
    factor_residual_H = sp.simplify(q_H**2 - sp.diff(q_H, w) - V_H)
    runtime_guard(_zero(factor_residual_H), "PASS_TRANSVERSE_FACTORIZATION",
                  {"A_dagger_A_minus_O_potential": factor_residual_H,
                   "A": "d/dw + 2*tanh(w/ell)/ell"})
    runtime_guard(_zero(O_f0), "PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE",
                  {"O_perp_f0": O_f0, "f0": f0})

    H_inner_weight = primitive("H_inner_weight", sp.Rational(2))
    N0 = sp.integrate(H_inner_weight*f0**2, (w, -sp.oo, sp.oo))
    expected_N0 = 8/(3*ell)
    runtime_guard(_zero(N0 - expected_N0) and _positive(N0),
                  "PASS_LOCALIZED_H_NORM",
                  {"N0": N0, "expected_card_value": expected_N0})

    ell_star = sp.Rational(1, 20)
    M4_star = primitive("M4_star", sp.Rational(3, 160))
    physical_H_norm = sp.simplify(M4_star*N0.subs(ell, ell_star))
    runtime_guard(_positive(physical_H_norm),
                  "PASS_LOCALIZED_H_PHYSICAL_NORM",
                  {"M4*N0": physical_H_norm, "M4*": M4_star})
    runtime_guard(_zero(physical_H_norm - Mh) and physical_H_norm == 1
                  and _positive(Mh),
                  "PASS_REDUCED_H_INERTIA",
                  {"N0*M4": physical_H_norm, "declared_Mh": Mh,
                   "positive": _positive(Mh)})

    # ------------------------------------------------------------------
    # Part B.1: derive the planar wall Hessian from the second variation of
    # V_chi.  Derive A_chi independently from the centered logistic order
    # parameter 2*r0-1; it is not defined from the Hessian residual.
    x, kappa_chi, r_symbol = sp.symbols(
        "x kappa_chi r_symbol", positive=True, real=True)
    wall_lambda_ratio = primitive("wall_lambda_ratio", sp.Rational(2))
    lambda_chi = wall_lambda_ratio*kappa_chi/ell**2
    V_chi = lambda_chi*r_symbol**2*(1-r_symbol)**2/4
    V_chi_second = sp.diff(V_chi, r_symbol, 2)
    r0 = 1/(1 + sp.exp(-x/ell))
    Lchi_potential_over_kappa = sp.simplify(
        V_chi_second.subs(r_symbol, r0)/kappa_chi)
    q_chi = sp.simplify((2*r0 - 1)/ell)
    Achi_dagger_Achi_potential = sp.simplify(q_chi**2 - sp.diff(q_chi, x))
    wall_factor_residual = sp.simplify(sp.together(
        (Lchi_potential_over_kappa - Achi_dagger_Achi_potential).rewrite(sp.exp)))
    r0_prime = sp.diff(r0, x)
    Achi_on_r0_prime = sp.simplify(sp.together(
        (sp.diff(r0_prime, x) + q_chi*r0_prime).rewrite(sp.exp)))
    Lchi_on_r0_prime = sp.simplify(sp.together((
        -sp.diff(r0_prime, x, 2)
        + Lchi_potential_over_kappa*r0_prime).rewrite(sp.exp)))
    runtime_guard(
        _zero(wall_factor_residual) and _zero(Achi_on_r0_prime)
        and _zero(Lchi_on_r0_prime),
        "PASS_PLANAR_WALL_FACTORIZATION",
        {"second_variation_over_kappa": Lchi_potential_over_kappa,
         "Achi_superpotential_from_kink": q_chi,
         "factor_residual": wall_factor_residual,
         "Achi_r0prime": Achi_on_r0_prime,
         "Lchi_r0prime": Lchi_on_r0_prime,
         "scope": "planar unconstrained kink only"})

    # ------------------------------------------------------------------
    # Part B.2: derive U'', K*, hbar*, and mu* from U and the section 7
    # witness.  The Fourier coefficient is derived separately from the
    # quantum-gradient energy term.
    rho_symbol = sp.symbols("rho_symbol", positive=True, real=True)
    k_wave = sp.symbols("k", real=True)
    rho0 = sp.Rational(1)
    m_bulk = sp.Rational(1)
    cs = sp.Rational(2)
    cs_squared = cs**2
    xi_Q = sp.Rational(1, 20)
    K_bulk = sp.simplify(m_bulk*cs_squared/(5*rho0**4))
    hbar_bulk = sp.simplify(sp.sqrt(2)*m_bulk*cs*xi_Q)
    U_denominator = primitive("U_denominator", sp.Rational(4))
    U_bulk = K_bulk*rho_symbol**5/U_denominator
    mu_from_U = sp.simplify(sp.diff(U_bulk, rho_symbol).subs(rho_symbol, rho0))
    U_second = sp.simplify(sp.diff(U_bulk, rho_symbol, 2).subs(rho_symbol, rho0))
    bulk_witness_ok = (
        K_bulk == sp.Rational(4, 5)
        and hbar_bulk == sp.sqrt(2)/10
        and mu_from_U == 1
        and U_second == 4
        and _positive(U_second)
    )
    runtime_guard(bool(bulk_witness_ok), "PASS_BULK_U_SECOND_DERIVATIVE",
                  {"U(rho)": U_bulk, "K*": K_bulk,
                   "hbar*": hbar_bulk, "mu_from_dU": mu_from_U,
                   "U_second_from_U": U_second})

    quantum_gradient_sign = primitive("quantum_gradient_sign", sp.Rational(1))
    # E_Q^(2)=(1/2)[hbar^2/(4*m*rho0)] k^2 |delta rho|^2.
    quantum_energy_prefactor = (
        quantum_gradient_sign*hbar_bulk**2/(8*m_bulk*rho0))
    fourier_k2_coefficient = sp.simplify(2*quantum_energy_prefactor)
    density_fourier_hessian = sp.simplify(
        U_second + fourier_k2_coefficient*k_wave**2)
    fourier_positive_for_all_real_k = (
        _positive(U_second) and _positive(fourier_k2_coefficient)
    )
    runtime_guard(bool(fourier_positive_for_all_real_k),
                  "PASS_BULK_FOURIER_DENSITY_HESSIAN",
                  {"H_rhorho(k)": density_fourier_hessian,
                   "constant_coefficient": U_second,
                   "k_squared_coefficient": fourier_k2_coefficient,
                   "domain": "all real k"})

    phase_quotient_norm = primitive("phase_quotient_norm", sp.Rational(1))
    phase_k2_coefficient = sp.simplify(hbar_bulk**2*rho0/m_bulk)
    phase_k0_coefficient = primitive("bulk_phase_mass", sp.Rational(0))
    phase_quotient_ok = (
        _positive(phase_k2_coefficient)
        and phase_k0_coefficient == 0
        and phase_quotient_norm == 1
    )
    runtime_guard(bool(phase_quotient_ok),
                  "PASS_PHASE_CONSTANT_MODE_QUOTIENT",
                  {"phase_Hessian": phase_k0_coefficient + phase_k2_coefficient*k_wave**2,
                   "k0": phase_k0_coefficient,
                   "integral_W_IR": phase_quotient_norm})

    # ------------------------------------------------------------------
    # Gate 4: operational no-pinning means no bulk h gap/Yukawa mass.  It
    # does not assert a globally unpinned constant in the presence of the
    # localized positive mouth Robin term.
    omega = sp.symbols("omega", real=True)
    reduced_h_mass = primitive("reduced_h_mass", sp.Rational(0))
    reduced_u_mass = sp.Rational(0)
    Qs = sp.Matrix([
        [Aeff*omega**2 - Beff*k_wave**2 - reduced_u_mass,
         -Chu*k_wave**2],
        [-Chu*k_wave**2,
         Mh*omega**2 - Kh*k_wave**2 - reduced_h_mass],
    ])
    Qs_at_origin = Qs.subs({omega: 0, k_wave: 0})
    runtime_guard(reduced_h_mass == 0 and _matrix_zero(Qs_at_origin),
                  "PASS_REDUCED_H_MASSLESSNESS",
                  {"k0_h2_coefficient": reduced_h_mass,
                   "Qs(0,0)": Qs_at_origin,
                   "operational_scope": "no bulk spectral gap/Yukawa mass"})

    # ------------------------------------------------------------------
    # Bare mouth-source prerequisite, including an independently integrated
    # normalization of the actual annular eta_i.
    mouth_a = sp.Rational(1)
    mouth_ell = sp.Rational(1, 20)
    eta_denominator = 4*sp.pi*((mouth_a + mouth_ell)**3 - mouth_a**3)
    eta_coefficient = sp.Rational(3)/eta_denominator
    eta_integral = sp.simplify(
        eta_coefficient*4*sp.pi*((mouth_a + mouth_ell)**3 - mouth_a**3)/3)
    J_m_star = primitive("J_m_star", sp.Rational(1, 20))
    # Section 3 separately declares the parent J_m, the parent zero-mode
    # projection f0(0), and the assembled reduced coefficient g_chih.
    f0_at_mouth = sp.simplify(f0.subs({w: 0, ell: mouth_ell}))
    g_chih = primitive("g_chih", sp.Rational(1))
    s_i = sp.symbols("s_i", nonzero=True, real=True)
    projected_mouth_source = sp.simplify(
        eta_integral*J_m_star*f0_at_mouth*s_i)
    mouth_ok = (
        eta_integral == 1 and f0_at_mouth == 1/mouth_ell
        and g_chih == J_m_star/mouth_ell
        and g_chih != 0
        and _zero(projected_mouth_source - g_chih*s_i)
    )
    runtime_guard(bool(mouth_ok), "PASS_BARE_MOUTH_SOURCE",
                  {"integral_eta": eta_integral, "f0(0)": f0_at_mouth,
                   "g_chih=J_m*f0(0)=J_m/ell": g_chih,
                   "integral_eta_g_s": projected_mouth_source})

    # ------------------------------------------------------------------
    # Part B.3: independently integrate the actual B_3, B_1, D_i, and R_0
    # normalizations before applying the controller identities.
    sigma, radial, Rret, shell_a = sp.symbols(
        "sigma radial Rret shell_a", positive=True, real=True)
    bump_power = primitive("bump_power", sp.Rational(2))
    B1_coefficient = sp.gamma(sp.Rational(1, 2) + 3)/(
        2*sp.pi**sp.Rational(1, 2)*sigma)
    B1_integral = sp.simplify(
        2*B1_coefficient*sp.integrate(
            (1-radial**2/sigma**2)**bump_power, (radial, 0, sigma)))
    B3_coefficient = sp.gamma(sp.Rational(3, 2) + 3)/(
        2*sp.pi**sp.Rational(3, 2)*sigma**3)
    B3_integral = sp.simplify(
        4*sp.pi*B3_coefficient*sp.integrate(
            radial**2*(1-radial**2/sigma**2)**bump_power,
            (radial, 0, sigma)))
    D_integral = sp.simplify(B3_integral*B1_integral)
    eta_return_integral = sp.simplify(
        (sp.Rational(3)/(4*sp.pi*((Rret+shell_a)**3-Rret**3)))
        * 4*sp.pi*((Rret+shell_a)**3-Rret**3)/3)
    R0_integral = sp.simplify(
        eta_return_integral*sp.Rational(1, 2)*(B1_integral+B1_integral))
    Gamma0 = sp.Rational(1, 1000)
    kernel_ok = (
        B1_integral == 1 and B3_integral == 1
        and D_integral == 1 and R0_integral == 1 and _positive(Gamma0)
    )
    runtime_guard(bool(kernel_ok), "PASS_DRAIN_KERNEL_NORMALIZATIONS",
                  {"integral_B1": B1_integral,
                   "integral_B3": B3_integral,
                   "integral_Di": D_integral,
                   "integral_eta_ret": eta_return_integral,
                   "integral_R0": R0_integral,
                   "Gamma0": Gamma0})

    n_bodies = sp.symbols("N_bodies", positive=True, integer=True)
    leakage_mass_gain = primitive("leakage_mass_gain", sp.Rational(1))
    mass_controller_integral = sp.simplify(
        -n_bodies*Gamma0*D_integral
        + leakage_mass_gain*n_bodies*Gamma0*R0_integral)
    runtime_guard(_zero(mass_controller_integral),
                  "PASS_DRAIN_MASS_CONTROLLER",
                  {"integral(S_drain+S_leakage)": mass_controller_integral})

    # These independent moments represent arbitrary spatially varying v(y):
    # Sum_i integral(D_i v) and integral(R_0 v) are not identified or replaced
    # by a constant velocity.
    VDx, VDy, VDz, VRx, VRy, VRz = sp.symbols(
        "VDx VDy VDz VRx VRy VRz", real=True)
    sum_D_v = sp.Matrix([VDx, VDy, VDz])
    return_v = sp.Matrix([VRx, VRy, VRz])
    momentum_drain = -m_bulk*Gamma0*sum_D_v
    momentum_leak0 = m_bulk*Gamma0*n_bodies*return_v
    I_ret_sign = primitive("I_ret_sign", sp.Rational(-1))
    I_ret = I_ret_sign*(momentum_drain + momentum_leak0)
    momentum_controller_integral = sp.simplify(
        momentum_drain + momentum_leak0 + R0_integral*I_ret)
    runtime_guard(_matrix_zero(momentum_controller_integral),
                  "PASS_DRAIN_MOMENTUM_CONTROLLER",
                  {"integral(f_total)": momentum_controller_integral,
                   "arbitrary_spatial_moments":
                       "Sum_i integral(D_i v) and integral(R0 v) independent"})

    # Likewise, these moments are independent functionals of an arbitrary
    # spatially varying carrier enthalpy e_c(y).
    sum_D_ec, return_ec = sp.symbols("ED_sum ER", real=True)
    energy_drain = -Gamma0*sum_D_ec
    energy_leak0 = Gamma0*n_bodies*return_ec
    P_ret_sign = primitive("P_ret_sign", sp.Rational(-1))
    P_ret = P_ret_sign*(energy_drain + energy_leak0)
    energy_controller_integral = sp.simplify(
        energy_drain + energy_leak0 + R0_integral*P_ret)
    runtime_guard(_zero(energy_controller_integral),
                  "PASS_DRAIN_ENERGY_CONTROLLER",
                  {"integral(w_total)": energy_controller_integral,
                   "arbitrary_spatial_moments":
                       "Sum_i integral(D_i e_c) and integral(R0 e_c) independent"})

    # These are global controller identities only.  They are not the gate-6
    # finite-volume PDE/control-surface residual (class 2/3 above).

    # ------------------------------------------------------------------
    # Conservative-Hessian symmetry as an integrability prerequisite.  The
    # action Hessian and the two independently transcribed displayed Euler
    # operators are compared; the solved pair-force integrability gate is not
    # claimed here.
    grad_u, grad_h = sp.symbols("grad_u grad_h", real=True)
    scalar_energy = (
        Beff*grad_u**2/2 + Kh*grad_h**2/2 + Chu*grad_u*grad_h)
    action_hessian = sp.hessian(scalar_energy, (grad_u, grad_h))
    u_equation_cross = primitive(
        "u_equation_cross", sp.Rational(1, 2))
    h_equation_cross = primitive(
        "h_equation_cross", sp.Rational(1, 2))
    displayed_euler_operator = sp.Matrix([
        [Beff, u_equation_cross], [h_equation_cross, Kh]])
    conservative_symmetry_ok = (
        displayed_euler_operator == action_hessian
        and displayed_euler_operator == displayed_euler_operator.T
    )
    runtime_guard(bool(conservative_symmetry_ok),
                  "PASS_CONSERVATIVE_HESSIAN_SYMMETRY",
                  {"action_Hessian": action_hessian,
                   "displayed_Euler_operator": displayed_euler_operator,
                   "scope": "static integrability prerequisite, not solved pair-force integrability"})

    # ------------------------------------------------------------------
    # Retained dimensional headline, with the Part C vacuous comparisons
    # replaced by independent dimensions of the balance LHS and source-law
    # primitives.
    L = sp.Matrix([1, 0, 0])
    T = sp.Matrix([0, 1, 0])
    Mass = sp.Matrix([0, 0, 1])
    E = 2*L - 2*T + Mass
    zero_dimension = sp.zeros(3, 1)
    H_dimension = -L
    h_dimension = zero_dimension
    u_dimension = L
    Zchi_dimension = Mass - 2*L
    kappa_chi_dimension = Mass - 2*T
    lambda_chi_dimension = Mass - 2*L - 2*T
    M4_dimension = Mass
    K4_dimension = E
    Aeff_dimension = Mass - 3*L
    Mh_dimension = Mass - L
    Beff_dimension = Mass - L - 2*T
    Kh_dimension = Mass + L - 2*T
    Chu_dimension = Mass - 2*T
    Km_dimension = E + 2*L
    Jm_dimension = E + L
    eta_dimension = -3*L
    rho_dimension = -4*L
    velocity_dimension = L - T
    energy_density_dimension = E - 4*L
    Gamma0_time_power = primitive("Gamma0_time_power", sp.Rational(-1))
    Gamma0_dimension = Gamma0_time_power*T
    Di_dimension = -4*L
    R0_dimension = -4*L
    S_from_primitives_dimension = Gamma0_dimension + Di_dimension
    e_c_dimension = E
    I_ret_dimension = Mass + L - 2*T
    P_ret_dimension = E - T

    dimensional_checks = {
        "wall_kin": Zchi_dimension - 2*T,
        "wall_grad": kappa_chi_dimension - 2*L,
        "wall_pot": lambda_chi_dimension,
        "H_kin": M4_dimension + 2*H_dimension - 2*T,
        "H_grad": K4_dimension + 2*(H_dimension - L),
        "H_pot": K4_dimension + 2*H_dimension - 2*L,
        "u_kin": Aeff_dimension + 2*(u_dimension - T),
        "u_stiff": Beff_dimension + 2*(u_dimension - L),
        "h_kin": Mh_dimension - 2*T,
        "h_stiff": Kh_dimension - 2*L,
        "mix": Chu_dimension + (u_dimension-L) + (h_dimension-L),
        "mouth_robin_integrand": eta_dimension + Km_dimension + 2*H_dimension,
        "mouth_source_integrand": eta_dimension + Jm_dimension + H_dimension,
        # Corrected: [dt rho] is independently compared with [Gamma0]+[Di].
        "mass_source": S_from_primitives_dimension,
        "momentum_source": Mass + S_from_primitives_dimension + velocity_dimension,
        # Corrected: [dt epsilon] is checked through both source decompositions.
        "energy_source_from_ec_plus_S": e_c_dimension + S_from_primitives_dimension,
        "energy_source_from_R0_plus_Pret": R0_dimension + P_ret_dimension,
    }
    dimensional_targets = {
        **{name: E-4*L for name in
           ("wall_kin", "wall_grad", "wall_pot", "H_kin", "H_grad", "H_pot")},
        **{name: E-3*L for name in
           ("u_kin", "u_stiff", "h_kin", "h_stiff", "mix",
            "mouth_robin_integrand", "mouth_source_integrand")},
        "mass_source": rho_dimension - T,
        "momentum_source": Mass + rho_dimension + velocity_dimension - T,
        "energy_source_from_ec_plus_S": energy_density_dimension - T,
        "energy_source_from_R0_plus_Pret": energy_density_dimension - T,
    }
    dimensional_residuals = {
        name: tuple(sp.simplify(dimensional_checks[name]-dimensional_targets[name]))
        for name in dimensional_checks
    }
    runtime_guard(all(value == (0, 0, 0)
                      for value in dimensional_residuals.values()),
                  "PASS_DIMENSIONAL_HOMOGENEITY",
                  {"residuals": dimensional_residuals,
                   "mass_comparison": "[dt rho] vs [Gamma0]+[Di]",
                   "energy_comparisons":
                       "[dt epsilon] vs [ec]+[S] and [R0]+[P_ret]"})

    print("SCOPE static algebraic prerequisites only; class-(2)/(3) assembled-solve and far-field claims remain unproved")


def run_mutation_self_test(script_path: str) -> int:
    """Launch one failing child process per predicate and emit the matrix."""
    matrix = []
    protocol_ok = True
    for predicate, (primitive_description, _mutation) in ABLATIONS.items():
        completed = subprocess.run(
            [sys.executable, script_path, "--mutation", predicate,
             "--no-disposition"],
            text=True, capture_output=True, check=False)
        first_failure = None
        for line in completed.stdout.splitlines():
            if line.startswith("FIRST_FAILURE="):
                first_failure = line.split("=", 1)[1]
                break
        row_ok = completed.returncode != 0 and first_failure == predicate
        protocol_ok = protocol_ok and row_ok
        matrix.append((predicate, primitive_description,
                       first_failure or "<none>", completed.returncode,
                       "PASS" if row_ok else "FAIL"))

    print("ABLATION_MATRIX")
    print("| predicate | mutated card-derived primitive | first-failing check ID | exit code | protocol |")
    print("|---|---|---|---:|---|")
    for predicate, primitive_description, first_failure, exit_code, status in matrix:
        print(f"| {predicate} | {primitive_description} | {first_failure} | {exit_code} | {status} |")
    if not protocol_ok:
        print("FAIL MUTATION_SELF_TEST: at least one mutation did not fail first at its own predicate")
        return 2
    print("PASS MUTATION_SELF_TEST all predicates fail first under a card-primitive mutation with nonzero exit")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true",
                        help="run all primitive-mutation child-process ablations")
    parser.add_argument("--mutation", choices=list(ABLATIONS),
                        help=argparse.SUPPRESS)
    parser.add_argument("--no-disposition", action="store_true",
                        help=argparse.SUPPRESS)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.no_disposition:
        print_disposition_table()

    mutation_values = ABLATIONS[args.mutation][1] if args.mutation else None
    run_suite(mutation_values)
    if args.mutation:
        # A mutation reaching here failed to ablate its predicate.
        print(f"MUTATION_DID_NOT_FAIL={args.mutation}")
        return 2
    if args.self_test:
        return run_mutation_self_test(os.path.abspath(__file__))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
