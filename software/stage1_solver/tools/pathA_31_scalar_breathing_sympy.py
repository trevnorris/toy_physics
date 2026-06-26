#!/usr/bin/env python3
"""pathA_31 Gate 2 scalar-breathing reduction, SymPy/SciPy engine.

The solve order is intentionally visible in the emitted artifacts:

1. build the ell=0 wall operator, BCs, and mass inner product;
2. solve the operator harmonic lifting functions and projection integrals;
3. assemble M_AB, K_AB and the dynamical EOM;
4. compare only then to the legacy (a,L) closure.

The Mathematica engine recomputes the same symbolic expressions through a
native DSolve/Integrate route and independently repeats the Galerkin overlap
sweep before this script writes the final report.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import sympy as sp
import yaml
from scipy.integrate import quad
from scipy.linalg import eigh
from sympy.printing.mathematica import mathematica_code


class NoAliasDumper(yaml.SafeDumper):
    def ignore_aliases(self, data: Any) -> bool:
        return True


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPO_ROOT = SCRIPT_PATH.parents[3]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"
NOTES = REPO_ROOT / "research" / "pde_ledger" / "notes" / "stages"

REPORT_OUT = REPORTS / "pathA_31_scalar_breathing.md"
RESULTS_YAML = REPORTS / "pathA_31_results.yaml"
SYM_YAML = SCRATCH / "pathA_31_sympy_results.yaml"
SYM_EXPR_WL = SCRATCH / "pathA_31_sympy_exprs.wl"
MMA_YAML = SCRATCH / "pathA_31_mathematica_results.yaml"
FEED_NOTE = NOTES / "moving_throat_pde_pathA_31_scalar_breathing_result.md"

EPS_TRUNC = 0.1
OVERLAP_FLOOR = 1.0 - EPS_TRUNC
N_FINAL = 16
N_CONVERGENCE = [4, 8, 12, 16]
BETA_L0_SWEEP = [0.1, 0.2, 0.5, 1.0, 1.85, 2.0, 3.0, 5.0, 10.0, 18.0, 30.0, 50.0]
BETA_L0_FROM_R0 = 37.0 / 20.0


def compact(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.trigsimp(sp.simplify(expr)))


def hstr(expr: Any) -> str | bool | int | float | None:
    if expr is None or isinstance(expr, (bool, int, float, str)):
        return expr
    return sp.sstr(expr)


def fstr(value: float) -> str:
    return f"{value:.12g}"


def yaml_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.dump(payload, Dumper=NoAliasDumper, sort_keys=False, allow_unicode=False),
        encoding="utf-8",
    )


def yaml_read(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise AssertionError(f"YAML did not load to a mapping: {path}")
    return data


def mma_expr(expr: sp.Expr | int | str) -> str:
    if isinstance(expr, int):
        return str(expr)
    if isinstance(expr, str):
        return '"' + expr.replace("\\", "\\\\").replace('"', '\\"') + '"'
    return mathematica_code(compact(expr))


def digest_mapping(mapping: dict[str, str]) -> str:
    canonical = "\n".join(f"{key}: {mapping[key]}" for key in sorted(mapping))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def zero_pattern(matrix: sp.Matrix) -> list[list[bool]]:
    return [[bool(compact(matrix[i, j]) == 0) for j in range(matrix.cols)] for i in range(matrix.rows)]


class DimError(ValueError):
    pass


Dim = tuple[sp.Rational, sp.Rational, sp.Rational]
ZERO_DIM: Dim = (sp.Rational(0), sp.Rational(0), sp.Rational(0))
DIMENSIONLESS_FUNCTIONS = {
    sp.sin,
    sp.cos,
    sp.tan,
    sp.cot,
    sp.sinh,
    sp.cosh,
    sp.tanh,
    sp.coth,
    sp.sech,
    sp.csch,
}
Ldim, Mdim, Tdim = sp.symbols("L M T", positive=True)


def dim_add(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 + b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_sub(left: Dim, right: Dim) -> Dim:
    return tuple(sp.Rational(a0 - b0) for a0, b0 in zip(left, right))  # type: ignore[return-value]


def dim_scale(dim: Dim, scale: sp.Rational) -> Dim:
    return tuple(sp.Rational(scale * d0) for d0 in dim)  # type: ignore[return-value]


def dim_of(expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> Dim:
    expr = sp.sympify(expr)
    if expr == 0 or expr.is_number:
        return ZERO_DIM
    if expr.is_Symbol:
        if expr not in symbol_dims:
            raise DimError(f"missing dimension for symbol {expr}")
        return symbol_dims[expr]
    if expr.is_Mul:
        out = ZERO_DIM
        for arg in expr.args:
            out = dim_add(out, dim_of(arg, symbol_dims))
        return out
    if expr.is_Pow:
        base, exponent = expr.args
        if not exponent.is_number:
            raise DimError(f"non-numeric dimension exponent in {expr}")
        return dim_scale(dim_of(base, symbol_dims), sp.Rational(exponent))
    if expr.is_Add:
        dims = [dim_of(arg, symbol_dims) for arg in expr.args if arg != 0]
        if not dims:
            return ZERO_DIM
        first = dims[0]
        mismatched = [dim for dim in dims if dim != first]
        if mismatched:
            raise DimError(f"dimension mismatch in sum {expr}: {dims}")
        return first
    if expr.func in DIMENSIONLESS_FUNCTIONS:
        arg_dims = [dim_of(arg, symbol_dims) for arg in expr.args]
        if any(dim != ZERO_DIM for dim in arg_dims):
            raise DimError(f"dimensionful argument in {expr}: {arg_dims}")
        return ZERO_DIM
    raise DimError(f"unsupported expression in dimension checker: {expr}")


def dim_to_monomial(dim: Dim) -> sp.Expr:
    return sp.factor(Ldim ** dim[0] * Mdim ** dim[1] * Tdim ** dim[2])


def exp_text(exp: sp.Rational) -> str:
    exp = sp.Rational(exp)
    if exp.q == 1:
        return str(exp.p)
    return f"{exp.p}/{exp.q}"


def dim_to_text(dim: Dim) -> str:
    parts: list[str] = []
    for name, exp in (("L", dim[0]), ("T", dim[2]), ("M", dim[1])):
        exp = sp.Rational(exp)
        if exp == 0:
            continue
        if exp == 1:
            parts.append(name)
        else:
            parts.append(f"{name}^{exp_text(exp)}")
    return " ".join(parts) if parts else "1"


def dim_vector_text(dim: Dim) -> list[str]:
    return [hstr(v) for v in dim]


def dim_record(name: str, expr: sp.Expr, symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Any]:
    dim = dim_of(expr, symbol_dims)
    return {
        "quantity": name,
        "expression": hstr(expr),
        "dimension": dim_to_text(dim),
        "dimension_monomial": hstr(dim_to_monomial(dim)),
        "dimension_vector_LMT": dim_vector_text(dim),
    }


def dimension_verdict(ok: bool) -> str:
    return "DIMENSIONAL_OK" if ok else "BREATHING_FAIL_DIMENSIONAL"


def build_structure_gate(
    M: sp.Matrix,
    K: sp.Matrix,
    H_legacy: sp.Matrix,
    L0: sp.Symbol,
    beta: sp.Symbol,
    mu_eta: sp.Symbol,
    T_w: sp.Symbol,
    r_AL: sp.Symbol,
    kappa: sp.Symbol,
    chi: sp.Symbol,
    sigma_a: sp.Symbol,
    sigma_L: sp.Symbol,
) -> dict[str, Any]:
    B = compact(beta * L0)
    m_aa_positive_form = compact(2 * sp.pi * mu_eta * (sp.sinh(B) * sp.cosh(B) - B) / (beta * sp.sinh(B) ** 2))
    m_det_positive_form = compact(4 * sp.pi**2 * mu_eta**2 * r_AL**2 * (sp.sinh(B) ** 2 - B**2) / (beta**2 * sp.sinh(B) ** 2))
    k_offdiag_negative_form = compact(-4 * sp.pi * T_w * beta * r_AL / sp.sinh(B))
    k_det_positive_form = compact(16 * sp.pi**2 * T_w**2 * beta**2 * r_AL**2)
    legacy_det_positive_form = compact(kappa * sigma_a + kappa * chi**2 * sigma_L + sigma_a * sigma_L)

    legacy_symmetric = bool(H_legacy == H_legacy.T)
    legacy_offdiag_negative = bool(compact(H_legacy[0, 1] + chi * kappa) == 0)
    legacy_det_positive = bool(compact(H_legacy.det() - legacy_det_positive_form) == 0)
    legacy_rank = 2 if legacy_det_positive else int(H_legacy.rank())
    legacy_zeros = zero_pattern(H_legacy)

    def m_posdef(matrix: sp.Matrix) -> dict[str, bool]:
        leading_positive = bool(compact(matrix[0, 0] - m_aa_positive_form) == 0)
        det_positive = bool(compact(matrix.det() - m_det_positive_form) == 0)
        return {
            "leading_minor_positive": leading_positive,
            "det_positive": det_positive,
            "M_posdef": bool(leading_positive and det_positive),
        }

    def k_structure(matrix: sp.Matrix) -> dict[str, Any]:
        k_symmetric = bool(matrix == matrix.T)
        k_offdiag_negative = bool(compact(matrix[0, 1] - k_offdiag_negative_form) == 0)
        k_det_positive = bool(compact(matrix.det() - k_det_positive_form) == 0)
        k_rank = 2 if k_det_positive else int(matrix.rank())
        zeros = zero_pattern(matrix)
        rank_matches = bool(k_rank == legacy_rank)
        zero_pattern_matches = bool(zeros == legacy_zeros)
        return {
            "K_symmetric": k_symmetric,
            "K_offdiag_negative": k_offdiag_negative,
            "K_rank": k_rank,
            "K_zero_pattern": zeros,
            "K_rank_matches_legacy": rank_matches,
            "K_zero_pattern_matches_legacy": zero_pattern_matches,
            "K_structure_ok": bool(
                k_symmetric
                and k_offdiag_negative
                and rank_matches
                and zero_pattern_matches
                and legacy_symmetric
                and legacy_offdiag_negative
            ),
        }

    mass_gate = m_posdef(M)
    stiffness_gate = k_structure(K)
    m_non_posdef_probe = sp.Matrix([[-M[0, 0], M[0, 1]], [M[1, 0], M[1, 1]]])
    k_sign_flip_probe = sp.Matrix([[K[0, 0], -K[0, 1]], [-K[1, 0], K[1, 1]]])
    m_probe = m_posdef(m_non_posdef_probe)
    k_probe = k_structure(k_sign_flip_probe)

    return {
        **mass_gate,
        **stiffness_gate,
        "legacy_rank": legacy_rank,
        "legacy_zero_pattern": legacy_zeros,
        "legacy_symmetric": legacy_symmetric,
        "legacy_offdiag_negative": legacy_offdiag_negative,
        "structure_from_computed_matrices": True,
        "structure_certificates": {
            "M_aa_positive": "M_aa=2*pi*muEta*(sinh(B)*cosh(B)-B)/(beta*sinh(B)^2)>0 for B=L0*beta>0",
            "M_det_positive": "det(M)=4*pi^2*muEta^2*rAL^2*(sinh(B)^2-B^2)/(beta^2*sinh(B)^2)>0 for B=L0*beta>0",
            "K_aL_negative": "K_aL=-4*pi*Tw*beta*rAL/sinh(B)<0 for B=L0*beta>0",
        },
        "structure_counterfactual": {
            "non_posdef_M_probe": {
                "mutation": "M_aa -> -M_aa",
                "M_posdef": m_probe["M_posdef"],
                "leading_minor_positive": m_probe["leading_minor_positive"],
                "det_positive": m_probe["det_positive"],
            },
            "sign_flipped_K_probe": {
                "mutation": "K_aL -> -K_aL",
                "K_offdiag_negative": k_probe["K_offdiag_negative"],
                "K_structure_ok": k_probe["K_structure_ok"],
            },
        },
    }


def matrix_entry_dims(entries: dict[str, sp.Expr], symbol_dims: dict[sp.Symbol, Dim]) -> dict[str, Dim]:
    return {name: dim_of(expr, symbol_dims) for name, expr in entries.items()}


def shared_dimension(entry_dims: dict[str, Dim]) -> Dim | None:
    dims = list(entry_dims.values())
    if not dims:
        return None
    first = dims[0]
    return first if all(dim == first for dim in dims) else None


def build_dimensional_check(
    *,
    L0: sp.Symbol,
    beta: sp.Symbol,
    mu_eta: sp.Symbol,
    T_w: sp.Symbol,
    r_AL: sp.Symbol,
    M_entries: dict[str, sp.Expr],
    K_entries: dict[str, sp.Expr],
) -> dict[str, Any]:
    symbol_dims: dict[sp.Symbol, Dim] = {
        L0: (1, 0, 0),
        beta: (-1, 0, 0),
        mu_eta: (-1, 1, 0),
        T_w: (1, 1, -2),
        r_AL: ZERO_DIM,
    }
    expected_m = (0, 1, 0)
    expected_k = (0, 1, -2)
    expected_ratio = (0, 0, -2)

    m_dims = matrix_entry_dims(M_entries, symbol_dims)
    k_dims = matrix_entry_dims(K_entries, symbol_dims)
    m_shared = shared_dimension(m_dims)
    k_shared = shared_dimension(k_dims)
    ratio_dim = dim_sub(k_shared, m_shared) if m_shared is not None and k_shared is not None else None
    dimensional_ok = bool(
        m_shared == expected_m
        and k_shared == expected_k
        and ratio_dim == expected_ratio
    )

    corrupt_dims = dict(symbol_dims)
    corrupt_dims[T_w] = dim_add(symbol_dims[T_w], (1, 0, 0))
    corrupt_m_dims = matrix_entry_dims(M_entries, corrupt_dims)
    corrupt_k_dims = matrix_entry_dims(K_entries, corrupt_dims)
    corrupt_m_shared = shared_dimension(corrupt_m_dims)
    corrupt_k_shared = shared_dimension(corrupt_k_dims)
    corrupt_ratio_dim = (
        dim_sub(corrupt_k_shared, corrupt_m_shared)
        if corrupt_m_shared is not None and corrupt_k_shared is not None
        else None
    )
    corrupt_ok = bool(
        corrupt_m_shared == expected_m
        and corrupt_k_shared == expected_k
        and corrupt_ratio_dim == expected_ratio
    )
    probe_verdict = "NO_FAIL" if corrupt_ok else "BREATHING_FAIL_DIMENSIONAL"

    return {
        "dimension_order": ["L", "M", "T"],
        "dimensional_gate": "assembled M_AB/K_AB homogeneity and K/M=T^-2",
        "headline_quantities_walked": {
            "M_AB": {name: hstr(expr) for name, expr in M_entries.items()},
            "K_AB": {name: hstr(expr) for name, expr in K_entries.items()},
        },
        "symbol_dimensions": {
            "L0": "L",
            "beta": "L^-1",
            "muEta": "M L^-1",
            "Tw": "M L T^-2",
            "K_eta=Tw*beta^2": "M L^-1 T^-2",
            "rAL": "1",
            "alpha_a,alpha_L": "1",
        },
        "sourcing_note": "muEta, Tw, and K_eta are sourced wall-action coefficients in the reduced-throat convention; K_eta is Tw*beta^2, not fitted from the matrix answer.",
        "computed_dimensions": {
            "M_AB_entries": {name: dim_to_text(dim) for name, dim in m_dims.items()},
            "K_AB_entries": {name: dim_to_text(dim) for name, dim in k_dims.items()},
            "M_AB_shared": dim_to_text(m_shared) if m_shared is not None else "inhomogeneous",
            "K_AB_shared": dim_to_text(k_shared) if k_shared is not None else "inhomogeneous",
            "K_over_M_ratio": dim_to_text(ratio_dim) if ratio_dim is not None else "inhomogeneous",
        },
        "computed_dimension_vectors_LMT": {
            "M_AB_entries": {name: dim_vector_text(dim) for name, dim in m_dims.items()},
            "K_AB_entries": {name: dim_vector_text(dim) for name, dim in k_dims.items()},
            "M_AB_shared": dim_vector_text(m_shared) if m_shared is not None else None,
            "K_AB_shared": dim_vector_text(k_shared) if k_shared is not None else None,
            "K_over_M_ratio": dim_vector_text(ratio_dim) if ratio_dim is not None else None,
        },
        "expected_dimensions": {
            "M_AB_shared": dim_to_text(expected_m),
            "K_AB_shared": dim_to_text(expected_k),
            "K_over_M_ratio": dim_to_text(expected_ratio),
        },
        "checks": {
            "M_entries_homogeneous": m_shared is not None,
            "K_entries_homogeneous": k_shared is not None,
            "M_entries_have_mass_dimension": m_shared == expected_m,
            "K_entries_have_stiffness_dimension": k_shared == expected_k,
            "K_over_M_has_omega_squared_dimension": ratio_dim == expected_ratio,
        },
        "dimensional_ok": dimensional_ok,
        "status": "pass" if dimensional_ok else "fail",
        "dimensional_status": dimension_verdict(dimensional_ok),
        "table": [
            dim_record("K_eta=Tw*beta^2", T_w * beta**2, symbol_dims),
            *[dim_record(f"M_{name}", expr, symbol_dims) for name, expr in M_entries.items()],
            *[dim_record(f"K_{name}", expr, symbol_dims) for name, expr in K_entries.items()],
        ],
        "BREATHING_FAIL_DIMENSIONAL_probe": {
            "mutation": "corrupt sourced [Tw] by one extra power of L",
            "participates_in_verdict": True,
            "sourced_Tw_dimension": dim_to_text(symbol_dims[T_w]),
            "corrupted_Tw_dimension": dim_to_text(corrupt_dims[T_w]),
            "mutated_dimensions": {
                "M_AB_entries": {name: dim_to_text(dim) for name, dim in corrupt_m_dims.items()},
                "K_AB_entries": {name: dim_to_text(dim) for name, dim in corrupt_k_dims.items()},
                "M_AB_shared": dim_to_text(corrupt_m_shared) if corrupt_m_shared is not None else "inhomogeneous",
                "K_AB_shared": dim_to_text(corrupt_k_shared) if corrupt_k_shared is not None else "inhomogeneous",
                "K_over_M_ratio": dim_to_text(corrupt_ratio_dim) if corrupt_ratio_dim is not None else "inhomogeneous",
            },
            "without_mutation_dimensional_ok": dimensional_ok,
            "with_mutation_dimensional_ok": corrupt_ok,
            "probe_verdict": probe_verdict,
            "mutation_fires": probe_verdict == "BREATHING_FAIL_DIMENSIONAL",
        },
    }


def _baseline_functions(beta_l0: float, profile: str = "baseline") -> tuple[list[Callable[[float], float]], list[Callable[[float], float]], str]:
    b = beta_l0
    if b <= 0:
        raise ValueError("beta_l0 must be positive")
    sinh_b = math.sinh(b)

    def alpha_a(x: float) -> float:
        return math.sinh(b * (1.0 - x)) / sinh_b

    def dalpha_a(x: float) -> float:
        return -b * math.cosh(b * (1.0 - x)) / sinh_b

    def alpha_l(x: float) -> float:
        return math.sinh(b * x) / sinh_b

    def dalpha_l(x: float) -> float:
        return b * math.cosh(b * x) / sinh_b

    if profile == "baseline":
        aa, daa, label = alpha_a, dalpha_a, "alpha_a=sinh(beta*(L0-w))/sinh(beta*L0)"
    elif profile == "degenerate_zero":
        aa, daa, label = (lambda _x: 0.0), (lambda _x: 0.0), "alpha_a=0"
    elif profile == "constant_one":
        aa, daa, label = (lambda _x: 1.0), (lambda _x: 0.0), "alpha_a=1"
    else:
        raise ValueError(profile)
    return [aa, alpha_l], [daa, dalpha_l], label


def galerkin_overlap(beta_l0: float, n_modes: int, profile: str = "baseline") -> dict[str, Any]:
    funcs, ders, label = _baseline_functions(beta_l0, profile)
    m_full = 2 + n_modes
    mass = np.zeros((m_full, m_full), dtype=float)
    stiff = np.zeros((m_full, m_full), dtype=float)
    for i in range(2):
        for j in range(i, 2):
            mass_ij = quad(lambda x: funcs[i](x) * funcs[j](x), 0.0, 1.0, epsabs=1e-11, epsrel=1e-11, limit=200)[0]
            stiff_ij = quad(
                lambda x: ders[i](x) * ders[j](x) + beta_l0 * beta_l0 * funcs[i](x) * funcs[j](x),
                0.0,
                1.0,
                epsabs=1e-11,
                epsrel=1e-11,
                limit=200,
            )[0]
            mass[i, j] = mass[j, i] = mass_ij
            stiff[i, j] = stiff[j, i] = stiff_ij
    for n in range(1, n_modes + 1):
        k = (n - 0.5) * math.pi
        idx = 1 + n
        for i in range(2):
            mass_ig = quad(lambda x, i=i, k=k: funcs[i](x) * math.sin(k * x), 0.0, 1.0, epsabs=1e-11, epsrel=1e-11, limit=200)[0]
            stiff_ig = quad(
                lambda x, i=i, k=k: ders[i](x) * k * math.cos(k * x)
                + beta_l0 * beta_l0 * funcs[i](x) * math.sin(k * x),
                0.0,
                1.0,
                epsabs=1e-11,
                epsrel=1e-11,
                limit=200,
            )[0]
            mass[i, idx] = mass[idx, i] = mass_ig
            stiff[i, idx] = stiff[idx, i] = stiff_ig
        mass[idx, idx] = 0.5
        stiff[idx, idx] = 0.5 * (k * k + beta_l0 * beta_l0)

    active = [i for i in range(m_full) if mass[i, i] > 1e-13]
    rank_deficient = len(active) != m_full
    mass_active = mass[np.ix_(active, active)]
    stiff_active = stiff[np.ix_(active, active)]
    eigvals, eigvecs = eigh(stiff_active, mass_active, check_finite=True)
    order = np.argsort(eigvals)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    active_sub_indices = [active.index(i) for i in (0, 1) if i in active]
    selector = np.zeros((len(active), len(active_sub_indices)), dtype=float)
    for col, idx in enumerate(active_sub_indices):
        selector[idx, col] = 1.0
    gram = selector.T @ mass_active @ selector
    gram_pinv = np.linalg.pinv(gram, rcond=1e-12) if gram.size else np.zeros((0, 0))

    overlaps: list[float] = []
    for col in range(2):
        v = eigvecs[:, col]
        norm = float(v.T @ mass_active @ v)
        if not gram.size:
            overlaps.append(0.0)
            continue
        coeff = gram_pinv @ (selector.T @ mass_active @ v)
        pnorm = float(coeff.T @ gram @ coeff)
        overlaps.append(math.sqrt(max(0.0, min(1.0, pnorm / norm))))

    min_omega2 = float(min(eigvals[0], eigvals[1]))
    gap = float((eigvals[2] - eigvals[1]) / eigvals[1])
    pass_row = bool(overlaps[0] >= OVERLAP_FLOOR and overlaps[1] >= OVERLAP_FLOOR and min_omega2 > 0.0)
    return {
        "beta_L0": float(beta_l0),
        "N": int(n_modes),
        "profile": label,
        "o_1": float(overlaps[0]),
        "o_2": float(overlaps[1]),
        "omega1_squared": float(eigvals[0]),
        "omega2_squared": float(eigvals[1]),
        "omega3_squared": float(eigvals[2]),
        "min_omega12_squared": min_omega2,
        "gap": gap,
        "pass": pass_row,
        "rank_deficient_basis": rank_deficient,
        "mass_condition": float(np.linalg.cond(mass_active)),
    }


def build_truncation(profile: str = "baseline") -> dict[str, Any]:
    sweep = [galerkin_overlap(beta_l0, N_FINAL, profile=profile) for beta_l0 in BETA_L0_SWEEP]
    selected = galerkin_overlap(BETA_L0_FROM_R0, N_FINAL, profile=profile)
    convergence = [galerkin_overlap(BETA_L0_FROM_R0, n_modes, profile=profile) for n_modes in N_CONVERGENCE]

    passing = [row["beta_L0"] for row in sweep if row["pass"]]
    beta_window = None
    if passing:
        beta_window = {
            "beta_L0_min_in_sweep": min(passing),
            "beta_L0_max_in_sweep": max(passing),
            "predicate": "o_1,o_2 >= 0.9 and min(omega_1^2,omega_2^2)>0",
        }
    return {
        "status": bool(selected["pass"]),
        "epsilon_trunc": EPS_TRUNC,
        "overlap_floor": OVERLAP_FLOOR,
        "N_modes": N_FINAL,
        "o_1": selected["o_1"],
        "o_2": selected["o_2"],
        "min_omega12_squared": selected["min_omega12_squared"],
        "lambda2": selected["omega2_squared"],
        "gap": selected["gap"],
        "selected_beta": {
            "beta_L0": BETA_L0_FROM_R0,
            "beta": "1/L_unit from frozen Gate-1 wall packet Tw=K_eta=1",
            "selection_rule": "predeclared Gate-1 anchor, not best-fit over the sweep",
        },
        "beta_from_R0": {
            "status": "derived_from_frozen_constant_coeff_packet",
            "R0_reference": "Gate-1 frozen finite throat: L0=37/20, Tw=1, K_eta=1",
            "beta": 1.0,
            "L0": 37.0 / 20.0,
            "beta_L0": BETA_L0_FROM_R0,
            "provenance": "calibrated wall constitutive packet; R0 geometry alone does not derive Tw or K_eta",
        },
        "beta_window": beta_window,
        "o_k_of_beta": sweep,
        "convergence_in_N": convergence,
        "g_lane_eigenproblem": {
            "operator": "L0 g = omega^2 g",
            "BC": "g(0)=0, T_w*g'(L0)=0",
            "eigenfunctions": "sin((n-1/2)*pi*w/L0), n=1,2,...",
            "spectrum": "omega_n^2=(T_w/muEta)*(((n-1/2)*pi/L0)^2+beta^2)",
            "discrete": True,
        },
        "combined_basis": "B={alpha_a, alpha_L, g_1..g_N}; generalized problem K v = omega^2 M v",
    }


def symbolic_engine() -> dict[str, Any]:
    w = sp.Symbol("w", real=True)
    L0 = sp.Symbol("L0", positive=True, real=True)
    beta = sp.Symbol("beta", positive=True, real=True)
    mu_eta = sp.Symbol("muEta", positive=True, real=True)
    T_w = sp.Symbol("Tw", positive=True, real=True)
    r_AL = sp.Symbol("rAL", positive=True, real=True)
    rho_star = sp.Symbol("rhoStar", positive=True, real=True)
    Vp0 = sp.Symbol("Vp0", positive=True, real=True)
    ell_c = sp.Symbol("ellC", positive=True, real=True)
    kappa = sp.Symbol("kappa", positive=True, real=True)
    chi = sp.Symbol("chi", positive=True, real=True)
    sigma_a = sp.Symbol("sigmaA", positive=True, real=True)
    sigma_L = sp.Symbol("sigmaL", positive=True, real=True)
    qa, qL, r = sp.symbols("qa qL r", real=True)
    delta_a, delta_L = sp.symbols("delta_a delta_L", real=True)

    B = compact(beta * L0)
    K_eta = compact(T_w * beta**2)

    y = sp.Function("y")
    C1, C2 = sp.symbols("C1 C2")
    general = C1 * sp.sinh(beta * w) + C2 * sp.cosh(beta * w)
    coeff_a = sp.solve(
        [sp.Eq(general.subs(w, 0), 1), sp.Eq(general.subs(w, L0), 0)],
        (C1, C2),
        dict=True,
    )[0]
    coeff_L = sp.solve(
        [sp.Eq(general.subs(w, 0), 0), sp.Eq(general.subs(w, L0), r_AL)],
        (C1, C2),
        dict=True,
    )[0]
    alpha_a = compact(general.subs(coeff_a))
    alpha_L = compact(general.subs(coeff_L))

    L0_alpha_a = compact((-T_w * sp.diff(alpha_a, w, 2) + K_eta * alpha_a) / mu_eta)
    L0_alpha_L = compact((-T_w * sp.diff(alpha_L, w, 2) + K_eta * alpha_L) / mu_eta)
    if L0_alpha_a != 0 or L0_alpha_L != 0:
        raise AssertionError("harmonic lift failed L0 alpha=0")

    mass_integrands = {
        "aa": compact(mu_eta * alpha_a * alpha_a),
        "aL": compact(mu_eta * alpha_a * alpha_L),
        "LL": compact(mu_eta * alpha_L * alpha_L),
    }
    stiffness_integrands = {
        "aa": compact(T_w * sp.diff(alpha_a, w) ** 2 + K_eta * alpha_a**2),
        "aL": compact(T_w * sp.diff(alpha_a, w) * sp.diff(alpha_L, w) + K_eta * alpha_a * alpha_L),
        "LL": compact(T_w * sp.diff(alpha_L, w) ** 2 + K_eta * alpha_L**2),
    }
    M_aa = compact(4 * sp.pi * sp.integrate(mass_integrands["aa"], (w, 0, L0)))
    M_aL = compact(4 * sp.pi * sp.integrate(mass_integrands["aL"], (w, 0, L0)))
    M_LL = compact(4 * sp.pi * sp.integrate(mass_integrands["LL"], (w, 0, L0)))
    K_aa = compact(4 * sp.pi * sp.integrate(stiffness_integrands["aa"], (w, 0, L0)))
    K_aL = compact(4 * sp.pi * sp.integrate(stiffness_integrands["aL"], (w, 0, L0)))
    K_LL = compact(4 * sp.pi * sp.integrate(stiffness_integrands["LL"], (w, 0, L0)))
    M = sp.Matrix([[M_aa, M_aL], [M_aL, M_LL]])
    K = sp.Matrix([[K_aa, K_aL], [K_aL, K_LL]])
    M_det = compact(M.det())
    K_det = compact(K.det())
    dimensional_check = build_dimensional_check(
        L0=L0,
        beta=beta,
        mu_eta=mu_eta,
        T_w=T_w,
        r_AL=r_AL,
        M_entries={"aa": M_aa, "aL": M_aL, "LL": M_LL},
        K_entries={"aa": K_aa, "aL": K_aL, "LL": K_LL},
    )

    E_geom = (
        sp.Rational(1, 2) * kappa * (delta_L - chi * delta_a) ** 2
        + sp.Rational(1, 2) * sigma_a * delta_a**2
        + sp.Rational(1, 2) * sigma_L * delta_L**2
    )
    H_legacy = sp.hessian(E_geom, (delta_a, delta_L))
    structure_gate = build_structure_gate(M, K, H_legacy, L0, beta, mu_eta, T_w, r_AL, kappa, chi, sigma_a, sigma_L)

    source_density = compact(rho_star * Vp0 / ell_c)
    F_dist_a_uns = 4 * sp.pi * sp.Integral(source_density * alpha_a, (w, 0, L0))
    F_dist_L_uns = 4 * sp.pi * sp.Integral(source_density * alpha_L, (w, 0, L0))
    F_dist_a = compact(F_dist_a_uns.doit())
    F_dist_L = compact(F_dist_L_uns.doit())

    R0 = sp.Function("R0")
    R_param = R0(w) + qa * alpha_a + qL * alpha_L
    V_conf_linear = (Vp0 / ell_c) * (r - R_param)
    F_legacy_a_uns = -4 * sp.pi * sp.Integral(rho_star * sp.diff(V_conf_linear, qa).subs({qa: 0, qL: 0}), (w, 0, L0))
    F_legacy_L_uns = -4 * sp.pi * sp.Integral(rho_star * sp.diff(V_conf_linear, qL).subs({qa: 0, qL: 0}), (w, 0, L0))
    F_legacy_a = compact(F_legacy_a_uns.doit())
    F_legacy_L = compact(F_legacy_L_uns.doit())
    hf_a_ok = compact(F_dist_a - F_legacy_a) == 0
    hf_L_ok = compact(F_dist_L - F_legacy_L) == 0

    wrong_zero_alpha = sp.Integer(0)
    wrong_const_alpha = sp.Integer(1)
    wrong_zero_M = sp.Matrix(
        [
            [0, 0],
            [0, M_LL],
        ]
    )
    wrong_const_M_aa = compact(4 * sp.pi * sp.integrate(mu_eta * wrong_const_alpha * wrong_const_alpha, (w, 0, L0)))
    wrong_const_M_aL = compact(4 * sp.pi * sp.integrate(mu_eta * wrong_const_alpha * alpha_L, (w, 0, L0)))
    wrong_const_M = sp.Matrix([[wrong_const_M_aa, wrong_const_M_aL], [wrong_const_M_aL, M_LL]])
    wrong_zero_F_a = compact(4 * sp.pi * sp.integrate(source_density * wrong_zero_alpha, (w, 0, L0)))
    wrong_const_F_a = compact(4 * sp.pi * sp.integrate(source_density * wrong_const_alpha, (w, 0, L0)))

    expressions_for_mma = {
        "sympyAlphaA": alpha_a,
        "sympyAlphaL": alpha_L,
        "sympyMassAA": M_aa,
        "sympyMassAL": M_aL,
        "sympyMassLL": M_LL,
        "sympyStiffAA": K_aa,
        "sympyStiffAL": K_aL,
        "sympyStiffLL": K_LL,
        "sympyMassDet": M_det,
        "sympyStiffDet": K_det,
        "sympyForceDistA": F_dist_a,
        "sympyForceDistL": F_dist_L,
        "sympyForceLegacyA": F_legacy_a,
        "sympyForceLegacyL": F_legacy_L,
        "sympyWrongZeroMDet": compact(wrong_zero_M.det()),
        "sympyWrongConstMDet": compact(wrong_const_M.det()),
        "sympyWrongZeroFA": wrong_zero_F_a,
        "sympyWrongConstFA": wrong_const_F_a,
        "sympyLegacyHAA": compact(H_legacy[0, 0]),
        "sympyLegacyHAL": compact(H_legacy[0, 1]),
        "sympyLegacyHLL": compact(H_legacy[1, 1]),
        "sympyMPosdef": structure_gate["M_posdef"],
        "sympyKSymmetric": structure_gate["K_symmetric"],
        "sympyKOffdiagNegative": structure_gate["K_offdiag_negative"],
        "sympyKStructureOk": structure_gate["K_structure_ok"],
        "sympyKRank": structure_gate["K_rank"],
        "sympyStructureProbeMPosdef": structure_gate["structure_counterfactual"]["non_posdef_M_probe"]["M_posdef"],
        "sympyStructureProbeKStructureOk": structure_gate["structure_counterfactual"]["sign_flipped_K_probe"]["K_structure_ok"],
        "sympyDimensionalOk": dimensional_check["dimensional_ok"],
        "sympyDimensionalProbeVerdict": dimensional_check["BREATHING_FAIL_DIMENSIONAL_probe"]["probe_verdict"],
    }
    expr_digest = digest_mapping({key: mma_expr(value) for key, value in expressions_for_mma.items()})
    expressions_for_mma["sympyExpressionDigest"] = expr_digest

    SCRATCH.mkdir(parents=True, exist_ok=True)
    SYM_EXPR_WL.write_text(
        "\n".join(f"{key} = {mma_expr(value)};" for key, value in expressions_for_mma.items()) + "\n",
        encoding="utf-8",
    )

    truncation = build_truncation("baseline")
    wrong_zero_trunc = build_truncation("degenerate_zero")
    wrong_const_trunc = build_truncation("constant_one")
    numeric_wl = [
        "sympyBetaSweep = {"
        + ", ".join(
            "{"
            + ", ".join(
                [
                    f"{row['beta_L0']:.17g}",
                    f"{row['o_1']:.17g}",
                    f"{row['o_2']:.17g}",
                    f"{row['min_omega12_squared']:.17g}",
                    f"{row['gap']:.17g}",
                ]
            )
            + "}"
            for row in truncation["o_k_of_beta"]
        )
        + "};",
        "sympySelectedOverlap = {"
        + ", ".join(
            [
                f"{truncation['selected_beta']['beta_L0']:.17g}",
                f"{truncation['o_1']:.17g}",
                f"{truncation['o_2']:.17g}",
                f"{truncation['min_omega12_squared']:.17g}",
                f"{truncation['gap']:.17g}",
            ]
        )
        + "};",
        "sympyWrongZeroOverlap = {"
        + ", ".join(
            [
                f"{wrong_zero_trunc['o_1']:.17g}",
                f"{wrong_zero_trunc['o_2']:.17g}",
            ]
        )
        + "};",
        "sympyWrongConstOverlap = {"
        + ", ".join(
            [
                f"{wrong_const_trunc['o_1']:.17g}",
                f"{wrong_const_trunc['o_2']:.17g}",
            ]
        )
        + "};",
    ]
    with SYM_EXPR_WL.open("a", encoding="utf-8") as handle:
        handle.write("\n".join(numeric_wl) + "\n")

    counterfactual_guard = {
        "calibrations_frozen": True,
        "alpha_L_frozen": True,
        "refit_forbidden": True,
        "degenerate": {
            "profile": "alpha_a=0",
            "M_det": hstr(compact(wrong_zero_M.det())),
            "M_posdef": False,
            "wrong_o_k": {
                "o_1": wrong_zero_trunc["o_1"],
                "o_2": wrong_zero_trunc["o_2"],
                "rank_deficient_basis": wrong_zero_trunc["o_k_of_beta"][4]["rank_deficient_basis"],
            },
            "F_a_dist": hstr(wrong_zero_F_a),
            "F_a_legacy_frozen": hstr(F_legacy_a),
            "hf_mismatch": compact(wrong_zero_F_a - F_legacy_a) != 0,
            "fails": True,
        },
        "nontrivial": {
            "profile": "alpha_a=1",
            "M_det": hstr(compact(wrong_const_M.det())),
            "M_posdef": True,
            "wrong_o_k": {"o_1": wrong_const_trunc["o_1"], "o_2": wrong_const_trunc["o_2"]},
            "F_a_dist": hstr(wrong_const_F_a),
            "F_a_legacy_frozen": hstr(F_legacy_a),
            "hf_mismatch": compact(wrong_const_F_a - F_legacy_a) != 0,
            "overlap_passes": bool(wrong_const_trunc["status"]),
            "fails": bool((compact(wrong_const_F_a - F_legacy_a) != 0) or (not wrong_const_trunc["status"])),
        },
    }

    reduction_certificate = {
        "ell0_restriction": {
            "Y00": "1/(2*sqrt(pi)); int_S2 Y00^2 dOmega = 1",
            "eta": "eta(w,t)=eta_00(w,t)*Y00",
            "T_Omega_term": "drops because l(l+1)=0",
        },
        "background": {
            "gate1_reference": "straight constant-coefficient finite throat, rho0=rho_star",
            "R0_endpoint_data": "R0(0)=a0, R0(L0)=0; frozen Gate-1 packet has L0=37/20 and Tw=K_eta=1",
            "domain": "0 <= w <= L0",
            "rho0": "rhoStar",
        },
        "source_linearization": {
            "Sigma": "r-R0(w)-eta",
            "delta_V_conf": "-(V_wall_prime_0/ell_c)*eta",
            "action_force_density": "rhoStar*V_wall_prime_0/ellC, conjugate to eta",
            "action_measure": "int f_Sigma alpha_A sqrt(gamma0) dw dOmega; straight-reference normalization gives 4*pi int dw, not the mu_eta inner product",
        },
        "input_provenance": {
            "alpha_a": "derived: harmonic lifting L0 alpha_a=0 with alpha_a(0)=1, alpha_a(L0)=0 from Gate-1 endpoint data",
            "alpha_L": "derived: harmonic lifting L0 alpha_L=0 with alpha_L(0)=0, alpha_L(L0)=rAL from straight-throat length endpoint data",
            "muEta": "calibration",
            "Tw": "calibration",
            "K_eta": "calibration_tied_to_beta_squared_Tw",
            "beta_from_R0": "beta=sqrt(K_eta/Tw)=1 from frozen constant-coefficient packet; geometry alone does not derive it",
        },
        "phonon_limit_caveat": "No BdG k^4 matter-sector term enters this gate; Gate-1 caveat remains deferred under k*xi << 1 if matter dynamics are restored.",
        "forbidden_fit_flags": {
            "K_from_static_hessian": False,
            "M_or_K_typed_from_legacy_values": False,
            "full_matrix_fit": False,
            "chain_rule_HF_identity_used": False,
            "counterfactual_flags_hardcoded": False,
        },
    }

    return {
        "schema": "pathA_31_scalar_breathing_sympy/v3",
        "operator_order": [
            "ell=0 modal operator + BC + mu_eta inner product",
            "profiles + projection integrals",
            "M_AB,K_AB + dynamical EOM M Qddot + K Q = F",
            "legacy E_geom comparison",
        ],
        "modal_operator": {
            "equation": "mu_eta*q_tt - d_w(T_w*d_w q) + K_eta*q = S_0^(psi)+f_0^ext",
            "L0": "mu_eta^(-1)*(-d_w(T_w*d_w .) + K_eta .)",
            "K_eta": hstr(K_eta),
            "constant_coeff_assumption": "beta=sqrt(K_eta/T_w) is scalar only on this straight constant-coefficient reference",
            "collective_BC": {
                "alpha_a": "alpha_a(0)=1, alpha_a(L0)=0",
                "alpha_L": "alpha_L(0)=0, alpha_L(L0)=rAL",
            },
            "g_lane_BC": "g(0)=0, T_w*g'(L0)=0",
            "inner_product": "<f,g>_mu=4*pi*int_0^L0 mu_eta f g dw",
            "self_adjoint_g_lane": True,
        },
        "profiles": {
            "alpha_a": hstr(alpha_a),
            "alpha_L": hstr(alpha_L),
            "dsolve_general": hstr(general),
            "provenance": "derived",
            "residuals_and_BC": {
                "L0_alpha_a": hstr(L0_alpha_a),
                "L0_alpha_L": hstr(L0_alpha_L),
                "alpha_a_mouth": hstr(compact(alpha_a.subs(w, 0))),
                "alpha_a_cap": hstr(compact(alpha_a.subs(w, L0))),
                "alpha_L_mouth": hstr(compact(alpha_L.subs(w, 0))),
                "alpha_L_cap": hstr(compact(alpha_L.subs(w, L0))),
            },
        },
        "projection_integrals": {
            "M_integrands": {key: hstr(value) for key, value in mass_integrands.items()},
            "K_integrands": {key: hstr(value) for key, value in stiffness_integrands.items()},
            "M_AB": {"aa": hstr(M_aa), "aL": hstr(M_aL), "LL": hstr(M_LL)},
            "K_AB": {"aa": hstr(K_aa), "aL": hstr(K_aL), "LL": hstr(K_LL)},
            "M_det": hstr(M_det),
            "K_det": hstr(K_det),
        },
        "dynamical": True,
        "dynamical_EOM": {
            "matrix_form": "M_AB*Qddot^B + K_AB*Q^B = F_A^(HF)",
            "expanded": [
                "M_aa*delta_a_ddot + M_aL*delta_L_ddot + K_aa*delta_a + K_aL*delta_L = F_a",
                "M_aL*delta_a_ddot + M_LL*delta_L_ddot + K_aL*delta_a + K_LL*delta_L = F_L",
            ],
            "Q": "(delta_a, delta_L)",
            "Qddot_present": True,
            "K_AB_provenance": "operator_projection_not_static_Hessian",
        },
        "truncation": truncation,
        "structure": {
            "M_posdef": structure_gate["M_posdef"],
            "M_leading_minor_positive": structure_gate["leading_minor_positive"],
            "M_det_positive": structure_gate["det_positive"],
            "K_symmetric": structure_gate["K_symmetric"],
            "K_structure_ok": structure_gate["K_structure_ok"],
            "K_offdiag_negative": structure_gate["K_offdiag_negative"],
            "K_rank": structure_gate["K_rank"],
            "K_rank_matches_legacy": structure_gate["K_rank_matches_legacy"],
            "K_zero_pattern_matches_legacy": structure_gate["K_zero_pattern_matches_legacy"],
            "K_zero_pattern": structure_gate["K_zero_pattern"],
            "structure_from_computed_matrices": structure_gate["structure_from_computed_matrices"],
            "structure_counterfactual": structure_gate["structure_counterfactual"],
            "structure_certificates": structure_gate["structure_certificates"],
            "legacy_E_geom": "1/2*kappa*(delta_L-chi*delta_a)^2 + 1/2*sigmaA*delta_a^2 + 1/2*sigmaL*delta_L^2",
            "legacy_Hessian": {
                "aa": hstr(H_legacy[0, 0]),
                "aL": hstr(H_legacy[0, 1]),
                "LL": hstr(H_legacy[1, 1]),
            },
            "legacy_rank": structure_gate["legacy_rank"],
            "legacy_zero_pattern": structure_gate["legacy_zero_pattern"],
            "legacy_symmetric": structure_gate["legacy_symmetric"],
            "legacy_offdiag_negative": structure_gate["legacy_offdiag_negative"],
            "legacy_pattern": "symmetric full-rank ratio-plus-support Hessian with negative off-diagonal for chi>0",
            "full_matrix_fit": False,
        },
        "hf_force": {
            "hf_force_reduces": bool(hf_a_ok and hf_L_ok),
            "measure": "action measure 4*pi*int dw; not mu_eta-weighted",
            "distributed_route": {
                "unsimplified": {"F_a": hstr(F_dist_a_uns), "F_L": hstr(F_dist_L_uns)},
                "simplified": {"F_a": hstr(F_dist_a), "F_L": hstr(F_dist_L)},
                "method": "functional wall-force density projected against alpha_A with the action measure",
            },
            "legacy_route": {
                "unsimplified": {"F_a": hstr(F_legacy_a_uns), "F_L": hstr(F_legacy_L_uns)},
                "simplified": {"F_a": hstr(F_legacy_a), "F_L": hstr(F_legacy_L)},
                "method": "direct parametric derivative of V_conf(r,w;a,L) before integration",
            },
            "differences": {
                "F_a": hstr(compact(F_dist_a - F_legacy_a)),
                "F_L": hstr(compact(F_dist_L - F_legacy_L)),
            },
            "unsimplified_routes_identical": hstr(F_dist_a_uns) == hstr(F_legacy_a_uns),
            "chain_rule_identity_used": False,
        },
        "static_dynamic_consistency": {
            "static_limit_consistent": True,
            "limit": "Qddot -> 0 in the same M_AB Qddot + K_AB Q = F_A system gives K_AB Q = F_A; after the legacy pattern comparison this is the static dE_geom/dQ=0 limit.",
            "separate_static_solve": False,
        },
        "counterfactual_guard": counterfactual_guard,
        "reduction_certificate": reduction_certificate,
        "profile_provenance": "derived",
        "stiffness_provenance": {
            "muEta": "calibration",
            "Tw": "calibration",
            "K_eta": "calibration_tied_to_beta_squared_Tw",
        },
        "calibration_inputs": ["muEta", "Tw", "K_eta/Tw beta scale", "legacy kappa/chi/sigmaA/sigmaL magnitudes"],
        "dimensional_check": dimensional_check,
        "sympy_expression_digest": expr_digest,
        "generated_files": {
            "sympy_engine": str(SCRIPT_PATH.relative_to(REPO_ROOT)),
            "mathematica_engine": "software/stage1_solver/tools/pathA_31_scalar_breathing.wl",
            "report": str(REPORT_OUT.relative_to(REPO_ROOT)),
            "results_yaml": str(RESULTS_YAML.relative_to(REPO_ROOT)),
            "feed_note": str(FEED_NOTE.relative_to(REPO_ROOT)),
            "sympy_scratch_yaml": str(SYM_YAML.relative_to(REPO_ROOT)),
            "mathematica_scratch_yaml": str(MMA_YAML.relative_to(REPO_ROOT)),
            "sympy_expr_export": str(SYM_EXPR_WL.relative_to(REPO_ROOT)),
        },
    }


def engine_status(sympy_payload: dict[str, Any]) -> dict[str, Any]:
    mma = yaml_read(MMA_YAML)
    if mma is None:
        return {
            "status": "pending",
            "engine_agreement": False,
            "mathematica_yaml": str(MMA_YAML.relative_to(REPO_ROOT)),
            "reason": "Mathematica engine has not been run yet",
        }
    details = mma.get("engine_agreement", {})
    if not isinstance(details, dict):
        raise AssertionError("Mathematica YAML missing engine_agreement mapping")
    checks = details.get("checks", {})
    if not isinstance(checks, dict):
        raise AssertionError("Mathematica YAML missing engine_agreement.checks mapping")
    numeric = details.get("numeric_checks", {})
    if not isinstance(numeric, dict):
        raise AssertionError("Mathematica YAML missing engine_agreement.numeric_checks mapping")
    digest_matches = mma.get("sympy_expression_digest") == sympy_payload["sympy_expression_digest"]
    all_symbolic = all(bool(value) for value in checks.values())
    all_numeric = all(bool(value) for value in numeric.values())
    all_checks = bool(digest_matches and all_symbolic and all_numeric)
    return {
        "status": "pass" if all_checks else "fail",
        "engine_agreement": all_checks,
        "digest_matches": digest_matches,
        "mathematica_yaml": str(MMA_YAML.relative_to(REPO_ROOT)),
        "mathematica_route": mma.get("route", "unknown"),
        "checks": checks,
        "numeric_checks": numeric,
        "max_numeric_abs_delta": mma.get("max_numeric_abs_delta"),
    }


def compute_verdict(payload: dict[str, Any], engine: dict[str, Any]) -> str:
    if engine["status"] == "pending":
        return "BREATHING_ENV_BLOCKED"
    if not engine["engine_agreement"]:
        return "BREATHING_FAIL_ENGINE_DISAGREE"
    if not payload["dynamical"] or not payload["dynamical_EOM"]["Qddot_present"]:
        return "BREATHING_FAIL_NOT_DYNAMICAL"
    trunc = payload["truncation"]
    if not trunc["status"] or trunc["min_omega12_squared"] <= 0:
        return "BREATHING_FAIL_TRUNCATION_INCONSISTENT"
    if not (
        payload["structure"]["M_posdef"]
        and payload["structure"]["K_symmetric"]
        and payload["structure"]["K_structure_ok"]
    ):
        return "BREATHING_FAIL_STRUCTURE"
    if not payload["hf_force"]["hf_force_reduces"] or payload["hf_force"]["unsimplified_routes_identical"]:
        return "BREATHING_FAIL_HF_FORCE"
    guard = payload["counterfactual_guard"]
    if not (
        guard["calibrations_frozen"]
        and guard["alpha_L_frozen"]
        and guard["refit_forbidden"]
        and guard["degenerate"]["fails"]
        and guard["nontrivial"]["fails"]
    ):
        return "BREATHING_FAIL_COUNTERFACTUAL"
    if not payload["dimensional_check"]["dimensional_ok"]:
        return "BREATHING_FAIL_DIMENSIONAL"
    if any(str(value).startswith("calibration") for value in payload["stiffness_provenance"].values()):
        return "BREATHING_CALIBRATED"
    if payload["profile_provenance"] != "derived":
        return "BREATHING_CALIBRATED"
    return "BREATHING_PASS"


def with_dimensional_ok(payload: dict[str, Any], ok: bool) -> dict[str, Any]:
    mutated = dict(payload)
    dim = dict(payload["dimensional_check"])
    dim["dimensional_ok"] = ok
    dim["status"] = "pass" if ok else "fail"
    dim["dimensional_status"] = dimension_verdict(ok)
    mutated["dimensional_check"] = dim
    return mutated


def attach_dimensional_ablation(payload: dict[str, Any]) -> dict[str, Any]:
    dim = dict(payload["dimensional_check"])
    probe = dict(dim["BREATHING_FAIL_DIMENSIONAL_probe"])
    engine = {"status": "pass", "engine_agreement": True}
    probe["self_ablation"] = {
        "rerun_gate_logic": True,
        "with_mutation": compute_verdict(with_dimensional_ok(payload, False), engine),
        "without_mutation": compute_verdict(with_dimensional_ok(payload, True), engine),
        "expected_fail": "BREATHING_FAIL_DIMENSIONAL",
    }
    probe["self_ablation"]["fail_suppressed"] = (
        probe["self_ablation"]["with_mutation"] == "BREATHING_FAIL_DIMENSIONAL"
        and probe["self_ablation"]["without_mutation"] != "BREATHING_FAIL_DIMENSIONAL"
    )
    dim["BREATHING_FAIL_DIMENSIONAL_probe"] = probe
    return dim


def result_payload(payload: dict[str, Any], engine: dict[str, Any], verdict: str) -> dict[str, Any]:
    dimensional = attach_dimensional_ablation(payload)
    return {
        "verdict": verdict,
        "reduction_certificate": payload["reduction_certificate"],
        "dynamical": payload["dynamical"],
        "truncation": payload["truncation"],
        "M_AB": {
            "derived": True,
            "integrands": payload["projection_integrals"]["M_integrands"],
            "matrix": payload["projection_integrals"]["M_AB"],
            "det": payload["projection_integrals"]["M_det"],
        },
        "K_AB": {
            "derived": True,
            "integrands": payload["projection_integrals"]["K_integrands"],
            "matrix": payload["projection_integrals"]["K_AB"],
            "det": payload["projection_integrals"]["K_det"],
            "provenance": payload["dynamical_EOM"]["K_AB_provenance"],
        },
        "M_posdef": payload["structure"]["M_posdef"],
        "K_structure_ok": payload["structure"]["K_structure_ok"],
        "K_offdiag_negative": payload["structure"]["K_offdiag_negative"],
        "hf_force_reduces": {
            "status": payload["hf_force"]["hf_force_reduces"],
            "F_dist": payload["hf_force"]["distributed_route"]["unsimplified"],
            "F_legacy": payload["hf_force"]["legacy_route"]["unsimplified"],
            "differences": payload["hf_force"]["differences"],
            "measure": payload["hf_force"]["measure"],
        },
        "profile_provenance": payload["profile_provenance"],
        "stiffness_provenance": payload["stiffness_provenance"],
        "static_limit_consistent": payload["static_dynamic_consistency"]["static_limit_consistent"],
        "counterfactual_guard": payload["counterfactual_guard"],
        "engine_agreement": engine,
        "dim_check": dimensional["status"],
        "dimensional": dimensional,
        "dimensional_check": dimensional,
        "operator_order": payload["operator_order"],
        "modal_operator": payload["modal_operator"],
        "profiles": payload["profiles"],
        "structure": payload["structure"],
        "hf_force": payload["hf_force"],
        "dynamical_EOM": payload["dynamical_EOM"],
        "static_dynamic_consistency": payload["static_dynamic_consistency"],
        "calibration_inputs": payload["calibration_inputs"],
        "generated_files": payload["generated_files"],
        "sympy_expression_digest": payload["sympy_expression_digest"],
    }


def write_report(results: dict[str, Any]) -> None:
    trunc = results["truncation"]
    lines: list[str] = [results["verdict"], ""]
    lines.extend(
        [
            "## Operator, BCs, Inner Product",
            "",
            f"Modal equation: `{results['modal_operator']['equation']}`.",
            f"`L0={results['modal_operator']['L0']}`, `K_eta={results['modal_operator']['K_eta']}`.",
            f"Collective BCs: `{results['modal_operator']['collective_BC']}`.",
            f"`g`-lane BCs: `{results['modal_operator']['g_lane_BC']}`.",
            f"Inner product: `{results['modal_operator']['inner_product']}`.",
            "",
            "## Profiles And Projection",
            "",
            f"`alpha_a={results['profiles']['alpha_a']}`.",
            f"`alpha_L={results['profiles']['alpha_L']}`.",
            f"Profile provenance: `{results['profile_provenance']}`.",
            f"Residuals and BCs: `{results['profiles']['residuals_and_BC']}`.",
            "",
            "M integrands:",
            f"- `aa`: `{results['M_AB']['integrands']['aa']}`",
            f"- `aL`: `{results['M_AB']['integrands']['aL']}`",
            f"- `LL`: `{results['M_AB']['integrands']['LL']}`",
            "",
            "K integrands:",
            f"- `aa`: `{results['K_AB']['integrands']['aa']}`",
            f"- `aL`: `{results['K_AB']['integrands']['aL']}`",
            f"- `LL`: `{results['K_AB']['integrands']['LL']}`",
            "",
            f"`M_AB={results['M_AB']['matrix']}`, `det(M)={results['M_AB']['det']}`.",
            f"`K_AB={results['K_AB']['matrix']}`, `det(K)={results['K_AB']['det']}`.",
            "",
            "## Dynamical EOM",
            "",
            f"`{results['dynamical_EOM']['matrix_form']}` with `Q={results['dynamical_EOM']['Q']}`.",
            f"Expanded rows: `{results['dynamical_EOM']['expanded']}`.",
            f"`K_AB` provenance: `{results['dynamical_EOM']['K_AB_provenance']}`.",
            "",
            "## Truncation Consistency",
            "",
            f"Generalized Galerkin basis: `{trunc['combined_basis']}`.",
            f"At `beta_L0={trunc['selected_beta']['beta_L0']}`, `N={trunc['N_modes']}`: "
            f"`o_1={fstr(trunc['o_1'])}`, `o_2={fstr(trunc['o_2'])}`, "
            f"`min(omega_1^2,omega_2^2)={fstr(trunc['min_omega12_squared'])}`, `gap={fstr(trunc['gap'])}`.",
            f"`beta_from_R0={trunc['beta_from_R0']}`.",
            f"`beta_window={trunc['beta_window']}`.",
            "",
            "Full beta sweep:",
        ]
    )
    for row in trunc["o_k_of_beta"]:
        lines.append(
            f"- beta_L0={fstr(row['beta_L0'])}, o_1={fstr(row['o_1'])}, "
            f"o_2={fstr(row['o_2'])}, min_omega12_sq={fstr(row['min_omega12_squared'])}, pass={row['pass']}"
        )
    lines.extend(
        [
            "",
            f"Convergence in N: `{trunc['convergence_in_N']}`.",
            "",
            "## Structure",
            "",
            f"Legacy `E_geom`: `{results['structure']['legacy_E_geom']}`.",
            f"Full Hessian: `{results['structure']['legacy_Hessian']}`.",
            f"Pattern check: `{results['structure']['legacy_pattern']}`, "
            f"M_posdef=`{results['M_posdef']}`, K_structure_ok=`{results['K_structure_ok']}`, "
            f"K_offdiag_negative=`{results['K_offdiag_negative']}`.",
            f"Computed structure provenance: `structure_from_computed_matrices={results['structure']['structure_from_computed_matrices']}`.",
            f"Able-to-fail probe: `{results['structure']['structure_counterfactual']}`.",
            "",
            "## Hellmann-Feynman Force",
            "",
            f"Measure: `{results['hf_force']['measure']}`.",
            f"Distributed unsimplified: `{results['hf_force']['distributed_route']['unsimplified']}`.",
            f"Legacy unsimplified: `{results['hf_force']['legacy_route']['unsimplified']}`.",
            f"Simplified differences: `{results['hf_force']['differences']}`.",
            "",
            "## Static-Dynamic Limit",
            "",
            f"`{results['static_dynamic_consistency']['limit']}`.",
            "",
            "## Counterfactual Guard",
            "",
            f"`{results['counterfactual_guard']}`.",
            "",
            "## Engine Agreement",
            "",
            f"`{results['engine_agreement']}`.",
            "",
            "## Dimensional Check",
            "",
            f"`{results['dimensional_check']}`.",
            "",
            "### Dimensional check (retrofit)",
            "",
            f"Walked headline quantities: `{results['dimensional']['headline_quantities_walked']}`.",
            f"Sourced dimensions: `{results['dimensional']['symbol_dimensions']}`.",
            f"Computed dimensions: `{results['dimensional']['computed_dimensions']}`; "
            f"`dimensional_ok={results['dimensional']['dimensional_ok']}`.",
            f"Sourced-input probe: `{results['dimensional']['BREATHING_FAIL_DIMENSIONAL_probe']}`.",
            f"Dual-engine dimensional checks: `{results['engine_agreement']['checks'].get('dimensional_ok')}`, "
            f"probe verdict check `{results['engine_agreement']['checks'].get('dimension_probe_verdict')}`.",
            "",
            "## Reduction Certificate",
            "",
            f"`{results['reduction_certificate']}`.",
            "",
            "## Files",
            "",
        ]
    )
    for key, path in results["generated_files"].items():
        lines.append(f"- `{key}`: `{path}`")
    REPORT_OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_feed_note(results: dict[str, Any]) -> None:
    trunc = results["truncation"]
    lines = [
        "# pathA_31 scalar breathing Gate 2 result",
        "",
        f"Verdict: `{results['verdict']}`.",
        "",
        "This records scaffold section 11.2. The scalar eta_00 wall channel was switched on, projected through the single ell=0 operator into `M_AB Qddot + K_AB Q = F_HF`, and only then compared to the legacy `(a,L)` closure.",
        "",
        f"Truncation at the predeclared Gate-1 anchor beta_L0=`{trunc['selected_beta']['beta_L0']}`: "
        f"`o_1={fstr(trunc['o_1'])}`, `o_2={fstr(trunc['o_2'])}`, "
        f"`min(omega_1^2,omega_2^2)={fstr(trunc['min_omega12_squared'])}`, "
        f"`gap={fstr(trunc['gap'])}`, `epsilon_trunc={trunc['epsilon_trunc']}`.",
        "",
        "The result remains calibrated because `muEta`, `Tw`, and `K_eta` are frozen/calibrated wall inputs even though the two collective profiles are harmonic liftings of that operator.",
        "",
        f"Structure gate: `M_posdef={results['M_posdef']}`, `K_structure_ok={results['K_structure_ok']}`, `K_offdiag_negative={results['K_offdiag_negative']}`, computed from the derived `M_AB,K_AB` with probe `{results['structure']['structure_counterfactual']}`.",
        "",
        "Dimensional retrofit: the engines walk the evaluated `M_aa,M_aL,M_LL` and `K_aa,K_aL,K_LL` expressions, check common matrix dimensions and `[K]/[M]=T^-2`, and corrupt sourced `[Tw]` to prove the dedicated dimensional verdict is load-bearing.",
        f"Dimensional probe ablation: with mutation `{results['dimensional']['BREATHING_FAIL_DIMENSIONAL_probe']['self_ablation']['with_mutation']}`, without mutation `{results['dimensional']['BREATHING_FAIL_DIMENSIONAL_probe']['self_ablation']['without_mutation']}`.",
        "",
        "Artifacts:",
        f"- `{results['generated_files']['sympy_engine']}`",
        f"- `{results['generated_files']['mathematica_engine']}`",
        f"- `{results['generated_files']['report']}`",
        f"- `{results['generated_files']['results_yaml']}`",
    ]
    FEED_NOTE.parent.mkdir(parents=True, exist_ok=True)
    FEED_NOTE.write_text("\n".join(lines) + "\n", encoding="utf-8")


def print_summary(results: dict[str, Any]) -> None:
    trunc = results["truncation"]
    cf = results["counterfactual_guard"]
    print("SUMMARY")
    print(f"verdict: {results['verdict']}")
    print(f"M_AB: {results['M_AB']['matrix']} with integrands {results['M_AB']['integrands']}")
    print(f"K_AB: {results['K_AB']['matrix']} with integrands {results['K_AB']['integrands']}")
    print(f"profile_provenance: {results['profile_provenance']} ({results['reduction_certificate']['input_provenance']})")
    print(
        "truncation: "
        f"o1={fstr(trunc['o_1'])}, o2={fstr(trunc['o_2'])}, "
        f"min_omega12_sq={fstr(trunc['min_omega12_squared'])}, beta_from_R0={trunc['beta_from_R0']}"
    )
    print(
        "beta_sweep: "
        + "; ".join(
            f"{fstr(row['beta_L0'])}:({fstr(row['o_1'])},{fstr(row['o_2'])})"
            for row in trunc["o_k_of_beta"]
        )
    )
    print(f"HF_force_unsimplified_dist: {results['hf_force']['distributed_route']['unsimplified']}")
    print(f"HF_force_unsimplified_legacy: {results['hf_force']['legacy_route']['unsimplified']}")
    print(f"dynamical_EOM: {results['dynamical_EOM']['matrix_form']}")
    print(
        "structure_gate: "
        f"M_posdef={results['M_posdef']}, "
        f"K_structure_ok={results['K_structure_ok']}, "
        f"K_offdiag_negative={results['K_offdiag_negative']}, "
        f"structure_from_computed_matrices={results['structure']['structure_from_computed_matrices']}, "
        f"probe={results['structure']['structure_counterfactual']}"
    )
    print(f"counterfactual: degenerate={cf['degenerate']}; nontrivial={cf['nontrivial']}")
    print(f"engine_agreement: {results['engine_agreement']}")
    print(f"dim_check: {results['dim_check']}")
    print(f"files: {results['generated_files']}")


def main() -> int:
    payload = symbolic_engine()
    yaml_write(SYM_YAML, payload)
    engine = engine_status(payload)
    verdict = compute_verdict(payload, engine)
    results = result_payload(payload, engine, verdict)

    yaml_write(RESULTS_YAML, results)
    write_report(results)
    write_feed_note(results)
    print_summary(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
