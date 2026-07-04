#!/usr/bin/env python3
"""pathA_39 Stage 3 operator-parity gate, SymPy engine.

This gate is a finite spurion-EFT selection-rule computation.  It does not
try to compute a literal reduced-slab dO/dV, because the imported pathA_38
background is w-only.  Instead it generates the declared O(V) canonical slot
manifest, filters Hermiticity, classifies P_w in two ways, and checks native
mode matrix-element witnesses.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
import yaml
from sympy.printing.mathematica import mathematica_code


SCRIPT_PATH = Path(__file__).resolve()
STAGE1_ROOT = SCRIPT_PATH.parents[1]
REPORTS = STAGE1_ROOT / "reports"
SCRATCH = STAGE1_ROOT / "_scratch"

PATH38_YAML = REPORTS / "pathA_38_results.yaml"
PATH39_STAGE01_REPORT = REPORTS / "pathA_39_scalar_admixture_screen.md"
PATH39_STAGE2_REPORT = REPORTS / "pathA_39_magnetic_force.md"

SYM_OUT = SCRATCH / "pathA_39_stage3_operator_parity_sympy.json"
WL_OUT = SCRATCH / "pathA_39_stage3_operator_parity_mathematica.json"
YAML_OUT = REPORTS / "pathA_39_stage3_operator_parity_results.yaml"
REPORT_OUT = REPORTS / "pathA_39_stage3_operator_parity.md"

SCHEMA = "pathA_39_stage3_operator_parity/v1"
WL_SCHEMA = "pathA_39_stage3_operator_parity_mathematica/v1"

VERDICT_CODES = {
    "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING": 1,
    "PASS_CONDITIONAL_ON_NO_OPERATOR_PARITY_MIXING": 2,
    "PASS_TO_O(V)_CONTROLLED_O(V2)_DEFERRED": 3,
    "WEAK_BREAKING_EPSILON_V": 4,
    "ACCIDENTAL_ZERO_SIM_DEFERRED": 5,
    "BASIS_INCOMPLETE": 6,
    "CLEAN_PARITY_PRESERVING": 7,
    "REJECT_ILLEGAL_INPUT": 8,
    "REJECT_UNDECLARED_SPURION": 9,
    "PARITY_TWO_WAYS_AGREE": 10,
    "CLEAN_STATIC_SPLIT": 11,
    "A1_IDS_RESTORED": 12,
}

FIELDS = ("eta", "n", "tau")
FIELD_SIGMA = {"eta": -1, "n": -1, "tau": 1}
X_FACTORS = ("iVdotk", "omegaVdotk")
X_ADJOINT = {"iVdotk": -1, "omegaVdotk": 1}
COEFF_TYPES = ("C", "B", "K")
W_MODULE = {
    "C": {"N_w": 0, "adjoint_sign": 1, "module": "C(w)"},
    "B": {"N_w": 1, "adjoint_sign": -1, "module": "A_B=1/2[B D_w + D_w B]"},
    "K": {"N_w": 2, "adjoint_sign": 1, "module": "-D_w K D_w"},
}
PF_SIGNS = {"even": 1, "odd": -1}
EXPECTED_RAW_COUNT_FULL = 2 * 2 * 9 * 3 * 2


@dataclass(frozen=True)
class MatrixDef:
    name: str
    kind: str
    left: str
    right: str
    entries: tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]

    @property
    def sigma_ab(self) -> int:
        return FIELD_SIGMA[self.left] * FIELD_SIGMA[self.right]

    @property
    def adjoint_sign(self) -> int:
        return 1 if self.kind == "symmetric" else -1

    @property
    def adjoint_class(self) -> str:
        return "self-adjoint" if self.kind == "symmetric" else "anti-self-adjoint"

    @property
    def matrix(self) -> sp.Matrix:
        return sp.Matrix(self.entries)


@dataclass(frozen=True)
class Slot:
    X: str
    a: int
    M: MatrixDef
    coeff_type: str
    p_f: str

    @property
    def id(self) -> str:
        return f"{self.X}|{self.a}|{self.M.name}|{self.coeff_type}|{self.p_f}"


class GrammarValidationError(ValueError):
    def __init__(self, verdict: str, reason: str, observation: dict[str, Any]):
        super().__init__(reason)
        self.verdict = verdict
        self.reason = reason
        self.observation = observation


def matrix_defs() -> list[MatrixDef]:
    z = ((0, 0, 0), (0, 0, 0), (0, 0, 0))

    def m(*entries: tuple[int, int, int]) -> tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]:
        rows = [list(row) for row in z]
        for i, j, value in entries:
            rows[i][j] = value
        return tuple(tuple(row) for row in rows)  # type: ignore[return-value]

    return [
        MatrixDef("etaeta", "symmetric", "eta", "eta", m((0, 0, 1))),
        MatrixDef("nn", "symmetric", "n", "n", m((1, 1, 1))),
        MatrixDef("tautau", "symmetric", "tau", "tau", m((2, 2, 1))),
        MatrixDef("etan", "symmetric", "eta", "n", m((0, 1, 1), (1, 0, 1))),
        MatrixDef("etatau", "symmetric", "eta", "tau", m((0, 2, 1), (2, 0, 1))),
        MatrixDef("ntau", "symmetric", "n", "tau", m((1, 2, 1), (2, 1, 1))),
        MatrixDef("[etan]", "antisymmetric", "eta", "n", m((0, 1, 1), (1, 0, -1))),
        MatrixDef("[etatau]", "antisymmetric", "eta", "tau", m((0, 2, 1), (2, 0, -1))),
        MatrixDef("[ntau]", "antisymmetric", "n", "tau", m((1, 2, 1), (2, 1, -1))),
    ]


def sign_str(value: int) -> str:
    return "+1" if value == 1 else "-1"


def verdict_code(label: str) -> int:
    return VERDICT_CODES[label]


def hstr(expr: Any) -> Any:
    if expr is None or isinstance(expr, (bool, int, str)):
        return expr
    if isinstance(expr, sp.MatrixBase):
        return [[hstr(entry) for entry in row] for row in expr.tolist()]
    return sp.sstr(sp.factor(sp.simplify(expr)))


def mma_expr(expr: sp.Expr | int) -> str:
    return mathematica_code(sp.simplify(expr))


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_json(data: Any) -> str:
    return sha256_text(json.dumps(data, sort_keys=True, separators=(",", ":")))


def load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise AssertionError(f"YAML did not load to mapping: {path}")
    return data


def import_banked_inputs() -> dict[str, Any]:
    p38 = load_yaml(PATH38_YAML)
    if not PATH39_STAGE01_REPORT.exists():
        raise FileNotFoundError(PATH39_STAGE01_REPORT)
    if not PATH39_STAGE2_REPORT.exists():
        raise FileNotFoundError(PATH39_STAGE2_REPORT)

    s01_text = PATH39_STAGE01_REPORT.read_text(encoding="utf-8")
    s2_text = PATH39_STAGE2_REPORT.read_text(encoding="utf-8")

    fields_match = p38["model"]["fields"] == list(FIELDS)
    sigma_import = [[int(entry) for entry in row] for row in p38["model"]["parity_matrix"]]
    sigma_match = sigma_import == [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    mode_names = [item["name"] for item in p38["transverse_mode_spectrum"]["modes"]]
    required_modes = ["translation_goldstone", "wall_shape_partner", "relative_lock_partner", "tau_source_tilt"]
    modes_match = all(name in mode_names for name in required_modes)
    engine_match = p38["engine_agreement"]["status"] == "ENGINE_AGREE"
    stage01_match = "FAIL_OBSERVABLE_SCALAR_ADMIXTURE" in s01_text and "ENGINE_AGREE" in s01_text
    stage2_match = "MAGNETIC_FORCE_DERIVED" in s2_text and "ENGINE_AGREE" in s2_text

    checks = {
        "pathA_38_fields": {"expected": list(FIELDS), "actual": p38["model"]["fields"], "match": fields_match},
        "pathA_38_Sigma": {
            "expected": "diag(-1,-1,+1)",
            "actual": sigma_import,
            "match": sigma_match,
        },
        "pathA_38_modes": {"expected": required_modes, "actual": mode_names, "match": modes_match},
        "pathA_38_engine": {"expected": "ENGINE_AGREE", "actual": p38["engine_agreement"]["status"], "match": engine_match},
        "pathA_39_stage0_plus_1": {
            "expected": "FAIL_OBSERVABLE_SCALAR_ADMIXTURE with ENGINE_AGREE",
            "actual": "report text",
            "match": stage01_match,
        },
        "pathA_39_stage2": {
            "expected": "MAGNETIC_FORCE_DERIVED with ENGINE_AGREE",
            "actual": "report text",
            "match": stage2_match,
        },
    }
    bad = [name for name, item in checks.items() if not item["match"]]
    if bad:
        raise AssertionError(f"banked import check failed: {bad}")

    return {
        "checks": checks,
        "imported": {
            "operator": "O=-partial_w(K_w partial_w)+M2(w), K_w=K_parallel=I_3",
            "Sigma": "diag(-1,-1,+1)",
            "f0": "[sech(w/ell)^2/ell, sech(w/ell)^2/ell, 0], P_w odd",
            "wall_partner": "[sech(w/ell)*tanh(w/ell), sech(w/ell)*tanh(w/ell), 0], P_w even",
            "tau": "[0,0,1], P_w even",
            "relative_lock": "[sech(w/ell)^2/ell,-sech(w/ell)^2/ell,0], P_w odd; imported but not targeted",
            "pathA_39_stage0_plus_1": "FAIL_OBSERVABLE_SCALAR_ADMIXTURE; source-level split not un-earned",
            "pathA_39_stage2": "MAGNETIC_FORCE_DERIVED; operator motion coefficient remains sim-deferred",
        },
        "declared": {
            "external_factors": list(X_FACTORS),
            "charge_spurion": "s is P_w odd; a in {0,1} at leading O(s^1)",
            "coefficients": "one free even or odd coefficient function per canonical slot",
            "forbidden_spurions": ["V_i n_i without a brane vector", "epsilon_ijk V_i k_j without an axial vector"],
        },
        "sim_deferred": [
            "realized nonlinear-throat O(V) coefficient",
            "finite-mouth coefficient profile and compactness",
            "contamination magnitude hierarchy",
        ],
    }


def generate_slots(a_values: Iterable[int] = (0, 1)) -> list[Slot]:
    mats = matrix_defs()
    return [
        Slot(X=X, a=a, M=M, coeff_type=coeff_type, p_f=p_f)
        for X in X_FACTORS
        for a in a_values
        for M in mats
        for coeff_type in COEFF_TYPES
        for p_f in ("even", "odd")
    ]


def concrete_coefficient_profile(p_f: str, z: sp.Symbol) -> sp.Expr:
    sech = lambda arg: 1 / sp.cosh(arg)
    if p_f == "even":
        return sech(z) ** 2
    if p_f == "odd":
        return sp.tanh(z) * sech(z) ** 2
    raise AssertionError(f"unknown coefficient parity: {p_f}")


def representative_derivative_tokens(coeff_type: str) -> tuple[sp.Symbol, ...]:
    if coeff_type == "C":
        return ()
    if coeff_type == "B":
        return (sp.Symbol("D_w_B"),)
    if coeff_type == "K":
        return (sp.Symbol("D_w_K_left"), sp.Symbol("D_w_K_right"))
    raise AssertionError(f"unknown coeff_type: {coeff_type}")


def product_expr(factors: Iterable[sp.Expr]) -> sp.Expr:
    total = sp.Integer(1)
    for factor in factors:
        total *= factor
    return total


def matrix_equal(lhs: sp.MatrixBase, rhs: sp.MatrixBase) -> bool:
    if lhs.shape != rhs.shape:
        return False
    return all(sp.simplify(lhs[i, j] - rhs[i, j]) == 0 for i in range(lhs.rows) for j in range(lhs.cols))


def parity_conjugation_sign(slot: Slot) -> int:
    z, s = sp.symbols("z s")
    sigma = sp.diag(-1, -1, 1)
    derivative_tokens = representative_derivative_tokens(slot.coeff_type)
    scalar = (s**slot.a) * concrete_coefficient_profile(slot.p_f, z) * product_expr(derivative_tokens)
    operator = scalar * slot.M.matrix

    substitutions = {z: -z, s: -s, **{token: -token for token in derivative_tokens}}
    conjugated_scalar = sp.simplify(scalar.subs(substitutions, simultaneous=True))
    conjugated = conjugated_scalar * (sigma * slot.M.matrix * sigma)

    if matrix_equal(conjugated, operator):
        return 1
    if matrix_equal(conjugated, -operator):
        return -1
    raise AssertionError(f"literal P O P^-1 is not a parity eigen-operator: {slot.id}")


def illegal_attempt_spec(illegal_input: str) -> dict[str, Any]:
    base: dict[str, Any] = {
        "X": "iVdotk",
        "a": 1,
        "M": "etaeta",
        "coeff_type": "C",
        "p_f": "even",
        "odd_spurions": ["s"],
    }
    if illegal_input == "extra_odd_spurion":
        return {
            **base,
            "attempted_structure": "iVdotk*s*q_odd_extra*C_even(w)*etaeta",
            "odd_spurions": ["s", "q_odd_extra"],
        }
    if illegal_input == "V_i n_i":
        return {
            **base,
            "attempted_structure": "V_i n_i*s*C_even(w)*etaeta",
            "X": "V_i n_i",
            "requires_declared_brane_vector": "d_i",
        }
    if illegal_input == "epsilon_ijk V_i k_j":
        return {
            **base,
            "attempted_structure": "epsilon_ijk V_i k_j*s*C_even(w)*etaeta",
            "X": "epsilon_ijk V_i k_j",
            "requires_declared_axial_spurion": "J_k",
        }
    raise AssertionError(f"unknown illegal control input: {illegal_input}")


def validate_operator_attempt(attempt: dict[str, Any]) -> Slot:
    declared_odd_spurions = {"s"}
    declared_brane_vectors: set[str] = set()
    declared_axial_spurions: set[str] = set()
    observations = {
        "declared_external_scalars": list(X_FACTORS),
        "declared_odd_spurions": sorted(declared_odd_spurions),
        "declared_brane_vectors": sorted(declared_brane_vectors),
        "declared_axial_spurions": sorted(declared_axial_spurions),
        "attempted_structure": attempt["attempted_structure"],
    }

    odd_spurions = set(attempt.get("odd_spurions", []))
    undeclared_odd = sorted(odd_spurions - declared_odd_spurions)
    if undeclared_odd:
        raise GrammarValidationError(
            "REJECT_ILLEGAL_INPUT",
            "undeclared P_w-odd spurion distinct from s",
            {**observations, "undeclared_odd_spurions": undeclared_odd},
        )

    required_brane_vector = attempt.get("requires_declared_brane_vector")
    if required_brane_vector and required_brane_vector not in declared_brane_vectors:
        raise GrammarValidationError(
            "REJECT_UNDECLARED_SPURION",
            f"{attempt['X']} requires an undeclared brane vector",
            {**observations, "missing_declared_brane_vector": required_brane_vector},
        )

    required_axial = attempt.get("requires_declared_axial_spurion")
    if required_axial and required_axial not in declared_axial_spurions:
        raise GrammarValidationError(
            "REJECT_UNDECLARED_SPURION",
            f"{attempt['X']} requires an undeclared axial spurion",
            {**observations, "missing_declared_axial_spurion": required_axial},
        )

    if attempt["X"] not in X_FACTORS:
        raise GrammarValidationError(
            "REJECT_ILLEGAL_INPUT",
            f"{attempt['X']} is not a declared external scalar",
            {**observations, "undeclared_external_scalar": attempt["X"]},
        )
    if attempt["a"] not in (0, 1):
        raise GrammarValidationError(
            "REJECT_ILLEGAL_INPUT",
            f"a={attempt['a']} is outside the leading O(s^1) grammar",
            {**observations, "invalid_a": attempt["a"]},
        )
    if attempt["coeff_type"] not in COEFF_TYPES:
        raise GrammarValidationError(
            "REJECT_ILLEGAL_INPUT",
            f"{attempt['coeff_type']} is not a declared coefficient module",
            {**observations, "invalid_coeff_type": attempt["coeff_type"]},
        )
    if attempt["p_f"] not in PF_SIGNS:
        raise GrammarValidationError(
            "REJECT_ILLEGAL_INPUT",
            f"{attempt['p_f']} is not a declared coefficient parity",
            {**observations, "invalid_p_f": attempt["p_f"]},
        )

    matrix_by_name = {matrix.name: matrix for matrix in matrix_defs()}
    if attempt["M"] not in matrix_by_name:
        raise GrammarValidationError(
            "REJECT_ILLEGAL_INPUT",
            f"{attempt['M']} is not in the declared field-matrix basis",
            {**observations, "invalid_matrix": attempt["M"]},
        )

    slot = Slot(
        X=attempt["X"],
        a=int(attempt["a"]),
        M=matrix_by_name[attempt["M"]],
        coeff_type=attempt["coeff_type"],
        p_f=attempt["p_f"],
    )
    generated_ids = {generated.id for generated in generate_slots((0, 1))}
    if slot.id not in generated_ids:
        raise GrammarValidationError(
            "REJECT_ILLEGAL_INPUT",
            f"{slot.id} was not generated by the canonical grammar",
            {**observations, "candidate_slot_id": slot.id},
        )
    return slot


def illegal_control_result(name: str, illegal_input: str) -> dict[str, Any]:
    attempt = illegal_attempt_spec(illegal_input)
    try:
        accepted_slot = validate_operator_attempt(attempt)
    except GrammarValidationError as exc:
        return {
            "name": name,
            "status": "FIRED",
            "verdict": exc.verdict,
            "verdict_code": verdict_code(exc.verdict),
            "reason": exc.reason,
            "attempted_structure": attempt["attempted_structure"],
            "validation_observation": "REFUSED_BY_DECLARED_GRAMMAR",
            "validator_observation": exc.observation,
        }
    return {
        "name": name,
        "status": "NOT_FIRED",
        "verdict": "ILLEGAL_INPUT_ACCEPTED",
        "reason": "declared grammar admitted an illegal control structure",
        "attempted_structure": attempt["attempted_structure"],
        "accepted_slot_id": accepted_slot.id,
    }


def classify_slot(slot: Slot) -> dict[str, Any]:
    x_sign = X_ADJOINT[slot.X]
    m_sign = slot.M.adjoint_sign
    w_sign = int(W_MODULE[slot.coeff_type]["adjoint_sign"])
    product = x_sign * m_sign * w_sign
    n_w = int(W_MODULE[slot.coeff_type]["N_w"])
    p_sign = PF_SIGNS[slot.p_f]
    parity_rule = ((-1) ** slot.a) * slot.M.sigma_ab * p_sign * ((-1) ** n_w)
    parity_conj = parity_conjugation_sign(slot)
    c0_sign = (-1) ** slot.a
    return {
        "id": slot.id,
        "X": slot.X,
        "a": slot.a,
        "M": slot.M.name,
        "coeff_type": slot.coeff_type,
        "p_f": slot.p_f,
        "sigma_Asigma_B": slot.M.sigma_ab,
        "N_w": n_w,
        "M_adjoint_class": slot.M.adjoint_class,
        "X_adjoint_sign": x_sign,
        "M_adjoint_sign": m_sign,
        "w_module_adjoint_sign": w_sign,
        "adjoint_product": product,
        "hermitian_retained": product == 1,
        "pi_P_rule": int(parity_rule),
        "pi_P_conjugation": int(parity_conj),
        "parity_methods_agree": int(parity_rule) == int(parity_conj),
        "P_w_invariant": int(parity_rule) == 1,
        "C0_sign": int(c0_sign),
        "ibp_representative": slot.id,
        "ibp_note": "singleton canonical operator class",
    }


def build_tables(a_values: Iterable[int] = (0, 1)) -> dict[str, Any]:
    slots = generate_slots(a_values)
    rows = [classify_slot(slot) for slot in slots]
    manifest_rows = [
        {
            "id": row["id"],
            "sigma_Asigma_B": row["sigma_Asigma_B"],
            "N_w": row["N_w"],
            "M_adjoint_class": row["M_adjoint_class"],
            "p_f": row["p_f"],
        }
        for row in rows
    ]
    adjoint_rows = [
        {
            "id": row["id"],
            "X_adjoint_sign": row["X_adjoint_sign"],
            "M_adjoint_sign": row["M_adjoint_sign"],
            "w_module_adjoint_sign": row["w_module_adjoint_sign"],
            "adjoint_product": row["adjoint_product"],
            "hermitian_retained": row["hermitian_retained"],
        }
        for row in rows
    ]
    hermitian_rows = [row for row in rows if row["hermitian_retained"]]
    parity_rows = [
        {
            "id": row["id"],
            "pi_P_rule": row["pi_P_rule"],
            "pi_P_conjugation": row["pi_P_conjugation"],
            "parity_methods_agree": row["parity_methods_agree"],
            "P_w_invariant": row["P_w_invariant"],
            "a_class": "combined-parity-mixing" if row["a"] == 1 and row["P_w_invariant"] else "parity-preserving",
        }
        for row in hermitian_rows
    ]
    physical_rows = [row for row in hermitian_rows if row["P_w_invariant"]]
    return {
        "slots": slots,
        "classified_rows": rows,
        "manifest_rows": manifest_rows,
        "adjoint_rows": adjoint_rows,
        "hermitian_rows": hermitian_rows,
        "parity_rows": parity_rows,
        "physical_rows": physical_rows,
    }


def module_on_vector(coeff_type: str, coeff: sp.Expr, vec: sp.Matrix, z: sp.Symbol, ell: sp.Symbol) -> sp.Matrix:
    def dw(expr: sp.Expr) -> sp.Expr:
        return sp.diff(expr, z) / ell

    if coeff_type == "C":
        return coeff * vec
    if coeff_type == "B":
        return coeff * vec.applyfunc(dw) + sp.Rational(1, 2) * dw(coeff) * vec
    if coeff_type == "K":
        return sp.Matrix([-dw(coeff * dw(entry)) for entry in vec])
    raise AssertionError(f"unknown coeff_type: {coeff_type}")


def integrate_w(expr: sp.Expr, z: sp.Symbol, ell: sp.Symbol) -> sp.Expr:
    """Integrate the fixed witness profiles over w.

    The witness expressions are finite sums of sech(z)^m tanh(z)^n after
    multiplying by dw=ell dz.  Direct hyperbolic integration is slow here, so
    use the beta-function identity
    int sech(z)^m tanh(z)^n dz = 0 for odd n, and
    Beta((n+1)/2, m/2) for even n.
    """
    s_sym, t_sym = sp.symbols("S_witness T_witness")
    raw = sp.trigsimp(sp.expand(expr * ell))
    raw = raw.xreplace({sp.tanh(z): t_sym, sp.cosh(z): 1 / s_sym, sp.sinh(z): t_sym / s_sym})
    raw = sp.expand(sp.powsimp(sp.cancel(raw), force=True))
    total = sp.Integer(0)
    for term in sp.Add.make_args(raw):
        if term == 0:
            continue
        powers = term.as_powers_dict()
        m_power = powers.get(s_sym, sp.Integer(0))
        n_power = powers.get(t_sym, sp.Integer(0))
        if not (m_power.is_integer and n_power.is_integer):
            raise AssertionError(f"non-integer witness powers in term {term}")
        m_int = int(m_power)
        n_int = int(n_power)
        coeff = sp.simplify(term / (s_sym**m_int * t_sym**n_int))
        if coeff.has(s_sym, t_sym):
            raise AssertionError(f"could not separate witness monomial {term} -> coeff {coeff}")
        if n_int % 2:
            continue
        if m_int <= 0:
            raise AssertionError(f"non-decaying witness monomial {term}")
        left = sp.Rational(n_int + 1, 2)
        right = sp.Rational(m_int, 2)
        total += coeff * sp.gamma(left) * sp.gamma(right) / sp.gamma(left + right)
    return sp.factor(sp.simplify(total))


def build_witnesses(physical_rows: list[dict[str, Any]], slots_by_id: dict[str, Slot]) -> dict[str, Any]:
    z = sp.symbols("z", real=True)
    ell = sp.symbols("ell", positive=True)
    sech = lambda arg: 1 / sp.cosh(arg)
    f0 = sp.Matrix([sech(z) ** 2 / ell, sech(z) ** 2 / ell, 0])
    wall = sp.Matrix([sech(z) * sp.tanh(z), sech(z) * sp.tanh(z), 0])
    tau = sp.Matrix([0, 0, 1])
    profiles = {
        "even": sech(z) ** 2,
        "odd": sp.tanh(z) * sech(z) ** 2,
    }
    profile_text = {
        "even": "sech(w/ell)^2",
        "odd": "tanh(w/ell)*sech(w/ell)^2",
    }

    pairs = {
        "f0_to_wall_partner": (f0, wall),
        "f0_to_tau": (f0, tau),
        "wall_partner_to_f0": (wall, f0),
    }
    a1_rows = [row for row in physical_rows if row["a"] == 1]
    witness_rows: list[dict[str, Any]] = []
    exprs: dict[str, sp.Expr] = {}

    for index, row in enumerate(a1_rows):
        slot = slots_by_id[row["id"]]
        coeff = profiles[slot.p_f]
        matrix = slot.M.matrix
        finite_values: dict[str, sp.Expr] = {}
        static_values: dict[str, sp.Expr] = {}
        generic_functionals: dict[str, str] = {}
        for pair_name, (bra, ket) in pairs.items():
            acted = matrix * module_on_vector(slot.coeff_type, coeff, ket, z, ell)
            overlap = integrate_w((bra.T * acted)[0], z, ell)
            finite_values[pair_name] = overlap
            static_value = sp.Integer(0) if slot.X == "omegaVdotk" else overlap
            static_values[pair_name] = static_value
            exprs[f"witness_{index:03d}_{pair_name}_finite"] = overlap
            exprs[f"witness_{index:03d}_{pair_name}_static"] = static_value
            generic_functionals[pair_name] = (
                f"{slot.X}*s^{slot.a}*Int[bra={pair_name.split('_to_')[0]}, "
                f"M={slot.M.name}, module={slot.coeff_type}, coeff_parity={slot.p_f}, "
                f"ket={pair_name.split('_to_')[-1]}]"
            )
        nonzero_pairs = [name for name, value in static_values.items() if sp.simplify(value) != 0]
        finite_nonzero_pairs = [name for name, value in finite_values.items() if sp.simplify(value) != 0]
        burden = slot.X == "iVdotk"
        if burden and nonzero_pairs:
            witness_status = "NONZERO_STATIC_WITNESS"
        elif burden:
            witness_status = "STRUCTURAL_ZERO_FOR_IMPORTED_TARGETS"
        else:
            witness_status = "OMEGA_STATIC_MODE_NO_WITNESS_BURDEN"
        witness_rows.append(
            {
                "slot_id": row["id"],
                "X": slot.X,
                "a": slot.a,
                "M": slot.M.name,
                "coeff_type": slot.coeff_type,
                "p_f": slot.p_f,
                "witness_profile": profile_text[slot.p_f],
                "functionals": generic_functionals,
                "finite_frequency_w_overlap_values": {key: hstr(value) for key, value in finite_values.items()},
                "static_omega0_values": {key: hstr(value) for key, value in static_values.items()},
                "static_nonzero_pairs": nonzero_pairs,
                "finite_frequency_nonzero_pairs": finite_nonzero_pairs,
                "witness_status": witness_status,
                "continuum_threshold_flag": "not_excited_by_native_discrete_witness",
            }
        )

    decisive = [
        row
        for row in witness_rows
        if row["X"] == "iVdotk" and row["witness_status"] == "NONZERO_STATIC_WITNESS"
    ]
    return {
        "a1_witness_rows": witness_rows,
        "decisive_static_witness_rows": decisive,
        "exprs": exprs,
    }


def neutral_overlap_control_expr(slots_by_id: dict[str, Slot]) -> dict[str, Any]:
    slot_id = "iVdotk|0|etaeta|B|odd"
    slot = slots_by_id[slot_id]
    z = sp.symbols("z", real=True)
    ell = sp.symbols("ell", positive=True)
    sech = lambda arg: 1 / sp.cosh(arg)
    f0 = sp.Matrix([sech(z) ** 2 / ell, sech(z) ** 2 / ell, 0])
    wall = sp.Matrix([sech(z) * sp.tanh(z), sech(z) * sp.tanh(z), 0])
    coeff = sp.tanh(z) * sech(z) ** 2
    acted = slot.M.matrix * module_on_vector(slot.coeff_type, coeff, wall, z, ell)
    value = integrate_w((f0.T * acted)[0], z, ell)
    return {
        "slot_id": slot_id,
        "profile": "tanh(w/ell)*sech(w/ell)^2",
        "overlap": value,
        "overlap_string": hstr(value),
    }


def decide_case(
    *,
    name: str,
    allowed_a: tuple[int, ...] = (0, 1),
    c0_symmetry: bool = False,
    delete_id: str | None = None,
    illegal_input: str | None = None,
    v_zero: bool = False,
    witness_by_id: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    if illegal_input is not None:
        return illegal_control_result(name, illegal_input)

    tables = build_tables(allowed_a)
    generated_ids = [row["id"] for row in tables["classified_rows"]]
    expected_ids = generated_ids.copy()
    if delete_id is not None:
        generated_ids = [slot_id for slot_id in generated_ids if slot_id != delete_id]

    missing = sorted(set(expected_ids) - set(generated_ids))
    extra = sorted(set(generated_ids) - set(expected_ids))
    if missing or extra:
        return {
            "name": name,
            "status": "FIRED",
            "verdict": "BASIS_INCOMPLETE",
            "expected_raw_count": len(expected_ids),
            "generated_raw_count": len(generated_ids),
            "missing_ids": missing,
            "extra_ids": extra,
        }

    hermitian = [row for row in tables["hermitian_rows"] if row["id"] in set(generated_ids)]
    parity_ok = [row for row in hermitian if row["P_w_invariant"]]
    c0_ok = [row for row in parity_ok if (not c0_symmetry or row["C0_sign"] == 1)]
    active = [] if v_zero else c0_ok
    a0 = [row for row in active if row["a"] == 0]
    a1 = [row for row in active if row["a"] == 1]
    nonzero_a1 = []
    if witness_by_id:
        nonzero_a1 = [
            row
            for row in a1
            if witness_by_id.get(row["id"], {}).get("witness_status") == "NONZERO_STATIC_WITNESS"
        ]

    if v_zero:
        verdict = "CLEAN_STATIC_SPLIT"
    elif c0_symmetry and not a1:
        verdict = "PASS_CONDITIONAL_ON_NO_OPERATOR_PARITY_MIXING"
    elif not a1:
        verdict = "CLEAN_PARITY_PRESERVING"
    elif nonzero_a1:
        verdict = "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING"
    else:
        verdict = "ACCIDENTAL_ZERO_SIM_DEFERRED"

    result = {
        "name": name,
        "status": "FIRED",
        "verdict": verdict,
        "verdict_code": verdict_code(verdict),
        "raw_count": len(generated_ids),
        "hermitian_count": len(hermitian),
        "physical_count": len(c0_ok),
        "active_operator_count": len(active),
        "a0_count": len(a0),
        "a1_count": len(a1),
        "a1_nonzero_static_witness_count": len(nonzero_a1),
        "a1_ids": [row["id"] for row in a1],
        "a1_nonzero_static_witness_ids": [row["id"] for row in nonzero_a1],
    }
    if c0_symmetry:
        result["c0_label_basis"] = (
            "observed_post_C0_filter_a1_count_eq_0"
            if len(a1) == 0
            else "observed_post_C0_filter_a1_count_nonzero"
        )
    return result


def build_symbolics() -> dict[str, Any]:
    tables = build_tables((0, 1))
    slots_by_id = {slot.id: slot for slot in tables["slots"]}
    witnesses = build_witnesses(tables["physical_rows"], slots_by_id)
    witness_by_id = {row["slot_id"]: row for row in witnesses["a1_witness_rows"]}
    neutral_control = neutral_overlap_control_expr(slots_by_id)

    parity_mismatches = [row for row in tables["hermitian_rows"] if not row["parity_methods_agree"]]
    if parity_mismatches:
        raise AssertionError(f"parity two-method mismatch: {[row['id'] for row in parity_mismatches]}")

    main = decide_case(name="charged_moving_throat", witness_by_id=witness_by_id)
    neutral = decide_case(name="neutral_moving_throat", allowed_a=(0,), witness_by_id=witness_by_id)
    c0 = decide_case(name="injected_C0", c0_symmetry=True, witness_by_id=witness_by_id)
    deleted = decide_case(
        name="deleted_operator",
        delete_id=tables["classified_rows"][0]["id"],
        witness_by_id=witness_by_id,
    )
    extra_odd = decide_case(name="injected_extra_odd_spurion", illegal_input="extra_odd_spurion")
    smuggled_vector = decide_case(name="smuggled_V_i_n_i", illegal_input="V_i n_i")
    smuggled_axial = decide_case(name="smuggled_epsilon_ijk_V_i_k_j", illegal_input="epsilon_ijk V_i k_j")
    v_zero = decide_case(name="V_equals_zero", v_zero=True, witness_by_id=witness_by_id)

    c0_removed_restores = {
        "name": "C0_removed_restores_a1",
        "status": "FIRED" if main["a1_count"] > 0 and c0["a1_count"] == 0 else "NOT_FIRED",
        "verdict": "A1_IDS_RESTORED" if main["a1_count"] > 0 and c0["a1_count"] == 0 else "RESTORE_FAILED",
        "without_C0_a1_count": main["a1_count"],
        "with_C0_a1_count": c0["a1_count"],
    }
    parity_two_ways = {
        "name": "parity_two_ways",
        "status": "FIRED" if len(parity_mismatches) == 0 else "NOT_FIRED",
        "verdict": "PARITY_TWO_WAYS_AGREE" if len(parity_mismatches) == 0 else "PARITY_TWO_WAYS_DISAGREE",
        "checked_retained_count": len(tables["hermitian_rows"]),
    }

    controls = {
        "neutral_moving_throat": {
            **neutral,
            "required_overlap_exercised": {
                "slot_id": neutral_control["slot_id"],
                "profile": neutral_control["profile"],
                "overlap": neutral_control["overlap_string"],
                "status": "FIRED" if sp.simplify(neutral_control["overlap"]) == 0 else "NOT_FIRED",
            },
        },
        "charged_moving_throat": main,
        "injected_C0": c0,
        "C0_removed_restores_a1": c0_removed_restores,
        "injected_extra_odd_spurion": extra_odd,
        "deleted_operator": deleted,
        "smuggled_V_i_n_i": smuggled_vector,
        "smuggled_epsilon_ijk_V_i_k_j": smuggled_axial,
        "parity_two_ways": parity_two_ways,
        "V_equals_zero": v_zero,
    }
    for name, item in controls.items():
        if item["status"] != "FIRED":
            raise AssertionError(f"control did not fire: {name}: {item}")

    physical = tables["physical_rows"]
    a0 = [row for row in physical if row["a"] == 0]
    a1 = [row for row in physical if row["a"] == 1]
    count_summary = {
        "raw_manifest_count": len(tables["classified_rows"]),
        "expected_raw_manifest_count": EXPECTED_RAW_COUNT_FULL,
        "hermitian_retained_pre_ibp": len(tables["hermitian_rows"]),
        "hermitian_retained_post_ibp": len({row["ibp_representative"] for row in tables["hermitian_rows"]}),
        "parity_retained_pre_ibp": len(physical),
        "parity_retained_post_ibp": len({row["ibp_representative"] for row in physical}),
        "a0_physical_count": len(a0),
        "a1_physical_count": len(a1),
        "a1_static_nonzero_witness_count": len(witnesses["decisive_static_witness_rows"]),
    }
    if count_summary["raw_manifest_count"] != EXPECTED_RAW_COUNT_FULL:
        raise AssertionError("raw manifest count mismatch")
    expected_landing_observation = {
        "expected": "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING",
        "computed": main["verdict"],
        "status": (
            "OBSERVED_EXPECTED_FAIL"
            if main["verdict"] == "FAIL_UNPROTECTED_OPERATOR_PARITY_MIXING"
            else "COMPUTED_NONFAIL_EMITTED"
        ),
    }

    exprs_for_agreement: dict[str, sp.Expr] = {
        "raw_manifest_count": sp.Integer(count_summary["raw_manifest_count"]),
        "hermitian_retained_post_ibp": sp.Integer(count_summary["hermitian_retained_post_ibp"]),
        "parity_retained_post_ibp": sp.Integer(count_summary["parity_retained_post_ibp"]),
        "a0_physical_count": sp.Integer(count_summary["a0_physical_count"]),
        "a1_physical_count": sp.Integer(count_summary["a1_physical_count"]),
        "a1_static_nonzero_witness_count": sp.Integer(count_summary["a1_static_nonzero_witness_count"]),
        "parity_mismatch_count": sp.Integer(len(parity_mismatches)),
        "neutral_control_overlap": neutral_control["overlap"],
        "main_verdict_code": sp.Integer(verdict_code(main["verdict"])),
        "neutral_verdict_code": sp.Integer(verdict_code(neutral["verdict"])),
        "c0_verdict_code": sp.Integer(verdict_code(c0["verdict"])),
        "deleted_verdict_code": sp.Integer(verdict_code(deleted["verdict"])),
        "extra_odd_verdict_code": sp.Integer(verdict_code(extra_odd["verdict"])),
        "smuggled_vector_verdict_code": sp.Integer(verdict_code(smuggled_vector["verdict"])),
        "smuggled_axial_verdict_code": sp.Integer(verdict_code(smuggled_axial["verdict"])),
        "v_zero_verdict_code": sp.Integer(verdict_code(v_zero["verdict"])),
    }
    exprs_for_agreement.update(witnesses["exprs"])

    return {
        "tables": tables,
        "count_summary": count_summary,
        "witnesses": witnesses,
        "neutral_control": neutral_control,
        "controls": controls,
        "main": main,
        "expected_landing_observation": expected_landing_observation,
        "exprs_for_agreement": exprs_for_agreement,
    }


def control_label(item: dict[str, Any]) -> str:
    return str(item.get("verdict", item.get("status", "UNKNOWN")))


def build_payload() -> dict[str, Any]:
    imports = import_banked_inputs()
    sym = build_symbolics()
    tables = sym["tables"]
    count_summary = sym["count_summary"]
    witnesses = sym["witnesses"]
    agreement_exprs = {key: mma_expr(value) for key, value in sym["exprs_for_agreement"].items()}
    expr_digest = sha256_json(agreement_exprs)

    physical_rows = tables["physical_rows"]
    a0_ids = [row["id"] for row in physical_rows if row["a"] == 0]
    a1_ids = [row["id"] for row in physical_rows if row["a"] == 1]
    agreement_payload = {
        "top_line_verdict": sym["main"]["verdict"],
        "main_verdict_code": verdict_code(sym["main"]["verdict"]),
        "counts": count_summary,
        "manifest_ids": [row["id"] for row in tables["classified_rows"]],
        "hermitian_ids": [row["id"] for row in tables["hermitian_rows"]],
        "physical_ids": [row["id"] for row in physical_rows],
        "a0_ids": a0_ids,
        "a1_ids": a1_ids,
        "adjoint_table": tables["adjoint_rows"],
        "parity_table": tables["parity_rows"],
        "witness_status_table": [
            {
                "slot_id": row["slot_id"],
                "witness_status": row["witness_status"],
                "static_nonzero_pairs": row["static_nonzero_pairs"],
                "finite_frequency_nonzero_pairs": row["finite_frequency_nonzero_pairs"],
            }
            for row in witnesses["a1_witness_rows"]
        ],
        "control_verdicts": {name: control_label(item) for name, item in sym["controls"].items()},
        "checked_expression_count": len(sym["exprs_for_agreement"]),
        "expr_digest": expr_digest,
    }

    return {
        "schema": SCHEMA,
        "engine": "sympy",
        "directive": "software/stage1_solver/directives/pathA_39_stage3_operator_parity.md",
        "top_line_verdict": sym["main"]["verdict"],
        "expected_landing_observation": sym["expected_landing_observation"],
        "engine_agreement": {
            "status": "PENDING_MATHEMATICA",
            "mathematica_exprs": agreement_exprs,
            "sympy_expression_digest": expr_digest,
        },
        "imports": imports,
        "grammar": {
            "fields": list(FIELDS),
            "field_sigma": FIELD_SIGMA,
            "external_factors": list(X_FACTORS),
            "a_values": [0, 1],
            "matrix_basis": [
                {
                    "name": M.name,
                    "kind": M.kind,
                    "sigma_Asigma_B": M.sigma_ab,
                    "adjoint_sign": M.adjoint_sign,
                }
                for M in matrix_defs()
            ],
            "coeff_types": W_MODULE,
            "p_f_values": PF_SIGNS,
            "expected_raw_count": EXPECTED_RAW_COUNT_FULL,
        },
        "manifest": {
            "summary": sym["count_summary"],
            "raw_table": tables["manifest_rows"],
            "raw_ids": [row["id"] for row in tables["classified_rows"]],
            "ibp_quotient": {
                "policy": "canonical C, A_B, and -D K D modules are already in slab-IBP normal form; no free-coefficient singleton slots merge",
                "raw_class_count": len(tables["classified_rows"]),
                "hermitian_class_count": sym["count_summary"]["hermitian_retained_post_ibp"],
                "physical_class_count": sym["count_summary"]["parity_retained_post_ibp"],
            },
        },
        "adjoint_table": tables["adjoint_rows"],
        "parity_table": tables["parity_rows"],
        "physical_slots": {
            "a0_parity_preserving_ids": a0_ids,
            "a1_combined_parity_mixing_ids": a1_ids,
        },
        "witnesses": {
            "a1_mixing_slots": witnesses["a1_witness_rows"],
            "decisive_static_iVdotk_witness_slots": witnesses["decisive_static_witness_rows"],
            "omega_static_note": "Imported native modes are omega=0; omegaVdotk rows are classified but carry no static witness burden.",
            "neutral_control_overlap": {
                "slot_id": sym["neutral_control"]["slot_id"],
                "profile": sym["neutral_control"]["profile"],
                "overlap": sym["neutral_control"]["overlap_string"],
            },
        },
        "controls": sym["controls"],
        "dimensional_firewall": {
            "allowed_external_scalars": list(X_FACTORS),
            "illegal_inputs_rejected": ["V_i n_i", "epsilon_ijk V_i k_j", "extra P_w-odd spurion distinct from s"],
            "status": "PASS",
        },
        "remediation": {
            "fix_1": "production path emits the computed verdict; the expected fail is recorded as an observation, not asserted",
            "fix_2": "illegal controls attempt admission through the declared grammar validator and derive rejection from refusal",
            "fix_3": "injected_C0 PASS_CONDITIONAL is based on observed post-C0-filter a=1 count equal to zero",
            "fix_4": "P O P^-1 check uses concrete coefficient profiles, explicit D_w tokens, s->-s, and Sigma M Sigma",
        },
        "agreement_payload": agreement_payload,
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def compare_with_mathematica(payload: dict[str, Any]) -> dict[str, Any]:
    if not WL_OUT.exists():
        raise FileNotFoundError(f"Mathematica payload missing: {WL_OUT}")
    other = json.loads(WL_OUT.read_text(encoding="utf-8"))
    if other.get("schema") != WL_SCHEMA:
        raise AssertionError(f"unexpected Mathematica schema: {other.get('schema')}")
    if other.get("status") != "OK":
        raise AssertionError(f"Mathematica engine did not complete cleanly: {other}")
    if other.get("sympy_expression_digest") != payload["engine_agreement"]["sympy_expression_digest"]:
        raise AssertionError("Mathematica compared a stale SymPy expression digest")
    if other.get("agreement_payload") != payload["agreement_payload"]:
        raise AssertionError(
            "ENGINE_DISAGREE\n"
            + json.dumps({"sympy": payload["agreement_payload"], "mathematica": other.get("agreement_payload")}, indent=2)
        )

    payload["engine_agreement"]["status"] = "ENGINE_AGREE"
    payload["engine_agreement"]["mathematica_payload"] = str(WL_OUT)
    payload["engine_agreement"]["sympy_payload"] = str(SYM_OUT)

    results = {
        "schema": SCHEMA,
        "verdict": payload["top_line_verdict"],
        "expected_landing_observation": payload["expected_landing_observation"],
        "headline": {
            "verdict": payload["top_line_verdict"],
            "engine_agreement": "ENGINE_AGREE",
            "raw_manifest_count": payload["manifest"]["summary"]["raw_manifest_count"],
            "hermitian_retained": payload["manifest"]["summary"]["hermitian_retained_post_ibp"],
            "parity_retained": payload["manifest"]["summary"]["parity_retained_post_ibp"],
            "a0_physical_count": payload["manifest"]["summary"]["a0_physical_count"],
            "a1_physical_count": payload["manifest"]["summary"]["a1_physical_count"],
            "a1_static_nonzero_witness_count": payload["manifest"]["summary"]["a1_static_nonzero_witness_count"],
        },
        "engine_agreement": {
            "status": "ENGINE_AGREE",
            "compared_payload": payload["agreement_payload"],
            "sympy_payload": str(SYM_OUT),
            "mathematica_payload": str(WL_OUT),
        },
        "imports": payload["imports"],
        "grammar": payload["grammar"],
        "manifest": payload["manifest"],
        "adjoint_table": payload["adjoint_table"],
        "parity_table": payload["parity_table"],
        "physical_slots": payload["physical_slots"],
        "witnesses": payload["witnesses"],
        "controls": payload["controls"],
        "dimensional_firewall": payload["dimensional_firewall"],
        "remediation": payload["remediation"],
        "honest_scope": {
            "derived": [
                "finite O(V^1) canonical slot manifest",
                "Hermiticity table under Convention A",
                "P_w parity by rule and by literal P O P^-1 representative conjugation",
                "native-mode symbolic witness overlaps for a=1 slots",
                "controls that fire in pass, fail, rejection, and basis-incomplete directions",
            ],
            "sim_deferred": payload["imports"]["sim_deferred"],
            "does_not_unearn": [
                "pathA_38 Coulomb localization and source-level odd/even split",
                "pathA_39 Stage 0+1 scalar-admixture floor",
                "pathA_39 Stage 2 magnetic-force derivation",
            ],
        },
    }
    YAML_OUT.write_text(yaml.safe_dump(results, sort_keys=False, width=140), encoding="utf-8")
    write_report(results)
    write_json(SYM_OUT, payload)
    return results


def table_bool(value: bool) -> str:
    return "yes" if value else "no"


def write_report(results: dict[str, Any]) -> None:
    summary = results["manifest"]["summary"]
    controls = results["controls"]
    imports = results["imports"]["imported"]
    declared = results["imports"]["declared"]
    sim_deferred = results["imports"]["sim_deferred"]

    lines: list[str] = [
        "# pathA_39 Stage 3 Operator Parity Under Motion",
        "",
        f"Computed headline: `{results['verdict']}` with dual-engine `ENGINE_AGREE`.",
        "",
        "Remediation note: fixes 1-4 are applied.  The production verdict is emitted from `decide_case`, illegal controls are validator-observed rejections, the `C0` label follows the observed post-filter `a=1` count, and the conjugation column uses a literal representative `P O P^-1` computation.",
        "",
        "The generated operator module contains allowed `a=1` O(V) slots.  The decisive static burden is carried by `iVdotk` slots, and at least one parity-respecting witness profile gives a nonzero imported-mode overlap.  The realized nonlinear-throat coefficient and magnitude remain sim-deferred.",
        "",
        "## Manifest Summary",
        "",
        "| quantity | count |",
        "|---|---:|",
        f"| raw manifest | `{summary['raw_manifest_count']}` |",
        f"| retained after Hermiticity, post-IBP | `{summary['hermitian_retained_post_ibp']}` |",
        f"| retained after P_w parity, post-IBP | `{summary['parity_retained_post_ibp']}` |",
        f"| `a=0` parity-preserving | `{summary['a0_physical_count']}` |",
        f"| `a=1` combined-parity-mixing | `{summary['a1_physical_count']}` |",
        f"| `a=1` static nonzero `iVdotk` witnesses | `{summary['a1_static_nonzero_witness_count']}` |",
        "",
        f"IBP quotient: {results['manifest']['ibp_quotient']['policy']}.  The quotient counts are unchanged: raw `{results['manifest']['ibp_quotient']['raw_class_count']}`, Hermitian `{results['manifest']['ibp_quotient']['hermitian_class_count']}`, physical `{results['manifest']['ibp_quotient']['physical_class_count']}`.",
        "",
        "## Full Raw Manifest",
        "",
        "| ID | sigma_Asigma_B | N_w | adjoint class | p_f |",
        "|---|---:|---:|---|---|",
    ]
    for row in results["manifest"]["raw_table"]:
        lines.append(
            f"| `{row['id']}` | `{sign_str(row['sigma_Asigma_B'])}` | `{row['N_w']}` | `{row['M_adjoint_class']}` | `{row['p_f']}` |"
        )

    lines.extend(
        [
            "",
            "## Adjoint Table",
            "",
            "| ID | X adj | M adj | w-module adj | product | retained |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in results["adjoint_table"]:
        lines.append(
            f"| `{row['id']}` | `{sign_str(row['X_adjoint_sign'])}` | `{sign_str(row['M_adjoint_sign'])}` | `{sign_str(row['w_module_adjoint_sign'])}` | `{sign_str(row['adjoint_product'])}` | `{table_bool(row['hermitian_retained'])}` |"
        )

    lines.extend(
        [
            "",
            "## Parity Two Ways",
            "",
            "The conjugation path builds a concrete representative operator for each retained slot, applies `z -> -z`, `D_w -> -D_w`, `s -> -s`, and `Sigma M Sigma`, then reads the sign by comparing the substituted operator to the original.  It does not reuse the rule's `p_f`, `N_w`, or `(-1)^a` factors.",
            "",
            "| Hermitian ID | rule pi_P | conjugation pi_P | invariant | class |",
            "|---|---:|---:|---:|---|",
        ]
    )
    for row in results["parity_table"]:
        lines.append(
            f"| `{row['id']}` | `{sign_str(row['pi_P_rule'])}` | `{sign_str(row['pi_P_conjugation'])}` | `{table_bool(row['P_w_invariant'])}` | `{row['a_class']}` |"
        )

    lines.extend(
        [
            "",
            "## a=1 Witnesses",
            "",
            "For each `a=1` physical slot, the table gives the allowed witness coefficient profile and the symbolic static `omega=0` native-mode overlaps.  `omegaVdotk` slots are still classified, but carry no static witness burden because the imported modes have `omega=0`.",
            "",
            "| slot | profile | <f0|O|wall_partner> | <f0|O|tau> | <wall_partner|O|f0> | status |",
            "|---|---|---:|---:|---:|---|",
        ]
    )
    for row in results["witnesses"]["a1_mixing_slots"]:
        vals = row["static_omega0_values"]
        lines.append(
            f"| `{row['slot_id']}` | `{row['witness_profile']}` | `{vals['f0_to_wall_partner']}` | `{vals['f0_to_tau']}` | `{vals['wall_partner_to_f0']}` | `{row['witness_status']}` |"
        )

    lines.extend(
        [
            "",
            "Decisive nonzero static witness IDs:",
        ]
    )
    for row in results["witnesses"]["decisive_static_iVdotk_witness_slots"]:
        lines.append(f"- `{row['slot_id']}` via {', '.join(row['static_nonzero_pairs'])}")

    lines.extend(
        [
            "",
            "## Controls",
            "",
            "| control | status | verdict/result |",
            "|---|---:|---|",
        ]
    )
    for name, item in controls.items():
        lines.append(f"| `{name}` | `{item['status']}` | `{control_label(item)}` |")
    neutral_overlap = controls["neutral_moving_throat"]["required_overlap_exercised"]
    lines.extend(
        [
            "",
            f"Neutral-throat overlap machinery exercised on `{neutral_overlap['slot_id']}` with profile `{neutral_overlap['profile']}`: overlap `{neutral_overlap['overlap']}`.",
            f"`injected_C0` verdict basis: `{controls['injected_C0'].get('c0_label_basis', 'not-applicable')}` with observed post-C0 `a=1` count `{controls['injected_C0'].get('a1_count', 'NA')}`.",
            "Illegal-input controls attempted admission through the declared grammar validator and recorded `REFUSED_BY_DECLARED_GRAMMAR` before emitting their reject verdicts.",
            "",
            "## Remediation",
            "",
            "| fix | status |",
            "|---|---|",
        ]
    )
    for key, value in results["remediation"].items():
        lines.append(f"| `{key}` | {value} |")
    lines.extend(
        [
            "",
            "## Provenance",
            "",
            "Imported:",
        ]
    )
    for key, value in imports.items():
        lines.append(f"- `{key}`: {value}")
    lines.append("")
    lines.append("Declared:")
    for key, value in declared.items():
        if isinstance(value, list):
            lines.append(f"- `{key}`: {'; '.join(value)}")
        else:
            lines.append(f"- `{key}`: {value}")
    lines.append("")
    lines.append("Sim-deferred:")
    for item in sim_deferred:
        lines.append(f"- {item}")

    lines.extend(
        [
            "",
            "## Honest Scope",
            "",
            "Stage 3 establishes allowed-and-unprotected O(V) operator parity mixing with a nonzero symbolic witness.  It does not compute the realized finite-mouth coefficient or contamination magnitude, and it does not undo the pathA_38 source-level Coulomb localization or the Stage 0+1/2 results.",
            "",
            "## Dual Engine",
            "",
            f"`ENGINE_AGREE` over `{results['engine_agreement']['compared_payload']['checked_expression_count']}` compared quantities, including manifest IDs/counts, adjoint rows, parity rows, a=0/a=1 split, witness overlaps, and control verdicts.",
            "",
            "Run commands:",
            "",
            "```text",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_stage3_operator_parity_sympy.py",
            "timeout 600 math -script software/stage1_solver/tools/pathA_39_stage3_operator_parity.wl",
            "timeout 600 python3 software/stage1_solver/tools/pathA_39_stage3_operator_parity_sympy.py --compare",
            "```",
        ]
    )
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare", action="store_true", help="compare with Mathematica payload and write YAML/report")
    args = parser.parse_args()

    payload = build_payload()
    write_json(SYM_OUT, payload)

    if args.compare:
        results = compare_with_mathematica(payload)
        print(json.dumps({"engine": "sympy", "status": "ENGINE_AGREE", "verdict": results["verdict"]}, sort_keys=True))
    else:
        print(json.dumps({"engine": "sympy", "status": "OK", "verdict": payload["top_line_verdict"], "json": str(SYM_OUT)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
