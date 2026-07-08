#!/usr/bin/env python3
"""Ledger stage007 SymPy audit: pathA_35 G0 shear-surface freeze.

Standalone, print-only, no arguments, no file output.  This reshapes the
pathA_35 G0 checker into a ledger citizen: byte-fidelity of the frozen action,
the full dimensional firewall, rank-computed flat-brane DOF, and the drift
count computed from an enumerated table.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
from typing import Any, Iterable

import sympy as sp


PASS_COUNT = 0
FAIL_COUNT = 0


class AuditFailure(AssertionError):
    pass


class MissingReport(AuditFailure):
    pass


class MissingFence(AuditFailure):
    pass


class AmbiguousFence(AuditFailure):
    pass


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
REPORT_ROOT = REPO_ROOT / "software" / "stage1_solver" / "reports"
T0_REPORT = REPORT_ROOT / "pathA_24_T0_freeze.md"
G0_REPORT = REPORT_ROOT / "pathA_35_G0_freeze.md"

EXPECTED_T0_HASH = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064"
EXPECTED_G0_HASH = "d9520d3819c3f718290f9d0be57138c07d5bf02d2237106478e17b6a1e389ac3"
EXPECTED_G0_SHORT = "d9520d3819c3"
EXPECTED_DRIFT_TOKEN = "SECOND_MEDIUM_DRIFT_AT_FREEZE(11)"
EXPECTED_FREEZE_TOKEN = f"T0_SHEAR_FROZEN({EXPECTED_G0_SHORT})"
DRIFT_DIM_TARGET_KEYS = {
    "rho_br": "rho_br",
    "mu_R": "mu_R",
    "lambda_Pu": "lambda_Pu",
    "Omega_w": "Omega_w",
    "g_ell(w)": "g_ell",
}


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


def expect_zero(name: str, residual: sp.Expr | int) -> None:
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


def expect_nonzero(name: str, residual: sp.Expr | int) -> None:
    assert_no_float(name, residual)
    clean = sp.simplify(residual)
    assert_no_float(name, clean)
    if clean != 0:
        _record_pass(f"PASS  {name} is nonzero as required (residual = {compact(clean)})")
        return
    _record_fail(f"FAIL  {name}: required nonzero residual vanished")
    raise AuditFailure(f"{name} unexpectedly had zero residual")


def expect_fail(name: str, residual: sp.Expr | int) -> None:
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


def dim_residual(actual: Dim, expected: Dim) -> sp.Expr:
    return sp.simplify(
        sum((have - want) ** 2 for have, want in zip(actual.components(), expected.components()))
    )


def expect_dim(name: str, actual: Dim, expected: Dim) -> None:
    expect_zero(name, dim_residual(actual, expected))


def homogeneity_residual(terms: dict[str, Dim]) -> sp.Expr:
    if not terms:
        raise AuditFailure("homogeneity check requires at least one term")
    dims = list(terms.values())
    reference = dims[0]
    return sp.simplify(sum(dim_residual(actual, reference) for actual in dims[1:]))


def matrix_rank(matrix: sp.Matrix) -> int:
    return int(matrix.rank())


@dataclass(frozen=True)
class ExtractedBlock:
    path: Path
    label: str
    block: bytes
    start_0: int
    end_0: int

    @property
    def start_1(self) -> int:
        return self.start_0 + 1

    @property
    def end_1(self) -> int:
        return self.end_0

    @property
    def length(self) -> int:
        return self.end_0 - self.start_0


def extract_fence_bytes(path: Path, label: str) -> ExtractedBlock:
    if not path.exists():
        raise MissingReport(f"missing report: {path}")

    opening = f"```{label}\n".encode("utf-8")
    closing = b"```\n"
    data = path.read_bytes()
    lines = data.splitlines(keepends=True)

    blocks: list[ExtractedBlock] = []
    in_block = False
    current: list[bytes] = []
    start = -1
    offset = 0

    for line in lines:
        if not in_block and line == opening:
            in_block = True
            current = []
            start = offset + len(line)
            offset += len(line)
            continue
        if in_block and line == closing:
            end = offset
            blocks.append(ExtractedBlock(path=path, label=label, block=b"".join(current), start_0=start, end_0=end))
            in_block = False
            current = []
            start = -1
            offset += len(line)
            continue
        if in_block:
            current.append(line)
        offset += len(line)

    if in_block:
        raise MissingFence(f"unterminated {label!r} fence in {path}")
    if not blocks:
        raise MissingFence(f"missing {label!r} fence in {path}")
    if len(blocks) > 1:
        raise AmbiguousFence(f"ambiguous {label!r} fence in {path}: found {len(blocks)} blocks")
    return blocks[0]


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_equality_residual(data: bytes, expected: str) -> sp.Integer:
    return sp.Integer(0) if sha256_hex(data) == expected else sp.Integer(1)


def run_freeze_fidelity() -> tuple[ExtractedBlock, ExtractedBlock]:
    subbanner("Freeze-fidelity byte audit")
    t0 = extract_fence_bytes(T0_REPORT, "freeze-action")
    g0 = extract_fence_bytes(G0_REPORT, "freeze-action")

    t0_hash = sha256_hex(t0.block)
    g0_hash = sha256_hex(g0.block)
    print(f"  T0 report: {t0.path}")
    print(f"  G0 report: {g0.path}")
    print(
        "  G0 freeze-action byte range (informative only): "
        f"0-based [{g0.start_0},{g0.end_0}), 1-based {g0.start_1}-{g0.end_1}, length {g0.length}"
    )
    print(
        "  T0 freeze-action byte range (informative only): "
        f"0-based [{t0.start_0},{t0.end_0}), 1-based {t0.start_1}-{t0.end_1}, length {t0.length}"
    )

    expect_bool("G0 short hash is prefix of frozen full SHA-256 constant", EXPECTED_G0_HASH.startswith(EXPECTED_G0_SHORT))
    expect_zero("T0 freeze-action SHA-256 matches frozen report", hash_equality_residual(t0.block, EXPECTED_T0_HASH))
    expect_zero("G0 freeze-action SHA-256 matches frozen report", hash_equality_residual(g0.block, EXPECTED_G0_HASH))
    expect_bool("byte-identical T0 freeze-action block is embedded inside G0 block", t0.block in g0.block)

    corrupted = bytearray(g0.block)
    corrupted[0] ^= 1
    expect_fail(
        "hash tooth: one-byte in-memory G0 corruption trips SHA-256 mismatch",
        hash_equality_residual(bytes(corrupted), EXPECTED_G0_HASH),
    )
    try:
        extract_fence_bytes(G0_REPORT, "nonexistent-freeze-action")
    except MissingFence:
        missing_fence_residual = sp.Integer(1)
    else:
        missing_fence_residual = sp.Integer(0)
    expect_fail("hash tooth: nonexistent fence tag trips extractor missing-fence path", missing_fence_residual)
    return t0, g0


def run_dimensional_firewall() -> dict[str, Dim]:
    subbanner("Dimensional firewall over the full frozen surface")
    bulk_lag = MASS / ((LENGTH**2) * (TIME**2))
    brane_lag = MASS / (LENGTH * (TIME**2))
    action_dim = MASS * (LENGTH**2) / TIME
    stress = bulk_lag
    eom_u_op = MASS / ((LENGTH**3) * (TIME**2))

    d_m = MASS
    d_hbar = action_dim
    d_rho = LENGTH**-4
    d_K = MASS * (LENGTH**18) * (TIME**-2)
    d_a = LENGTH
    d_u = LENGTH
    d_P = DIMENSIONLESS
    d_w = LENGTH
    d_ell_g = LENGTH
    d_g = LENGTH**-1
    d_grad = LENGTH**-1
    d_dt = TIME**-1
    d_dt_measure = TIME
    d_d4x = LENGTH**4
    d_v = LENGTH / TIME
    d_k = LENGTH**-1
    d_omega = TIME**-1
    d_rho_br = MASS * (LENGTH**-3)
    d_mu_R = brane_lag
    d_lambda_Pu = brane_lag
    d_Omega_w = d_omega

    d_cs2 = d_K * (d_rho**4) / d_m
    expect_dim("kept GNLS c_s^2(rho)=5 K rho^4/m", d_cs2, (LENGTH**2) * (TIME**-2))
    expect_dim("kept GNLS U(rho)=(K/4)rho^5", d_K * (d_rho**5), bulk_lag)
    expect_dim(
        "kept GNLS quantum pressure (hbar^2/(8 m rho))(partial_i rho)^2",
        (d_hbar**2) / (d_m * d_rho) * ((d_grad * d_rho) ** 2),
        bulk_lag,
    )
    expect_dim("kept GNLS bulk kinetic stress scale m rho v_i v_j", d_m * d_rho * (d_v**2), stress)

    d_DtP = d_dt
    d_gradP = d_grad * d_P
    expect_dim("kept T0 P^i dimensionless", d_P, DIMENSIONLESS)
    expect_dim("kept T0 L_pol kinetic term", d_m * d_rho * (d_a**2) * (d_DtP**2), bulk_lag)
    expect_dim("kept T0 L_pol Frank term", d_m * d_rho * d_cs2 * (d_a**2) * (d_gradP**2), bulk_lag)
    expect_dim("kept T0 L_pol radial term", d_m * d_rho * d_cs2, bulk_lag)
    expect_zero(
        "kept T0 L_pol all terms homogeneous",
        homogeneity_residual(
            {
                "P_kinetic": d_m * d_rho * (d_a**2) * (d_DtP**2),
                "P_Frank": d_m * d_rho * d_cs2 * (d_a**2) * (d_gradP**2),
                "P_radial": d_m * d_rho * d_cs2,
            }
        ),
    )
    expect_dim("kept T0 couple-stress inertia m rho a^2", d_m * d_rho * (d_a**2), MASS * (LENGTH**-2))
    expect_dim(
        "kept T0 couple-stress stiffness m rho c_s^2 a^2",
        d_m * d_rho * d_cs2 * (d_a**2),
        MASS * (TIME**-2),
    )
    expect_dim("kept T0 bulk radial scale m rho c_s^2", d_m * d_rho * d_cs2, bulk_lag)

    expect_dim("profile ratio w/ell_g is dimensionless", d_w / d_ell_g, DIMENSIONLESS)
    expect_dim("profile g_ell(w)=exp(-(w/ell_g)^2)/(sqrt(pi) ell_g)", d_g, LENGTH**-1)
    expect_dim("profile measure dw g_ell(w) is dimensionless", LENGTH * d_g, DIMENSIONLESS)
    w = sp.Symbol("w", real=True)
    ell_g = sp.Symbol("ell_g", positive=True)
    g_expr = sp.exp(-(w / ell_g) ** 2) / (sp.sqrt(sp.pi) * ell_g)
    normalization = sp.integrate(g_expr, (w, -sp.oo, sp.oo))
    expect_zero("derived Gaussian normalization integral int g_ell(w) dw = 1", sp.simplify(normalization - 1))

    d_dtu = d_dt * d_u
    d_curlu = d_grad * d_u
    d_uw = d_u
    d_dtuw = d_dt * d_uw
    d_varpi = d_P
    expect_dim("target [u^a]=L", d_u, LENGTH)
    expect_dim("target [u_w]=L", d_uw, LENGTH)
    expect_dim("target [rho_br]=M L^-3", d_rho_br, MASS * (LENGTH**-3))
    expect_dim("target [mu_R]=M L^-1 T^-2", d_mu_R, brane_lag)
    expect_dim("target [lambda_Pu]=M L^-1 T^-2", d_lambda_Pu, brane_lag)
    expect_dim("target [Omega_w]=T^-1", d_Omega_w, TIME**-1)
    expect_dim("target [g_ell]=L^-1", d_g, LENGTH**-1)
    expect_dim("target [ell_g]=L", d_ell_g, LENGTH)

    mac_kin = d_rho_br * (d_dtu**2)
    mac_curl = d_mu_R * (d_curlu**2)
    expect_dim("L_Mac kinetic rho_br (partial_t u^a)^2", mac_kin, brane_lag)
    expect_dim("L_Mac MacCullagh curl mu_R Omega_u^a Omega_u^a", mac_curl, brane_lag)
    expect_zero("L_Mac homogeneous", homogeneity_residual({"kinetic": mac_kin, "curl": mac_curl}))
    pu_term = d_lambda_Pu * d_varpi * d_curlu
    expect_dim("L_Pu parity-repaired lambda_Pu varpi_a Omega_u^a", pu_term, brane_lag)
    uw_kin = d_rho_br * (d_dtuw**2)
    uw_gap = d_rho_br * (d_Omega_w**2) * (d_uw**2)
    expect_dim("L_uw kinetic rho_br (partial_t u_w)^2", uw_kin, brane_lag)
    expect_dim("L_uw bare gap rho_br Omega_w^2 u_w^2", uw_gap, brane_lag)
    expect_zero("L_uw homogeneous", homogeneity_residual({"kinetic": uw_kin, "gap": uw_gap}))

    for name, dim in {
        "g_ell L_Mac kinetic": mac_kin,
        "g_ell L_Mac curl": mac_curl,
        "g_ell L_Pu": pu_term,
        "g_ell L_uw kinetic": uw_kin,
        "g_ell L_uw gap": uw_gap,
    }.items():
        expect_dim(f"brane bulk representation {name}", d_g * dim, bulk_lag)

    expect_dim("action measure int dt d^4X bulk_lag", d_dt_measure * d_d4x * bulk_lag, action_dim)
    expect_dim("action measure int dt d^4X g_ell(w) L_brane", d_dt_measure * d_d4x * d_g * brane_lag, action_dim)

    t_wa = d_m * d_rho * d_v * d_v
    slope = d_grad * d_uw
    stress_slope = stress * slope
    expect_dim("projected traction T_wa=m rho v_w v_a", t_wa, stress)
    expect_dim("projected traction partial_b u_w is dimensionless", slope, DIMENSIONLESS)
    expect_dim("projected traction slope mixing", stress_slope, stress)
    expect_zero("full projected traction T_na homogeneous", homogeneity_residual({"T_wa": t_wa, "slope": stress_slope}))

    expect_fail("dim ablation drop_m_from_T_wa", dim_residual(d_rho * d_v * d_v, stress))
    expect_fail("dim ablation MacCullagh_without_curl", dim_residual(d_mu_R * (d_u**2), brane_lag))

    op_rho = d_rho_br * (d_omega**2)
    op_mu = d_mu_R * (d_k**2)
    expect_dim("linearization O_u rho_br omega^2 term", op_rho, eom_u_op)
    expect_dim("linearization O_u mu_R k^2 term", op_mu, eom_u_op)
    expect_zero("linearization O_u homogeneous", homogeneity_residual({"rho_br omega^2": op_rho, "mu_R k^2": op_mu}))
    c_gamma_sq = d_mu_R / d_rho_br
    expect_dim("target [c_gamma^2=mu_R/rho_br]=L^2 T^-2", c_gamma_sq, (LENGTH**2) * (TIME**-2))
    expect_dim("linearization omega_T^2=c_gamma^2 k^2", c_gamma_sq * (d_k**2), TIME**-2)
    expect_dim("linearization omega_P^2=c_s^2 k^2", d_cs2 * (d_k**2), TIME**-2)
    expect_dim("linearization P radial gap 2 c_s^2/a^2", d_cs2 / (d_a**2), TIME**-2)
    expect_zero(
        "linearization omega_radial^2 homogeneous",
        homogeneity_residual({"c_s^2 k^2": d_cs2 * (d_k**2), "2 c_s^2/a^2": d_cs2 / (d_a**2)}),
    )
    expect_dim("linearization omega_uw,bare^2=Omega_w^2", d_Omega_w**2, TIME**-2)
    expect_dim("linearization Fourier P-u monomial lambda_Pu P k u", d_lambda_Pu * d_P * d_k * d_u, brane_lag)

    print("  parity record: direct P_parallel.Omega_u is parity-ODD and excluded.")
    print("  parity record: w_hat.(P_parallel x Omega_u) is parity-EVEN using the conceded normal and is frozen.")

    return {
        "P": d_P,
        "u": d_u,
        "u_w": d_uw,
        "rho_br": d_rho_br,
        "mu_R": d_mu_R,
        "lambda_Pu": d_lambda_Pu,
        "Omega_w": d_Omega_w,
        "g_ell": d_g,
        "ell_g": d_ell_g,
        "c_gamma_sq": c_gamma_sq,
        "mu_R_4D": MASS * (LENGTH**-2) * (TIME**-2),
    }


def compute_flat_brane_dof(
    *,
    u_active: tuple[int, int, int] = (1, 1, 1),
    p_kinetic_present: bool = True,
    p_frank_present: bool = True,
    p_radial_present: bool = True,
    u_w_kinetic_present: bool = True,
    u_w_gap_present: bool = True,
    phi_present: bool = False,
) -> dict[str, Any]:
    """Count G0 flat-brane field-space DOF from frozen quadratic terms only."""

    eye_u = sp.eye(3)
    active_u = sp.diag(*u_active)
    curl_projector = sp.diag(0, 1, 1)
    u_kinetic_form = active_u * eye_u * active_u
    u_curl_form = active_u * curl_projector * active_u
    u_kinetic_rank = matrix_rank(u_kinetic_form)
    u_curl_rank = matrix_rank(u_curl_form)
    u_curl_nullity = u_kinetic_rank - u_curl_rank
    if u_curl_nullity < 0:
        raise AuditFailure("u curl rank exceeds active u kinetic rank")

    eye_p = sp.eye(4)
    zero_p = sp.zeros(4, 4)
    tangent_p = sp.diag(1, 1, 1, 0)
    radial_p = sp.diag(0, 0, 0, 1)
    p_kinetic_form = eye_p if p_kinetic_present else zero_p
    p_frank_form = eye_p if p_frank_present else zero_p
    p_radial_hessian = radial_p if p_radial_present else zero_p
    p_tangent_kinetic_rank = matrix_rank(tangent_p * p_kinetic_form * tangent_p)
    p_tangent_frank_rank = matrix_rank(tangent_p * p_frank_form * tangent_p)
    p_radial_kinetic_rank = matrix_rank(radial_p * p_kinetic_form * radial_p)
    p_radial_hessian_rank = matrix_rank(radial_p * p_radial_hessian * radial_p)

    u_w_kinetic_rank = matrix_rank(sp.Matrix([[1 if u_w_kinetic_present else 0]]))
    u_w_gap_rank = matrix_rank(sp.Matrix([[1 if u_w_gap_present else 0]]))
    phi_kinetic_rank = matrix_rank(sp.Matrix([[1 if phi_present else 0]]))

    computed_counts = {
        "u_transverse": u_curl_rank,
        "u_longitudinal": u_curl_nullity,
        "P_spin_wave": min(p_tangent_kinetic_rank, p_tangent_frank_rank),
        "P_soft_spin_radial": min(p_radial_kinetic_rank, p_radial_hessian_rank),
        "u_w": min(u_w_kinetic_rank, u_w_gap_rank),
        "phi": phi_kinetic_rank,
    }
    computed_total = sum(computed_counts.values())
    return {
        "counts": computed_counts,
        "total": computed_total,
        "rank_bookkeeping": {
            "u_kinetic_rank": u_kinetic_rank,
            "u_curl_rank": u_curl_rank,
            "u_curl_nullity_within_active_kinetic_space": u_curl_nullity,
            "P_tangent_kinetic_rank": p_tangent_kinetic_rank,
            "P_tangent_Frank_rank": p_tangent_frank_rank,
            "P_radial_kinetic_rank": p_radial_kinetic_rank,
            "P_radial_soft_spin_hessian_rank": p_radial_hessian_rank,
            "u_w_kinetic_rank": u_w_kinetic_rank,
            "u_w_gap_rank": u_w_gap_rank,
            "phi_kinetic_rank": phi_kinetic_rank,
        },
    }


def run_flat_brane_dof() -> dict[str, Any]:
    subbanner("Flat-brane linear DOF rank computation")
    k = sp.Symbol("k", nonzero=True)
    curl_stiffness = sp.Matrix([[0, 0, 0], [0, k**2, 0], [0, 0, k**2]])
    curl_rank = matrix_rank(curl_stiffness)
    curl_nullity = 3 - curl_rank
    expect_zero("curl stiffness rank transverse to k_a is 2", curl_rank - 2)
    expect_zero("curl stiffness nullity is one k-parallel u component", curl_nullity - 1)

    computed = compute_flat_brane_dof()
    reported_counts = dict(computed["counts"])
    expect_zero("reported u^a curl-transverse count is rank-computed", reported_counts["u_transverse"] - curl_rank)
    expect_zero("reported u^a kinetic-minus-curl count is rank-computed", reported_counts["u_longitudinal"] - curl_nullity)
    expect_zero("flat-brane total DOF rank-computes to frozen value 8", computed["total"] - 8)

    expected_breakdown = {
        "u_transverse": 2,
        "u_longitudinal": 1,
        "P_spin_wave": 3,
        "P_soft_spin_radial": 1,
        "u_w": 1,
        "phi": 0,
    }
    for name, expected in expected_breakdown.items():
        expect_zero(f"G0.5 breakdown {name} reported->computed", reported_counts[name] - expected)
    print("  structural-postulate-6 fact: C5 phi is absent in the active baseline, so phi contributes 0 DOF.")
    print(f"  rank bookkeeping: {computed['rank_bookkeeping']}")
    print(f"  computed G0.5 counts: {reported_counts}; total={computed['total']}")

    mutations = [
        ("drop_u_w_gap_term", compute_flat_brane_dof(u_w_gap_present=False)),
        ("drop_P_soft_spin_radial_term", compute_flat_brane_dof(p_radial_present=False)),
        ("zero_u_longitudinal_component", compute_flat_brane_dof(u_active=(0, 1, 1))),
    ]
    for name, mutated in mutations:
        expect_zero(f"DOF ablation {name} computes total 7", mutated["total"] - 7)
        expect_fail(f"DOF ablation {name} changes computed total away from 8", mutated["total"] - computed["total"])
    return computed


@dataclass(frozen=True)
class DriftEntry:
    name: str
    category: str
    dim: Dim | None
    note: str


STRUCTURAL_POSTULATES = (
    "imposed `ŵ` axis + `w=0` surface (conceded-wall)",
    "`uᵃ` same-medium surface collective, tangentially free-slip (`u̇ᵃ ≠ vᵃ`)",
    "T0 `Pⁱ` reused as the Cosserat micro-rotation reservoir (0 new DOF)",
    "baseline `Pⁱ` spin-wave status = `massless` (alternates `gapped`/`slaved-rigid` named-inactive)",
    "the `ŵ`-dependent parity-EVEN P–u operator re-admits the ε-contracted/chiral class excluded by T0 and REQUIRES the conceded axis `ŵ` (a structural-postulate cost, not a free choice; the direct `P_∥·Ω_u` is parity-ODD, excluded)",
    "no C5 `φ` analog / no longitudinal constraint",
)

FROZEN_CONSTANT_NAMES = ("rho_br", "mu_R", "lambda_Pu", "Omega_w")
FROZEN_FUNCTION_NAMES = ("g_ell(w)",)
ANTI_ABSORPTION_NAMES = frozenset({"rho_B0", "chi_c", "C_hu"})
VALID_CATEGORIES = frozenset({"constant", "function", "structural_postulate"})


def drift_expected_n() -> int:
    match = re.fullmatch(r"SECOND_MEDIUM_DRIFT_AT_FREEZE\((\d+)\)", EXPECTED_DRIFT_TOKEN)
    if match is None:
        raise AuditFailure(f"malformed expected drift token: {EXPECTED_DRIFT_TOKEN}")
    return int(match.group(1))


def build_drift_table() -> tuple[DriftEntry, ...]:
    return (
        DriftEntry("rho_br", "constant", MASS * (LENGTH**-3), "surface inertia"),
        DriftEntry("mu_R", "constant", MASS * (LENGTH**-1) * (TIME**-2), "MacCullagh modulus"),
        DriftEntry("lambda_Pu", "constant", MASS * (LENGTH**-1) * (TIME**-2), "parity-repaired P-u coupling"),
        DriftEntry("Omega_w", "constant", TIME**-1, "bare u_w gap scale"),
        DriftEntry("g_ell(w)", "function", LENGTH**-1, "fixed Gaussian shape, ONE width knob ell_g; no free-form profile"),
        *(DriftEntry(postulate, "structural_postulate", None, "verbatim structural postulate") for postulate in STRUCTURAL_POSTULATES),
    )


def category_counts(table: Iterable[DriftEntry]) -> Counter[str]:
    return Counter(entry.category for entry in table)


def set_residual(actual: set[str], expected: set[str]) -> sp.Integer:
    return sp.Integer(0) if actual == expected else sp.Integer(1)


def drift_table_residual(table: tuple[DriftEntry, ...], *, verdict_n_delta: int = 0) -> sp.Expr:
    counts = category_counts(table)
    names = {entry.name for entry in table}
    residual = sp.Integer(0)
    residual += sp.Integer(sum(1 for entry in table if entry.category not in VALID_CATEGORIES))
    residual += (counts["constant"] - len(FROZEN_CONSTANT_NAMES)) ** 2
    residual += (counts["function"] - len(FROZEN_FUNCTION_NAMES)) ** 2
    residual += (counts["structural_postulate"] - len(STRUCTURAL_POSTULATES)) ** 2
    residual += set_residual({entry.name for entry in table if entry.category == "constant"}, set(FROZEN_CONSTANT_NAMES))
    residual += set_residual({entry.name for entry in table if entry.category == "function"}, set(FROZEN_FUNCTION_NAMES))
    residual += sp.Integer(0 if names.isdisjoint(ANTI_ABSORPTION_NAMES) else 1)
    n = sum(counts[category] for category in VALID_CATEGORIES)
    residual += (n - drift_expected_n()) ** 2
    verdict = f"SECOND_MEDIUM_DRIFT_AT_FREEZE({n + verdict_n_delta})"
    residual += sp.Integer(0 if verdict == EXPECTED_DRIFT_TOKEN else 1)
    return sp.simplify(residual)


def drift_table_dim_residual(table: tuple[DriftEntry, ...], dims: dict[str, Dim]) -> sp.Expr:
    residual = sp.Integer(0)
    for entry in table:
        if entry.category not in {"constant", "function"}:
            continue
        target_key = DRIFT_DIM_TARGET_KEYS.get(entry.name)
        if target_key is None or entry.dim is None or target_key not in dims:
            residual += sp.Integer(1)
            continue
        residual += dim_residual(entry.dim, dims[target_key])
    return sp.simplify(residual)


def run_computed_drift_ledger(dims: dict[str, Dim]) -> tuple[DriftEntry, ...]:
    subbanner("Computed drift ledger enumeration")
    table = build_drift_table()
    counts = category_counts(table)
    constant_subcount = counts["constant"]
    function_subcount = counts["function"]
    structural_subcount = counts["structural_postulate"]
    n = sum(counts[category] for category in VALID_CATEGORIES)
    verdict = f"SECOND_MEDIUM_DRIFT_AT_FREEZE({n})"

    print("  Enumerated drift members:")
    for entry in table:
        dim_text = "" if entry.dim is None else f", [{entry.dim}]"
        print(f"    - {entry.name}: {entry.category}{dim_text}; {entry.note}")
    expect_zero("constant subcount computed from enumeration", constant_subcount - len(FROZEN_CONSTANT_NAMES))
    expect_zero("function subcount computed from enumeration", function_subcount - len(FROZEN_FUNCTION_NAMES))
    expect_zero("structural-postulate subcount computed from enumeration", structural_subcount - len(STRUCTURAL_POSTULATES))
    expect_zero("independent new input count n computed from subcounts", n - drift_expected_n())
    expect_zero("verdict string built from computed n equals frozen token", 0 if verdict == EXPECTED_DRIFT_TOKEN else 1)
    expect_zero("drift table anti-absorption guard excludes Part-VI trio", drift_table_residual(table))
    expect_zero("drift table Dim fields match dimensional-firewall targets", drift_table_dim_residual(table, dims))
    print(f"  computed n = {constant_subcount}+{function_subcount}+{structural_subcount} = {n}")
    print(f"  computed verdict = {verdict}")

    new_fields = (
        ("u^a", ("u_x", "u_y", "u_z")),
        ("u_w", ("u_w",)),
    )
    new_field_subcount = sum(len(components) for _, components in new_fields)
    expect_zero("new-field subcount computed separately from field list", new_field_subcount - 4)
    expect_bool("new-field names are kept out of the 11-input drift table", {"u^a", "u_w"}.isdisjoint({entry.name for entry in table}))
    print(f"  new-field subcount = {new_field_subcount} from u^a (3) + u_w (1), separate from the 11-input drift count.")

    kept_t0_couple_stress = (
        ("m rho a^2", False),
        ("m rho c_s^2 a^2", False),
        ("m rho c_s^2", False),
    )
    t0_new_count = sum(1 for _, is_new in kept_t0_couple_stress if is_new)
    expect_zero("T0 couple-stress coefficients contribute 0 new entries", t0_new_count)
    print("  T0 couple-stress coefficients are kept-not-new: m rho a^2; m rho c_s^2 a^2; m rho c_s^2.")

    dropped = table[:-1]
    expect_fail("enumeration tooth: drop one entry gives n=10 and trips drift validation", drift_table_residual(dropped))
    miscategorized = tuple(
        DriftEntry(entry.name, "structural_postulate", entry.dim, entry.note) if entry.name == "Omega_w" else entry
        for entry in table
    )
    expect_fail("enumeration tooth: miscategorize Omega_w trips subcount assertions", drift_table_residual(miscategorized))
    injected = table + (DriftEntry("rho_B0", "constant", MASS * (LENGTH**-3), "forbidden Part-VI injection"),)
    expect_fail("enumeration tooth: inject rho_B0 trips anti-absorption guard", drift_table_residual(injected))
    expect_fail("enumeration tooth: corrupt computed n before verdict assembly trips token equality", drift_table_residual(table, verdict_n_delta=1))
    dim_corrupted = tuple(
        DriftEntry(entry.name, entry.category, entry.dim * LENGTH, entry.note)
        if entry.name == "rho_br" and entry.dim is not None
        else entry
        for entry in table
    )
    expect_fail("enumeration tooth: corrupt rho_br table Dim trips firewall consistency", drift_table_dim_residual(dim_corrupted, dims))
    expect_zero("baseline drift table remains valid after copy-mutation teeth", drift_table_residual(table))
    return table


def structural_distinctness_failure_residual(mu_r: Dim, mu_r_4d: Dim) -> sp.Integer:
    return sp.Integer(1) if dim_residual(mu_r, mu_r_4d) == 0 else sp.Integer(0)


def run_mu_r_firewall(dims: dict[str, Dim]) -> None:
    subbanner("mu_R notational firewall and R17 pending edge")
    mu_r = dims["mu_R"]
    mu_r_4d = dims["mu_R_4D"]
    print(f"  [mu_R] 3D brane modulus = {mu_r}")
    print(f"  [mu_R^(4)] 4D shear-stiffness density = {mu_r_4d}")
    expect_nonzero("[mu_R] != [mu_R^(4)] as exponent triples", dim_residual(mu_r, mu_r_4d))
    expect_dim("R17 dim consistency [mu_R^(4)]*L = [mu_R]", mu_r_4d * LENGTH, mu_r)
    print("  R17 status: PENDING (mu_R = int chi_B mu_R^(4) dw; deferred nonlinear throat/projection).")
    expect_fail("mu_R firewall tooth: forced equality trips distinctness inequality", structural_distinctness_failure_residual(mu_r, mu_r))


def print_structural_postulates() -> None:
    subbanner("Six structural postulates printed verbatim")
    for index, postulate in enumerate(STRUCTURAL_POSTULATES, start=1):
        print(f"  {index}. {postulate}")


def print_provenance_blocks() -> None:
    subbanner("Erratum, supersession, methodology, Gate-L scope, and carried debt")
    print("  ERRATUM (2026-07-04): the SECOND_MEDIUM_DRIFT_AT_FREEZE(11) count is NOT inflated by a rho_br overcount; this freeze's count STANDS.")
    print("  pathA_25 varrho_br[rho] belongs to the CLOSED density-smectic candidate, FAIL_NOT_CODIM1, OUT_OF_ACTIVE_NG5.")
    print("  This rho_br/mu_R is genuine postulated shear-surface inertia/modulus with registered-pending pathA_40 Route-A reduction.")
    print("  Corroboration token: NO_OVERCOUNT_ROUTE_A_PENDING.")
    print("  Honest cross-sector drift {rho_B0, chi_c, C_hu} is a Part-VI item, not absorbed into the 11.")
    print("")
    print("  Supersession fact 1: stage006 chi_B order-field wall superseded fixed-shape g_ell(w) as the MATERIAL-STATE closure.")
    print("  Supersession fact 2: G0 REMAINS the light-sector CONSTITUTIVE freeze; stage003 consumes L_Mac as-is.")
    print("")
    print("  Methodology: postulating an ingredient is allowed; postulating an outcome is not.")
    print("  Methodology: late ingredient = AD_HOC_RESCUE -> fresh G0; every knob counted; >=2 new inputs => drift reported plainly.")
    print("  Methodology: g(w) admitted on locality/minimality grounds ONLY, target-blind.")
    print("  Methodology anti-impose: G0 freezes TERMS, not gate answers; no bounded-below, traction, or longitudinal-is-gauge claims.")
    print("  Methodology: a clean all-pass is suspicious, so able-to-fail teeth are live.")
    print("")
    print("  Scope: G0 freeze only; no Gate L verdict computed.")
    print("  classification_guard: counts only; no gate verdict, no boundedness, no gauge, no leak claim.")
    print("  Gate-L exposure names are prose provenance only: FAIL_HIDDEN_PROPAGATING_MODE; FAIL_GYROSTAT_NO_CLOSURE; FAIL_NOT_BOUNDED_BELOW; linked FAIL_COUPLE_STRESS_NOGO chain remains able to fire.")
    print("  Gate-L artifacts are not imported; exposure strings are not used as computed predicates.")
    print("")
    print("  Route-A carried debt: {rho_br, mu_R} reduction = Route-A PENDING (R10), ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT, free-unreduced brane constants on the deferred nonlinear throat.")
    print("  Labeled postulated constants registered as knobs: {lambda_Pu, Omega_w, ell_g}.")
    print("")
    print("  Downstream consumers: ledger_stage003 consumes mu_R, rho_br, c_gamma^2=mu_R/rho_br, and frozen L_Mac.")
    print("  Downstream consumers: ledger_stage006 cites c_gamma^2=mu_R/rho_br plus the chi_B/G0 supersession relationship.")


def print_verdict_labels() -> None:
    print("")
    print("Verdict labels:")
    print("  ledger earned-label (NOT a source verdict token): G0_FREEZE_FIDELITY_DOF_DIM_VERIFIED  (freeze-action SHA-256 re-verified incl. T0 embedding; flat-brane DOF=8 rank-computed; dimensional firewall over every frozen term; drift count COMPUTED from the enumerated table)")
    print(f"  source top-line verdicts: {EXPECTED_FREEZE_TOKEN} + {EXPECTED_DRIFT_TOKEN}")
    print("  honest landing: the 11 = POSTULATED/CALIBRATED freeze inputs (4 constants {rho_br, mu_R, lambda_Pu, Omega_w} + 1 function g_l(w; l_g) + 6 structural postulates), enumerated + computed — the freeze freezes TERMS, not gate answers")
    print("  erratum (2026-07-04): the 11 STANDS (no rho_br overcount; pathA_25 varrho_br = CLOSED candidate, OUT_OF_ACTIVE_NG5); honest cross-sector drift {rho_B0, chi_c, C_hu} = Part-VI (pathA_41), NOT absorbed here [NO_OVERCOUNT_ROUTE_A_PENDING]")
    print("  supersession: stage006 chi_B wall = the MATERIAL-STATE closure (supersedes fixed-shape g_l(w) as material wall); G0 REMAINS the light-sector CONSTITUTIVE freeze (stage003 consumes L_Mac as-is)")
    print("  notational firewall: mu_R (3D brane, M L^-1 T^-2) != mu_R_4D (4D density, M L^-2 T^-2); related only by PENDING R17 projection")
    print("  Gate-L: EXCLUDED — no gate verdict computed or imported; exposure names printed as provenance only")
    print("  carried: Route-A reduction PENDING (pathA_40, R10) for {rho_br, mu_R}; new postulated constants {lambda_Pu, Omega_w, l_g} registered as knobs")


def main() -> None:
    banner("ledger_stage007_shear_surface_g0_freeze SymPy audit")
    run_freeze_fidelity()
    dims = run_dimensional_firewall()
    run_flat_brane_dof()
    run_computed_drift_ledger(dims)
    run_mu_r_firewall(dims)
    print_structural_postulates()
    print_provenance_blocks()
    print_verdict_labels()
    print("")
    print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
    print("OVERALL PASS: SymPy verified ledger_stage007 G0 freeze fidelity, dimensions, DOF, and computed drift")


if __name__ == "__main__":
    try:
        main()
    except AuditFailure as exc:
        print("")
        print(f"PASS tally: {PASS_COUNT}; FAIL tally: {FAIL_COUNT}")
        print(f"OVERALL FAIL: SymPy stage007 audit did not close ({exc})")
        raise SystemExit(1)
