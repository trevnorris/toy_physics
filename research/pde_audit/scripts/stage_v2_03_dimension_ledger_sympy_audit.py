#!/usr/bin/python3
"""
Stage V2-03: Dimension ledger and normalization audit.

This script implements a small exact dimension-vector engine and uses SymPy for
series identities in the grouped-P2 / outgoing quadrupole bridge.

Base dimensions are
    M  = mass
    L  = length
    T  = time
    Q  = electric charge
    O  = abstract reduced wall/operator unit

The O axis is intentionally included because the reduced wall operator D0 is a
mechanical/kernel unit whose absolute physical dimension depends on the chosen
wall coordinate normalization. Ratios such as N0/D0 should cancel O.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Dict, List, Tuple
import json
import sys

# Make the project venv's SymPy visible even when this script is run with
# /usr/bin/python3, which avoids the slow site-startup issue in this container.
sys.path.append('/opt/pyvenv/lib/python3.13/site-packages')
try:
    import sympy as sp
except Exception:  # pragma: no cover - fallback for external users without sympy
    sp = None

AXES = ("M", "L", "T", "Q", "O")


def F(x) -> Fraction:
    return Fraction(x)


@dataclass(frozen=True)
class Dim:
    powers: Tuple[Fraction, ...]

    def __post_init__(self):
        if len(self.powers) != len(AXES):
            raise ValueError(f"Dim must have {len(AXES)} axes")

    @staticmethod
    def of(vals: Iterable[object]) -> "Dim":
        return Dim(tuple(F(v) for v in vals))

    def __mul__(self, other: "Dim") -> "Dim":
        return Dim(tuple(a + b for a, b in zip(self.powers, other.powers)))

    def __truediv__(self, other: "Dim") -> "Dim":
        return Dim(tuple(a - b for a, b in zip(self.powers, other.powers)))

    def __pow__(self, n: object) -> "Dim":
        n = F(n)
        return Dim(tuple(n * a for a in self.powers))

    def sqrt(self) -> "Dim":
        return self ** Fraction(1, 2)

    def is_dimensionless(self) -> bool:
        return all(p == 0 for p in self.powers)

    def pretty(self) -> str:
        chunks = []
        for axis, p in zip(AXES, self.powers):
            if p == 0:
                continue
            if p == 1:
                chunks.append(axis)
            else:
                chunks.append(f"{axis}^{p}")
        return "1" if not chunks else " ".join(chunks)


ONE = Dim.of([0, 0, 0, 0, 0])
MASS = Dim.of([1, 0, 0, 0, 0])
LENGTH = Dim.of([0, 1, 0, 0, 0])
TIME = Dim.of([0, 0, 1, 0, 0])
CHARGE = Dim.of([0, 0, 0, 1, 0])
OPER = Dim.of([0, 0, 0, 0, 1])

# Physical constants / primitive units.
G = LENGTH**3 / (MASS * TIME**2)
c = LENGTH / TIME
c_s = c
radius_a = LENGTH
omega = ONE / TIME
hbar = MASS * LENGTH**2 / TIME
energy = MASS * LENGTH**2 / TIME**2

# EM zero-mode bookkeeping.
Z_int = LENGTH  # if Z(w) is dimensionless and dw has length.
mu0_eff_4d = MASS * LENGTH / (CHARGE**2)
mu0_5d = mu0_eff_4d * Z_int
q_eff = CHARGE
q_star = q_eff * Z_int.sqrt()

# Gravitational quadrupole/response target dimensions.
P0_target = G * c_s**5 / (radius_a**5 * c**5)
Lambda0 = P0_target
Gamma_GR = G / c**5
G5_out = radius_a**5 / c_s**5
A_out = radius_a**2 / c_s**2
B_out = radius_a**4 / c_s**4

# Abstract reduced operator moments.
D0 = OPER
D2 = OPER * TIME**2
D4 = OPER * TIME**4
P0 = P0_target
P2 = P0 * TIME**2
P4 = P0 * TIME**4
N0 = P0 * D0
N2 = N0 * TIME**2
N4 = N0 * TIME**4
u2 = TIME**2
u4 = TIME**4

# Raw mechanical Maxwell/mixed prototype dimensions with canonical internal
# coordinates. This reveals the raw-vs-normalized bridge requirement.
# For a canonical oscillator coordinate U, 1/2 Udot^2 has energy dimension.
U_can = (energy * TIME**2).sqrt()
q_wall = LENGTH  # representative wall displacement coordinate.
D_mech = energy / (q_wall**2)
g_can = energy / (q_wall * U_can)
Omega2 = ONE / TIME**2
R_mix = Omega2
Delta_mix = Omega2 * Omega2
P_mix_numer = Omega2 * g_can  # Omega_U^2 g_W + R g_U
N0_raw_mix = (P_mix_numer**2) / (Delta_mix**2)
P0_raw_mix = N0_raw_mix / D_mech
bridge_scale = P0_target / P0_raw_mix


def assert_same(name: str, lhs: Dim, rhs: Dim, records: List[Dict]) -> None:
    ok = lhs == rhs
    records.append({"name": name, "lhs": lhs.pretty(), "rhs": rhs.pretty(), "status": "PASS" if ok else "FAIL"})
    if not ok:
        raise AssertionError(f"{name}: {lhs.pretty()} != {rhs.pretty()}")


def assert_dimless(name: str, dim: Dim, records: List[Dict]) -> None:
    ok = dim.is_dimensionless()
    records.append({"name": name, "dim": dim.pretty(), "status": "PASS" if ok else "FAIL"})
    if not ok:
        raise AssertionError(f"{name}: expected dimensionless, got {dim.pretty()}")


def sympy_identity_checks() -> List[Dict]:
    if sp is None:
        return [{"name": "SymPy import", "status": "SKIP", "detail": "SymPy unavailable"}]

    w = sp.symbols('omega')
    D0s, D2s, D4s, N0s, N2s, N4s = sp.symbols('D0 D2 D4 N0 N2 N4', nonzero=True)
    D = D0s + D2s*w**2 + D4s*w**4

    # Normalized response Y = D0/D.
    Y = sp.series(D0s / D, w, 0, 6).removeO().expand()
    u2_expr = sp.simplify(Y.coeff(w, 2))
    u4_expr = sp.simplify(Y.coeff(w, 4))
    expected_u2 = -D2s/D0s
    expected_u4 = (D2s**2 - D0s*D4s)/D0s**2

    # Outgoing prefactor Pref = D0 N / D^2.
    N = N0s + N2s*w**2 + N4s*w**4
    Pref = sp.series(D0s * N / D**2, w, 0, 6).removeO().expand()
    P0_expr = sp.simplify(Pref.coeff(w, 0))
    P2_expr = sp.simplify(Pref.coeff(w, 2))
    P4_expr = sp.simplify(Pref.coeff(w, 4))
    expected_P0 = N0s/D0s
    expected_P2 = (D0s*N2s - 2*D2s*N0s)/D0s**2
    expected_P4 = (D0s**2*N4s - 2*D0s*(D2s*N2s + D4s*N0s) + 3*D2s**2*N0s)/D0s**3

    # Stage-4 one-lane transfer expansion.
    Delta, S2, P, gW = sp.symbols('Delta S2 P gW', nonzero=True)
    Nmix = ((P - gW*w**2)**2 / (Delta - S2*w**2 + w**4)**2)
    Nmix_series = sp.series(Nmix, w, 0, 6).removeO().expand()
    expected_Nmix0 = P**2/Delta**2
    expected_Nmix2 = 2*P*(P*S2 - Delta*gW)/Delta**3
    expected_Nmix4 = (Delta**2*gW**2 - 2*Delta*P**2 - 4*Delta*P*S2*gW + 3*P**2*S2**2)/Delta**4

    checks = [
        ("SymPy response u2 coefficient", u2_expr, expected_u2),
        ("SymPy response u4 coefficient", u4_expr, expected_u4),
        ("SymPy prefactor P0 coefficient", P0_expr, expected_P0),
        ("SymPy prefactor P2 coefficient", P2_expr, expected_P2),
        ("SymPy prefactor P4 coefficient", P4_expr, expected_P4),
        ("SymPy mixed transfer N0 coefficient", sp.simplify(Nmix_series.coeff(w,0)), expected_Nmix0),
        ("SymPy mixed transfer N2 coefficient", sp.simplify(Nmix_series.coeff(w,2)), expected_Nmix2),
        ("SymPy mixed transfer N4 coefficient", sp.simplify(Nmix_series.coeff(w,4)), expected_Nmix4),
    ]
    out = []
    for name, got, exp in checks:
        ok = sp.simplify(got - exp) == 0
        out.append({"name": name, "status": "PASS" if ok else "FAIL", "got": str(got), "expected": str(exp)})
        if not ok:
            raise AssertionError(f"{name}: {got} != {exp}")
    return out


def main() -> None:
    records: List[Dict] = []

    # 1. Core target packets.
    assert_same("P0_target and Lambda0 have same units", P0_target, Lambda0, records)
    assert_same("Gamma_GR equals P0_target * a^5 / c_s^5", Gamma_GR, P0_target * G5_out, records)
    assert_same("Kbar0_target same as P0_target", P0_target, G * c_s**5 / (radius_a**5 * c**5), records)

    # 2. Outgoing l=2 fingerprint dimensionlessness.
    assert_dimless("A_out * omega^2", A_out * omega**2, records)
    assert_dimless("B_out * omega^4", B_out * omega**4, records)
    assert_dimless("G5_out * omega^5", G5_out * omega**5, records)

    # 3. Operator / response conversion formulas.
    assert_same("u2 = -D2/D0", D2 / D0, u2, records)
    assert_same("u4 = (D2^2-D0D4)/D0^2", (D2**2) / (D0**2), u4, records)
    assert_same("D0*D4 term in u4", (D0 * D4) / (D0**2), u4, records)
    assert_same("one-pole identity u4 ~ u2^2", u4, u2**2, records)

    # 4. Prefactor moments P0, P2, P4 from N/D formulas.
    assert_same("P0 = N0/D0", N0 / D0, P0, records)
    assert_same("P2 formula term D0*N2/D0^2", (D0 * N2) / (D0**2), P2, records)
    assert_same("P2 formula term D2*N0/D0^2", (D2 * N0) / (D0**2), P2, records)
    assert_same("P4 formula term D0^2*N4/D0^3", (D0**2 * N4) / (D0**3), P4, records)
    assert_same("P4 formula term D0*D2*N2/D0^3", (D0 * D2 * N2) / (D0**3), P4, records)
    assert_same("P4 formula term D0*D4*N0/D0^3", (D0 * D4 * N0) / (D0**3), P4, records)
    assert_same("P4 formula term D2^2*N0/D0^3", (D2**2 * N0) / (D0**3), P4, records)

    # 5. Constant-prefactor branch conditions.
    assert_same("constant-prefactor N2 condition", N2, D2 * N0 / D0, records)
    assert_same("constant-prefactor N4 first term", N4, (D0 * D2 * N2) / (D0**2), records)
    assert_same("constant-prefactor N4 second term", N4, (D0 * D4 * N0) / (D0**2), records)
    assert_same("constant-prefactor N4 third term", N4, (D2**2 * N0) / (D0**2), records)

    # 6. Full-bundle isotropic one-pole surface.
    # D0*(B4+Z4)=3*(M+B2+Z2)^2 with B4+Z4 ~ D4 and M+B2+Z2 ~ D2.
    assert_same("isotropic one-pole surface dimensions", D0 * D4, D2**2, records)

    # 7. EM zero-mode ledger.
    assert_same("mu0_5d/Z_int gives mu0_eff_4d", mu0_5d / Z_int, mu0_eff_4d, records)
    assert_same("q_star/sqrt(Z_int) gives q_eff", q_star / Z_int.sqrt(), q_eff, records)

    # 8. Support-wave mass theorem ledger.
    omega_support = c_s / radius_a
    assert_same("support frequency c_s/a", omega_support, omega, records)
    assert_same("hbar*omega_support is energy", hbar * omega_support, energy, records)

    # 9. Raw mechanical mixed model: exact but not yet gravitationally normalized.
    assert_same("raw canonical mixed N0 has D_mech*T^2", N0_raw_mix, D_mech * TIME**2, records)
    assert_same("raw canonical mixed P0 has T^2", P0_raw_mix, TIME**2, records)
    assert_same("bridge scale maps raw P0 to gravitational P0", bridge_scale * P0_raw_mix, P0_target, records)

    sympy_records = sympy_identity_checks()

    summary = {
        "stage": "V2-03 dimension ledger and normalization audit",
        "sympy_version": None if sp is None else sp.__version__,
        "dimension_checks": len(records),
        "dimension_passes": sum(1 for r in records if r["status"] == "PASS"),
        "sympy_identity_checks": len(sympy_records),
        "sympy_identity_passes": sum(1 for r in sympy_records if r["status"] == "PASS"),
        "P0_target_dim": P0_target.pretty(),
        "Gamma_GR_dim": Gamma_GR.pretty(),
        "D0_dim": D0.pretty(),
        "N0_normalized_dim": N0.pretty(),
        "raw_mixed_P0_dim": P0_raw_mix.pretty(),
        "required_raw_to_grav_bridge_scale_dim": bridge_scale.pretty(),
        "verdict": (
            "PASS for published target identities and response-moment algebra; "
            "WARN that the raw canonical Maxwell/mixed transfer factor is a mechanical T^2 "
            "object and must be multiplied by an explicit port/source normalization scale "
            "to become the gravitational P0 used in the 2.5PN/4PN target."
        ),
    }

    print("=== Stage V2-03 Dimension Ledger ===")
    for r in records:
        if "lhs" in r:
            print(f"[{r['status']}] {r['name']}: {r['lhs']} == {r['rhs']}")
        else:
            print(f"[{r['status']}] {r['name']}: {r['dim']}")
    print("\n=== SymPy identity checks ===")
    for r in sympy_records:
        print(f"[{r['status']}] {r['name']}")
    print("\n=== Summary ===")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
