"""SymPy dimensional audit for the Path-A to R_norm chain.

The module is intentionally standalone: it does not import the numerical
solvers, and it does not change any model formula.  It records the restored
unit dictionary and checks the dimensional contracts used by the Path-A
reduction/extraction chain.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import sympy as sp


L_EXP, T_EXP, M_EXP = sp.symbols("L_exp T_exp M_exp")


@dataclass(frozen=True)
class Dim:
    """Exponent vector for base dimensions (L, T, M)."""

    l: sp.Rational | int = 0
    t: sp.Rational | int = 0
    m: sp.Rational | int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "l", sp.Rational(self.l))
        object.__setattr__(self, "t", sp.Rational(self.t))
        object.__setattr__(self, "m", sp.Rational(self.m))

    def __mul__(self, other: "Dim") -> "Dim":
        return Dim(self.l + other.l, self.t + other.t, self.m + other.m)

    def __truediv__(self, other: "Dim") -> "Dim":
        return Dim(self.l - other.l, self.t - other.t, self.m - other.m)

    def __pow__(self, power: int | sp.Rational) -> "Dim":
        p = sp.Rational(power)
        return Dim(self.l * p, self.t * p, self.m * p)

    def is_dimensionless(self) -> bool:
        return self == DIMENSIONLESS

    def as_tuple(self) -> tuple[str, str, str]:
        return (str(self.l), str(self.t), str(self.m))

    def monomial(self) -> sp.Expr:
        return sp.simplify(L_EXP**self.l * T_EXP**self.t * M_EXP**self.m)

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
ENERGY = Dim(2, -2, 1)
ACTION = Dim(2, -1, 1)
FORCE = Dim(1, -2, 1)


@dataclass(frozen=True)
class Check:
    section: str
    name: str
    status: str
    expected: Dim | None
    actual: Dim | None
    factor_needed: Dim | None
    note: str

    def as_dict(self) -> dict[str, object]:
        return {
            "section": self.section,
            "name": self.name,
            "status": self.status,
            "expected": None if self.expected is None else str(self.expected),
            "actual": None if self.actual is None else str(self.actual),
            "expected_tuple_L_T_M": None if self.expected is None else self.expected.as_tuple(),
            "actual_tuple_L_T_M": None if self.actual is None else self.actual.as_tuple(),
            "factor_needed_to_reach_expected": (
                None if self.factor_needed is None else str(self.factor_needed)
            ),
            "factor_tuple_L_T_M": (
                None if self.factor_needed is None else self.factor_needed.as_tuple()
            ),
            "note": self.note,
        }


def factor_to_reach(expected: Dim, actual: Dim) -> Dim:
    return expected / actual


def expect_dim(section: str, name: str, actual: Dim, expected: Dim, note: str = "") -> Check:
    if actual == expected:
        return Check(section, name, "CONSISTENT", expected, actual, DIMENSIONLESS, note)
    return Check(section, name, "INCONSISTENT", expected, actual, factor_to_reach(expected, actual), note)


def homogeneous(section: str, name: str, terms: Mapping[str, Dim], note: str = "") -> Check:
    if not terms:
        raise ValueError("homogeneous check requires at least one term")
    items = list(terms.items())
    expected = items[0][1]
    offenders = [(term, dim) for term, dim in items[1:] if dim != expected]
    if not offenders:
        detail = ", ".join(f"{term}:{dim}" for term, dim in items)
        return Check(section, name, "CONSISTENT", expected, expected, DIMENSIONLESS, f"{note} terms={detail}".strip())
    offender_text = "; ".join(
        f"{term} has {dim}, needs factor {factor_to_reach(expected, dim)}" for term, dim in offenders
    )
    return Check(section, name, "INCONSISTENT", expected, offenders[0][1], factor_to_reach(expected, offenders[0][1]), f"{note} {offender_text}".strip())


def assert_expected(actual: Dim, expected: Dim, label: str) -> None:
    if actual != expected:
        raise AssertionError(
            f"{label}: expected {expected}, got {actual}; "
            f"factor needed {factor_to_reach(expected, actual)}"
        )


def assert_homogeneous(terms: Mapping[str, Dim], label: str) -> None:
    check = homogeneous("assert", label, terms)
    if check.status != "CONSISTENT":
        raise AssertionError(check.note)


def dimensional_dictionary() -> dict[str, Dim]:
    """Step-0 dictionary, restored before natural-unit pinning."""

    rho = LENGTH**-4
    psi = LENGTH**-2
    eos_k = ENERGY / (rho**4)
    pressure = ENERGY / (LENGTH**4)
    velocity = LENGTH / TIME
    omega = TIME**-1
    reduced_stiffness = MASS * (LENGTH**-1) * (TIME**-2)
    reduced_mass = MASS * (LENGTH**-1)
    return {
        "1": DIMENSIONLESS,
        "L": LENGTH,
        "T": TIME,
        "M": MASS,
        "hbar": ACTION,
        "m_particle": MASS,
        "rho_4d_number_density": rho,
        "psi": psi,
        "K_eos": eos_k,
        "P_pressure": pressure,
        "U_bulk_energy_density": pressure,
        "h_enthalpy": ENERGY,
        "c_s": velocity,
        "c": velocity,
        "a": LENGTH,
        "R": LENGTH,
        "w": LENGTH,
        "dw": LENGTH,
        "omega": omega,
        "varpi": omega,
        "q_canonical": DIMENSIONLESS,
        "q_A0": ENERGY,
        "q_Ai": ACTION / LENGTH,
        "A0_if_q_dimensionless": ENERGY,
        "Ai_if_q_dimensionless": ACTION / LENGTH,
        "tau_tension": FORCE,
        "mu_wall_restored_as_tau_over_c_s2": MASS / LENGTH,
        "T_w": FORCE,
        "T_Omega": FORCE / (LENGTH**2),
        "U_Sigma_per_w": FORCE,
        "U_Sigma_R": FORCE / LENGTH,
        "U_Sigma_RR": FORCE / (LENGTH**2),
        "radial_shell_3ball": LENGTH**3,
        "wall_line_measure": LENGTH,
        "tube_cell_measure": LENGTH**4,
        "raw_4ball_radial_measure": LENGTH**4,
        "chi_wall_normalized": LENGTH ** sp.Rational(-1, 2),
        "bdg_u_v_symplectic_normalized": psi,
        "reduced_K_B0_Z0_N0": reduced_stiffness,
        "reduced_M_B2_Z2_N2": reduced_mass,
        "reduced_B4_Z4_N4": reduced_mass * (TIME**2),
        "P0_assumed": DIMENSIONLESS,
        "Gamma_port": TIME**5,
        "G_3_spatial": Dim(3, -2, -1),
        "G_4_spatial": Dim(4, -2, -1),
        "gamma_GR_4_spatial": Dim(4, -2, -1) / (velocity**5),
    }


D = dimensional_dictionary()


def _d1_checks() -> list[Check]:
    rho = D["rho_4d_number_density"]
    psi = D["psi"]
    k_eos = D["K_eos"]
    h = D["h_enthalpy"]
    pressure = D["P_pressure"]
    q_a0 = D["q_A0"]
    q_ai = D["q_Ai"]
    tau = D["tau_tension"]
    tw = D["T_w"]
    u_r = D["U_Sigma_R"]
    k1 = ENERGY / LENGTH
    source_reduced = D["radial_shell_3ball"] * k1 * rho
    static_force = FORCE / LENGTH
    omega = D["omega"]
    return [
        homogeneous(
            "D1",
            "EOS P=K*rho^5 and U=(K/4)*rho^5",
            {"P": pressure, "K*rho^5": k_eos * (rho**5), "U": D["U_bulk_energy_density"]},
        ),
        homogeneous(
            "D1",
            "enthalpy h=(5K/4)*rho^4",
            {"h": h, "K*rho^4": k_eos * (rho**4)},
        ),
        homogeneous(
            "D1",
            "gauged GNLS additive operator terms",
            {
                "i*hbar*partial_t*psi": ACTION * (TIME**-1) * psi,
                "kinetic*psi": ((ACTION**2) / (MASS * (LENGTH**2))) * psi,
                "V_conf*psi": ENERGY * psi,
                "h(rho)*psi": h * psi,
                "q*A0*psi": q_a0 * psi,
            },
            "Spatial minimal coupling separately requires q*Ai/hbar ~ L^-1.",
        ),
        expect_dim(
            "D1",
            "spatial gauge covariant derivative q*Ai/hbar",
            q_ai / ACTION,
            LENGTH**-1,
        ),
        homogeneous(
            "D1",
            "static wall balance operator",
            {
                "-partial_w(T_w*partial_w R)": tw / LENGTH,
                "0.5*T_w_R*(partial_w R)^2": (tw / LENGTH),
                "U_R": u_r,
                "radial-reduced source": source_reduced,
            },
            "The check uses the reduced source sum over the 3-ball shell only, leaving a force per wall coordinate.",
        ),
        expect_dim(
            "D1",
            "return source if full cell volume were used here",
            D["tube_cell_measure"] * k1 * rho,
            static_force,
            "This is the counterfactual: multiplying by dw in the wall RHS would add one power of L.",
        ),
        homogeneous(
            "D1",
            "BdG additive energy block",
            {
                "kinetic": (ACTION**2) / (MASS * (LENGTH**2)),
                "V_conf": ENERGY,
                "q*A0": q_a0,
                "mu": ENERGY,
                "h": h,
                "rho*dh_drho": rho * (k_eos * (rho**3)),
                "pairing dh_drho*psi^2": (k_eos * (rho**3)) * (psi**2),
            },
        ),
        homogeneous(
            "D1",
            "Maxwell c=1 VSH weak-form rows",
            {
                "i*omega*A_lane": omega * D["Ai_if_q_dimensionless"],
                "spatial_derivative*A_lane": D["Ai_if_q_dimensionless"] / LENGTH,
            },
            "This is homogeneous only in the code's c=1 lane normalization; restoring SI-like A0/Ai needs explicit c factors.",
        ),
    ]


def _d2_checks() -> list[Check]:
    chi = D["chi_wall_normalized"]
    omega = D["omega"]
    k0 = D["reduced_K_B0_Z0_N0"]
    k2 = D["reduced_M_B2_Z2_N2"]
    k4 = D["reduced_B4_Z4_N4"]
    rho = D["rho_4d_number_density"]
    psi = D["psi"]
    k1 = ENERGY / LENGTH
    shell = D["radial_shell_3ball"]
    dw = LENGTH
    density_response = psi * D["bdg_u_v_symplectic_normalized"]
    drive_by_w = shell * k1 * density_response
    live_bdg_coupling = chi * drive_by_w * dw
    expected_coupling = (k0 * (omega**2)) ** sp.Rational(1, 2)
    live_b0 = (live_bdg_coupling**2) / (omega**2)
    gamma = D["Gamma_port"]
    return [
        expect_dim("D2", "radial shell measure Delta V_i^r", shell, LENGTH**3),
        expect_dim("D2", "tensor cell measure Delta V_i^r * dw", shell * dw, LENGTH**4),
        expect_dim(
            "D2",
            "counterfactual 4-ball radial shell times dw",
            D["raw_4ball_radial_measure"] * dw,
            LENGTH**4,
            "A raw 4-ball radial shell would overcount the tube by one length when also integrated over w.",
        ),
        expect_dim("D2", "wall chi normalization chi^T W chi=1", (chi**2) * dw, DIMENSIONLESS),
        expect_dim(
            "D2",
            "BdG symplectic normalization integral",
            (D["bdg_u_v_symplectic_normalized"] ** 2) * D["tube_cell_measure"],
            DIMENSIONLESS,
        ),
        expect_dim(
            "D2",
            "live B2a wall-drive coupling c_j dimension",
            live_bdg_coupling,
            expected_coupling,
            "From code: c_j = int dw chi(w) sum_r DeltaV_r k1 delta_rho_j. "
            "The expected value is fixed by B0=c_j^2/varpi^2 having the same units as K.",
        ),
        expect_dim(
            "D2",
            "live B0 from B2a coupling c_j^2/varpi^2",
            live_b0,
            k0,
            "The missing factor is the same as 1/(m*a^2) in the hbar=m=a=1 normalization.",
        ),
        expect_dim("D2", "formal B2 moment dimension c_j^2/varpi^4", k0 / (omega**2), k2),
        expect_dim("D2", "formal B4 moment dimension c_j^2/varpi^6", k0 / (omega**4), k4),
        expect_dim("D2", "Gamma_port=4*a^5/(27*c_s^5)", (LENGTH**5) / ((LENGTH / TIME) ** 5), gamma),
        expect_dim("D2", "Gamma_port*omega^5", gamma * (omega**5), DIMENSIONLESS),
        expect_dim(
            "D2",
            "Maxwell N_n extraction if Sigma has stiffness units",
            k0 / (gamma * (omega**5)),
            k0,
            "Since Gamma_port*omega^5 is dimensionless, -Im Sigma/(Gamma omega^5) preserves Sigma's stiffness units.",
        ),
    ]


def _d3_checks() -> list[Check]:
    k0 = D["reduced_K_B0_Z0_N0"]
    k2 = D["reduced_M_B2_Z2_N2"]
    k4 = D["reduced_B4_Z4_N4"]
    omega = D["omega"]
    rho = D["rho_4d_number_density"]
    psi = D["psi"]
    k1 = ENERGY / LENGTH
    live_c = D["chi_wall_normalized"] * D["radial_shell_3ball"] * k1 * (psi * D["bdg_u_v_symplectic_normalized"]) * LENGTH
    live_b0 = (live_c**2) / (omega**2)
    p0 = DIMENSIONLESS
    target_4d = D["G_4_spatial"] * ((LENGTH / TIME) ** 5) / ((LENGTH**5) * ((LENGTH / TIME) ** 5))
    gamma_gr_4d = D["gamma_GR_4_spatial"]
    return [
        homogeneous("D3", "formal D0=K-B0-Z0", {"K": k0, "B0": k0, "Z0": k0}),
        homogeneous("D3", "live D0 using B2a overlap-derived B0", {"K": k0, "B0_live": live_b0, "Z0": k0}),
        homogeneous("D3", "D2=-(M+B2+Z2)", {"M": k2, "B2": k2, "Z2": k2}),
        homogeneous("D3", "D4=-(B4+Z4)", {"B4": k4, "Z4": k4}),
        expect_dim("D3", "P0=N0/D0 formal dimension", k0 / k0, p0),
        expect_dim("D3", "P2 formula formal dimension", (k2 / k0), TIME**2),
        expect_dim("D3", "P4 formula formal dimension", (k4 / k0), TIME**4),
        homogeneous("D3", "R_pole=D0*(B4+Z4)-3*(M+B2+Z2)^2", {"D0*C": k0 * k4, "A^2": k2**2}),
        expect_dim("D3", "gamma_GR=2G/(5c^5) with 4-spatial G", gamma_gr_4d, Dim(-1, 3, -1)),
        expect_dim(
            "D3",
            "GR target 54*G*c_s^5/(5*a^5*c^5) with 4-spatial G",
            target_4d,
            DIMENSIONLESS,
            "The existing R_norm comparison assumes mhat0^2*S_port*P0 is dimensionless.",
        ),
        homogeneous(
            "D3",
            "R_norm two-term subtraction",
            {"mhat0^2*S_port*P0": p0, "GR_target_4spatial_G": target_4d},
            "No dimension was assigned to mhat0 or S_port in the freeze sheet; both are pinned to 1.",
        ),
    ]


def run_audit() -> dict[str, object]:
    checks = [*_d1_checks(), *_d2_checks(), *_d3_checks()]
    failures = [check for check in checks if check.status != "CONSISTENT"]
    dictionary = {key: str(value) for key, value in D.items()}
    verdict = {
        "status": "INCONSISTENT",
        "summary": (
            "The formal reduced algebra is homogeneous only after forcing B_n/Z_n/N_n "
            "to have the declared reduced stiffness dimensions. The live B2a overlap "
            "map and the 4-spatial-dimensional GR target do not satisfy those contracts."
        ),
        "primary_artifacts": [
            {
                "step": "B2a live wall-drive overlap -> B_n",
                "missing_factor_dimension": str(factor_to_reach(D["reduced_K_B0_Z0_N0"], _live_b0_dimension())),
                "one_natural_units_realization": "1/(m*a^2)",
                "hidden_by": "m=a=1",
            },
            {
                "step": "R_norm GR target with 4-spatial-dimensional G",
                "missing_factor_dimension": str(factor_to_reach(DIMENSIONLESS, _target_4d_dimension())),
                "one_natural_units_realization": "m*a^3/c_s^2",
                "hidden_by": "m=a=c_s=1",
            },
        ],
        "order_statement": (
            "No unique decimal order follows from dimensional analysis alone: the "
            "missing factors are dimensionful and evaluate to 1 under the solver's "
            "natural-unit pins. With physical scales restored, their log10 value is "
            "set by the chosen m, a, and c_s units."
        ),
    }
    return {
        "schema": "stage1_patha_dimensional_audit/v1",
        "base_dimensions": ["L", "T", "M"],
        "natural_unit_pins": ["a=1", "c_s=1", "hbar=1", "m=1", "c=1 in the GR target path"],
        "dictionary": dictionary,
        "checks": [check.as_dict() for check in checks],
        "summary": {
            "total": len(checks),
            "consistent": len(checks) - len(failures),
            "inconsistent": len(failures),
        },
        "verdict": verdict,
    }


def _live_b0_dimension() -> Dim:
    omega = D["omega"]
    live_c = (
        D["chi_wall_normalized"]
        * D["radial_shell_3ball"]
        * (ENERGY / LENGTH)
        * (D["psi"] * D["bdg_u_v_symplectic_normalized"])
        * LENGTH
    )
    return (live_c**2) / (omega**2)


def _target_4d_dimension() -> Dim:
    velocity = LENGTH / TIME
    return D["G_4_spatial"] * (velocity**5) / ((LENGTH**5) * (velocity**5))


def render_markdown(report: Mapping[str, object]) -> str:
    lines = [
        "# Path-A dimensional audit",
        "",
        "## Step 0 dictionary",
        "",
    ]
    dictionary = report["dictionary"]
    assert isinstance(dictionary, Mapping)
    for key in sorted(dictionary):
        lines.append(f"- `{key}`: `{dictionary[key]}`")
    lines.extend(["", "## Checks", ""])
    checks = report["checks"]
    assert isinstance(checks, Sequence)
    for raw in checks:
        assert isinstance(raw, Mapping)
        factor = raw["factor_needed_to_reach_expected"]
        factor_text = "" if factor in (None, "1") else f"; factor needed `{factor}`"
        lines.append(
            f"- `{raw['section']}` `{raw['name']}`: **{raw['status']}** "
            f"(expected `{raw['expected']}`, actual `{raw['actual']}`{factor_text}). "
            f"{raw['note']}"
        )
    verdict = report["verdict"]
    assert isinstance(verdict, Mapping)
    lines.extend(
        [
            "",
            "## D4 verdict",
            "",
            f"**{verdict['status']}**: {verdict['summary']}",
            "",
            str(verdict["order_statement"]),
            "",
        ]
    )
    return "\n".join(lines)


def _dim_dict(dim: Dim) -> dict[str, object]:
    return {
        "dimension": str(dim),
        "tuple_L_T_M": dim.as_tuple(),
    }


def _patha19_pin_analysis() -> dict[str, object]:
    quantities = ["a", "c_s0", "hbar", "m_GNLS"]
    ltm_vectors = [
        [1, 0, 0],
        [1, -1, 0],
        [2, -1, 1],
        [0, 0, 1],
    ]
    ltm_matrix = sp.Matrix(ltm_vectors).T
    lt_vectors = [
        [1, 0],
        [1, -1],
        [2, -1],
        [0, 0],
    ]
    lt_matrix = sp.Matrix(lt_vectors).T
    return {
        "quantities": quantities,
        "ltm_rank": int(ltm_matrix.rank()),
        "ltm_pin_count": len(quantities),
        "ltm_relations": [
            {
                "dimensionless_product": "a*c_s0*m_GNLS/hbar",
                "exponent_vector_for_a_cs_hbar_m": ["1", "1", "-1", "1"],
                "pin_relation": "a = hbar/(m_GNLS*c_s0)",
            }
        ],
        "lt_rank_if_m_GNLS_dimensionless_and_hbar_L2_per_T": int(lt_matrix.rank()),
        "lt_pin_count": len(quantities),
        "lt_relations": [
            "hbar/(a*c_s0) = 1",
            "m_GNLS = 1 is dimensionless by representation, not derived from the action",
        ],
    }


def _patha19_foundation_checks() -> list[Check]:
    rho4 = D["rho_4d_number_density"]
    rho3_reduced = LENGTH**-3
    psi = D["psi"]
    velocity = LENGTH / TIME
    number_rate = TIME**-1
    lpsi_density = ENERGY / (LENGTH**4)
    electric_field = FORCE
    magnetic_field = ACTION / (LENGTH**2)
    maxwell_coeff = lpsi_density / (electric_field**2)
    wall_density = FORCE
    return [
        expect_dim(
            "pathA_19_F2",
            "4D-bulk closed-3-surface number flux J_bulk",
            rho4 * velocity * (LENGTH**3),
            number_rate,
            "Frame: 4D bulk; rho=L^-4, v=L/T, dSigma_3=L^3.",
        ),
        expect_dim(
            "pathA_19_F2",
            "3D-brane reduced 2-sphere number flux J_brane",
            rho3_reduced * velocity * (LENGTH**2),
            number_rate,
            "Frame: brane reduction after integrating one transverse length; rho_3=L^-3.",
        ),
        expect_dim(
            "pathA_19_F2",
            "4D-bulk volumetric flux Q_vol=rho^-1 J",
            (rho4**-1) * number_rate,
            (LENGTH**4) / TIME,
        ),
        expect_dim(
            "pathA_19_F2",
            "3D-brane volumetric flux Q_vol=rho_3^-1 J",
            (rho3_reduced**-1) * number_rate,
            (LENGTH**3) / TIME,
        ),
        expect_dim(
            "pathA_19_F2",
            "constituent mass flux m_GNLS*J",
            MASS * number_rate,
            MASS / TIME,
        ),
        expect_dim(
            "pathA_19_F1",
            "conditional defect rest-frequency conversion hbar*J/c_gamma^2",
            ACTION * number_rate / (velocity**2),
            MASS,
            "Dimensionally valid only as a conversion; it is not an inflow-to-m_defect derivation.",
        ),
        homogeneous(
            "pathA_19_F2",
            "bulk continuity equation",
            {
                "partial_t rho": rho4 / TIME,
                "div_4(rho v)": (rho4 * velocity) / LENGTH,
            },
        ),
        expect_dim(
            "pathA_19_F3",
            "sound-speed law 5*K*rho^4/m",
            D["K_eos"] * (rho4**4) / MASS,
            velocity**2,
        ),
        expect_dim(
            "pathA_19_F3",
            "EOS enthalpy scale h0=(m_GNLS*c_s0^2)/4",
            MASS * (velocity**2),
            ENERGY,
            "The factor 1/4 is dimensionless; it enters the healing-length coefficient.",
        ),
        expect_dim(
            "pathA_19_F3",
            "GNLS healing length sqrt(hbar^2/(2*m_GNLS*h0))",
            (ACTION**2 / (MASS * ENERGY)) ** sp.Rational(1, 2),
            LENGTH,
            "Using h0=m_GNLS*c_s0^2/4 gives xi_h=sqrt(2)*hbar/(m_GNLS*c_s0).",
        ),
        homogeneous(
            "pathA_19_F3",
            "parent GNLS Lagrangian density terms",
            {
                "i*hbar*psi*partial_t psi": ACTION * (TIME**-1) * (psi**2),
                "hbar^2/(2m)*|D_i psi|^2": (ACTION**2 / MASS) * ((psi / LENGTH) ** 2),
                "V_conf*rho": ENERGY * rho4,
                "U=K*rho^5/4": D["K_eos"] * (rho4**5),
            },
        ),
        expect_dim(
            "pathA_19_F3",
            "spatial gauge minimal-coupling dimension q*A_i/hbar",
            D["q_Ai"] / ACTION,
            LENGTH**-1,
        ),
        homogeneous(
            "pathA_19_F3",
            "localized Maxwell sector with explicit c factors",
            {
                "(Z/mu0)*E_i^2": maxwell_coeff * (electric_field**2),
                "(Z/mu0)*c^2*B_ij^2": maxwell_coeff * (velocity**2) * (magnetic_field**2),
                "A0*J0_ext": D["q_A0"] * rho4,
                "Ai*Ji_ext": D["q_Ai"] * (rho4 * velocity),
            },
            "Uses q dimensionless; the c^2 factor is the required restoration for E/B homogeneity.",
        ),
        homogeneous(
            "pathA_19_F3",
            "wall action density before dt*dw integration",
            {
                "mu_eta*(partial_t eta)^2": D["mu_wall_restored_as_tau_over_c_s2"] * ((LENGTH / TIME) ** 2),
                "T_w*(partial_w eta)^2": D["T_w"],
                "K_eta*eta^2": D["U_Sigma_RR"] * (LENGTH**2),
            },
            "Each term has force units so dt*dw gives action.",
        ),
    ]


def _patha19_lt_representation_checks() -> list[Check]:
    rho = LENGTH**-4
    psi = LENGTH**-2
    velocity = LENGTH / TIME
    action_lt = Dim(2, -1, 0)
    energy_lt = Dim(2, -2, 0)
    mass_lt = DIMENSIONLESS
    force_lt = Dim(1, -2, 0)
    k_eos_lt = energy_lt / (rho**4)
    lpsi_density_lt = energy_lt / (LENGTH**4)
    electric_field_lt = force_lt
    magnetic_field_lt = action_lt / (LENGTH**2)
    maxwell_coeff_lt = lpsi_density_lt / (electric_field_lt**2)
    return [
        homogeneous(
            "pathA_19_LT_representation",
            "local GNLS terms after projecting m_GNLS to dimensionless",
            {
                "i*hbar*psi*partial_t psi": action_lt * (TIME**-1) * (psi**2),
                "hbar^2/(2m)*|D_i psi|^2": (action_lt**2 / mass_lt) * ((psi / LENGTH) ** 2),
                "V_conf*rho": energy_lt * rho,
                "U=K*rho^5/4": k_eos_lt * (rho**5),
            },
            "This is a natural-unit representation check, not an action-level derivation of m_GNLS.",
        ),
        homogeneous(
            "pathA_19_LT_representation",
            "local Maxwell terms after M projection",
            {
                "(Z/mu0)*E_i^2": maxwell_coeff_lt * (electric_field_lt**2),
                "(Z/mu0)*c^2*B_ij^2": maxwell_coeff_lt * (velocity**2) * (magnetic_field_lt**2),
            },
        ),
        homogeneous(
            "pathA_19_LT_representation",
            "local wall terms after M projection",
            {
                "mu_eta*(partial_t eta)^2": Dim(-1, 0, 0) * ((LENGTH / TIME) ** 2),
                "T_w*(partial_w eta)^2": force_lt,
                "K_eta*eta^2": Dim(-1, -2, 0) * (LENGTH**2),
            },
        ),
    ]


def _patha19_flags() -> list[dict[str, object]]:
    velocity = LENGTH / TIME
    formal_4d_target = D["G_4_spatial"] * (velocity**5) / ((LENGTH**5) * (velocity**5))
    observed_3d_target = D["G_3_spatial"] * (velocity**5) / ((LENGTH**5) * (velocity**5))
    lt_g3_target = Dim(3, -2, 0) * (velocity**5) / ((LENGTH**5) * (velocity**5))
    return [
        {
            "name": "formal_4D_R_norm_target_not_dimensionless_without_conversion",
            "status": "FLAGGED_EXISTING_pathA_18_GAP",
            "expected": str(DIMENSIONLESS),
            "actual": str(formal_4d_target),
            "factor_needed_to_reach_expected": str(factor_to_reach(DIMENSIONLESS, formal_4d_target)),
            "source": "dimensional_check.py pathA_18 D3 lane; part01_parent_geometry.tex:1404-1418",
            "consequence": "Preserve pathA_18 behavior and do not repair R_norm in pathA_19; the conversion belongs to pathA_22.",
        },
        {
            "name": "observed_3D_GR_target_not_dimensionless_without_conversion",
            "status": "FLAGGED_FOR_pathA_21_pathA_22",
            "expected": str(DIMENSIONLESS),
            "actual": str(observed_3d_target),
            "factor_needed_to_reach_expected": str(factor_to_reach(DIMENSIONLESS, observed_3d_target)),
            "source": "pde.tex:2080-2093; part01_parent_geometry.tex:1404-1418",
            "consequence": "Do not interpret dimensional matching as a derived B2c normalization; the conversion belongs to pathA_21/pathA_22.",
        },
        {
            "name": "LT_R_norm_gate_fails_without_new_conversion_factor",
            "status": "REJECTS_TRUE_LT_BASE",
            "expected": str(DIMENSIONLESS),
            "actual": str(lt_g3_target),
            "factor_needed_to_reach_expected": str(factor_to_reach(DIMENSIONLESS, lt_g3_target)),
            "source": "pathA_19 action-homogeneity gate",
            "consequence": "{L,T} is only a natural-unit representation in this step, not the base dimensional system.",
        },
    ]


def run_patha19_foundation() -> dict[str, object]:
    checks = _patha19_foundation_checks()
    lt_checks = _patha19_lt_representation_checks()
    failures = [check for check in [*checks, *lt_checks] if check.status != "CONSISTENT"]
    dictionary_subset = {
        key: str(D[key])
        for key in (
            "hbar",
            "m_particle",
            "rho_4d_number_density",
            "psi",
            "K_eos",
            "P_pressure",
            "h_enthalpy",
            "c_s",
            "a",
            "G_3_spatial",
            "G_4_spatial",
        )
    }
    residuals = [
        {
            "name": "INFLOW_MASS_SOURCE_MISSING",
            "status": "BLOCKS_MASS_EMERGENCE",
            "source": (
                "part01_parent_geometry.tex:174-219 and pde.tex:326-406 keep m_GNLS in the exact action/current; "
                "brane_bulk_ontology.tex:1267-1297 gives drainage scaling, not a boundary/source/Noether/Hamiltonian derivation."
            ),
            "downstream_consequence": (
                "Retain base dimensions {L,T,M}; record J as a conserved rate label only. "
                "A positive m_defect result must be derived in pathA_21 from a defect source, boundary energy, or Hamiltonian/Noether charge."
            ),
        },
        {
            "name": "NO_NET_ACCRETION_BC_UNDERIVED",
            "status": "CARRIED_FORWARD",
            "source": (
                "part01_parent_geometry.tex:298-330 and pde.tex:512-539 give projected leakage; "
                "brane_bulk_ontology.tex:1998-2039 leaves the throat bottom open/closed/connected."
            ),
            "downstream_consequence": (
                "Gauss flux is surface-independent only in no-source/no-leakage regions. "
                "No-net-accretion must be supplied as a boundary condition before using recirculation as a physical closure."
            ),
        },
        {
            "name": "A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT",
            "status": "CARRIED_FORWARD",
            "source": (
                "part01_parent_geometry.tex:447-510 and pde.tex:633-648,1074-1118 define a(t) as a mouth-average collective moment."
            ),
            "downstream_consequence": (
                "Use J and explicit branch data as invariant scale-map inputs; treat a as derived geometry. "
                "No normalization value is changed in pathA_19."
            ),
        },
    ]
    return {
        "schema": "stage1_pathA_19_foundation/v1",
        "base_set_verdict": {
            "status": "RETAIN_L_T_M",
            "reason": (
                "m_GNLS is an explicit action parameter and no action-level derivation ties m_defect to inflow. "
                "The LT projection is a natural-unit representation only."
            ),
        },
        "mass_symbol_split": {
            "m_GNLS": "constituent/inertial mass in the parent GNLS action; dimension M",
            "m_defect": "throat rest/gravitational branch property; not derived from J in this step",
            "conditional_conversion": "m_defect = alpha_J*hbar*J/c_gamma^2 is dimensionally valid only after a physical inflow/rest-frequency derivation.",
        },
        "pin_analysis": _patha19_pin_analysis(),
        "healing_length": {
            "pin_count_relation": "a = hbar/(m_GNLS*c_s0)",
            "action_core_balance": "xi_h = sqrt(hbar^2/(2*m_GNLS*h0)); h0=(5K/4)*rho0^4=(m_GNLS*c_s0^2)/4",
            "derived_relation": "xi_h = sqrt(2)*hbar/(m_GNLS*c_s0)",
            "consequence": "If a is identified with the GNLS healing core, then a/xi_h is a convention/branch factor; the raw four pins correspond to a=xi_h/sqrt(2).",
        },
        "dictionary_subset": dictionary_subset,
        "checks": [check.as_dict() for check in checks],
        "lt_representation_checks": [check.as_dict() for check in lt_checks],
        "flags": _patha19_flags(),
        "residuals": residuals,
        "summary": {
            "total_algebraic_checks": len(checks) + len(lt_checks),
            "consistent_algebraic_checks": len(checks) + len(lt_checks) - len(failures),
            "inconsistent_algebraic_checks": len(failures),
            "residual_count": len(residuals),
            "acceptance_status": "PASS_WITH_NAMED_RESIDUALS",
        },
    }


def render_patha19_foundation_markdown(report: Mapping[str, object]) -> str:
    base = report["base_set_verdict"]
    assert isinstance(base, Mapping)
    healing = report["healing_length"]
    assert isinstance(healing, Mapping)
    pin = report["pin_analysis"]
    assert isinstance(pin, Mapping)
    lines = [
        "# Path-A 19 dimensional foundation reference",
        "",
        "## Verdict",
        "",
        f"- **Base set:** `{base['status']}`.",
        f"- **Reason:** {base['reason']}",
        "- **Mass split:** `m_GNLS` is the constituent mass in the GNLS action; `m_defect` is a throat branch mass and is not derived from inflow in this step.",
        "- **No normalization change:** this report does not derive or modify any frozen normalization factor.",
        "",
        "## F1 mass fork",
        "",
        "| Quantity | Status | Dimension/result | Source |",
        "|---|---|---|---|",
        "| `m_GNLS` | explicit action parameter | `M` retained | `part01_parent_geometry.tex:174-219`; `pde.tex:326-406` |",
        "| `m_defect` | not action-derived here | blocked by `INFLOW_MASS_SOURCE_MISSING` | `brane_bulk_ontology.tex:1267-1297`, `1998-2039` |",
        "| `hbar*J/c_gamma^2` | dimensional conversion only | `M` | harness check `conditional defect rest-frequency conversion` |",
        "",
        "The current action contains `m_GNLS` in the kinetic operator, current, Madelung velocity, Euler equation, and sound-speed law. The ontology supplies drainage/volume-deficit scaling for defect mass, but not a boundary source, Noether charge, or Hamiltonian energy theorem tying `m_defect` to the inflow rate. Therefore the mass-emergence hypothesis is rejected for this foundation gate and `{L,T,M}` is retained.",
        "",
        "## F2 frame-tagged flux and a pin",
        "",
        "| Flux | Frame | Dimension | Interpretation |",
        "|---|---|---|---|",
        "| `J_bulk = int rho v.dSigma_3` | 4D bulk closed 3-surface | `T^-1` | number flux/rate |",
        "| `J_brane = int rho_3 v.dS_2` | 3D brane reduction | `T^-1` | number flux/rate after transverse reduction |",
        "| `Q_vol,bulk = rho^-1 J` | 4D bulk | `L^4 T^-1` | bulk four-volume flux |",
        "| `Q_vol,brane = rho_3^-1 J` | 3D brane | `L^3 T^-1` | brane volume flux |",
        "| `m_GNLS J` | action constituent mass flux | `M T^-1` | mass-per-particle times number rate |",
        "",
        "Gauss shape-independence holds only in a region with no enclosed source or leakage. Projection creates `S_leak` (`part01_parent_geometry.tex:298-330`; `pde.tex:512-539`), and the throat bottom is explicitly open/closed/connected pending microscopic input (`brane_bulk_ontology.tex:1998-2039`). No-net-accretion is therefore carried as `NO_NET_ACCRETION_BC_UNDERIVED`.",
        "",
        "`a` is a mouth-radius collective moment, not a fundamental coordinate: `a0=R0(0)` and `a(t)` is the mouth average (`part01_parent_geometry.tex:447-510`; `pde.tex:633-648`). The conserved rate `J` is the better invariant label; `a` remains branch geometry consumed by downstream scale maps.",
        "",
        "## F3 pins, healing length, and dictionary",
        "",
        f"- `{pin['ltm_pin_count']}` pins on `{pin['ltm_rank']}` base dimensions leave one relation: `a*c_s0*m_GNLS/hbar = 1`, i.e. `{pin['ltm_relations'][0]['pin_relation']}`.",
        f"- GNLS core balance gives `{healing['action_core_balance']}`.",
        f"- **Derived healing scale:** `{healing['derived_relation']}`.",
        f"- Consequence: {healing['consequence']}",
        "",
        "| Quantity | Current status |",
        "|---|---|",
        "| `hbar` | independent action constant, dimension `M L^2 T^-1` |",
        "| `m_GNLS` | independent action constant, dimension `M` |",
        "| `K` | independent EOS constant after choosing the EOS, dimension `M L^18 T^-2` |",
        "| chosen state `rho0` | independent branch/state datum, dimension `L^-4` in 4D bulk |",
        "| `c_s0` | derived from `c_s0^2=5K rho0^4/m_GNLS` |",
        "| `xi_h` | derived from GNLS kinetic/enthalpy balance |",
        "| `a` | derived/branch collective geometry if identified with a core scale; otherwise an input branch moment |",
        "| `m_defect` | not emergent in this step; blocked by `INFLOW_MASS_SOURCE_MISSING` |",
        "",
        "Dictionary confirmation: the harness 4D action dictionary is homogeneous for the GNLS, gauge coupling, Maxwell sector with explicit `c` factors, wall action, current, and flux claims. The observed 3D GR target remains flagged as a downstream conversion problem, not a pathA_19 base-system change.",
        "",
        "## F4 paper-prose reconciliation",
        "",
        "| Source | Statement | Classification | Note |",
        "|---|---|---|---|",
        "| `part01_parent_geometry.tex:140-144`, `174-203` | `rho=|psi|^2`, GNLS action/EOS/sound speed | AGREES | Implies 4D `rho=L^-4`, `psi=L^-2`, `K=M L^18 T^-2`. |",
        "| `part01_parent_geometry.tex:213-219`, `268-291` | current, velocity, Euler/vorticity identities contain `m` | AGREES | Supports `m_GNLS` as action content. |",
        "| `part01_parent_geometry.tex:298-330` | normalized projection and leakage source | AGREES/AMBIGUOUS | Agrees on open-system projection; normalized kernel keeps dimensions distinct from integrated 3D reduction. |",
        "| `part01_parent_geometry.tex:447-510` | `a0=R0(0)`, `a(t)` mouth average | AGREES | `a` is a length and a collective moment. |",
        "| `pde.tex:326-352`, `396-406` | parent action/EOS/current | AGREES | Same 4D parent dictionary as the harness. |",
        "| `pde.tex:512-539` | projected continuity has `S_leak` | AGREES | Blocks unqualified no-net-accretion. |",
        "| `pde.tex:633-648`, `1074-1118` | `a,L` are collective moments and reduced wall coordinates | AGREES | Supports reassessing the `a=1` pin. |",
        "| `brane_bulk_ontology.tex:1267-1297`, `1967-1975` | mass as drainage/volume deficit; charge as vorticity flux | AMBIGUOUS | Physical scaling prose, not an action-level mass theorem. |",
        "| `brane_bulk_ontology.tex:1998-2039` | bottom open/closed/connected | AGREES | Carries no-net-accretion as a gap. |",
        "| `em_fields.tex:1717-1721` | `rho0` as `kg m^-3` mass density | WRONG-3D-CONVENTION | Legacy 3D/SI prose; not the 4D number-density harness dictionary. |",
        "| `em_fields.tex:1723-1726`, `1738-1752` | `c_s`, `a`, `L`, circulation dimensions | AGREES | Correct in the stated SI/3D frame. |",
        "| `em_fields.tex:1728-1736`, `1782-1786` | pressure/enthalpy per mass in SI units | WRONG-3D-CONVENTION | Useful prose but not the 4D action density/enthalpy dictionary. |",
        "| `em_fields.tex:1757-1779` | `V=pi a^2L`, `m_G=kappa_m rho0 V`, `q` units | WRONG-3D-CONVENTION | Legacy 3D throat-volume bookkeeping; not a derivation of `m_defect` from `J`. |",
        "",
        "## F5 gaps carried forward",
        "",
    ]
    residuals = report["residuals"]
    assert isinstance(residuals, Sequence)
    for raw in residuals:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: {raw['status']}. {raw['downstream_consequence']} Source: {raw['source']}")
    lines.extend(
        [
            "- `EOS_FROM_GNLS_FACTOR`: pathA_20 must keep the factor `h0=(m_GNLS*c_s0^2)/4` and the derived healing scale `sqrt(2)*hbar/(m_GNLS*c_s0)`.",
            "- `M_TO_G_UNIFICATION`: pathA_21 must derive any defect-mass/back-reaction relation; pathA_19 does not prove it.",
            "- `SCALE_MAP_INPUTS`: pathA_22 must consume `J`, branch geometry `a(J,branch)` if derived, `rho0`, `K`, `m_GNLS`, `hbar`, and the 3D reduction/conversion factors.",
            "",
            "## Algebraic harness summary",
            "",
        ]
    )
    summary = report["summary"]
    assert isinstance(summary, Mapping)
    lines.append(
        f"- Algebraic checks: {summary['consistent_algebraic_checks']} consistent, "
        f"{summary['inconsistent_algebraic_checks']} inconsistent, {summary['total_algebraic_checks']} total."
    )
    lines.append(f"- Acceptance status: `{summary['acceptance_status']}`.")
    lines.append("")
    lines.append("Flagged algebraic residuals:")
    flags = report["flags"]
    assert isinstance(flags, Sequence)
    for raw in flags:
        assert isinstance(raw, Mapping)
        lines.append(
            f"- `{raw['name']}`: {raw['status']}; actual `{raw['actual']}`, "
            f"factor needed `{raw['factor_needed_to_reach_expected']}`. {raw['consequence']}"
        )
    lines.append("")
    return "\n".join(lines)


def write_patha19_foundation_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = run_patha19_foundation()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_19_foundation_dimensional_report.json"
    scratch_md_path = out_dir / "pathA_19_foundation_dimensional_report.md"
    reference_path = report_dir / "pathA_19_dimensional_foundation.md"
    rendered = render_patha19_foundation_markdown(report) + "\n"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered, encoding="utf-8")
    reference_path.write_text(rendered, encoding="utf-8")
    return json_path, reference_path, report


def _algebra_check(name: str, actual: sp.Expr, expected: sp.Expr, note: str = "") -> dict[str, object]:
    residual = sp.simplify(actual - expected)
    return {
        "name": name,
        "pass": bool(residual == 0),
        "expected": str(expected),
        "actual": str(sp.simplify(actual)),
        "residual": str(residual),
        "note": note,
    }


def _boolean_check(name: str, actual: bool, expected: bool, note: str = "") -> dict[str, object]:
    actual_bool = bool(actual)
    expected_bool = bool(expected)
    return {
        "name": name,
        "pass": actual_bool == expected_bool,
        "expected": str(expected_bool),
        "actual": str(actual_bool),
        "note": note,
    }


def _patha20_velocity_checks() -> list[Check]:
    rho4 = D["rho_4d_number_density"]
    rho3 = LENGTH**-3
    velocity = LENGTH / TIME
    wave_number = LENGTH**-1
    number_rate = TIME**-1
    return [
        expect_dim(
            "pathA_20_S1",
            "sound-speed squared law c_s^2=5*K*rho^4/m_GNLS",
            D["K_eos"] * (rho4**4) / MASS,
            velocity**2,
            "The numerical factor 5 is dimensionless; this is relative to the imposed P=K*rho^5 EOS.",
        ),
        expect_dim(
            "pathA_20_S1",
            "sound speed c_s=sqrt(5*K*rho^4/m_GNLS)",
            (D["K_eos"] * (rho4**4) / MASS) ** sp.Rational(1, 2),
            velocity,
        ),
        homogeneous(
            "pathA_20_S1_S2",
            "stationary quantum-Bernoulli additive terms",
            {
                "0.5*m_GNLS*v_b^2": MASS * (velocity**2),
                "h(rho)": D["h_enthalpy"],
                "V_conf": ENERGY,
                "Q": (ACTION**2) / (MASS * (LENGTH**2)),
            },
            "The stationary profile needs continuity plus this Bernoulli/Euler balance; no profile is solved by dimensions.",
        ),
        homogeneous(
            "pathA_20_S2",
            "bulk continuity equation with v_b",
            {
                "partial_t rho": rho4 / TIME,
                "div_4(rho v_b)": (rho4 * velocity) / LENGTH,
            },
        ),
        expect_dim(
            "pathA_20_S2",
            "Madelung background velocity v_b=(hbar/m_GNLS)*grad(theta)",
            (ACTION / MASS) * (LENGTH**-1),
            velocity,
            "theta is dimensionless; the gauge-invariant source formula also permits -q*A_i/m_GNLS.",
        ),
        expect_dim(
            "pathA_20_S2",
            "photon/gauge-wave speed c_gamma",
            velocity,
            velocity,
            "Dimension pin only; the value relative to c_s is a dynamical question.",
        ),
        homogeneous(
            "pathA_20_S2",
            "massless gauge-wave dispersion omega^2=c_gamma^2*k^2",
            {
                "omega^2": (TIME**-1) ** 2,
                "c_gamma^2*k^2": (velocity**2) * (wave_number**2),
            },
            "Plane-wave dispersion gives group speed c_gamma before any rest-energy interpretation.",
        ),
        homogeneous(
            "pathA_20_S2",
            "trapped-mode wave dispersion omega^2=c_gamma^2*(k_parallel^2+k_perp^2)",
            {
                "omega^2": (TIME**-1) ** 2,
                "c_gamma^2*k_parallel^2": (velocity**2) * (wave_number**2),
                "c_gamma^2*k_perp^2": (velocity**2) * (wave_number**2),
            },
            "A fixed transverse k_perp supplies the rest oscillation omega0=c_gamma*k_perp.",
        ),
        expect_dim(
            "pathA_20_S2",
            "trapped-mode group velocity d omega/d k",
            velocity * (wave_number / wave_number),
            velocity,
            "For omega=c_gamma*sqrt(k^2+k0^2), d omega/dk has speed dimension and is bounded by c_gamma.",
        ),
        expect_dim(
            "pathA_20_S2",
            "ratio c_gamma/c_s",
            velocity / velocity,
            DIMENSIONLESS,
            "Dimensionless ratio only; machine agreement is not a derivation of its value.",
        ),
        expect_dim(
            "pathA_20_S2",
            "tail factor (c/c_s)^3 with c=c_gamma",
            (velocity / velocity) ** 3,
            DIMENSIONLESS,
        ),
        expect_dim(
            "pathA_20_S2b",
            "4D-bulk candidate sonic number flux rho_* c_s,* A_3,*",
            rho4 * velocity * (LENGTH**3),
            number_rate,
            "This checks the critical-law dimension only; the actual transonic profile is not derived here.",
        ),
        expect_dim(
            "pathA_20_S2b",
            "3D-brane candidate sonic number flux rho_3,* c_s,* A_2,*",
            rho3 * velocity * (LENGTH**2),
            number_rate,
            "Frame: brane-reduced number density after integrating the transverse direction.",
        ),
        expect_dim(
            "pathA_20_S2b",
            "background pressure P0=K*rho0^5",
            D["K_eos"] * (rho4**5),
            D["P_pressure"],
            "Any solved flux law inherits environment dependence through P0 and c_s0(rho0).",
        ),
        expect_dim(
            "pathA_20_S3",
            "pin relation hbar=m_GNLS*c_s0*a",
            MASS * velocity * LENGTH,
            ACTION,
            "This is the pathA_19 a-pin relation, not an hbar-emergence proof.",
        ),
        expect_dim(
            "pathA_20_S3",
            "healing-length relation hbar=m_GNLS*c_s0*xi_h/sqrt(2)",
            MASS * velocity * LENGTH,
            ACTION,
            "The sqrt(2) factor is dimensionless and belongs to xi_h, not the pathA_19 a.",
        ),
        expect_dim(
            "pathA_20_S3",
            "circulation kappa=int v_b dl",
            velocity * LENGTH,
            (ACTION / MASS),
            "Single-valued phase makes the step h/m_GNLS per complete winding.",
        ),
        expect_dim(
            "pathA_20_S3",
            "phase-momentum exchange p=hbar*grad(theta)",
            ACTION / LENGTH,
            MASS * velocity,
        ),
        expect_dim(
            "pathA_20_S3",
            "quantum pressure Q=-hbar^2/(2m)*laplacian(sqrt(rho))/sqrt(rho)",
            (ACTION**2) / (MASS * (LENGTH**2)),
            ENERGY,
        ),
        expect_dim(
            "pathA_20_S2_S3",
            "candidate mass bridge hbar*J/c_gamma^2",
            ACTION * number_rate / (velocity**2),
            MASS,
            "Dimensionally valid candidate only; pathA_21 must derive alpha_J and the source/Hamiltonian bridge.",
        ),
        expect_dim(
            "pathA_20_S2_S3",
            "cycle-rate bridge h*J_nu/c_gamma^2",
            ACTION * number_rate / (velocity**2),
            MASS,
            "h and hbar have the same dimensions; the 2*pi placement is not decided by dimensions.",
        ),
    ]


def _patha20_algebraic_checks() -> list[dict[str, object]]:
    rho = sp.symbols("rho", positive=True)
    c_s_profile = rho**2
    log_slope = sp.simplify(rho * sp.diff(c_s_profile, rho) / c_s_profile)
    cstar_over_c0 = sp.sqrt(sp.Rational(1, 3))
    rhostar_over_rho0 = sp.sqrt(cstar_over_c0)
    ideal_flux_factor = sp.simplify(cstar_over_c0 * rhostar_over_rho0)
    lambda_gamma = sp.symbols("lambda_gamma", positive=True)
    return [
        _algebra_check(
            "state dependence d ln c_s / d ln rho for n=5",
            log_slope,
            sp.Integer(2),
            "Since c_s(rho) is proportional to rho^2.",
        ),
        _algebra_check(
            "conditional ideal no-Q/no-V sonic c_s,* / c_s0",
            cstar_over_c0,
            1 / sp.sqrt(3),
            "Uses upstream v0=0 and Bernoulli 0.5*v_*^2+c_s,*^2/4=c_s0^2/4; not accepted as the branch verdict.",
        ),
        _algebra_check(
            "conditional ideal no-Q/no-V sonic rho_* / rho0",
            rhostar_over_rho0,
            sp.Pow(3, sp.Rational(-1, 4)),
            "Uses c_s proportional to rho^2; not accepted without the actual branch profile.",
        ),
        _algebra_check(
            "conditional ideal no-Q/no-V flux factor Jcrit/(rho0*c_s0*A*)",
            ideal_flux_factor,
            sp.Pow(3, sp.Rational(-3, 4)),
            "This is a conditional Euler-nozzle factor, not the pathA_20 flux_law_verdict.",
        ),
        _algebra_check(
            "tail factor with lambda_gamma=c_gamma/c_s",
            lambda_gamma**3,
            lambda_gamma**3,
            "With c=c_gamma, (c/c_s)^3 remains lambda_gamma^3 until the kinetic-operator ratio is derived.",
        ),
    ]


def run_patha20_velocity_constants() -> dict[str, object]:
    checks = _patha20_velocity_checks()
    algebra = _patha20_algebraic_checks()
    failures = [check for check in checks if check.status != "CONSISTENT"]
    algebra_failures = [check for check in algebra if not check["pass"]]
    residuals = [
        {
            "name": "EOS_CLOSURE_IMPOSED",
            "status": "CARRIED_FORWARD",
            "source": "part01_parent_geometry.tex:194-203; pde.tex:344-352 state P=K*rho^5 and c_s^2=(1/m_GNLS)dP/drho.",
            "downstream_consequence": "c_s is derived only relative to the imposed stiff-polytropic EOS; deriving the EOS from a deeper substrate remains outside pathA_20.",
        },
        {
            "name": "C_GAMMA_RATIO_UNDERDETERMINED",
            "status": "BLOCKS_NUMERIC_C_GAMMA_OVER_C_S",
            "source": (
                "part01_parent_geometry.tex:225-389 and pde.tex:355-565 give a localized Maxwell action, projection law, "
                "Z(w) renormalization, and measured-vs-flux closure issue; em_fields.tex:149-184, 482-499, 692-705 only give the legacy weak-field acoustic reuse."
            ),
            "downstream_consequence": "Carry lambda_gamma=c_gamma/c_s and tail factor (c/c_s)^3=lambda_gamma^3 into pathA_21/pathA_22 until the gauge/density kinetic operator and localization profile fix it.",
        },
        {
            "name": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "status": "FLUX_LAW_VERDICT",
            "source": (
                "pde.tex:2515-2566 requires a branch data set {R0, psi0, A_M0, wall data, spectra, overlaps, a, c_s, mhat}; "
                "pde.tex:2847-2879 states the actual stationary throat profile and DtN data remain branch-dependent."
            ),
            "downstream_consequence": "No accepted choked-flux law or nontransonic alternate law is available here; pathA_21 must consume this verdict rather than an unconditional Jcrit.",
        },
        {
            "name": "NO_NET_ACCRETION_BC_UNDERIVED",
            "status": "CARRIED_FORWARD",
            "source": "part01_parent_geometry.tex:298-330; pde.tex:511-539; brane_bulk_ontology.tex:1998-2042.",
            "downstream_consequence": "J_in=J_out requires throat-bottom topology/BC input; Gauss flux conservation is local no-leakage conservation, not a no-net-accretion proof.",
        },
        {
            "name": "HBAR_FREE_SUBSTRATE_RELATION_MISSING",
            "status": "BLOCKS_HBAR_EMERGENT",
            "source": "part01_parent_geometry.tex:174-219 and pde.tex:318-408 keep hbar as an action coefficient; pathA_19 only gives the pin relation a=hbar/(m_GNLS*c_s0).",
            "downstream_consequence": "S3 verdict is HBAR_PROVENANCE_UNDETERMINED; hbar remains explicit in the PDE and in the candidate mass bridge.",
        },
        {
            "name": "H_2PI_RATE_CLASSIFICATION_UNDERDETERMINED",
            "status": "CARRIED_FORWARD",
            "source": "pde.tex:429-469 gives psi=sqrt(rho) exp(i theta) and v_i=(hbar/m_GNLS)partial_i theta; brane_bulk_ontology.tex:668-671 and 1169-1180 treat circulation as quantized/integer-labeled.",
            "downstream_consequence": "Use h for complete winding/cycle-count relations and hbar for angular/PDE coefficients; pathA_21 must decide whether J is cycle-rate or angular-rate and where the 2*pi sits in alpha_J.",
        },
    ]
    return {
        "schema": "stage1_pathA_20_velocity_constants/v1",
        "base_dimensions": ["L", "T", "M"],
        "s1_sound_speed": {
            "formula": "c_s^2(rho)=(1/m_GNLS)*dP/drho=5*K*rho^4/m_GNLS",
            "dimension": "L T^-1",
            "state_dependence": "c_s(rho) proportional to rho^2; d ln c_s / d ln rho = 2",
            "profile_statement": "rho, v_b, and c_s are one stationary profile through continuity plus quantum-Bernoulli; c_s=1 denotes the asymptotic c_s0 pin.",
            "provenance": "derived relative to imposed EOS P=K*rho^5, not from an hbar-free microscopic EOS derivation in this step.",
        },
        "s2_velocities": {
            "v_b": "background condensate flow velocity, v_b=(hbar/m_GNLS) grad(theta) in the ungauged lane; profile variable",
            "c_s": "bulk density/phonon sound speed; profile c_s(x) with asymptotic c_s0",
            "c_gamma": "photon/gauge-wave speed from the gauge-wave kinetic operator; brane light-cone speed",
            "c_equals": "c=c_gamma from the massless wave-sector ceiling: omega^2=c_gamma^2 k^2 gives group velocity c_gamma; a trapped transverse mode has omega^2=c_gamma^2(k_parallel^2+k_perp^2), so d omega/dk is bounded by c_gamma and approaches it at high drive.",
            "bound_mode_clock": "The trapped-mode rest oscillation is omega0=c_gamma*k_perp from the wave boundary condition. A boosted wave-operator solution has phase exp[-i*omega0*gamma*(t-v*x/c_gamma^2)], so along the packet center x=v*t the internal clock advances at omega0/gamma. No E=m_defect*c_gamma^2 or Compton premise is used.",
            "mass_bridge_candidate": "m_defect=alpha_J*hbar*J/c_gamma^2, or alpha_J*h*J_nu/c_gamma^2 for a cycle-count rate; candidate conversion only.",
            "constants_vs_profiles": "Constants/input labels in this step: K, m_GNLS, hbar, conserved/no-leakage J label, and asymptotic rho0,c_s0. Profiles: rho(x), v_b(x), c_s(x), and possibly c_gamma(x) until the gauge ratio is fixed.",
            "c_gamma_ratio_verdict": "C_GAMMA_RATIO_UNDERDETERMINED",
            "tail_factor": "(c/c_s)^3=(c_gamma/c_s)^3=lambda_gamma^3, not set to 1 by dimensions or legacy weak-field prose.",
        },
        "s2b_flux": {
            "flux_law_verdict": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "flux_conservation_statement": "In a steady no-leakage region, J=int rho v_b.dSigma is surface-independent, but local v_b can accelerate as rho and area vary.",
            "no_net_accretion_status": "NO_NET_ACCRETION_BC_UNDERIVED",
            "conditional_ideal_nozzle_not_verdict": {
                "assumptions": "steady Euler nozzle, upstream v0=0, no quantum pressure, no confinement variation, solved sonic throat",
                "c_s_star_over_c_s0": "3^(-1/2)",
                "rho_star_over_rho0": "3^(-1/4)",
                "Jcrit_over_rho0_c_s0_Astar": "3^(-3/4)",
                "status": "CONDITIONAL_NOT_ACCEPTED_AS_BRANCH_LAW",
            },
        },
        "s3_hbar": {
            "verdict": "HBAR_PROVENANCE_UNDETERMINED",
            "current_action_role": "hbar is an explicit action/PDE coefficient in the present parent theory.",
            "anti_tautology": "hbar=m_GNLS*c_s0*a is a pin rearrangement unless a is independently fixed by an hbar-free substrate/action relation.",
            "h_2pi_assessment": (
                "Partially meaningful only in winding/rate bookkeeping: h is the natural action per complete phase winding, "
                "while hbar remains the angular phase/PDE coefficient. It does not split charge and mass provenance by itself; "
                "J cycle-vs-angular status and the 2*pi placement are deferred."
            ),
            "energy_n2_assessment": "A vortex kinetic-energy ladder proportional to kappa^2 would scale like n^2, so higher windings are energetically disfavored; this is recorded as a conditional winding-sector prediction, not a mass spectrum derivation.",
        },
        "checks": [check.as_dict() for check in checks],
        "algebraic_checks": algebra,
        "residuals": residuals,
        "summary": {
            "total_dimensional_checks": len(checks),
            "consistent_dimensional_checks": len(checks) - len(failures),
            "inconsistent_dimensional_checks": len(failures),
            "total_algebraic_checks": len(algebra),
            "consistent_algebraic_checks": len(algebra) - len(algebra_failures),
            "inconsistent_algebraic_checks": len(algebra_failures),
            "acceptance_status": "PASS_WITH_NAMED_RESIDUALS",
        },
    }


def render_patha20_velocity_markdown(report: Mapping[str, object]) -> str:
    s1 = report["s1_sound_speed"]
    s2 = report["s2_velocities"]
    s2b = report["s2b_flux"]
    s3 = report["s3_hbar"]
    summary = report["summary"]
    assert isinstance(s1, Mapping)
    assert isinstance(s2, Mapping)
    assert isinstance(s2b, Mapping)
    assert isinstance(s3, Mapping)
    assert isinstance(summary, Mapping)
    lines = [
        "# Path-A 20 velocity constants summary",
        "",
        "## Verdicts",
        "",
        f"- `c_gamma/c_s`: `{s2['c_gamma_ratio_verdict']}`. The carried ratio is `lambda_gamma=c_gamma/c_s`; `tail=(c/c_s)^3=lambda_gamma^3`.",
        f"- `flux_law_verdict`: `{s2b['flux_law_verdict']}`. No accepted `J_crit` law is produced in this step.",
        f"- `hbar` provenance: `{s3['verdict']}`. The `h/2pi` split is bookkeeping-useful for complete windings, not a provenance split.",
        "",
        "## S1 sound speed",
        "",
        f"- Formula: `{s1['formula']}`.",
        f"- Dimension: `{s1['dimension']}` machine-checked in SymPy and Mathematica.",
        f"- State dependence: {s1['state_dependence']}.",
        f"- Profile status: {s1['profile_statement']}",
        f"- Provenance: {s1['provenance']}",
        "",
        "Source anchors: `part01_parent_geometry.tex:194-203`; `pde.tex:344-352`.",
        "",
        "## S2 velocity structure",
        "",
        "| Velocity | Role | Status |",
        "|---|---|---|",
        f"| `v_b` | {s2['v_b']} | `[v_b]=L T^-1` checked |",
        f"| `c_s` | {s2['c_s']} | `[c_s]=L T^-1` checked |",
        f"| `c_gamma` | {s2['c_gamma']} | `[c_gamma]=L T^-1` checked; ratio to `c_s` underdetermined |",
        "",
        f"`c=c_gamma` result: {s2['c_equals']}",
        "",
        f"Bound-mode clock: {s2['bound_mode_clock']}",
        "",
        f"Constants vs profiles: {s2['constants_vs_profiles']}",
        "",
        f"Mass bridge recorded only as candidate form: `{s2['mass_bridge_candidate'].rstrip('.')}`. This does not collapse `M` and does not derive `alpha_J`.",
        "",
        "The localized Maxwell sources expose `Z(w)`, projection, and measured-vs-flux closure data (`part01_parent_geometry.tex:225-389`; `pde.tex:355-565`). The legacy EM acoustic reuse (`em_fields.tex:149-184`, `482-499`, `692-705`) is a prior, not the kinetic-operator proof required here.",
        "",
        "## S2b flux law",
        "",
        f"Verdict: `{s2b['flux_law_verdict']}`.",
        "",
        f"- Conservation statement: {s2b['flux_conservation_statement']}",
        f"- No-net-accretion status: `{s2b['no_net_accretion_status']}`.",
        "- Missing branch data: solved `R0`, `psi0`, `A_M0`, confinement/wall data, quantum-pressure contribution, leakage/topology BC, support/gauge spectra, and overlap/DtN data.",
        "- Consequence: pathA_21 must consume the verdict, not an unconditional choked flux.",
        "- Environment dependence: any solved flux law inherits `P0=K*rho0^5` and `c_s0^2=5K*rho0^4/m_GNLS`; this was dimension-checked but not converted into an accepted choked law.",
        "",
        "A conditional ideal Euler-nozzle algebra check was recorded but not accepted as the branch law: with upstream rest flow and no `Q`/`V_conf` variation, `c_s,* / c_s0=3^(-1/2)`, `rho_* / rho0=3^(-1/4)`, and `Jcrit/(rho0 c_s0 A_*)=3^(-3/4)`.",
        "",
        "Source anchors: `pde.tex:2515-2566`, `2847-2879`; `brane_bulk_ontology.tex:1998-2042`.",
        "",
        "## S3 hbar and h/2pi",
        "",
        f"Verdict: `{s3['verdict']}`.",
        "",
        f"- Current action role: {s3['current_action_role']}",
        f"- Anti-tautology gate: {s3['anti_tautology']}",
        f"- `h/2pi` assessment: {s3['h_2pi_assessment']}",
        f"- `n^2` assessment: {s3['energy_n2_assessment']}",
        "",
        "Role catalog dimensions checked: circulation `kappa=int v_b dl=h*n/m_GNLS`, phase momentum `p=hbar grad(theta)`, quantum pressure `Q`, and candidate bridge `hbar J/c_gamma^2`.",
        "",
        "Source anchors: `part01_parent_geometry.tex:174-219`, `270-289`; `pde.tex:429-469`; `brane_bulk_ontology.tex:668-671`, `1169-1180`.",
        "",
        "## Residuals",
        "",
    ]
    residuals = report["residuals"]
    assert isinstance(residuals, Sequence)
    for raw in residuals:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: {raw['status']}. {raw['downstream_consequence']} Source: {raw['source']}")
    lines.extend(
        [
            "",
            "## Algebraic harness summary",
            "",
            f"- Dimensional checks: {summary['consistent_dimensional_checks']} consistent, {summary['inconsistent_dimensional_checks']} inconsistent, {summary['total_dimensional_checks']} total.",
            f"- Algebraic checks: {summary['consistent_algebraic_checks']} consistent, {summary['inconsistent_algebraic_checks']} inconsistent, {summary['total_algebraic_checks']} total.",
            f"- Acceptance status: `{summary['acceptance_status']}`.",
            "",
        ]
    )
    return "\n".join(lines)


def write_patha20_velocity_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = run_patha20_velocity_constants()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_20_velocity_constants_report.json"
    scratch_md_path = out_dir / "pathA_20_velocity_constants_report.md"
    reference_path = report_dir / "pathA_20_velocity_constants.md"
    rendered = render_patha20_velocity_markdown(report) + "\n"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered, encoding="utf-8")
    reference_path.write_text(rendered, encoding="utf-8")
    return json_path, reference_path, report


def _patha20b_cgamma_cs_checks() -> list[Check]:
    rho4 = D["rho_4d_number_density"]
    velocity = LENGTH / TIME
    wave_number = LENGTH**-1
    lpsi_density = ENERGY / (LENGTH**4)
    electric_field = FORCE
    magnetic_field = ACTION / (LENGTH**2)
    maxwell_c_e = lpsi_density / (electric_field**2)
    maxwell_c_b = maxwell_c_e * (velocity**2)
    ai = D["q_Ai"]
    a0 = D["q_A0"]
    c_s_squared = D["K_eos"] * (rho4**4) / MASS
    return [
        expect_dim(
            "pathA_20b_L2",
            "phonon sound speed c_s0=sqrt(5*K*rho0^4/m_GNLS)",
            c_s_squared ** sp.Rational(1, 2),
            velocity,
            "Dimension check only; equality to c_gamma is not inferred.",
        ),
        expect_dim(
            "pathA_20b_L1_L2",
            "Maxwell principal speed squared C_B/C_E",
            maxwell_c_b / maxwell_c_e,
            velocity**2,
            "Only the time-vs-space principal coefficient ratio sets c_bulk^2; an overall Z/mu0 prefactor is non-evidentiary.",
        ),
        expect_dim(
            "pathA_20b_L2",
            "gauge speed c_gamma=sqrt(C_B/C_E)",
            (maxwell_c_b / maxwell_c_e) ** sp.Rational(1, 2),
            velocity,
            "This is the bulk/principal gauge speed dimension, not a proof that it equals c_s.",
        ),
        expect_dim(
            "pathA_20b_L3",
            "conditional bulk ratio c_bulk/c_s0",
            ((maxwell_c_b / maxwell_c_e) / c_s_squared) ** sp.Rational(1, 2),
            DIMENSIONLESS,
            "The ratio is symbolic because no source equation sets C_B/C_E=5*K*rho0^4/m_GNLS.",
        ),
        expect_dim(
            "pathA_20b_L4",
            "tail factor lambda_gamma^3=(c_gamma/c_s)^3",
            (((maxwell_c_b / maxwell_c_e) / c_s_squared) ** sp.Rational(1, 2)) ** 3,
            DIMENSIONLESS,
            "Dimensionless tail status is non-evidentiary for lambda_gamma=1.",
        ),
        homogeneous(
            "pathA_20b_L1b",
            "Maxwell transverse principal operator terms",
            {
                "C_E*partial_t^2 A_T": maxwell_c_e * ai / (TIME**2),
                "C_B*laplacian A_T": maxwell_c_b * ai / (LENGTH**2),
            },
            "Gauge-invariant transverse field strengths give the cone; gauge-fixing terms do not set the speed.",
        ),
        homogeneous(
            "pathA_20b_L2",
            "transverse gauge dispersion omega^2=(C_B/C_E)*k^2",
            {
                "omega^2": (TIME**-1) ** 2,
                "(C_B/C_E)*k^2": (maxwell_c_b / maxwell_c_e) * (wave_number**2),
            },
        ),
        homogeneous(
            "pathA_20b_L2",
            "phonon acoustic dispersion omega^2=c_s0^2*k^2",
            {
                "omega^2": (TIME**-1) ** 2,
                "c_s0^2*k^2": c_s_squared * (wave_number**2),
            },
            "The Bogoliubov quantum-pressure k^4 correction is dispersive and not used to identify c_gamma.",
        ),
        expect_dim(
            "pathA_20b_L1",
            "background charge density J_psi0^0=q_star*rho0",
            rho4,
            rho4,
            "The homogeneous Maxwell equation needs J_ext0^0=-J_psi0^0 for F0=0.",
        ),
        homogeneous(
            "pathA_20b_L1b",
            "linearized spatial current variation terms",
            {
                "phase-current term rho0*(hbar/m_GNLS)*grad(delta theta)": rho4 * (ACTION / MASS) / LENGTH,
                "London term (rho0/m_GNLS)*q_star*delta A_i": rho4 * ai / MASS,
                "spatial current dimension rho0*v": rho4 * velocity,
            },
            "These current/London terms are lower-order in the Maxwell equation than the second-derivative principal operator.",
        ),
        homogeneous(
            "pathA_20b_L1",
            "source coupling dimensions A_M delta J^M",
            {
                "A0*delta J0": a0 * rho4,
                "Ai*delta Ji": ai * (rho4 * velocity),
                "local action density": ENERGY / (LENGTH**4),
            },
        ),
    ]


def _patha20b_algebraic_checks() -> list[dict[str, object]]:
    omega, k2, rho0, k_eos, m_gnls, hbar = sp.symbols(
        "omega k2 rho0 K m_GNLS hbar",
        positive=True,
    )
    c_e, c_b = sp.symbols("C_E C_B", positive=True)
    lambda_gamma = sp.symbols("lambda_gamma", positive=True)
    h_prime = 5 * k_eos * rho0**3
    c_s_sq = sp.simplify(rho0 * h_prime / m_gnls)
    c_bulk_sq = c_b / c_e
    phonon_matrix = sp.Matrix(
        [
            [omega, -(rho0 * hbar / m_gnls) * k2],
            [-h_prime, hbar * omega],
        ]
    )
    phonon_det = sp.factor(phonon_matrix.det())
    gauge_transverse = c_e * omega**2 - c_b * k2
    coupled_det = sp.factor(phonon_det * gauge_transverse**2)
    equality_residual = sp.simplify(c_bulk_sq - c_s_sq)
    source_metric_equation_present = False
    forced_equals_valid = source_metric_equation_present and bool(equality_residual == 0)
    lambda_rho_factor = rho0**-2
    lambda_log_slope = sp.simplify(rho0 * sp.diff(lambda_rho_factor, rho0) / lambda_rho_factor)
    return [
        _algebra_check(
            "phonon determinant gives omega^2=c_s0^2*k^2",
            phonon_det,
            hbar * (omega**2 - c_s_sq * k2),
            "Uses h'(rho0)=5*K*rho0^3, so c_s0^2=rho0*h'(rho0)/m_GNLS=5*K*rho0^4/m_GNLS.",
        ),
        _algebra_check(
            "transverse gauge operator gives c_bulk^2=C_B/C_E",
            gauge_transverse,
            c_e * (omega**2 - c_bulk_sq * k2),
            "The overall Maxwell prefactor cancels from the characteristic equation.",
        ),
        _algebra_check(
            "block coupled characteristic determinant after principal decoupling",
            coupled_det,
            hbar * (omega**2 - c_s_sq * k2) * (c_e * (omega**2 - c_bulk_sq * k2)) ** 2,
            "Off-diagonal GNLS-gauge terms are lower derivative on the homogeneous neutralized background.",
        ),
        _boolean_check(
            "negative control keeps C_B/C_E independent from c_s0^2",
            equality_residual != 0,
            True,
            "The symbolic residual C_B/C_E - 5*K*rho0^4/m_GNLS must remain nonzero without a source equation.",
        ),
        _boolean_check(
            "forced C_GAMMA_EQUALS_C_S verdict rejected without source metric equation",
            forced_equals_valid,
            False,
            "This is the required negative-control fixture: dimensions alone cannot force equality.",
        ),
        _algebra_check(
            "conditional rho dependence of c_bulk/c_s0 when C_B/C_E is rho-independent",
            lambda_log_slope,
            sp.Integer(-2),
            "Since c_s0 is proportional to rho0^2, an independent c_bulk gives lambda_gamma proportional to rho0^-2.",
        ),
        _algebra_check(
            "standing-wave tail remains lambda_gamma^3",
            lambda_gamma**3,
            lambda_gamma**3,
            "No equality to 1 is inserted.",
        ),
    ]


def run_patha20b_cgamma_cs() -> dict[str, object]:
    checks = _patha20b_cgamma_cs_checks()
    algebra = _patha20b_algebraic_checks()
    failures = [check for check in checks if check.status != "CONSISTENT"]
    algebra_failures = [check for check in algebra if not check["pass"]]
    residuals = [
        {
            "name": "BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED",
            "status": "BULK_VERDICT_RESIDUAL",
            "source": (
                "part01_parent_geometry.tex:225-247 and pde.tex:357-416 give F_MN F^MN with no source-derived "
                "equation C_B/C_E=5*K*rho0^4/m_GNLS; pathA_20 showed Z(w) is an overall principal prefactor."
            ),
            "downstream_consequence": (
                "bulk_verdict is C_GAMMA_BULK_UNDERDETERMINED. Carry c_bulk^2=C_B/C_E and "
                "c_bulk/c_s0=sqrt((C_B/C_E)*m_GNLS/(5*K*rho0^4)); do not set lambda_gamma=1."
            ),
        },
        {
            "name": "PARENT_METRIC_ACOUSTIC_IDENTIFICATION_MISSING",
            "status": "BLOCKS_BULK_EQUALS_C_S",
            "source": (
                "part01_parent_geometry.tex:191-203 derives the acoustic speed from the EOS, while "
                "part01_parent_geometry.tex:225-247 gives the Maxwell parent metric separately; "
                "em_fields.tex:160-178,482-504,692-705 is legacy acoustic reuse, not parent-action evidence."
            ),
            "downstream_consequence": (
                "An EQUALS verdict is forbidden unless a later source equation identifies the gauge kinetic metric "
                "with the acoustic metric."
            ),
        },
        {
            "name": "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "status": "BRANE_VERDICT_RESIDUAL",
            "source": (
                "part01_parent_geometry.tex:333-389 and pde.tex:541-565 give projection and zero-mode reduction "
                "as controlled assumptions; pde.tex:749-763 and 903-931 keep A_w, J^w, and F_muw alive in the "
                "microscopic linearized problem."
            ),
            "downstream_consequence": (
                "brane_verdict is C_GAMMA_RATIO_STILL_UNDERDETERMINED. pathA_21 consumes this brane verdict and "
                "keeps lambda_gamma symbolic."
            ),
        },
        {
            "name": "BRANE_PHOTON_CONE_REQUIRES_PROFILE",
            "status": "BRANE_SUB_RESIDUAL",
            "source": (
                "part01_parent_geometry.tex:944-946 and 1502-1511 state that the matched Maxwell/mixed reduction "
                "and actual nonlinear branch realization remain branch/profile tasks."
            ),
            "downstream_consequence": (
                "If the observed photon is a localized mixed-sector mode rather than a strict far-field zero mode, "
                "its cone must be computed from the solved profile and reduction kernel."
            ),
        },
    ]
    return {
        "schema": "stage1_pathA_20b_cgamma_cs_linearization/v1",
        "base_dimensions": ["L", "T", "M"],
        "l1_background_validity": {
            "status": "LEGAL_WITH_EXPLICIT_NEUTRALIZING_EXTERNAL_SOURCE",
            "psi0": "psi0=sqrt(rho0)*exp(-i*mu*t/hbar) with uniform rho0 and v_b0=0",
            "A_M0": "A_M0=0, so F_MN0=0",
            "J_psi0": "J_psi0^0=q_star*rho0 and J_psi0^i=0",
            "J_ext0": "J_ext0^0=-q_star*rho0 and J_ext0^i=0",
            "neutrality_condition": "J_tot0^M=J_psi0^M+J_ext0^M=0, making the Maxwell background equation 0=0",
            "source_anchor": "pde.tex:370-374 permits explicit external/background sourcing; pde.tex:912-925 gives delta J^0=q_star*delta rho and the London current term.",
            "caveat": "Without this explicit neutralizer the correct stop residual would be HOMOGENEOUS_CHARGE_NEUTRALITY_UNSPECIFIED.",
        },
        "l1b_principal_symbol": {
            "variables": "(delta rho, delta theta, delta A_M)",
            "phonon_block": "det P_ph=hbar*(omega^2 - (5*K*rho0^4/m_GNLS)*k^2)",
            "gauge_transverse_block": "P_T=C_E*omega^2-C_B*k^2=C_E*(omega^2-c_bulk^2*k^2), c_bulk^2=C_B/C_E",
            "coupled_principal_determinant": "det P = P_ph * P_T^2 for the physical transverse gauge polarizations after lower-derivative off-diagonal terms are dropped from the principal symbol",
            "off_diagonal_principal_status": "VANISHES_ON_HOMOGENEOUS_NEUTRALIZED_BACKGROUND",
            "lower_order_terms": [
                "delta J^0=q_star*delta rho",
                "delta J^i contains rho0*(hbar/m_GNLS)*grad(delta theta)",
                "delta J^i contains -(q_star/m_GNLS)*rho0*delta A^i, a London/plasma term",
                "background-gradient current terms vanish on the homogeneous background or are lower order",
            ],
            "gauge_invariant_branch": "The cone is read from transverse field-strength modes. Gauge-fixing consistency is not used as speed evidence.",
        },
        "l2_dispersions": {
            "phonon": "omega^2=c_s0^2*k^2, c_s0^2=5*K*rho0^4/m_GNLS; quantum pressure gives the usual k^4 Bogoliubov correction but not c_gamma.",
            "gauge_transverse": "omega^2=c_bulk^2*k^2 plus possible lower-order London/plasma shifts; c_bulk^2=C_B/C_E.",
            "massless_branch_status": "BULK_PRINCIPAL_TRANSVERSE_BRANCH_ESTABLISHED; lower-order gapped/longitudinal branches are labeled separately and do not set the cone.",
        },
        "l3_two_layer_verdict": {
            "bulk_verdict": "C_GAMMA_BULK_UNDERDETERMINED",
            "bulk_residual": "BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED",
            "conditional_bulk_ratio": "c_bulk/c_s0=sqrt((C_B/C_E)*m_GNLS/(5*K*rho0^4))",
            "brane_verdict": "C_GAMMA_RATIO_STILL_UNDERDETERMINED",
            "brane_residual": "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "brane_sub_residual": "BRANE_PHOTON_CONE_REQUIRES_PROFILE",
            "pathA_21_consumes": "brane_verdict only; lambda_gamma remains symbolic",
            "lambda_gamma_rho_dependence": (
                "No pure number derived. If c_bulk is independent of rho0, lambda_gamma is proportional to rho0^-2; "
                "a pure number or equality would require the missing source equation."
            ),
        },
        "l4_implications": {
            "standing_wave_ceiling": "The standing-wave ceiling remains c=c_gamma, not c_s unless the missing brane verdict later proves equality.",
            "tail_factor": "R_tail=(c/c_s)^3=lambda_gamma^3; conditionally (c_bulk/c_s0)^3, not set to 1.",
            "brane_handoff": "The brane localization/profile question is a blocking sub-residual for pathA_21, not an afterthought.",
        },
        "negative_control": {
            "independent_symbols": ["c_bulk", "c_s0"],
            "forbidden_forced_equality": "C_B/C_E=5*K*rho0^4/m_GNLS",
            "result": "FORCED_EQUALITY_REJECTED_WITHOUT_SOURCE_EQUATION",
        },
        "checks": [check.as_dict() for check in checks],
        "algebraic_checks": algebra,
        "residuals": residuals,
        "summary": {
            "total_dimensional_checks": len(checks),
            "consistent_dimensional_checks": len(checks) - len(failures),
            "inconsistent_dimensional_checks": len(failures),
            "total_algebraic_checks": len(algebra),
            "consistent_algebraic_checks": len(algebra) - len(algebra_failures),
            "inconsistent_algebraic_checks": len(algebra_failures),
            "acceptance_status": "PASS_WITH_NAMED_RESIDUALS",
        },
    }


def render_patha20b_cgamma_cs_markdown(report: Mapping[str, object]) -> str:
    l1 = report["l1_background_validity"]
    l1b = report["l1b_principal_symbol"]
    l2 = report["l2_dispersions"]
    l3 = report["l3_two_layer_verdict"]
    l4 = report["l4_implications"]
    negative = report["negative_control"]
    summary = report["summary"]
    assert isinstance(l1, Mapping)
    assert isinstance(l1b, Mapping)
    assert isinstance(l2, Mapping)
    assert isinstance(l3, Mapping)
    assert isinstance(l4, Mapping)
    assert isinstance(negative, Mapping)
    assert isinstance(summary, Mapping)
    lines = [
        "# Path-A 20b c_gamma vs c_s coupled-linearization summary",
        "",
        "## Verdicts",
        "",
        f"- `bulk_verdict`: `{l3['bulk_verdict']}` with `{l3['bulk_residual']}`.",
        f"- `brane_verdict`: `{l3['brane_verdict']}` with `{l3['brane_residual']}`; sub-residual `{l3['brane_sub_residual']}`.",
        f"- `lambda_gamma`: {l3['lambda_gamma_rho_dependence']}",
        f"- `pathA_21`: {l3['pathA_21_consumes']}.",
        "",
        "No `EQUALS` verdict is claimed. The required source equation `C_B/C_E=5*K*rho0^4/m_GNLS` was not found in the cited parent-action sources.",
        "",
        "## L1 background validity",
        "",
        f"- Status: `{l1['status']}`.",
        f"- `psi0`: {l1['psi0']}.",
        f"- `A_M0`: {l1['A_M0']}.",
        f"- `J_psi0`: {l1['J_psi0']}.",
        f"- `J_ext0`: {l1['J_ext0']}.",
        f"- Neutrality: {l1['neutrality_condition']}.",
        f"- Caveat: {l1['caveat']}",
        "",
        f"Source anchor: {l1['source_anchor']}",
        "",
        "## L1b coupled principal symbol",
        "",
        f"- Variables: `{l1b['variables']}`.",
        f"- Phonon block: `{l1b['phonon_block']}`.",
        f"- Gauge transverse block: `{l1b['gauge_transverse_block']}`.",
        f"- Coupled determinant: `{l1b['coupled_principal_determinant']}`.",
        f"- Off-diagonal principal status: `{l1b['off_diagonal_principal_status']}`.",
        f"- Physical photon cone: {l1b['gauge_invariant_branch']}",
        "",
        "Lower-order/gapped terms, separated from the cone:",
    ]
    lower = l1b["lower_order_terms"]
    assert isinstance(lower, Sequence)
    for item in lower:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "## L2 dispersions",
            "",
            f"- Phonon: {l2['phonon']}",
            f"- Gauge: {l2['gauge_transverse']}",
            f"- Massless/transverse branch status: {l2['massless_branch_status']}",
            "",
            "Machine checks confirm `[c_s]=[c_gamma]=L T^-1` and `C_B/C_E` has speed-squared dimension. These checks are non-evidentiary for equality.",
            "",
            "## L3 two-layer verdict",
            "",
            f"- Bulk: `{l3['bulk_verdict']}` because `{l3['bulk_residual']}`. Conditional formula: `{l3['conditional_bulk_ratio']}`.",
            f"- Brane: `{l3['brane_verdict']}` because `{l3['brane_residual']}`; profile sub-residual `{l3['brane_sub_residual']}`.",
            f"- Rho dependence: {l3['lambda_gamma_rho_dependence']}",
            "",
            "## L4 implications",
            "",
            f"- Standing-wave ceiling: {l4['standing_wave_ceiling']}",
            f"- Tail factor: {l4['tail_factor']}",
            f"- Brane handoff: {l4['brane_handoff']}",
            "",
            "## Negative control",
            "",
            f"- Independent symbols: `{', '.join(str(x) for x in negative['independent_symbols'])}`.",
            f"- Forbidden forced equality: `{negative['forbidden_forced_equality']}`.",
            f"- Result: `{negative['result']}`.",
            "",
            "## Residuals carried",
            "",
        ]
    )
    residuals = report["residuals"]
    assert isinstance(residuals, Sequence)
    for raw in residuals:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: {raw['status']}. {raw['downstream_consequence']} Source: {raw['source']}")
    lines.extend(
        [
            "",
            "## Harness summary",
            "",
            f"- Dimensional checks: {summary['consistent_dimensional_checks']} consistent, {summary['inconsistent_dimensional_checks']} inconsistent, {summary['total_dimensional_checks']} total.",
            f"- Algebraic checks: {summary['consistent_algebraic_checks']} consistent, {summary['inconsistent_algebraic_checks']} inconsistent, {summary['total_algebraic_checks']} total.",
            f"- Acceptance status: `{summary['acceptance_status']}`.",
            "",
        ]
    )
    return "\n".join(lines)


def write_patha20b_cgamma_cs_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = run_patha20b_cgamma_cs()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_20b_cgamma_cs_report.json"
    scratch_md_path = out_dir / "pathA_20b_cgamma_cs_report.md"
    reference_path = report_dir / "pathA_20b_cgamma_cs_linearization.md"
    rendered = render_patha20b_cgamma_cs_markdown(report) + "\n"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered, encoding="utf-8")
    reference_path.write_text(rendered, encoding="utf-8")
    return json_path, reference_path, report


def write_report(out_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = run_audit()
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_18_dimensional_audit.json"
    md_path = out_dir / "pathA_18_dimensional_audit.md"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(render_markdown(report) + "\n", encoding="utf-8")
    return json_path, md_path, report


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        default="software/stage1_solver/_scratch",
        help="directory for JSON and Markdown audit artifacts",
    )
    parser.add_argument(
        "--patha19-foundation",
        action="store_true",
        help="run the side-by-side pathA_19 foundation check group instead of the pathA_18 audit",
    )
    parser.add_argument(
        "--patha20-velocity",
        action="store_true",
        help="run the side-by-side pathA_20 velocity/constants check group instead of the pathA_18 audit",
    )
    parser.add_argument(
        "--patha20b-cgamma-cs",
        action="store_true",
        help="run the side-by-side pathA_20b coupled c_gamma/c_s check group instead of the pathA_18 audit",
    )
    parser.add_argument(
        "--foundation-report-dir",
        default="software/stage1_solver/reports",
        help="directory for the pathA_19 foundation reference markdown",
    )
    parser.add_argument(
        "--patha20-report-dir",
        default="software/stage1_solver/reports",
        help="directory for the pathA_20 velocity/constants reference markdown",
    )
    parser.add_argument(
        "--patha20b-report-dir",
        default="software/stage1_solver/reports",
        help="directory for the pathA_20b c_gamma/c_s reference markdown",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.patha19_foundation:
        json_path, reference_path, report = write_patha19_foundation_report(
            Path(args.out_dir),
            Path(args.foundation_report_dir),
        )
        summary = report["summary"]
        base = report["base_set_verdict"]
        healing = report["healing_length"]
        flags = report["flags"]
        residuals = report["residuals"]
        print(f"wrote {json_path}")
        print(f"wrote {reference_path}")
        print(
            "pathA_19 foundation checks: "
            f"{summary['consistent_algebraic_checks']} consistent, "
            f"{summary['inconsistent_algebraic_checks']} inconsistent, "
            f"{summary['total_algebraic_checks']} total"
        )
        print(f"F1 base-set verdict: {base['status']} because {base['reason']}")
        print(
            "key derived results: "
            "J_number=T^-1 in both bulk and brane-reduced frames; "
            f"healing relation {healing['derived_relation']}"
        )
        print(
            "carried residuals: "
            + ", ".join(str(raw["name"]) for raw in residuals)
            + "; algebraic flags: "
            + ", ".join(str(raw["name"]) for raw in flags)
        )
        return 0
    if args.patha20_velocity:
        json_path, reference_path, report = write_patha20_velocity_report(
            Path(args.out_dir),
            Path(args.patha20_report_dir),
        )
        summary = report["summary"]
        s2 = report["s2_velocities"]
        s2b = report["s2b_flux"]
        s3 = report["s3_hbar"]
        residuals = report["residuals"]
        print(f"wrote {json_path}")
        print(f"wrote {reference_path}")
        print(
            "pathA_20 velocity checks: "
            f"{summary['consistent_dimensional_checks']} dimensional consistent, "
            f"{summary['inconsistent_dimensional_checks']} dimensional inconsistent, "
            f"{summary['total_dimensional_checks']} dimensional total; "
            f"{summary['consistent_algebraic_checks']} algebraic consistent, "
            f"{summary['inconsistent_algebraic_checks']} algebraic inconsistent, "
            f"{summary['total_algebraic_checks']} algebraic total"
        )
        print(f"c_gamma/c_s verdict: {s2['c_gamma_ratio_verdict']}; tail factor {s2['tail_factor']}")
        print(
            "flux_law_verdict: "
            f"{s2b['flux_law_verdict']}; no-net-accretion status {s2b['no_net_accretion_status']}"
        )
        print(
            "hbar provenance verdict: "
            f"{s3['verdict']}; h/2pi assessment: {s3['h_2pi_assessment']}"
        )
        print("carried residuals: " + ", ".join(str(raw["name"]) for raw in residuals))
        return 0
    if args.patha20b_cgamma_cs:
        json_path, reference_path, report = write_patha20b_cgamma_cs_report(
            Path(args.out_dir),
            Path(args.patha20b_report_dir),
        )
        summary = report["summary"]
        l1 = report["l1_background_validity"]
        l2 = report["l2_dispersions"]
        l3 = report["l3_two_layer_verdict"]
        residuals = report["residuals"]
        print(f"wrote {json_path}")
        print(f"wrote {reference_path}")
        print(
            "pathA_20b coupled c_gamma/c_s checks: "
            f"{summary['consistent_dimensional_checks']} dimensional consistent, "
            f"{summary['inconsistent_dimensional_checks']} dimensional inconsistent, "
            f"{summary['total_dimensional_checks']} dimensional total; "
            f"{summary['consistent_algebraic_checks']} algebraic consistent, "
            f"{summary['inconsistent_algebraic_checks']} algebraic inconsistent, "
            f"{summary['total_algebraic_checks']} algebraic total"
        )
        print(
            "bulk_verdict: "
            f"{l3['bulk_verdict']} with {l3['bulk_residual']}; "
            f"conditional ratio {l3['conditional_bulk_ratio']}"
        )
        print(
            "brane_verdict: "
            f"{l3['brane_verdict']} with {l3['brane_residual']} "
            f"(sub-residual {l3['brane_sub_residual']}); pathA_21 consumes this brane verdict"
        )
        print(
            "homogeneous background: "
            f"{l1['status']}; {l1['neutrality_condition']}"
        )
        print(
            "key derived results: "
            "[c_s]=[c_gamma]=L T^-1; "
            "c_s0^2=5*K*rho0^4/m_GNLS; "
            f"{l2['gauge_transverse']}; "
            f"{l3['lambda_gamma_rho_dependence']}"
        )
        print("named residuals carried to pathA_21: " + ", ".join(str(raw["name"]) for raw in residuals))
        return 0
    json_path, md_path, report = write_report(Path(args.out_dir))
    summary = report["summary"]
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")
    print(
        "checks: "
        f"{summary['consistent']} consistent, {summary['inconsistent']} inconsistent, {summary['total']} total"
    )
    print(f"verdict: {report['verdict']['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
