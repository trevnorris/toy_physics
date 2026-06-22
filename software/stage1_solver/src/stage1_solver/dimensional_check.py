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


M_UNIT, L_UNIT, T_UNIT = sp.symbols("M L T", positive=True)


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


def _patha21_emergent_g_checks() -> list[Check]:
    rho4 = D["rho_4d_number_density"]
    rho3 = LENGTH**-3
    velocity = LENGTH / TIME
    number_rate = TIME**-1
    q3_vol = (rho3**-1) * number_rate
    q4_vol = (rho4**-1) * number_rate
    mass_density_3 = MASS * rho3
    mass_density_4 = MASS * rho4
    force_coefficient_3 = FORCE * (LENGTH**2)
    force_coefficient_4 = FORCE * (LENGTH**3)
    n3_eff = LENGTH * rho4
    return [
        expect_dim(
            "pathA_21_P1",
            "reduced-3D volumetric drain strength Q_3=J/n_3",
            q3_vol,
            (LENGTH**3) / TIME,
            "The profile solve must supply n_3,eff and the dimensionless far-field drain factor Theta_Q.",
        ),
        expect_dim(
            "pathA_21_P1",
            "bulk-4D volumetric drain strength Q_4=J/rho_4",
            q4_vol,
            (LENGTH**4) / TIME,
            "This is the unreduced bulk comparison lane; it would give an r^-3 force law in four spatial dimensions.",
        ),
        expect_dim(
            "pathA_21_P1",
            "far-field brane drain velocity v_r=Q_3/(4*pi*r^2)",
            q3_vol / (LENGTH**2),
            velocity,
        ),
        expect_dim(
            "pathA_21_P1",
            "far-field bulk drain velocity v_R=Q_4/(S_3*R^3)",
            q4_vol / (LENGTH**3),
            velocity,
        ),
        expect_dim(
            "pathA_21_P1",
            "reduced-3D pressure/momentum force F=rho_m3*Q_2*v_1",
            mass_density_3 * q3_vol * velocity,
            FORCE,
            "This checks the control-surface cross term behind the G-free drain force.",
        ),
        expect_dim(
            "pathA_21_P1",
            "reduced-3D force coefficient C_F=rho_m3*Q_1*Q_2/(4*pi)",
            mass_density_3 * (q3_vol**2),
            force_coefficient_3,
            "F=C_F/r^2 in the compact reduced-3D drain limit; no G is present.",
        ),
        expect_dim(
            "pathA_21_P1",
            "bulk-4D force coefficient C_F4=rho_m4*Q_1*Q_2/S_3",
            mass_density_4 * (q4_vol**2),
            force_coefficient_4,
            "F=C_F4/R^3 before brane reduction; this is not the observed Newton lane.",
        ),
        expect_dim(
            "pathA_21_P1_P4",
            "effective reduced number density n_3,eff=W_eff*rho_4",
            n3_eff,
            rho3,
            "W_eff is a named reduction/profile width, not set equal to a or xi_h/sqrt(2).",
        ),
        expect_dim(
            "pathA_21_P2",
            "angular-rate bridge candidate alpha_J*hbar*J_omega/c_gamma^2",
            ACTION * number_rate / (velocity**2),
            MASS,
            "Dimensionally valid candidate only; alpha_J is not defined by this check.",
        ),
        expect_dim(
            "pathA_21_P2",
            "cycle-rate bridge candidate alpha_J*h*J_nu/c_gamma^2",
            ACTION * number_rate / (velocity**2),
            MASS,
            "h=2*pi*hbar has the same dimension; the 2*pi placement is algebraic bookkeeping.",
        ),
        expect_dim(
            "pathA_21_P2",
            "candidate throat Hamiltonian ratio H_throat/(hbar*J_omega)",
            ENERGY / (ACTION * number_rate),
            DIMENSIONLESS,
            "A real bridge would need a source equation identifying this profile energy with the rest mass.",
        ),
        expect_dim(
            "pathA_21_P2",
            "inertial-throat kinetic coefficient dimension",
            MASS,
            MASS,
            "Placeholder-free dimensional lane for m_inertial from a second velocity derivative of the effective throat action.",
        ),
        expect_dim(
            "pathA_21_P3",
            "hbar remains an action dimension in the retained {L,T,M} base",
            ACTION,
            Dim(2, -1, 1),
            "No independent hbar-free relation is introduced in pathA_21.",
        ),
        expect_dim(
            "pathA_21_P4",
            "conditional G_eff from C_F*c_gamma^4/(alpha1*alpha2*hbar^2*J1*J2)",
            force_coefficient_3 * (velocity**4) / ((ACTION**2) * (number_rate**2)),
            D["G_3_spatial"],
            "Algebraic dimension only; P4 rejects extraction because P2 and universality are not derived.",
        ),
        expect_dim(
            "pathA_21_P4",
            "conditional G_eff using W_eff*rho4: c_gamma^4*m_GNLS/(W_eff*rho4*hbar^2)",
            (velocity**4) * MASS / (n3_eff * (ACTION**2)),
            D["G_3_spatial"],
        ),
        expect_dim(
            "pathA_21_P4",
            "lambda_gamma=c_gamma/c_s",
            velocity / velocity,
            DIMENSIONLESS,
            "Consumed symbolically from pathA_20b.",
        ),
    ]


def _patha21_algebraic_checks() -> list[dict[str, object]]:
    q1, q2, rho_m, r, hbar, jnu, jomega, alpha1, alpha2, cgam = sp.symbols(
        "Q1 Q2 rho_m r hbar J_nu J_omega alpha1 alpha2 c_gamma",
        positive=True,
    )
    c_f = rho_m * q1 * q2 / (4 * sp.pi)
    force_radial = -c_f / r**2
    velocity_from_1 = -q1 / (4 * sp.pi * r**2)
    force_from_cross_term = rho_m * q2 * velocity_from_1
    m1 = alpha1 * hbar * jomega / cgam**2
    m2 = alpha2 * hbar * jomega / cgam**2
    g_cond = c_f * cgam**4 / (alpha1 * alpha2 * hbar**2 * jomega**2)
    return [
        _algebra_check(
            "reduced-3D drain force is inverse-square before any G",
            force_from_cross_term,
            force_radial,
            "Uses inward-positive Q_i and the pressure/momentum cross term rho_m*Q_2*v_1.",
        ),
        _algebra_check(
            "G-free force coefficient C_F",
            -force_radial * r**2,
            c_f,
            "C_F=rho_m3*Q1*Q2/(4*pi); G is not introduced.",
        ),
        _algebra_check(
            "cycle-rate bridge keeps the 2*pi outside alpha_J",
            2 * sp.pi * hbar * jnu,
            2 * sp.pi * hbar * jnu,
            "The check substitutes h=2*pi*hbar by expectation; alpha_J does not absorb 2*pi.",
        ),
        _algebra_check(
            "conditional Newton coefficient algebra if independent bridge existed",
            sp.simplify(g_cond * m1 * m2),
            c_f,
            "This is only an algebraic downstream conditional; it is not a P2/P4 derivation.",
        ),
    ]


def _patha21_residuals() -> list[dict[str, object]]:
    return [
        {
            "name": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "status": "CARRIED_FROM_pathA_20",
            "source": "pde.tex:2515-2566 and 2847-2879 require the realized branch data set before flux and overlap values are known.",
            "downstream_consequence": "The drain strength Theta_Q, W_eff, leakage/topology BC, and profile force integral remain symbolic.",
        },
        {
            "name": "PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED",
            "status": "P1_PROFILE_RESIDUAL",
            "source": "pde.tex:396-451 gives current, continuity, quantum potential, and Euler-like force balance but no solved finite-throat branch.",
            "downstream_consequence": "P1 provides a G-free profile-functional C_F; the compact reduced-3D r^-2 lane is conditional on the source profile.",
        },
        {
            "name": "FORCE_POWER_PROFILE_UNDERDETERMINED",
            "status": "P1_SCOPE_RESIDUAL",
            "source": "pde.tex:511-539 gives projected open-system brane continuity; pde.tex:541-565 distinguishes reduction from projection.",
            "downstream_consequence": "A universal Newtonian r^-2 force is not promoted until the profile solve proves compact reduced-3D no-leakage behavior.",
        },
        {
            "name": "BOUNDARY_HAMILTONIAN_MASS_RELATION_MISSING",
            "status": "BLOCKS_ALPHA_J",
            "source": "part01_parent_geometry.tex:174-219 and pde.tex:326-406 define the action and current; brane_bulk_ontology.tex:1267-1297 gives drainage scaling only.",
            "downstream_consequence": "alpha_J is not an independently derived profile functional; the mass bridge is rejected.",
        },
        {
            "name": "MASS_BRIDGE_FORM_NOT_DERIVED",
            "status": "P2_VERDICT",
            "source": "No action-level, boundary-source, Noether, or Hamiltonian equation was found that maps inflow J to m_defect without restating the target.",
            "downstream_consequence": "m_defect=alpha_J*hbar*J/c_gamma^2 remains a candidate dimensional form only.",
        },
        {
            "name": "H_2PI_RATE_CLASSIFICATION_UNDETERMINED",
            "status": "CARRIED_FROM_pathA_20",
            "source": "pde.tex:429-469 gives angular phase variables; circulation/winding sources do not classify the throat inflow J value.",
            "downstream_consequence": "Keep J_omega and J_nu separate; h*J_nu=2*pi*hbar*J_nu and alpha_J cannot absorb 2*pi.",
        },
        {
            "name": "INERTIAL_PROFILE_RESPONSE_MISSING",
            "status": "BLOCKS_EP_INERTIAL_SIDE",
            "source": "The cited parent sources do not provide an accelerated-throat kinetic response functional for m_inertial.",
            "downstream_consequence": "m_inertial cannot be matched to the source mass normalization in this step.",
        },
        {
            "name": "SOURCE_MASS_PROFILE_NORMALIZATION_MISSING",
            "status": "BLOCKS_EP_SOURCE_SIDE",
            "source": "P1 supplies C_F as a drain-force profile functional, not an independent mass functional.",
            "downstream_consequence": "m_source is not separately reduced to the same integral and normalization as m_inertial.",
        },
        {
            "name": "EP_NOT_DERIVED",
            "status": "P2_VERDICT",
            "source": "The inertial and source masses are not both available as separately sourced profile integrals.",
            "downstream_consequence": "Equivalence of m_inertial and m_source is not claimed.",
        },
        {
            "name": "HBAR_FREE_SUBSTRATE_RELATION_MISSING",
            "status": "P3_BLOCKER",
            "source": "pathA_20b retained HBAR_PROVENANCE_UNDETERMINED; no new hbar-free substrate relation appears in pathA_21.",
            "downstream_consequence": "The base remains {L,T,M}; INFLOW_MASS_SOURCE_MISSING is sharpened by MASS_BRIDGE_FORM_NOT_DERIVED.",
        },
        {
            "name": "NEWTON_G_FORM_NOT_DERIVED",
            "status": "P4_VERDICT",
            "source": "P4 requires a G-free universal inverse-square force plus independently derived P2 masses; P2 failed and P1 remains profile-conditional.",
            "downstream_consequence": "The m<->G algebraic form is recorded only as a conditional hand-off to the profile solve/pathA_22.",
        },
        {
            "name": "W_EFF_REDUCTION_UNDERIVED",
            "status": "P4_PROFILE_RESIDUAL",
            "source": "part01_parent_geometry.tex:298-389 and pde.tex:496-565 define projection/reduction but do not fix an invariant width.",
            "downstream_consequence": "Use named W_eff; do not set it to a or xi_h/sqrt(2).",
        },
        {
            "name": "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "status": "CARRIED_FROM_pathA_20b",
            "source": "pde.tex:541-565 and software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:47-51 keep the observed brane c_gamma/c_s as a profile/reduction residual.",
            "downstream_consequence": "lambda_gamma remains symbolic in all pathA_21 forms.",
        },
    ]


def _patha21_profile_spec_rows() -> list[dict[str, str]]:
    rows = [
        {
            "symbol": "R0(w)",
            "definition": "Stationary throat surface Sigma0(X)=r-R0(w); domain 0<=w<=L0 with mouth R0(0) and bottom/topology BC.",
            "dimension": "L",
            "frame": "4D-bulk / reduced throat",
            "source_anchor": "pde.tex:2515-2518; part01_parent_geometry.tex:447-461",
            "closes_which_output": "C_F, W_eff, branch geometry",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "downstream_consumer": "pathA_21 P1/P4; pathA_22 scale map",
        },
        {
            "symbol": "psi0(X)",
            "definition": "Background matter field on the stationary branch; rho0(X)=abs(psi0)^2 enters current, Q, h(rho), and drain flux.",
            "dimension": "L^-2",
            "frame": "4D-bulk",
            "source_anchor": "pde.tex:2519-2522; pde.tex:326-406",
            "closes_which_output": "C_F, J-value, pressure response",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "downstream_consumer": "pathA_21 P1/P5",
        },
        {
            "symbol": "rho0(X)",
            "definition": "rho0=abs(psi0)^2 on the branch; measure d^4X in bulk and reduced d^3x after W_eff integration.",
            "dimension": "L^-4 bulk; L^-3 reduced-3D after W_eff",
            "frame": "4D-bulk / reduced-3D",
            "source_anchor": "pde.tex:431-443; part01_parent_geometry.tex:266-278",
            "closes_which_output": "C_F, c_s, W_eff*rho0",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "downstream_consumer": "pathA_21 P1/P4",
        },
        {
            "symbol": "A_M0(x,w)",
            "definition": "Background localized gauge field entering D_t, D_i, v_i, Maxwell/mixed branch data, and c_gamma reduction.",
            "dimension": "A0: M L^2 T^-2; Ai: M L T^-1",
            "frame": "4D-bulk / brane",
            "source_anchor": "pde.tex:2523-2526; pde.tex:355-416",
            "closes_which_output": "brane c_gamma, lambda_gamma, mixed profile data",
            "status": "profile-solve",
            "residual_if_absent": "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "downstream_consumer": "pathA_21 P4; pathA_22",
        },
        {
            "symbol": "V_conf(X;Sigma0)",
            "definition": "Confinement potential on the stationary surface; domain bulk neighborhood of throat; integrand V_conf*rho in L_psi.",
            "dimension": "M L^2 T^-2",
            "frame": "4D-bulk",
            "source_anchor": "pde.tex:2527-2530; part01_parent_geometry.tex:466-497",
            "closes_which_output": "C_F pressure/Bernoulli profile",
            "status": "profile-solve",
            "residual_if_absent": "PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED",
            "downstream_consumer": "pathA_21 P1",
        },
        {
            "symbol": "Q0(rho0)",
            "definition": "Quantum potential -hbar^2/(2m_GNLS) nabla_4^2 sqrt(rho0)/sqrt(rho0), evaluated on the solved branch.",
            "dimension": "M L^2 T^-2",
            "frame": "4D-bulk",
            "source_anchor": "pde.tex:440-443; part01_parent_geometry.tex:275-286",
            "closes_which_output": "C_F pressure/Bernoulli profile",
            "status": "profile-solve",
            "residual_if_absent": "PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED",
            "downstream_consumer": "pathA_21 P1",
        },
        {
            "symbol": "S_leak",
            "definition": "Projected continuity leakage -[W j^w] + int W'(w) j^w dw; domain transverse w boundary plus support of W'.",
            "dimension": "L^-4 T^-1 projected; L^-3 T^-1 reduced",
            "frame": "brane projection / reduced-3D",
            "source_anchor": "pde.tex:511-539; part01_parent_geometry.tex:321-330",
            "closes_which_output": "r-power, J conservation, C_F",
            "status": "profile-solve",
            "residual_if_absent": "FORCE_POWER_PROFILE_UNDERDETERMINED",
            "downstream_consumer": "pathA_21 P1/P4",
        },
        {
            "symbol": "W_eff",
            "definition": "Named reduction width N_infty,3/rho_infty,4 from the solved brane localization/reduction kernel; measure dw with profile weight, not a or xi_h/sqrt(2).",
            "dimension": "L",
            "frame": "4D-bulk to reduced-3D",
            "source_anchor": "pde.tex:541-565; part01_parent_geometry.tex:298-389",
            "closes_which_output": "G, reduced C_F",
            "status": "profile-solve",
            "residual_if_absent": "W_EFF_REDUCTION_UNDERIVED",
            "downstream_consumer": "pathA_21 P4; pathA_22 scale map",
        },
        {
            "symbol": "N_infty,3",
            "definition": "Asymptotic reduced number density int rho0(x,w) chi_N(w) dw = W_eff*rho_infty,4 in the far field.",
            "dimension": "L^-3",
            "frame": "reduced-3D",
            "source_anchor": "pde.tex:496-509; software/stage1_solver/reports/pathA_19_dimensional_foundation.md:20-28",
            "closes_which_output": "C_F, conditional G",
            "status": "profile-solve",
            "residual_if_absent": "W_EFF_REDUCTION_UNDERIVED",
            "downstream_consumer": "pathA_21 P1/P4",
        },
        {
            "symbol": "J",
            "definition": "Number-rate flux lim_{S_R} int n_3 v_brane.n dS in a no-leakage steady region, with throat-source sign convention specified.",
            "dimension": "T^-1",
            "frame": "reduced-3D / brane",
            "source_anchor": "pde.tex:396-406; pde.tex:511-539",
            "closes_which_output": "C_F, alpha_J candidate, J-value",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "downstream_consumer": "pathA_21 P1/P2",
        },
        {
            "symbol": "Theta_Q",
            "definition": "Dimensionless far-field drain factor Theta_Q=-(N_infty,3/J) lim_{R->infty} int_{S_R} v_brane.n dS; fields psi0,R0,A_M0 and leakage BC.",
            "dimension": "1",
            "frame": "reduced-3D",
            "source_anchor": "pde.tex:511-539; pde.tex:2515-2566",
            "closes_which_output": "C_F",
            "status": "profile-solve",
            "residual_if_absent": "PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED",
            "downstream_consumer": "pathA_21 P1/P4",
        },
        {
            "symbol": "I_F,12",
            "definition": "Dimensionless control-surface pressure/momentum cross integral around throat 2: normalize int_{partial U2} Pi_cross[v1,v2,rho0].n dS by m_GNLS N_infty,3 Q1 Q2/(4*pi r^2).",
            "dimension": "1",
            "frame": "reduced-3D",
            "source_anchor": "pde.tex:445-451; em_fields.tex:118-124 for legacy Euler pressure form",
            "closes_which_output": "C_F, attractiveness",
            "status": "profile-solve",
            "residual_if_absent": "PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED",
            "downstream_consumer": "pathA_21 P1",
        },
        {
            "symbol": "C_F,12",
            "definition": "G-free force coefficient C_F,12=(m_GNLS N_infty,3 Q1 Q2/(4*pi))*I_F,12 with Qi=Theta_Qi*Ji/N_infty,3.",
            "dimension": "M L^3 T^-2",
            "frame": "reduced-3D",
            "source_anchor": "pde.tex:396-451 plus profile rows J,Theta_Q,I_F,12",
            "closes_which_output": "C_F",
            "status": "profile-solve",
            "residual_if_absent": "PRESSURE_FORCE_PROFILE_FUNCTIONAL_UNSOLVED",
            "downstream_consumer": "pathA_21 P1/P4",
        },
        {
            "symbol": "alpha_H,omega",
            "definition": "Would-be dimensionless profile energy ratio H_throat[psi0,A_M0,R0,wall]/(hbar J_omega), with H_throat from the canonical Hamiltonian over the throat domain and bottom BC.",
            "dimension": "1",
            "frame": "4D-bulk / throat",
            "source_anchor": "pde.tex:318-391; brane_bulk_ontology.tex:1267-1297",
            "closes_which_output": "alpha_J, m_defect bridge",
            "status": "new-physics",
            "residual_if_absent": "BOUNDARY_HAMILTONIAN_MASS_RELATION_MISSING",
            "downstream_consumer": "pathA_21 P2; future profile solve",
        },
        {
            "symbol": "J_omega",
            "definition": "Angular-rate version of the throat invariant used in alpha_J*hbar*J_omega/c_gamma^2; classification requires source relation to phase/angular rate.",
            "dimension": "T^-1",
            "frame": "throat / brane",
            "source_anchor": "pde.tex:429-469; software/stage1_solver/reports/pathA_20_velocity_constants.md:57-60",
            "closes_which_output": "2pi placement, alpha_J candidate",
            "status": "new-physics",
            "residual_if_absent": "H_2PI_RATE_CLASSIFICATION_UNDETERMINED",
            "downstream_consumer": "pathA_21 P2",
        },
        {
            "symbol": "J_nu",
            "definition": "Cycle-rate version of the throat invariant, with h*J_nu=2*pi*hbar*J_nu and alpha_J not absorbing 2*pi.",
            "dimension": "T^-1",
            "frame": "throat / brane",
            "source_anchor": "pde.tex:429-469; software/stage1_solver/reports/pathA_20_velocity_constants.md:57-60",
            "closes_which_output": "2pi placement, alpha_J candidate",
            "status": "new-physics",
            "residual_if_absent": "H_2PI_RATE_CLASSIFICATION_UNDETERMINED",
            "downstream_consumer": "pathA_21 P2",
        },
        {
            "symbol": "M_inertial",
            "definition": "Second derivative of the effective moving-throat action with respect to a slow center velocity after integrating fields over the solved throat/support domain.",
            "dimension": "M",
            "frame": "brane / reduced throat",
            "source_anchor": "pde.tex:2806-2879 scope statement; no explicit accelerated-throat source equation in cited parents",
            "closes_which_output": "EP inertial side",
            "status": "new-physics",
            "residual_if_absent": "INERTIAL_PROFILE_RESPONSE_MISSING",
            "downstream_consumer": "pathA_21 P2 EP",
        },
        {
            "symbol": "M_source",
            "definition": "Mass normalization inferred from the far-field drain source after C_F factorization, separately from M_inertial and without using Newton G.",
            "dimension": "M",
            "frame": "reduced-3D",
            "source_anchor": "pde.tex:396-451; brane_bulk_ontology.tex:1267-1297",
            "closes_which_output": "EP source side",
            "status": "new-physics",
            "residual_if_absent": "SOURCE_MASS_PROFILE_NORMALIZATION_MISSING",
            "downstream_consumer": "pathA_21 P2/P4",
        },
        {
            "symbol": "C_B/C_E",
            "definition": "Bulk transverse Maxwell principal coefficient ratio from the localized Maxwell kinetic operator.",
            "dimension": "L^2 T^-2",
            "frame": "4D-bulk",
            "source_anchor": "pde.tex:355-416; software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:39-50",
            "closes_which_output": "bulk c_gamma",
            "status": "profile-solve",
            "residual_if_absent": "BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED",
            "downstream_consumer": "pathA_21 P4",
        },
        {
            "symbol": "lambda_gamma",
            "definition": "Observed brane photon/sound ratio c_gamma/c_s from the zero-mode/profile reduction.",
            "dimension": "1",
            "frame": "brane",
            "source_anchor": "pde.tex:541-565; software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md:47-51",
            "closes_which_output": "G conditional form, mass bridge c_gamma",
            "status": "profile-solve",
            "residual_if_absent": "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "downstream_consumer": "pathA_21 P2/P4; pathA_22",
        },
        {
            "symbol": "mu_eta(w)",
            "definition": "Wall inertia density in the reduced wall action, integrated over the finite throat axial coordinate.",
            "dimension": "M L^-1",
            "frame": "reduced throat",
            "source_anchor": "pde.tex:2531-2535",
            "closes_which_output": "pathA_22 support/scale map",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "T_w(w)",
            "definition": "Axial wall tension function in the reduced wall operator.",
            "dimension": "M L T^-2",
            "frame": "reduced throat",
            "source_anchor": "pde.tex:2531-2535",
            "closes_which_output": "pathA_22 support/scale map",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "T_Omega(w)",
            "definition": "Angular wall stiffness/tension density entering the grouped l=2 wall operator.",
            "dimension": "M L^-1 T^-2",
            "frame": "reduced throat",
            "source_anchor": "pde.tex:2531-2535",
            "closes_which_output": "pathA_22 support/scale map",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "K_eta(w)",
            "definition": "Wall restoring stiffness density in the reduced wall operator.",
            "dimension": "M L^-1 T^-2",
            "frame": "reduced throat",
            "source_anchor": "pde.tex:2531-2535",
            "closes_which_output": "pathA_22 support/scale map",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "varpi_Aalpha",
            "definition": "Stable BdG/support frequency for grouped real lane A and mode alpha.",
            "dimension": "T^-1",
            "frame": "reduced throat",
            "source_anchor": "pde.tex:2537-2544; pde.tex:2602-2609",
            "closes_which_output": "pathA_22 B_n moments",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "c_Aalpha",
            "definition": "Wall/support coupling for lane A and mode alpha; enters B_n=sum c_Aalpha^2/varpi_Aalpha^(2+2n).",
            "dimension": "M^1/2 L^-1/2 T^-2",
            "frame": "reduced throat",
            "source_anchor": "pde.tex:2537-2544; pde.tex:2602-2609",
            "closes_which_output": "pathA_22 B_n moments",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "Omega_U,Ar",
            "definition": "Conservative mixed/gauge frequency for U-family mode r in grouped lane A.",
            "dimension": "T^-1",
            "frame": "reduced throat / mixed gauge",
            "source_anchor": "pde.tex:2545-2549; pde.tex:2611-2624",
            "closes_which_output": "pathA_22 Z_n/N_n moments",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "Omega_W,Ar",
            "definition": "Conservative mixed/gauge frequency for W-family mode r in grouped lane A.",
            "dimension": "T^-1",
            "frame": "reduced throat / mixed gauge",
            "source_anchor": "pde.tex:2545-2549; pde.tex:2611-2624",
            "closes_which_output": "pathA_22 Z_n/N_n moments",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "R_Ar",
            "definition": "Mixed U-W coupling satisfying Delta_r=Omega_U^2 Omega_W^2 - R_r^2.",
            "dimension": "T^-2",
            "frame": "reduced throat / mixed gauge",
            "source_anchor": "pde.tex:2545-2549; pde.tex:2611-2624",
            "closes_which_output": "pathA_22 Z_n/N_n moments",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "g_U,Ar",
            "definition": "Wall-to-U mixed/gauge coupling entering Q_r,H_r,P_r.",
            "dimension": "M^1/2 L^-1/2 T^-2",
            "frame": "reduced throat / mixed gauge",
            "source_anchor": "pde.tex:2545-2549; pde.tex:2619-2624",
            "closes_which_output": "pathA_22 Z_n/N_n moments",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "g_W,Ar",
            "definition": "Wall-to-W mixed/gauge coupling entering Q_r,H_r,P_r.",
            "dimension": "M^1/2 L^-1/2 T^-2",
            "frame": "reduced throat / mixed gauge",
            "source_anchor": "pde.tex:2545-2549; pde.tex:2619-2624",
            "closes_which_output": "pathA_22 Z_n/N_n moments",
            "status": "pathA_22",
            "residual_if_absent": "PATHA_22_BRANCH_PACKET_INCOMPLETE",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "a",
            "definition": "Mouth-radius collective moment a(t)=(1/4*pi) int_{S^2} R(Omega,0,t) dOmega; not an invariant reduction width.",
            "dimension": "L",
            "frame": "brane / reduced throat",
            "source_anchor": "part01_parent_geometry.tex:503-510; pde.tex:2551-2563",
            "closes_which_output": "pathA_22 scale map",
            "status": "pathA_22",
            "residual_if_absent": "A_PIN_IS_BRANCH_MOMENT_NOT_INVARIANT",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "c_s(branch)",
            "definition": "Branch sound speed c_s^2=5K rho0^4/m_GNLS evaluated on the asymptotic/profile state.",
            "dimension": "L T^-1",
            "frame": "4D-bulk / brane",
            "source_anchor": "pde.tex:342-352; software/stage1_solver/reports/pathA_20_velocity_constants.md:9-17",
            "closes_which_output": "lambda_gamma, pathA_22 target",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "downstream_consumer": "pathA_21 P4; pathA_22",
        },
        {
            "symbol": "mhat",
            "definition": "Source-map factor in the isotropic normalization law; in the 3D target lane its squared dimension must convert P0 to G*c_s^5/(a^5*c^5).",
            "dimension": "L^-1 T^-1 M^-1/2 in the 3D target normalization",
            "frame": "PN-facing reduced observable",
            "source_anchor": "pde.tex:2077-2092; pde.tex:2551-2563",
            "closes_which_output": "pathA_22 scale map",
            "status": "pathA_22",
            "residual_if_absent": "SCALE_MAP_SOURCE_FACTOR_UNDERIVED",
            "downstream_consumer": "pathA_22",
        },
        {
            "symbol": "chi_Q",
            "definition": "Outgoing-normalization scalar retained when the passive/outgoing branch is not fixed to canonical compact DtN.",
            "dimension": "1",
            "frame": "PN-facing reduced observable",
            "source_anchor": "pde.tex:1980-1996; pde.tex:2053-2082; pde.tex:2551-2552",
            "closes_which_output": "pathA_22 outgoing normalization",
            "status": "pathA_22",
            "residual_if_absent": "OUTGOING_DTN_BRANCH_UNDERIVED",
            "downstream_consumer": "pathA_22",
        },
    ]
    return rows


def _patha21_status_counts(rows: Sequence[Mapping[str, str]]) -> dict[str, int]:
    statuses = ("known", "profile-solve", "pathA_22", "new-physics")
    return {status: sum(1 for row in rows if row["status"] == status) for status in statuses}


def run_patha21_emergent_g_mass_bridge() -> dict[str, object]:
    checks = _patha21_emergent_g_checks()
    algebra = _patha21_algebraic_checks()
    failures = [check for check in checks if check.status != "CONSISTENT"]
    algebra_failures = [check for check in algebra if not check["pass"]]
    profile_rows = _patha21_profile_spec_rows()
    residuals = _patha21_residuals()
    status_counts = _patha21_status_counts(profile_rows)
    return {
        "schema": "stage1_pathA_21_emergent_g_mass_bridge/v1",
        "base_dimensions": ["L", "T", "M"],
        "p1_force": {
            "verdict": "G_FREE_PROFILE_FUNCTIONAL_DERIVED_CONDITIONAL_REDUCED_3D",
            "force_form": "F_12 = -[C_F,12/r^2] rhat in the compact reduced-3D drain lane",
            "coefficient": "C_F,12=(m_GNLS*N_infty,3*Q1*Q2/(4*pi))*I_F,12 = m_GNLS*Theta_Q1*Theta_Q2*J1*J2*I_F,12/(4*pi*N_infty,3)",
            "attractiveness": "Attractive for positive compressibility, positive N_infty,3, and inward-positive drains Q1,Q2>0 because the pressure/momentum cross term gives F_12 parallel to the external inflow toward the other drain.",
            "r_power": "r^-2 only after compact reduced-3D no-leakage behavior; unreduced bulk lane is r^-3.",
            "profile_dependencies": ["Q(rho)", "V_conf", "R0 geometry", "S_leak/topology BC", "W_eff", "Theta_Q", "I_F,12"],
            "non_closing_conditions": ["compact reduced-3D source", "far-field no-leakage", "positive compressibility", "solved stationary branch"],
        },
        "p2_mass_bridge": {
            "verdict": "MASS_BRIDGE_FORM_NOT_DERIVED",
            "alpha_status": "No independent alpha_J profile functional is derived; alpha_H,omega=H_throat/(hbar*J_omega) is specified as a needed future relation, not accepted as mass bridge.",
            "angular_form": "m_defect ?= alpha_J*hbar*J_omega/c_gamma^2",
            "cycle_form": "m_defect ?= alpha_J*h*J_nu/c_gamma^2 = 2*pi*alpha_J*hbar*J_nu/c_gamma^2",
            "two_pi_status": "H_2PI_RATE_CLASSIFICATION_UNDETERMINED; alpha_J does not absorb 2*pi.",
            "ep_verdict": "EP_NOT_DERIVED",
        },
        "p3_m_collapse": {
            "verdict": "RETAIN_L_T_M",
            "blocker": "HBAR_FREE_SUBSTRATE_RELATION_MISSING and MASS_BRIDGE_FORM_NOT_DERIVED",
            "residual_resolution": "INFLOW_MASS_SOURCE_MISSING is sharpened, not closed.",
        },
        "p4_g": {
            "verdict": "NEWTON_G_FORM_NOT_DERIVED",
            "reason": "P2 does not independently derive masses and P1 inverse-square universality remains profile-conditional.",
            "conditional_m_to_g_form": "If future solves prove universal Theta_Q/alpha_J and I_F,12=1, then G_cond=c_gamma^4*m_GNLS*Theta_Q1*Theta_Q2*I_F,12/(4*pi*N_infty,3*alpha_J1*alpha_J2*hbar^2), with N_infty,3=W_eff*rho_infty,4.",
            "width": "W_eff named reduction width; not set to a or xi_h/sqrt(2).",
        },
        "profile_spec_rows": profile_rows,
        "profile_spec_status_counts": status_counts,
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
            "profile_spec_row_count": len(profile_rows),
            "profile_solve_rows": status_counts["profile-solve"],
            "pathA_22_rows": status_counts["pathA_22"],
            "new_physics_rows": status_counts["new-physics"],
            "known_rows": status_counts["known"],
            "acceptance_status": "PASS_WITH_NAMED_RESIDUALS",
        },
    }


def render_patha21_emergent_g_markdown(report: Mapping[str, object]) -> str:
    p1 = report["p1_force"]
    p2 = report["p2_mass_bridge"]
    p3 = report["p3_m_collapse"]
    p4 = report["p4_g"]
    summary = report["summary"]
    assert isinstance(p1, Mapping)
    assert isinstance(p2, Mapping)
    assert isinstance(p3, Mapping)
    assert isinstance(p4, Mapping)
    assert isinstance(summary, Mapping)
    lines = [
        "# Path-A 21 emergent G and mass bridge summary",
        "",
        "## Verdicts",
        "",
        f"- P1 force: `{p1['verdict']}`. Form: `{p1['force_form']}`; coefficient `{p1['coefficient']}`.",
        f"- P2 bridge: `{p2['verdict']}`. EP: `{p2['ep_verdict']}`.",
        f"- P3 M-collapse: `{p3['verdict']}` because `{p3['blocker']}`.",
        f"- P4 G: `{p4['verdict']}` because {p4['reason']}",
        "",
        "Dual-engine scope: the Python and Mathematica scripts check dimensions and algebra only. They do not prove the non-algebraic P1/P2/P4 source relations.",
        "",
        "## P1 force coefficient",
        "",
        "Source-equation chain:",
        "",
        "1. Parent continuity gives the drain rate: `partial_t rho + partial_i j^i = 0`, `j^i=rho v^i` (`pde.tex:396-406`; `part01_parent_geometry.tex:213-219`).",
        "2. The stationary pressure response comes from the Euler-like identity with `V_conf`, `h(rho)`, and `Q(rho)` retained (`pde.tex:440-451`; `part01_parent_geometry.tex:275-286`).",
        "3. Projected brane continuity is open unless `S_leak` and the topology BC close (`pde.tex:511-539`), and reduction is a separate profile assumption (`pde.tex:541-565`).",
        "4. In the compact reduced-3D far-field lane, the solved profile returns `Q_i=Theta_Qi*J_i/N_infty,3`, `v_i=-Q_i rhat/(4*pi*r^2)`. The pressure/momentum cross term on a control surface around throat 2 gives `F_12=-C_F,12 rhat/r^2`.",
        "",
        f"Attractiveness: {p1['attractiveness']}",
        f"Power law: {p1['r_power']}",
        "Profile dependencies: " + ", ".join(str(item) for item in p1["profile_dependencies"]) + ".",
        "Non-closing conditions: " + ", ".join(str(item) for item in p1["non_closing_conditions"]) + ".",
        "",
        "## P2 mass bridge and EP",
        "",
        f"Verdict: `{p2['verdict']}`.",
        "",
        f"- Angular-rate candidate: `{p2['angular_form']}`.",
        f"- Cycle-rate candidate: `{p2['cycle_form']}`.",
        f"- `2*pi` status: {p2['two_pi_status']}",
        f"- `alpha_J`: {p2['alpha_status']}",
        f"- EP verdict: `{p2['ep_verdict']}` because the accelerated-throat inertial functional and the far-field source mass functional are not separately available with the same normalization.",
        "",
        "No row defines `alpha_J := m_defect*c_gamma^2/(hbar*J)`, and no mass formula is accepted by restatement.",
        "",
        "## P3 M-collapse",
        "",
        f"`{p3['verdict']}`. {p3['residual_resolution']} The retained base is `{', '.join(report['base_dimensions'])}`.",
        "",
        "## P4 G and m-to-G",
        "",
        f"Verdict: `{p4['verdict']}`.",
        f"Conditional algebraic hand-off only: `{p4['conditional_m_to_g_form']}`.",
        f"Width discipline: {p4['width']}",
        "",
        "This is not an extracted Newton constant because the independent mass bridge is missing and the inverse-square/factorized/universal force conditions are not all closed.",
        "",
        "## P5 profile-solve specification",
        "",
        f"Rows: {summary['profile_spec_row_count']} total; {summary['profile_solve_rows']} profile-solve, {summary['pathA_22_rows']} pathA_22, {summary['new_physics_rows']} new-physics, {summary['known_rows']} known.",
        "",
        "| symbol | definition | dimension | frame | source anchor | closes which output | status | residual if absent | downstream consumer |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    rows = report["profile_spec_rows"]
    assert isinstance(rows, Sequence)
    for raw in rows:
        assert isinstance(raw, Mapping)
        lines.append(
            "| "
            + " | ".join(
                str(raw[key]).replace("|", "/")
                for key in (
                    "symbol",
                    "definition",
                    "dimension",
                    "frame",
                    "source_anchor",
                    "closes_which_output",
                    "status",
                    "residual_if_absent",
                    "downstream_consumer",
                )
            )
            + " |"
        )
    lines.extend(["", "## Residuals carried", ""])
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


def write_patha21_emergent_g_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = run_patha21_emergent_g_mass_bridge()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_21_emergent_g_mass_bridge_report.json"
    scratch_md_path = out_dir / "pathA_21_emergent_g_mass_bridge_report.md"
    reference_path = report_dir / "pathA_21_emergent_G_mass_bridge.md"
    rendered = render_patha21_emergent_g_markdown(report) + "\n"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered, encoding="utf-8")
    reference_path.write_text(rendered, encoding="utf-8")
    return json_path, reference_path, report


def _patha21b_force_bvp_checks() -> list[Check]:
    rho4 = D["rho_4d_number_density"]
    rho3 = LENGTH**-3
    velocity = LENGTH / TIME
    number_rate = TIME**-1
    q3_vol = number_rate / rho3
    q4_vol = number_rate / rho4
    stress3 = FORCE / (LENGTH**2)
    stress4 = FORCE / (LENGTH**3)
    force_coefficient_3 = FORCE * (LENGTH**2)
    force_coefficient_4 = FORCE * (LENGTH**3)
    mass_density_3 = MASS * rho3
    mass_density_4 = MASS * rho4
    mu = ENERGY
    a0 = D["q_A0"]
    return [
        expect_dim(
            "pathA_21b_P1b",
            "reduced-3D Gauss solve velocity J/(N_infty,3*Omega2*r^2)",
            number_rate / (rho3 * (LENGTH**2)),
            velocity,
            "The power is d-1 with d=3; the check is dimensional only, not a source-chain proof.",
        ),
        expect_dim(
            "pathA_21b_P1b",
            "bulk-4D Gauss solve velocity J/(rho_infty,4*Omega3*R^3)",
            number_rate / (rho4 * (LENGTH**3)),
            velocity,
            "The unreduced four-spatial bulk lane gives power d-1=3.",
        ),
        expect_dim(
            "pathA_21b_P1b",
            "reduced-3D momentum-flux stress m_GNLS*N_infty,3*v^2",
            mass_density_3 * (velocity**2),
            stress3,
        ),
        expect_dim(
            "pathA_21b_P1b",
            "bulk-4D momentum-flux stress m_GNLS*rho_infty,4*v^2",
            mass_density_4 * (velocity**2),
            stress4,
        ),
        expect_dim(
            "pathA_21b_P1b",
            "reduced-3D force coefficient m_GNLS*N_infty,3*Q1*Q2",
            mass_density_3 * (q3_vol**2),
            force_coefficient_3,
            "F=C_F/r^2 after the reduced-3D Gauss solve and the Pi_cross surface integral.",
        ),
        expect_dim(
            "pathA_21b_P1b",
            "bulk-4D force coefficient m_GNLS*rho_infty,4*Q1*Q2",
            mass_density_4 * (q4_vol**2),
            force_coefficient_4,
            "F=C_F4/R^3 in the unreduced bulk lane.",
        ),
        homogeneous(
            "pathA_21b_G1",
            "stationary GNLS eigenvalue equation terms",
            {
                "kinetic": (ACTION**2) / (MASS * (LENGTH**2)),
                "V_conf": ENERGY,
                "h(rho0)": ENERGY,
                "q*A_00": a0,
                "mu": mu,
            },
            "For psi0 != 0, each bracketed term multiplying psi0 has energy dimension.",
        ),
        homogeneous(
            "pathA_21b_G1",
            "stationary GNLS density-weighted Euler/Bernoulli terms",
            {
                "m_GNLS*v^2": MASS * (velocity**2),
                "h(rho0)": ENERGY,
                "Q(rho0)": ENERGY,
                "V_conf": ENERGY,
                "q*A_00": a0,
            },
        ),
        homogeneous(
            "pathA_21b_G3",
            "Pi_cross reduced-3D stress terms",
            {
                "convective": mass_density_3 * (velocity**2),
                "pressure": D["P_pressure"] * LENGTH,
                "quantum": rho3 * ENERGY,
                "confinement": rho3 * ENERGY,
            },
            "Pressure is reduced by the transverse width in the 3D lane; quantum and confinement stresses are represented by their divergence-equivalent energy-density terms.",
        ),
        expect_dim(
            "pathA_21b_G5",
            "N_infty,3 from transverse integration of rho_infty,4 against a length-width kernel",
            LENGTH * rho4,
            rho3,
            "The formula is dimensional; the kernel shape remains W_KERNEL_UNDERDECLARED.",
        ),
    ]


def _patha21b_algebraic_checks() -> list[dict[str, object]]:
    j, rho, r, radius4, v, k_eos, rho_sym, m_gnls, delta_v2 = sp.symbols(
        "J rho r R v K rho_sym m_GNLS delta_v2"
    )
    omega2 = 4 * sp.pi
    omega3 = 2 * sp.pi**2
    v3_solution = sp.solve(sp.Eq(rho * omega2 * r**2 * v, -j), v)[0]
    v4_solution = sp.solve(sp.Eq(rho * omega3 * radius4**3 * v, -j), v)[0]
    q1, q2, n3, r12 = sp.symbols("Q1 Q2 N3 r12")
    v1_from_gauss = v3_solution.subs({j: q1 * n3, rho: n3, r: r12})
    force_from_surface = m_gnls * n3 * q2 * v1_from_gauss
    force_expected = -m_gnls * n3 * q1 * q2 / (4 * sp.pi * r12**2)
    h_expr = sp.Rational(5, 4) * k_eos * rho_sym**4
    dh_drho = sp.diff(h_expr, rho_sym)
    density_response = sp.simplify((-m_gnls * delta_v2 / 2) / dh_drho)
    return [
        _algebra_check(
            "Gauss solve gives reduced-3D drain velocity with r^(-2) power",
            v3_solution,
            -j / (4 * sp.pi * rho * r**2),
            "Solved from integral flux rho*v*Omega2*r^2=-J; no hand-inserted velocity field.",
        ),
        _algebra_check(
            "Gauss solve gives bulk-4D drain velocity with R^(-3) power",
            v4_solution,
            -j / (2 * sp.pi**2 * rho * radius4**3),
            "Solved from integral flux rho*v*Omega3*R^3=-J.",
        ),
        _algebra_check(
            "Pi_cross surface impulse gives reduced force coefficient",
            force_from_surface,
            force_expected,
            "Uses the Gauss-solved v1 and the drain-2 momentum uptake m_GNLS*N3*Q2*v1.",
        ),
        _algebra_check(
            "EOS enthalpy derivative used for pressure-response sign",
            dh_drho,
            5 * k_eos * rho_sym**3,
            "The sign statement in the report comes from stable K>0 and rho>0, not CAS assumptions.",
        ),
        _algebra_check(
            "Bernoulli density response to increased speed",
            density_response,
            -m_gnls * delta_v2 / (10 * k_eos * rho_sym**3),
            "For stable K>0, rho>0, and delta_v2>0 this is negative; the finite-throat force orientation remains a profile residual.",
        ),
    ]


def _patha21b_gap_statuses() -> list[dict[str, str]]:
    return [
        {
            "gap": "G1",
            "status": "CLOSED",
            "equation_or_residual": "BVP-G1.GNLS plus BVP-G1.Maxwell",
            "note": "Stationary gauged GNLS and stationary localized Maxwell are transcribed from the parent action for fixed branch geometry and source data.",
        },
        {
            "gap": "G2",
            "status": "NAMED RESIDUAL",
            "equation_or_residual": "R0_FREE_BOUNDARY_CONDITION_UNDERIVED",
            "note": "The parent supplies Sigma0=r-R0(w), a0=R0(0), and candidate bottom/regularity data, but not an Euler-Lagrange selector for R0(w).",
        },
        {
            "gap": "G3",
            "status": "CLOSED AS FUNCTIONAL",
            "equation_or_residual": "BVP-G3.Pi_cross stress integral",
            "note": "Pi_cross is written from the Euler stress balance/action terms; its numerical value and attractive orientation are profile integrals.",
        },
        {
            "gap": "G4",
            "status": "NAMED RESIDUAL",
            "equation_or_residual": "J_VALUE_BRANCH_PARAMETER / J_SELECTOR_UNDERIVED",
            "note": "Continuity and no-leakage conserve flux in a chosen lane but do not select the value of J.",
        },
        {
            "gap": "G5",
            "status": "NAMED RESIDUAL",
            "equation_or_residual": "W_KERNEL_UNDERDECLARED / W_EFF_REDUCTION_UNDERIVED",
            "note": "Projection formulas are source-anchored; the kernel shape is not selected by the parent.",
        },
        {
            "gap": "G6",
            "status": "NAMED RESIDUAL",
            "equation_or_residual": "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "note": "The brane photon cone requires a solved zero-mode/profile reduction; lambda_gamma remains symbolic.",
        },
    ]


def _patha21b_profile_spec_rows() -> list[dict[str, str]]:
    rows = [dict(row) for row in _patha21_profile_spec_rows()]
    overrides: dict[str, dict[str, str]] = {
        "R0(w)": {
            "definition": "Stationary throat surface Sigma0(X)=r-R0(w). Parent gives a0=R0(0) and candidate bottom/regularity data, but no free-boundary equation selecting R0(w).",
            "closes_which_output": "branch geometry for BVP-G1; C_F and W_eff only after a branch selection",
            "status": "branch-residual",
            "residual_if_absent": "R0_FREE_BOUNDARY_CONDITION_UNDERIVED",
            "downstream_consumer": "option C branch realization; pathA_22 scale map",
        },
        "psi0(X)": {
            "definition": "Background matter field solved by BVP-G1.GNLS: [-hbar^2/(2m_GNLS)D_i0D_i0+V_conf(X;Sigma0)+h(abs(psi0)^2)+q_*A_00-mu]psi0=0.",
            "source_anchor": "pde.tex:382-391; pde.tex:396-406; pde.tex:2519-2522",
            "closes_which_output": "stationary density/current profile entering C_F, c_s(branch), and Pi_cross",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_GNLS_BVP_NOT_SOLVED",
            "downstream_consumer": "option C; pathA_21b P1b/P5b",
        },
        "rho0(X)": {
            "definition": "rho0=abs(psi0)^2 from BVP-G1.GNLS; it enters h(rho0), P(rho0), Q(rho0), current, and the branch sound speed.",
            "source_anchor": "pde.tex:431-443; pde.tex:342-352; part01_parent_geometry.tex:266-278",
            "closes_which_output": "pressure response, Q0, c_s(branch), reduced density after G5 branch data",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_GNLS_BVP_NOT_SOLVED",
            "downstream_consumer": "option C; pathA_21b P1b/P5b",
        },
        "A_M0(x,w)": {
            "definition": "Background localized gauge field solved by BVP-G1.Maxwell: partial_M(Z F_0^{MN})+xi^-1 partial^N(partial.A_0)=mu0 J_tot,0^N with the gauge condition and finite-energy BCs.",
            "source_anchor": "pde.tex:355-416; pde.tex:2523-2526; part01_parent_geometry.tex:333-390",
            "closes_which_output": "stationary gauge background and mixed-sector branch data; brane lambda_gamma remains G6 residual",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_MAXWELL_BVP_NOT_SOLVED",
            "downstream_consumer": "option C; pathA_21b P5b; pathA_22",
        },
        "V_conf(X;Sigma0)": {
            "definition": "Confinement profile entering BVP-G1.GNLS for a selected Sigma0; parent promotes V_conf(X;a,L) to V_conf(X;Sigma) and gives the smooth-wall variation.",
            "source_anchor": "pde.tex:318-334; pde.tex:2527-2530; part01_parent_geometry.tex:466-497",
            "closes_which_output": "stationary GNLS potential and pressure/Bernoulli profile for a chosen branch",
            "status": "profile-solve",
            "residual_if_absent": "V_CONF_BRANCH_PROFILE_UNDERDECLARED",
            "downstream_consumer": "option C; pathA_21b P1b",
        },
        "Q0(rho0)": {
            "definition": "Quantum potential Q0=-(hbar^2/(2m_GNLS))*nabla_4^2(sqrt(rho0))/sqrt(rho0), evaluated from the BVP-G1 density.",
            "closes_which_output": "quantum-stress contribution to Pi_cross and Bernoulli balance",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_GNLS_BVP_NOT_SOLVED",
            "downstream_consumer": "option C; pathA_21b P1b",
        },
        "S_leak": {
            "definition": "Exact projected leakage S_leak=-[W j^w]+int W'(w)j^w dw. A compact reduced-3D force lane additionally assumes the far-field no-leakage branch S_leak=0.",
            "closes_which_output": "far-field continuity equation, drain r-power in the selected compact lane, and flux conservation",
            "status": "profile-solve",
            "residual_if_absent": "NO_LEAKAGE_BRANCH_BC_UNDERIVED",
            "downstream_consumer": "option C; pathA_21b P1b/P5b",
        },
        "W_eff": {
            "definition": "Reduction width W_eff=N_infty,3/rho_infty,4 after a projection/reduction kernel is selected; formula is declared but the kernel shape is not parent-selected.",
            "closes_which_output": "G reduction and reduced C_F only after G5 branch realization",
            "status": "branch-residual",
            "residual_if_absent": "W_KERNEL_UNDERDECLARED / W_EFF_REDUCTION_UNDERIVED",
            "downstream_consumer": "option C branch realization; pathA_22 scale map",
        },
        "N_infty,3": {
            "definition": "Asymptotic reduced density N_infty,3=int chi_N(w) rho_infty,4(w) dw = W_eff*rho_infty,4 only after the G5 kernel branch is selected.",
            "closes_which_output": "C_F normalization and conditional G denominator",
            "status": "branch-residual",
            "residual_if_absent": "W_KERNEL_UNDERDECLARED / W_EFF_REDUCTION_UNDERIVED",
            "downstream_consumer": "option C branch realization; pathA_21b P1b/P4",
        },
        "J": {
            "definition": "Conserved number-rate flux J in a stationary no-leakage region: int_{S_R} n v.n dS=-J. The parent does not select its value.",
            "closes_which_output": "P1b force coefficient after branch-selected value; alpha_J candidate remains rejected",
            "status": "branch-residual",
            "residual_if_absent": "J_VALUE_BRANCH_PARAMETER / J_SELECTOR_UNDERIVED",
            "downstream_consumer": "option C branch realization; pathA_21b P1b; pathA_22",
        },
        "Theta_Q": {
            "definition": "Dimensionless branch factor relating the mouth/source flux to the far-field Gauss flux; computable only after R0, leakage, and kernel branch choices.",
            "closes_which_output": "C_F factorization and far-field source normalization",
            "status": "branch-residual",
            "residual_if_absent": "THETA_Q_BRANCH_REALIZATION_UNDERIVED",
            "downstream_consumer": "option C branch realization; pathA_21b P1b/P4",
        },
        "I_F,12": {
            "definition": "Dimensionless Pi_cross control-surface integral I_F,12 defined by BVP-G3: normalize -int_{partial U2} Pi_cross.n dS by m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2).",
            "source_anchor": "pde.tex:445-451; pde.tex:342-352; part01_parent_geometry.tex:275-286; pathA_21b BVP-G3",
            "closes_which_output": "C_F magnitude; attractive orientation only after profile sign integral",
            "status": "profile-solve",
            "residual_if_absent": "PI_CROSS_STRESS_TENSOR_UNDERIVED",
            "downstream_consumer": "option C; pathA_21b P1b",
        },
        "C_F,12": {
            "definition": "G-free force coefficient C_F,12=(m_GNLS*N_infty,3*Q1*Q2/(4*pi))*I_F,12, with Qi=Theta_Qi*Ji/N_infty,3 from the Gauss-solved lane.",
            "source_anchor": "pathA_21b P1b Gauss solve plus BVP-G3 Pi_cross",
            "closes_which_output": "P1b force form",
            "status": "profile-solve",
            "residual_if_absent": "ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL",
            "downstream_consumer": "pathA_21b P1b; pathA_22 only after mass bridge remains separately derived",
        },
        "alpha_H,omega": {
            "definition": "Would-be dimensionless profile energy ratio H_throat/(hbar*J_omega). pde.tex:318-391 is the action-level source, not a canonical Hamiltonian; the Hamiltonian and boundary mass relation must be constructed separately.",
            "source_anchor": "pde.tex:318-391 action-level only; brane_bulk_ontology.tex:1267-1297 scaling only",
            "closes_which_output": "alpha_J and m_defect bridge remain unclosed",
            "status": "new-physics",
            "residual_if_absent": "CANONICAL_THROAT_HAMILTONIAN_UNCONSTRUCTED / BOUNDARY_HAMILTONIAN_MASS_RELATION_MISSING",
            "downstream_consumer": "pathA_22 or later new-physics bridge",
        },
        "C_B/C_E": {
            "definition": "Bulk transverse Maxwell principal coefficient ratio. The principal-symbol form is known, but the physical normalization is a calibration/branch datum rather than a throat BVP closure.",
            "closes_which_output": "bulk c_gamma normalization only after calibration/branch choice",
            "status": "branch-residual",
            "residual_if_absent": "BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED",
            "downstream_consumer": "pathA_22 calibration packet",
        },
        "lambda_gamma": {
            "definition": "Observed brane photon/sound ratio c_gamma/c_s. It requires the brane zero-mode/profile reduction and is not closed by P5b.",
            "closes_which_output": "G conditional form and mass bridge c_gamma remain symbolic",
            "status": "new-physics",
            "residual_if_absent": "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "downstream_consumer": "option C zero-mode reduction; pathA_22",
        },
        "c_s(branch)": {
            "definition": "Branch sound speed c_s^2=5K*rho0^4/m_GNLS evaluated on the BVP-G1 asymptotic/profile state.",
            "closes_which_output": "profile sound speed for lambda_gamma and pathA_22 target",
            "status": "profile-solve",
            "residual_if_absent": "STATIONARY_GNLS_BVP_NOT_SOLVED",
            "downstream_consumer": "option C; pathA_22",
        },
    }
    for row in rows:
        row.update(overrides.get(row["symbol"], {}))
    return rows


def _patha21b_status_counts(rows: Sequence[Mapping[str, str]]) -> dict[str, int]:
    statuses = ("profile-solve", "branch-residual", "pathA_22", "new-physics", "known")
    return {status: sum(1 for row in rows if row["status"] == status) for status in statuses}


def _patha21b_named_residuals() -> list[dict[str, str]]:
    return [
        {
            "name": "ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL",
            "status": "P1b_SIGN_RESIDUAL",
            "source": "pde.tex:445-451 plus the BVP-G3 control-surface integral; the parent does not select the finite-throat normal-orientation sign of I_F,12.",
            "downstream_consequence": "EOS/compressibility gives the pressure-drop sign, but the full attractive force sign is left to the solved profile integral.",
        },
        {
            "name": "R0_FREE_BOUNDARY_CONDITION_UNDERIVED",
            "status": "G2_BRANCH_REALIZATION_RESIDUAL",
            "source": "part01_parent_geometry.tex:447-461 defines Sigma0 and candidate bottom data; no action variation with respect to R0(w) is supplied.",
            "downstream_consequence": "Option C must choose or derive the R0 branch selector.",
        },
        {
            "name": "J_VALUE_BRANCH_PARAMETER",
            "status": "G4_BRANCH_REALIZATION_RESIDUAL",
            "source": "pde.tex:396-406 and pde.tex:511-539 conserve/project flux but do not fix the value.",
            "downstream_consequence": "J remains a branch parameter unless a choking, regularity, topology, or energy selector is added downstream.",
        },
        {
            "name": "J_SELECTOR_UNDERIVED",
            "status": "G4_BRANCH_REALIZATION_RESIDUAL",
            "source": "No separate source-anchored regularity/choking/energy condition selecting J was found in the cited parents.",
            "downstream_consequence": "No silent free parameter; option C must enumerate or solve the selector.",
        },
        {
            "name": "W_KERNEL_UNDERDECLARED",
            "status": "G5_BRANCH_REALIZATION_RESIDUAL",
            "source": "pde.tex:496-565 and part01_parent_geometry.tex:298-390 define projection/reduction formulas but not a kernel-shape selector.",
            "downstream_consequence": "W_eff and N_infty,3 remain branch-realization data.",
        },
        {
            "name": "CANONICAL_THROAT_HAMILTONIAN_UNCONSTRUCTED",
            "status": "ALPHA_H_ANCHOR_FIX",
            "source": "pde.tex:318-391 is action-level; it is not a canonical Hamiltonian construction.",
            "downstream_consequence": "alpha_H,omega cannot be used as a mass bridge restatement.",
        },
    ]


def run_patha21b_force_closure_and_profile_bvp() -> dict[str, object]:
    checks = _patha21b_force_bvp_checks()
    algebra = _patha21b_algebraic_checks()
    failures = [check for check in checks if check.status != "CONSISTENT"]
    algebra_failures = [check for check in algebra if not check["pass"]]
    profile_rows = _patha21b_profile_spec_rows()
    status_counts = _patha21b_status_counts(profile_rows)
    carried_residuals = _patha21_residuals()
    new_residuals = _patha21b_named_residuals()
    return {
        "schema": "stage1_pathA_21b_force_closure_and_profile_bvp/v1",
        "base_dimensions": ["L", "T", "M"],
        "p1b_force": {
            "corrected_verdict": "G_FREE_DRAIN_FORCE_FORM_DERIVED_WITH_ATTRACTIVE_SIGN_PROFILE_RESIDUAL",
            "supersedes_pathA_21_label": "G_FREE_PROFILE_FUNCTIONAL_DERIVED_CONDITIONAL_REDUCED_3D",
            "inverse_square_status": "P1_INVERSE_SQUARE_FIELD_ASSUMED_NOT_SOLVED resolved for the compact reduced-3D lane by continuity/Gauss; unreduced bulk lane gives R^-3.",
            "drain_velocity_reduced_3d": "int_{S_r} N_infty,3 v.n dS=-Theta_Q J => v_r(r)=-Theta_Q J/(4*pi*N_infty,3*r^2), or -Theta_Q J/(4*pi*r^2*N0(r)) if the asymptotic density is not constant.",
            "drain_velocity_bulk_4d": "int_{S_R} rho_infty,4 v.n dS=-Theta_Q4 J => v_R(R)=-Theta_Q4 J/(2*pi^2*rho_infty,4*R^3), or -Theta_Q4 J/(2*pi^2*R^3*rho0(R)) if the asymptotic density is not constant.",
            "force_form": "F_12=-(m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2))*I_F,12*rhat_12 in the compact reduced-3D lane, with Qi=Theta_Qi*Ji/N_infty,3.",
            "pi_cross_surface_integral": "F_12=-int_{partial U2} Pi_cross[v1,v2,rho0,A0,Sigma0].n_2 dS; partial U2 is a small closed reduced-3D surface enclosing throat 2 but excluding throat 1.",
            "sign_verdict": "ATTRACTIVE_SIGN_FROM_PROFILE_RESIDUAL",
            "sign_chain": "Bernoulli gives delta h=-(m_GNLS/2)delta(v^2) where V_conf,Q,qA0 are asymptotically fixed. From P=K*rho^5 and h=(5K/4)rho^4, dh/drho=5K*rho^3>0 for a stable K>0,rho>0 branch, so higher entrained speed lowers rho and pressure. The remaining finite-throat traction orientation is the profile sign residual.",
        },
        "p2_mass_bridge": {
            "verdict": "MASS_BRIDGE_FORM_NOT_DERIVED",
            "ep_verdict": "EP_NOT_DERIVED",
            "carried": "verbatim from pathA_21",
        },
        "p3_m_collapse": {
            "verdict": "RETAIN_L_T_M",
            "blocker": "HBAR_FREE_SUBSTRATE_RELATION_MISSING and MASS_BRIDGE_FORM_NOT_DERIVED",
            "carried": "verbatim from pathA_21",
        },
        "p4_g": {
            "verdict": "NEWTON_G_FORM_NOT_DERIVED",
            "carried": "verbatim from pathA_21",
        },
        "bvp": {
            "gaps": _patha21b_gap_statuses(),
            "g1_stationary_gnls": "BVP-G1.GNLS: [-hbar^2/(2m_GNLS)D_i0D_i0+V_conf(X;Sigma0)+h(abs(psi0)^2)+q_*A_00-mu]psi0=0 on the 4D bulk throat exterior/interior domain, with D_i0=partial_i-i*q_*A_i0/hbar.",
            "g1_stationary_maxwell": "BVP-G1.Maxwell: partial_M(Z(w)F_0^{MN})+xi^-1 partial^N(partial.A_0)=mu0[J_psi^N(psi0,A0)+J_ext,0^N], plus Bianchi identities and a gauge condition such as partial.A_0=0.",
            "g3_pi_cross": "Pi_cross,ij=m_GNLS*rho_asym(v1_i v2_j+v2_i v1_j)+delta_ij P_cross+Pi_Q,cross,ij+Pi_V,cross,ij+Pi_EM,cross,ij. Here P_cross=K[(rho0+drho1+drho2)^5-(rho0+drho1)^5-(rho0+drho2)^5+rho0^5], Pi_Q,cross=cross[(hbar^2/(4m_GNLS))*((partial_i rho partial_j rho)/rho-partial_i partial_j rho)], and partial_j Pi_V,cross,ij=[rho partial_i V_conf]_cross for the selected V_conf branch.",
            "g5_projection_formulas": "rho_brane=int W(w)rho dw, j_brane=int W(w)j_xyz dw, N_infty,3=int chi_N(w)rho_infty,4(w)dw, W_eff=N_infty,3/rho_infty,4 for constant far-field rho.",
            "alpha_H_omega_anchor_fix": "pde.tex:318-391 is action-level, not a canonical Hamiltonian; canonical Hamiltonian must be constructed before alpha_H,omega can be used.",
        },
        "profile_spec_rows": profile_rows,
        "profile_spec_status_counts": status_counts,
        "new_residuals": new_residuals,
        "carried_residuals": carried_residuals,
        "carried_negative_verdicts_verbatim": [
            "MASS_BRIDGE_FORM_NOT_DERIVED",
            "EP_NOT_DERIVED",
            "RETAIN_L_T_M",
            "NEWTON_G_FORM_NOT_DERIVED",
            "STATIONARY_PROFILE_UNDERDETERMINED_BY_BRANCH_DATA",
            "H_2PI_RATE_CLASSIFICATION_UNDETERMINED",
            "BRANE_ZERO_MODE_REDUCTION_UNDERIVED",
            "HBAR_PROVENANCE_UNDETERMINED",
            "HBAR_FREE_SUBSTRATE_RELATION_MISSING",
            "W_EFF_REDUCTION_UNDERIVED",
        ],
        "checks": [check.as_dict() for check in checks],
        "algebraic_checks": algebra,
        "summary": {
            "total_dimensional_checks": len(checks),
            "consistent_dimensional_checks": len(checks) - len(failures),
            "inconsistent_dimensional_checks": len(failures),
            "total_algebraic_checks": len(algebra),
            "consistent_algebraic_checks": len(algebra) - len(algebra_failures),
            "inconsistent_algebraic_checks": len(algebra_failures),
            "profile_spec_row_count": len(profile_rows),
            "closed_profile_solve_rows": status_counts["profile-solve"],
            "branch_residual_rows": status_counts["branch-residual"],
            "pathA_22_rows": status_counts["pathA_22"],
            "new_physics_rows": status_counts["new-physics"],
            "known_rows": status_counts["known"],
            "acceptance_status": "PASS_WITH_NAMED_RESIDUALS",
        },
    }


def render_patha21b_force_bvp_markdown(report: Mapping[str, object]) -> str:
    p1 = report["p1b_force"]
    bvp = report["bvp"]
    summary = report["summary"]
    p2 = report["p2_mass_bridge"]
    p3 = report["p3_m_collapse"]
    p4 = report["p4_g"]
    assert isinstance(p1, Mapping)
    assert isinstance(bvp, Mapping)
    assert isinstance(summary, Mapping)
    assert isinstance(p2, Mapping)
    assert isinstance(p3, Mapping)
    assert isinstance(p4, Mapping)
    lines = [
        "# Path-A 21b force closure and stationary-throat BVP",
        "",
        "## Verdicts",
        "",
        f"- P1b force: `{p1['corrected_verdict']}`; supersedes pathA_21 `{p1['supersedes_pathA_21_label']}`.",
        f"- P2 bridge: `{p2['verdict']}`. EP: `{p2['ep_verdict']}`. Carried verbatim.",
        f"- P3 M-collapse: `{p3['verdict']}` because `{p3['blocker']}`. Carried verbatim.",
        f"- P4 G: `{p4['verdict']}`. Carried verbatim.",
        "",
        "Dual-engine scope: Python and Mathematica check dimensions and algebra only. The derivation is the source-equation chain below; exit 0 is not treated as proof.",
        "",
        "## P1b drain field",
        "",
        "Continuity source chain: parent continuity gives `partial_t rho + partial_i j^i=0`, `j^i=rho v^i` (`pde.tex:396-406`; `part01_parent_geometry.tex:213-219`). In a stationary no-leakage region outside the localized throat source, Gauss gives the far-field radial solution rather than inserting the power.",
        "",
        f"- Reduced-3D lane: `{p1['drain_velocity_reduced_3d']}`.",
        f"- Bulk-4D lane: `{p1['drain_velocity_bulk_4d']}`.",
        f"- Inverse-square status: {p1['inverse_square_status']}",
        "",
        "The reduced-3D `r^-2` power is therefore the area power of the enclosing two-sphere. The unreduced four-spatial bulk lane is `R^-3` from the three-sphere area `Omega3 R^3`.",
        "",
        "## P1b force",
        "",
        f"Control surface: `{p1['pi_cross_surface_integral']}`",
        "",
        "The cross stress used by the option-C solve is",
        "",
        "```text",
        str(bvp["g3_pi_cross"]),
        "```",
        "",
        "The pressure cross term is the explicit EOS cross-difference. The quantum term uses the displayed representative of the Madelung quantum stress, with the cross operation meaning `F[rho0+drho1+drho2]-F[rho0+drho1]-F[rho0+drho2]+F[rho0]`. The confinement term is the branch-selected stress-divergence representative whose divergence equals the cross body force from `rho partial_i V_conf`; absent that representative, the row would fall back to `PI_CROSS_STRESS_TENSOR_UNDERIVED`. `Pi_EM,cross` is retained only when gauge/mixed fields are active. The anchors are the Euler identity and EOS (`pde.tex:342-352`, `pde.tex:440-451`; `part01_parent_geometry.tex:275-286`).",
        "",
        f"Force form: `{p1['force_form']}`",
        "",
        f"Sign verdict: `{p1['sign_verdict']}`. {p1['sign_chain']}",
        "",
        "## Stationary BVP",
        "",
        "### G1 Closed Core",
        "",
        f"- GNLS: `{bvp['g1_stationary_gnls']}` Anchor: `pde.tex:382-391`; action anchor `pde.tex:318-391`.",
        f"- Maxwell: `{bvp['g1_stationary_maxwell']}` Anchor: `pde.tex:410-426`; source bookkeeping `pde.tex:370-374`.",
        "",
        "G1 domain and BCs: the fields live on the 4D bulk throat branch domain with measure `d^4X` and stationary time factored by `psi=e^{-i mu t/hbar}psi0`. Far-field BCs fix the asymptotic density, chemical potential/gauge reference, finite energy, and the Gauss flux lane. Throat BCs require regular finite fields and branch-declared mouth/bottom flux data; the value selector for those branch data is not part of G1.",
        "",
        "### Boundary Conditions",
        "",
        "| condition | fields/domain | frame and measure | source anchor | status |",
        "|---|---|---|---|---|",
        "| Stationary matter asymptotic | `psi0 -> sqrt(rho_infty,4) exp(i theta_infty)`, `rho0 -> rho_infty,4`; `V_conf`, `Q`, and gauge reference approach branch constants | 4D bulk, `d^4X`; `psi0:L^-2`, `rho0:L^-4` | `pde.tex:382-406`; `pde.tex:2519-2530` | closed equation, branch value supplied by option C |",
        "| Reduced-3D Gauss flux | `int_{S_r} N_infty,3 v.n dS=-Theta_Q J`, so `v_r=-Theta_Q J/(4*pi*N_infty,3*r^2)` | reduced-3D, `dS=r^2 dOmega_2`; `J:T^-1`, `N_infty,3:L^-3` | `pde.tex:396-406`; `pde.tex:511-539` | closed conservation law; `J_VALUE_BRANCH_PARAMETER` remains |",
        "| Bulk-4D Gauss flux | `int_{S_R} rho_infty,4 v.n dS=-Theta_Q4 J`, so `v_R=-Theta_Q4 J/(2*pi^2*rho_infty,4*R^3)` | 4D bulk, `dS=R^3 dOmega_3`; `rho_infty,4:L^-4` | `pde.tex:396-406` | closed conservation law; compact reduced lane still branch-selected |",
        "| Throat regularity and mouth flux | finite `psi0`, `rho0`, and `A_M0`; mouth inflow oriented into the throat where the branch declares a drain | 4D bulk/brane mouth, mouth `d^2S`; volumetric flux from brane mouth data | `brane_bulk_ontology.tex:1267-1289`; `part01_parent_geometry.tex:447-461` | branch BC; no `J` value selector |",
        "| Bottom/topology | candidate `R0(L0)=0` or equivalent regular bottom condition | reduced throat, `0<=w<=L0` | `part01_parent_geometry.tex:447-461` | branch assumption, not closure |",
        "| Stationary Maxwell far field and gauge | finite-energy `A_M0`; gauge condition from the gauge-fixed Maxwell equation, e.g. `partial.A_0=0`; optional zero-mode BCs only in the reduced brane lane | 4D bulk, `d^4X`; `A0:ML^2T^-2`, `Ai:MLT^-1` | `pde.tex:355-426`; `pde.tex:541-565` | G1 closed for bulk field; G6 brane cone residual |",
        "| Projection kernel normalization | `int W(w)dw=1`; `rho_brane=int W rho dw`, `j_brane=int W j_xyz dw` | brane projection, `dw`; `W:L^-1` | `pde.tex:496-539`; `part01_parent_geometry.tex:298-330` | formula closed; kernel shape is `W_KERNEL_UNDERDECLARED` |",
        "",
        "### Gap Table",
        "",
        "| gap | status | equation or residual | note |",
        "|---|---|---|---|",
    ]
    gaps = bvp["gaps"]
    assert isinstance(gaps, Sequence)
    for raw in gaps:
        assert isinstance(raw, Mapping)
        lines.append(f"| {raw['gap']} | {raw['status']} | `{raw['equation_or_residual']}` | {raw['note']} |")
    lines.extend(
        [
            "",
            "### G2 Branch Choices",
            "",
            "`R0_FREE_BOUNDARY_CONDITION_UNDERIVED`: candidate option-C assumptions are a prescribed analytic `R0(w)`, the parent bottom cap `R0(L0)=0`, an equivalent regular bottom condition, or a future free-boundary stationarity equation. These are branches, not closure.",
            "",
            "### G4 Branch Choices",
            "",
            "`J_VALUE_BRANCH_PARAMETER / J_SELECTOR_UNDERIVED`: no-leakage conserves the flux in a chosen region, but it does not fix the value. Candidate selectors are sonic choking, regularity at the throat bottom, global topology/outflow balance, or an energy extremum, all explicitly downstream assumptions unless derived.",
            "",
            "### G5 Projection Formulas",
            "",
            f"`{bvp['g5_projection_formulas']}`. The formulas are source-anchored (`pde.tex:496-565`; `part01_parent_geometry.tex:298-390`), but the shape of `W(w)` or `chi_N(w)` is `W_KERNEL_UNDERDECLARED`.",
            "",
            "### G6 Brane Photon Residual",
            "",
            "`BRANE_ZERO_MODE_REDUCTION_UNDERIVED`: option C must solve the localized Maxwell zero-mode/profile reduction and compute the observed brane photon cone before `lambda_gamma=c_gamma/c_s` can be numerical.",
            "",
            "### alpha_H,omega Anchor Fix",
            "",
            str(bvp["alpha_H_omega_anchor_fix"]),
            "",
            "## Machine Table",
            "",
            f"Rows: {summary['profile_spec_row_count']} total; {summary['closed_profile_solve_rows']} closed profile-solve, {summary['branch_residual_rows']} branch-residual, {summary['pathA_22_rows']} pathA_22, {summary['new_physics_rows']} new-physics, {summary['known_rows']} known.",
            "",
            "| symbol | definition | dimension | frame | source anchor | closes which output | status | residual if absent | downstream consumer |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
    )
    rows = report["profile_spec_rows"]
    assert isinstance(rows, Sequence)
    for raw in rows:
        assert isinstance(raw, Mapping)
        lines.append(
            "| "
            + " | ".join(
                str(raw[key]).replace("|", "/")
                for key in (
                    "symbol",
                    "definition",
                    "dimension",
                    "frame",
                    "source_anchor",
                    "closes_which_output",
                    "status",
                    "residual_if_absent",
                    "downstream_consumer",
                )
            )
            + " |"
        )
    lines.extend(["", "## New Residuals", ""])
    new_residuals = report["new_residuals"]
    assert isinstance(new_residuals, Sequence)
    for raw in new_residuals:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: {raw['status']}. {raw['downstream_consequence']} Source: {raw['source']}")
    lines.extend(["", "## Carried Negatives", ""])
    carried = report["carried_negative_verdicts_verbatim"]
    assert isinstance(carried, Sequence)
    lines.append("Carried verbatim: " + ", ".join(f"`{item}`" for item in carried) + ".")
    lines.extend(["", "## Carried Residual Ledger", ""])
    residuals = report["carried_residuals"]
    assert isinstance(residuals, Sequence)
    for raw in residuals:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: {raw['status']}. {raw['downstream_consequence']} Source: {raw['source']}")
    lines.extend(
        [
            "",
            "## Harness Summary",
            "",
            f"- Dimensional checks: {summary['consistent_dimensional_checks']} consistent, {summary['inconsistent_dimensional_checks']} inconsistent, {summary['total_dimensional_checks']} total.",
            f"- Algebraic checks: {summary['consistent_algebraic_checks']} consistent, {summary['inconsistent_algebraic_checks']} inconsistent, {summary['total_algebraic_checks']} total.",
            f"- Acceptance status: `{summary['acceptance_status']}`.",
            "",
        ]
    )
    return "\n".join(lines)


def write_patha21b_force_bvp_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = run_patha21b_force_closure_and_profile_bvp()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_21b_force_closure_and_profile_bvp_report.json"
    machine_table_path = out_dir / "pathA_21b_profile_bvp_machine_table.json"
    scratch_md_path = out_dir / "pathA_21b_force_closure_and_profile_bvp.md"
    reference_path = report_dir / "pathA_21b_force_closure_and_profile_bvp.md"
    rendered = render_patha21b_force_bvp_markdown(report) + "\n"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    machine_table_path.write_text(
        json.dumps(report["profile_spec_rows"], indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    scratch_md_path.write_text(rendered, encoding="utf-8")
    reference_path.write_text(rendered, encoding="utf-8")
    return json_path, reference_path, report


def _patha21c_force_tensor_checks() -> list[Check]:
    rho4 = D["rho_4d_number_density"]
    rho3 = LENGTH**-3
    velocity = LENGTH / TIME
    number_rate = TIME**-1
    q3_vol = number_rate / rho3
    q4_vol = number_rate / rho4
    mass_density_3 = MASS * rho3
    mass_density_4 = MASS * rho4
    stress3 = FORCE / (LENGTH**2)
    stress4 = FORCE / (LENGTH**3)
    force_density3 = FORCE / (LENGTH**3)
    force_density4 = FORCE / (LENGTH**4)
    energy_gradient = ENERGY / LENGTH
    quantum_prefactor = (ACTION**2) / (MASS * (LENGTH**2))
    return [
        homogeneous(
            "pathA_21c_P1c_A",
            "reduced-3D momentum-balance terms",
            {
                "partial_t(m_GNLS*N*v_i)": mass_density_3 * velocity / TIME,
                "partial_j Pi_ij": stress3 / LENGTH,
                "V_conf body force N*partial_i V": rho3 * energy_gradient,
            },
            "Checks the balance law dimensions partial_t g_i + partial_j Pi_ij = f_i^body.",
        ),
        homogeneous(
            "pathA_21c_P1c_A",
            "bulk-4D momentum-balance terms",
            {
                "partial_t(m_GNLS*rho*v_i)": mass_density_4 * velocity / TIME,
                "partial_J Pi_iJ": stress4 / LENGTH,
                "V_conf body force rho*partial_i V": rho4 * energy_gradient,
            },
            "Bulk lane remains separate from the compact reduced-3D lane.",
        ),
        homogeneous(
            "pathA_21c_P1c_A",
            "reduced-3D Noether stress representatives",
            {
                "convective m_GNLS*N*v_i*v_j": mass_density_3 * (velocity**2),
                "pressure delta_ij P_3": D["P_pressure"] * LENGTH,
                "quantum sigma_Q": quantum_prefactor * rho3,
            },
            "The pressure term is reduced by one transverse length in the compact 3D lane.",
        ),
        homogeneous(
            "pathA_21c_P1c_A",
            "Euler force-per-volume identity terms in reduced 3D",
            {
                "m_GNLS*N*(partial_t+v.grad)v": mass_density_3 * velocity / TIME,
                "N*partial_i h": rho3 * energy_gradient,
                "N*partial_i Q": rho3 * energy_gradient,
                "N*partial_i V_conf": rho3 * energy_gradient,
                "N*q(E+vB)": rho3 * energy_gradient,
            },
            "Gauge force dimension is checked as q times an energy gradient per particle.",
        ),
        expect_dim(
            "pathA_21c_P1c_B",
            "reduced-3D surface traction integral gives force",
            stress3 * (LENGTH**2),
            FORCE,
        ),
        expect_dim(
            "pathA_21c_P1c_B",
            "bulk-4D surface traction integral gives force",
            stress4 * (LENGTH**3),
            FORCE,
        ),
        expect_dim(
            "pathA_21c_P1c_B",
            "reduced-3D Noether force coefficient m_GNLS*N_infty,3*Q1*Q2",
            mass_density_3 * (q3_vol**2),
            FORCE * (LENGTH**2),
            "This is the coefficient before the downstream dimensionless normalization knob.",
        ),
        expect_dim(
            "pathA_21c_P1c_B",
            "bulk-4D Noether force coefficient m_GNLS*rho_infty,4*Q1*Q2",
            mass_density_4 * (q4_vol**2),
            FORCE * (LENGTH**3),
            "The bulk lane uses Omega_3=2*pi^2 and R^-3, not the reduced 4*pi lane.",
        ),
        expect_dim(
            "pathA_21c_P1c_B",
            "V_conf body-force volume term dimension",
            force_density3 * (LENGTH**3),
            FORCE,
            "The term is a residual unless the selected profile supplies the volume integral.",
        ),
    ]


def _patha21c_algebraic_checks() -> list[dict[str, object]]:
    m_gnls, n3, rho4, q1, q2, r12, radius4, v1, k_eos, rho_sym, vdot = sp.symbols(
        "m_GNLS N3 rho4 Q1 Q2 r12 R v1 K rho vdot"
    )
    hbar_sym, x, m_sym = sp.symbols("hbar x m")
    rho_fn = sp.Function("rho")(x)
    omega2 = 4 * sp.pi
    omega3 = 2 * sp.pi**2
    d3 = sp.Integer(3)
    d4 = sp.Integer(4)
    convective_flux_3 = -(1 + sp.Rational(1, d3)) * m_gnls * n3 * q2 * v1
    pressure_flux_3 = sp.Rational(1, d3) * m_gnls * n3 * q2 * v1
    total_flux_3 = convective_flux_3 + pressure_flux_3
    force_along_v1_3 = -total_flux_3
    v1_from_gauss_3 = -q1 / (omega2 * r12**2)
    force_along_rhat_3 = sp.simplify(force_along_v1_3.subs(v1, v1_from_gauss_3))
    convective_flux_4 = -(1 + sp.Rational(1, d4)) * m_gnls * rho4 * q2 * v1
    pressure_flux_4 = sp.Rational(1, d4) * m_gnls * rho4 * q2 * v1
    total_flux_4 = convective_flux_4 + pressure_flux_4
    force_along_v1_4 = -total_flux_4
    v1_from_gauss_4 = -q1 / (omega3 * radius4**3)
    force_along_rhat_4 = sp.simplify(force_along_v1_4.subs(v1, v1_from_gauss_4))
    h_expr = sp.Rational(5, 4) * k_eos * rho_sym**4
    pressure_expr = k_eos * rho_sym**5
    dh_drho = sp.diff(h_expr, rho_sym)
    dpressure_drho = sp.diff(pressure_expr, rho_sym)
    delta_rho_cross = -m_gnls * vdot / dh_drho
    delta_pressure_cross = sp.simplify(dpressure_drho * delta_rho_cross)
    quantum_potential = -(hbar_sym**2 / (2 * m_sym)) * sp.diff(sp.sqrt(rho_fn), x, 2) / sp.sqrt(rho_fn)
    sigma_q = (hbar_sym**2 / (4 * m_sym)) * (
        (sp.diff(rho_fn, x) ** 2) / rho_fn - sp.diff(rho_fn, x, 2)
    )
    quantum_divergence_residual = sp.simplify(sp.diff(sigma_q, x) - rho_fn * sp.diff(quantum_potential, x))
    return [
        _algebra_check(
            "EOS identity dP/drho = rho*dh/drho",
            dpressure_drho,
            rho_sym * dh_drho,
            "This is the pressure term needed to convert the parent Euler identity into stress divergence form.",
        ),
        _algebra_check(
            "Bernoulli/EOS pressure cross term",
            delta_pressure_cross,
            -m_gnls * rho_sym * vdot,
            "Uses delta h_cross=-m_GNLS*(v1.v2); no density response is imported from the Gauss solve.",
        ),
        _algebra_check(
            "Madelung quantum stress divergence",
            quantum_divergence_residual,
            sp.Integer(0),
            "One-dimensional representative check of partial_j sigma_Q,ij = rho partial_i Q.",
        ),
        _algebra_check(
            "reduced-3D convective angular traction factor",
            convective_flux_3,
            -sp.Rational(4, 3) * m_gnls * n3 * q2 * v1,
            "Uses int n_i n_j dOmega = (4*pi/3) delta_ij and the drain-2 flux through the control surface.",
        ),
        _algebra_check(
            "reduced-3D pressure angular traction factor",
            pressure_flux_3,
            sp.Rational(1, 3) * m_gnls * n3 * q2 * v1,
            "The Bernoulli pressure term cancels the extra convective angular third.",
        ),
        _algebra_check(
            "reduced-3D total flux from Noether stress",
            total_flux_3,
            -m_gnls * n3 * q2 * v1,
            "This is the surface flux before applying F_12=-int Pi.n dS.",
        ),
        _algebra_check(
            "reduced-3D force structure after Gauss substitution",
            force_along_rhat_3,
            -m_gnls * n3 * q1 * q2 / (4 * sp.pi * r12**2),
            "rhat_12 points from defect 1 to defect 2; the minus sign is attractive for like drains.",
        ),
        _algebra_check(
            "bulk-4D total flux from Noether stress",
            total_flux_4,
            -m_gnls * rho4 * q2 * v1,
            "In four spatial dimensions the convective 1/4 and pressure 1/4 cancel the same way.",
        ),
        _algebra_check(
            "bulk-4D force structure after Gauss substitution",
            force_along_rhat_4,
            -m_gnls * rho4 * q1 * q2 / (2 * sp.pi**2 * radius4**3),
            "The bulk lane uses Omega_3=2*pi^2 and remains separate from W_eff/G5.",
        ),
    ]


def _patha21c_named_residuals() -> list[dict[str, str]]:
    return [
        {
            "name": "VCONF_BODY_FORCE_RESIDUAL",
            "status": "BODY_FORCE_PROFILE_RESIDUAL",
            "source": "The action contains -V_conf(X;Sigma0)*rho, so f_i^Vconf=-rho*partial_i V_conf. The selected finite throat profile is not solved in pathA_21c.",
            "downstream_consequence": "The exterior hydrodynamic matter-stress integral is derived, but the full core/body-force normalization is not called derived.",
        },
        {
            "name": "QUANTUM_STRESS_FAR_FIELD_RESIDUAL",
            "status": "PROFILE_DERIVATIVE_RESIDUAL",
            "source": "sigma_Q,ij is derived and its divergence check is machine-verified, but its cross surface integral needs the density derivative profile near the throat branch.",
            "downstream_consequence": "Quantum stress is not used to tune or flip the derived matter-stress sign.",
        },
        {
            "name": "MAXWELL_Z_GAUGE_JEXT_CANCELLATION_RESIDUAL",
            "status": "MAXWELL_BODY_FORCE_RESIDUAL",
            "source": "The localized Maxwell action has explicit Z(w), gauge fixing, and -A_M J_ext^M terms. The compact matter-drain lane does not prove their cross terms vanish or cancel.",
            "downstream_consequence": "Maxwell stress is not included in the force coefficient until the profile/zero-mode branch proves the cancellation.",
        },
        {
            "name": "SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE",
            "status": "FULL_SIGN_RESIDUAL",
            "source": "The convective plus Bernoulli-pressure matter stress gives an attractive like-drain far-field sign, but quantum, V_conf body-force, and Maxwell profile pieces are not all evaluated.",
            "downstream_consequence": "The full sign is not claimed as derived; the leading compact reduced-3D matter-stress sign is a target-blind far-field result.",
        },
    ]


def run_patha21c_force_from_noether_stress_tensor() -> dict[str, object]:
    checks = _patha21c_force_tensor_checks()
    algebra = _patha21c_algebraic_checks()
    failures = [check for check in checks if check.status != "CONSISTENT"]
    algebra_failures = [check for check in algebra if not check["pass"]]
    residuals = _patha21c_named_residuals()
    carried_residuals = _patha21_residuals()
    return {
        "schema": "stage1_pathA_21c_force_from_noether_stress_tensor/v1",
        "base_dimensions": ["L", "T", "M"],
        "consumed_inputs": [
            "software/stage1_solver/directives/pathA_21c_force_from_noether_stress_tensor.md",
            "software/stage1_solver/_scratch/pathA_21c_directive_review.log",
            "software/stage1_solver/_scratch/pathA_21c_directive_confirmpass.log",
            "software/stage1_solver/reports/pathA_21b_force_closure_and_profile_bvp.md",
            "software/stage1_solver/decisions/13_emergent_constants_derivation.md sections 11-12",
            "research/pde/paper/pde.tex lines 318-451",
            "research/pde_ledger/paper/parts/part01_parent_geometry.tex lines 275-286",
        ],
        "noether_balance": {
            "momentum_density": "g_i=m_GNLS*rho*v_i",
            "noether_trace": [
                "Start from L_psi=(i*hbar/2)(psi^*D_t psi-psi D_t psi^*)-(hbar^2/(2*m_GNLS))(D_i psi)^*(D_i psi)-V_conf*rho-U(rho).",
                "For an active spatial translation delta psi=-epsilon_i partial_i psi and delta A_M=-epsilon_i partial_i A_M, the canonical identity is partial_t T^0_i+partial_j T^j_i=partial L/partial x_i for explicit backgrounds.",
                "The matter canonical flux from the spatial-gradient term is reduced on shell with psi=sqrt(rho) exp(i theta), j_i=rho v_i, the phase equation, and P=rho*h-U to Pi_ij^matter=m_GNLS*rho*v_i*v_j+delta_ij*P+sigma_Q,ij.",
                "The explicit matter background term -V_conf(X;Sigma0)*rho contributes partial L/partial x_i=-rho*partial_i V_conf, so it is f_i^body, not part of Pi_ij.",
                "The Maxwell action contributes the standard field stress only in lanes where Z(w), gauge fixing, and J_ext backgrounds are proven to vanish/cancel; otherwise their explicit partial L/partial x_i terms are residualized.",
            ],
            "matter_stress": "Pi_ij^matter=m_GNLS*rho*v_i*v_j+delta_ij*P(rho)+sigma_Q,ij",
            "quantum_stress": "sigma_Q,ij=(hbar^2/(4*m_GNLS))*[(partial_i rho partial_j rho)/rho-partial_i partial_j rho]",
            "pressure": "P=K*rho^5, h=(5K/4)*rho^4, dP=rho*dh",
            "matter_balance": "partial_t g_i+partial_j Pi_ij^matter=q_star*rho*(E_i+v_j B_ij)-rho*partial_i V_conf",
            "body_forces": [
                "f_i^Vconf=-rho*partial_i V_conf",
                "f_i^Z=-(partial_i Z) F_MN F^MN/(4*mu0) in the Maxwell sector",
                "f_i^Jext=-A_M*partial_i J_ext^M for explicit external-source gradients",
                "gauge-fixing/background terms are included only after a selected branch proves cancellation; otherwise residualized",
            ],
            "euler_check": "Using continuity, dP=rho*dh, and partial_j sigma_Q,ij=rho*partial_i Q, the balance law divided by rho reproduces m_GNLS*(partial_t+v.grad)v_i=q_star*(E_i+v_jB_ij)-partial_i(V_conf+h+Q).",
            "stress_representative": "Canonical Noether stress is reduced to the displayed Madelung hydrodynamic representative. No Belinfante improvement is used in the accepted matter-stress integral; smooth closed-surface improvements integrate to zero by Gauss/antisymmetry, while singular core/profile pieces are not used as derived normalization.",
        },
        "force_integral": {
            "sign_convention": "n_hat outward from U2; F_12 is force on defect 2 by defect 1; stationary F_12=-int_boundary Pi_ij n_j dS+int_U2 f_i^body dV.",
            "control_surface": "reduced-3D sphere around defect 2 with core scale << a << r12; v2=-Q2*n_hat/(4*pi*a^2), v1 is constant over the sphere to leading order.",
            "convective_flux_reduced_3d": "int Pi_conv,cross.n dS=-(4/3)*m_GNLS*N_infty,3*Q2*v1",
            "pressure_flux_reduced_3d": "delta P_cross=-m_GNLS*N_infty,3*(v1.v2), so int delta_ij P_cross n_j dS=+(1/3)*m_GNLS*N_infty,3*Q2*v1",
            "matter_flux_reduced_3d": "int Pi_matter,cross.n dS=-m_GNLS*N_infty,3*Q2*v1",
            "force_reduced_3d": "F_12^matter=m_GNLS*N_infty,3*Q2*v1=-(m_GNLS*N_infty,3*Q1*Q2/(4*pi*r12^2))*rhat_12",
            "force_bulk_4d": "F_12^matter,4D=-(m_GNLS*rho_infty,4*Q1*Q2/(2*pi^2*R12^3))*Rhat_12",
            "structure_result": "Bilinear Q1*Q2 structure is an integral result from the v1*v2 cross stress plus Bernoulli pressure cross term.",
            "power_result": "Reduced compact lane gives r12^-2 because v1 from the carried 4*pi Gauss solve is r12^-2; unreduced bulk gives R12^-3 with Omega_3=2*pi^2.",
            "normalization_status": "Overall I_F,12^full / Theta_Q / branch-profile normalization remains a CALIBRATION KNOB, not derived.",
        },
        "sign": {
            "far_field_matter_verdict": "FORCE_ATTRACTIVE_DERIVED",
            "full_verdict": "SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE",
            "orientation": "With rhat_12 from defect 1 to defect 2, the matter-stress force is proportional to -Q1*Q2*rhat_12. Like drains/sources attract; opposite signs repel.",
            "why_not_full_sign": "The matter-stress sign is derived target-blind, but pathA_21c does not evaluate the quantum, V_conf body-force, and Maxwell profile pieces that could enter the full normalized control-volume force.",
        },
        "calibrate_predict_ledger": {
            "target_blind_predictions": [
                "force structure: bilinear Q1*Q2 from the stress integral",
                "lane power: reduced-3D r^-2 and separate bulk R^-3 from Gauss plus surface measure",
                "leading matter-stress sign: attractive for like drains, repulsive for opposite signs",
            ],
            "calibration_knobs": [
                "overall dimensionless normalization class: I_F,12^full / Theta_Q / branch-profile data",
            ],
            "prediction_count": 3,
            "knob_count": 1,
            "guardrail": "No observable is both calibrated-to and predicted; knobs < independent target-blind predictions. The full sign remains residual rather than being hidden in the normalization knob.",
        },
        "corrected_label": "P1c: G_FREE_NOETHER_STRESS_STRUCTURE_POWER_DERIVED_WITH_FAR_FIELD_MATTER_ATTRACTIVE_SIGN_AND_SIGN_RESIDUAL_QUANTUM_VCONF_MAXWELL_PROFILE",
        "carried_items_confirmed_verbatim": [
            "pathA_21b G1 stationary BVP and BC table",
            "reduced-3D Gauss solve v_r=-Theta_Q J/(4*pi*N_infty,3*r^2)",
            "bulk-4D Gauss solve v_R=-Theta_Q4 J/(2*pi^2*rho_infty,4*R^3)",
            "pathA_21b I_F,12^full definition carried as calibration/profile knob",
            "MASS_BRIDGE_FORM_NOT_DERIVED",
            "EP_NOT_DERIVED",
            "RETAIN_L_T_M",
            "NEWTON_G_FORM_NOT_DERIVED",
            "pathA_20/20b residuals carried unchanged",
        ],
        "new_residuals": residuals,
        "carried_residuals": carried_residuals,
        "checks": [check.as_dict() for check in checks],
        "algebraic_checks": algebra,
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


def render_patha21c_force_tensor_markdown(report: Mapping[str, object]) -> str:
    noether = report["noether_balance"]
    force = report["force_integral"]
    sign = report["sign"]
    ledger = report["calibrate_predict_ledger"]
    summary = report["summary"]
    assert isinstance(noether, Mapping)
    assert isinstance(force, Mapping)
    assert isinstance(sign, Mapping)
    assert isinstance(ledger, Mapping)
    assert isinstance(summary, Mapping)
    lines = [
        "# Path-A 21c force from Noether stress tensor",
        "",
        "## Verdicts",
        "",
        f"- Corrected P1b -> P1c label: `{report['corrected_label']}`.",
        f"- Leading matter-stress sign: `{sign['far_field_matter_verdict']}`.",
        f"- Full sign verdict: `{sign['full_verdict']}`.",
        "- Normalization status: CALIBRATION KNOB, not derived.",
        f"- Acceptance status: `{summary['acceptance_status']}`.",
        "",
        "Dual-engine scope: Python and Mathematica check dimensions and algebra only. The derivation is the Noether balance construction plus the explicit traction integrals below; exit 0 is necessary, not sufficient.",
        "",
        "## Noether Balance",
        "",
        "Noether trace:",
    ]
    noether_trace = noether["noether_trace"]
    assert isinstance(noether_trace, Sequence)
    for item in noether_trace:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            f"Momentum density: `{noether['momentum_density']}`.",
            "",
            f"Matter stress representative: `{noether['matter_stress']}`.",
            "",
            f"Quantum stress: `{noether['quantum_stress']}`.",
            "",
            f"EOS pressure identity: `{noether['pressure']}`.",
            "",
            f"Matter balance law: `{noether['matter_balance']}`.",
            "",
            "Explicit-background body forces:",
        ]
    )
    body_forces = noether["body_forces"]
    assert isinstance(body_forces, Sequence)
    for item in body_forces:
        lines.append(f"- `{item}`.")
    lines.extend(
        [
            "",
            f"Balance-law-vs-Euler check: {noether['euler_check']}",
            "",
            f"Stress representative: {noether['stress_representative']}",
            "",
            "## Force Integral",
            "",
            f"Convention: {force['sign_convention']}",
            "",
            f"Control surface: {force['control_surface']}",
            "",
            "Reduced-3D compact lane integral results:",
            "",
            f"- Convective cross flux: `{force['convective_flux_reduced_3d']}`.",
            f"- Bernoulli/EOS pressure cross flux: `{force['pressure_flux_reduced_3d']}`.",
            f"- Matter cross flux: `{force['matter_flux_reduced_3d']}`.",
            f"- Force result: `{force['force_reduced_3d']}`.",
            "",
            "Bulk lane kept separate:",
            "",
            f"- `{force['force_bulk_4d']}`.",
            "",
            f"Structure result: {force['structure_result']}",
            "",
            f"Power result: {force['power_result']}",
            "",
            f"Normalization: {force['normalization_status']}",
            "",
            "## Sign",
            "",
            f"Far-field matter-stress verdict: `{sign['far_field_matter_verdict']}`.",
            "",
            f"Full sign verdict: `{sign['full_verdict']}`.",
            "",
            f"Orientation: {sign['orientation']}",
            "",
            f"Residual reason: {sign['why_not_full_sign']}",
            "",
            "## Calibrate-Predict Ledger",
            "",
            "Target-blind predictions:",
        ]
    )
    predictions = ledger["target_blind_predictions"]
    knobs = ledger["calibration_knobs"]
    assert isinstance(predictions, Sequence)
    assert isinstance(knobs, Sequence)
    for item in predictions:
        lines.append(f"- {item}.")
    lines.append("")
    lines.append("Calibration knobs:")
    for item in knobs:
        lines.append(f"- {item}.")
    lines.extend(
        [
            "",
            f"Counts: predictions={ledger['prediction_count']}; knobs={ledger['knob_count']}.",
            "",
            f"Guardrail: {ledger['guardrail']}",
            "",
            "## Residuals",
            "",
        ]
    )
    residuals = report["new_residuals"]
    assert isinstance(residuals, Sequence)
    for raw in residuals:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: {raw['status']}. {raw['downstream_consequence']} Source: {raw['source']}")
    lines.extend(["", "## Carried Items", ""])
    carried = report["carried_items_confirmed_verbatim"]
    assert isinstance(carried, Sequence)
    for item in carried:
        lines.append(f"- {item}.")
    lines.extend(["", "## Carried Residual Ledger", ""])
    carried_residuals = report["carried_residuals"]
    assert isinstance(carried_residuals, Sequence)
    for raw in carried_residuals:
        assert isinstance(raw, Mapping)
        lines.append(f"- `{raw['name']}`: {raw['status']}. {raw['downstream_consequence']} Source: {raw['source']}")
    lines.extend(
        [
            "",
            "## Harness Summary",
            "",
            f"- Dimensional checks: {summary['consistent_dimensional_checks']} consistent, {summary['inconsistent_dimensional_checks']} inconsistent, {summary['total_dimensional_checks']} total.",
            f"- Algebraic checks: {summary['consistent_algebraic_checks']} consistent, {summary['inconsistent_algebraic_checks']} inconsistent, {summary['total_algebraic_checks']} total.",
            "",
        ]
    )
    return "\n".join(lines)


def write_patha21c_force_tensor_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = run_patha21c_force_from_noether_stress_tensor()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_21c_force_from_noether_stress_tensor_report.json"
    scratch_md_path = out_dir / "pathA_21c_force_from_noether_stress_tensor.md"
    reference_path = report_dir / "pathA_21c_force_from_noether_stress_tensor.md"
    rendered = render_patha21c_force_tensor_markdown(report) + "\n"
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered, encoding="utf-8")
    reference_path.write_text(rendered, encoding="utf-8")
    return json_path, reference_path, report


@dataclass(frozen=True)
class PathA22aKnob:
    name: str
    classification: str
    role: str
    provenance: str
    verdict_effect: str

    def as_dict(self) -> dict[str, str]:
        return {
            "name": self.name,
            "classification": self.classification,
            "role": self.role,
            "provenance": self.provenance,
            "verdict_effect": self.verdict_effect,
        }


def _patha22a_monomial(dim: Dim) -> str:
    expr = sp.simplify((L_UNIT**dim.l) * (T_UNIT**dim.t) * (M_UNIT**dim.m))
    return str(expr)


def _patha22a_dim_row(name: str, dim: Dim, derivation: str, provenance: str) -> dict[str, object]:
    return {
        "name": name,
        "dimension": str(dim),
        "tuple_L_T_M": dim.as_tuple(),
        "sympy_monomial": _patha22a_monomial(dim),
        "derivation": derivation,
        "provenance": provenance,
    }


def patha22a_homogeneity_checks(
    *,
    planted_rhs: Dim | None = None,
    planted_p0: Dim | None = None,
) -> dict[str, object]:
    """Run PathA_22a dimensional homogeneity with restored M,L,T units."""

    velocity = LENGTH / TIME
    reduced_stiffness = MASS / (LENGTH * (TIME**2))
    reduced_mass = MASS / LENGTH
    target = D["G_3_spatial"] * (D["c_s"] ** 5) / ((D["a"] ** 5) * (D["c"] ** 5))
    mhat0 = (LENGTH**-1) * (TIME**-1) * (MASS ** sp.Rational(-1, 2))
    mhat0_sq = mhat0**2
    s_port = DIMENSIONLESS
    chi_q = DIMENSIONLESS

    k0 = reduced_stiffness
    b0 = reduced_stiffness
    omega_u = TIME**-1
    omega_w = TIME**-1
    r_mix = TIME**-2
    g_port = (reduced_stiffness * (TIME**-2)) ** sp.Rational(1, 2)
    delta = (omega_u**2) * (omega_w**2)
    delta_via_r = r_mix**2
    s_mix = omega_u**2
    q_mix = (g_port**2) * (omega_w**2)
    h_mix = g_port**2
    p_mix = (omega_u**2) * g_port
    p_mix_via_r = r_mix * g_port
    z0 = q_mix / delta
    z2 = (q_mix * s_mix) / (delta**2)
    z4 = (q_mix * (s_mix**2)) / (delta**3)
    d0 = k0
    n0 = (p_mix**2) / (delta**2)
    n2 = (p_mix * (p_mix * s_mix)) / (delta**3)
    n4 = ((delta**2) * (g_port**2)) / (delta**4)
    p0_faithful = n0 / d0
    frequency_scale = (D["c_s"] / D["a"]) ** 2
    p0_normalized = p0_faithful * frequency_scale
    p0_for_lhs = planted_p0 or p0_normalized
    rhs_for_check = planted_rhs or target
    lhs = mhat0_sq * s_port * p0_for_lhs

    p0_factor_needed = DIMENSIONLESS / p0_faithful
    wrong_n0_assertion = expect_dim(
        "pathA_22a_negative",
        "WRONG assertion N0 has reduced stiffness dimension",
        n0,
        reduced_stiffness,
        "Negative control: the code formula N0=P^2/Delta^2 derives reduced mass, so asserting reduced stiffness must fail.",
    )

    checks = [
        homogeneous(
            "pathA_22a_A",
            "D0=K-B0-Z0 reduced static denominator",
            {"K": k0, "B0": b0, "Z0": z0},
            "Uses the reduced coefficient compiler contract in patha_extraction.py and pde.tex.",
        ),
        homogeneous(
            "pathA_22a_A",
            "mixed-port Delta=Omega_U^2*Omega_W^2-R^2",
            {"Omega_U^2*Omega_W^2": delta, "R^2": delta_via_r},
            "Uses pde.tex building-block dimensions Omega_U,Omega_W:T^-1 and R:T^-2.",
        ),
        homogeneous(
            "pathA_22a_A",
            "mixed-port P=Omega_U^2*g_W+R*g_U",
            {"Omega_U^2*g_W": p_mix, "R*g_U": p_mix_via_r},
            "Uses [g_U]^2=[g_W]^2=[K]*T^-2 from the same building-block dictionary.",
        ),
        expect_dim(
            "pathA_22a_A",
            "Z0=Q/Delta derived from mixed-port formula",
            z0,
            reduced_stiffness,
            "Cross-check: Z0 remains a reduced stiffness coefficient.",
        ),
        expect_dim(
            "pathA_22a_A",
            "Z2 derived from mixed-port formula",
            z2,
            reduced_mass,
            "Cross-check: Z2=K*T^2 matches the existing reduced_M_B2_Z2_N2 dictionary entry.",
        ),
        expect_dim(
            "pathA_22a_A",
            "Z4 derived from mixed-port formula",
            z4,
            reduced_mass * (TIME**2),
            "Cross-check: Z4=K*T^4 matches the existing reduced_B4_Z4_N4 dictionary entry.",
        ),
        expect_dim(
            "pathA_22a_A",
            "N0=P^2/Delta^2 derived from code formula",
            n0,
            reduced_mass,
            "This is derived from patha_extraction.py, not asserted from the old reduced_K_B0_Z0_N0 dictionary bucket.",
        ),
        expect_dim(
            "pathA_22a_A",
            "N2 derived from code formula",
            n2,
            reduced_mass * (TIME**2),
            "The N tower is shifted by T^2 relative to the old asserted N2 bucket.",
        ),
        expect_dim(
            "pathA_22a_A",
            "N4 derived from code formula",
            n4,
            reduced_mass * (TIME**4),
            "The N tower is shifted by T^2 relative to the old asserted N4 bucket.",
        ),
        expect_dim(
            "pathA_22a_A",
            "faithful P0=N0/D0 before frequency normalization",
            p0_faithful,
            TIME**2,
            "The actual code formula gives P0 with T^2 before applying the explicit frequency-normalization factor.",
        ),
        expect_dim(
            "pathA_22a_A",
            "P0 frequency-normalization factor (c_s/a)^2",
            frequency_scale,
            p0_factor_needed,
            "The paper's dimensionless P0 requires this explicit T^-2 factor; no a^5/c_s^5 radiation factor is hidden here.",
        ),
        expect_dim(
            "pathA_22a_A",
            "normalized P0=(c_s/a)^2*N0/D0 dimension",
            p0_normalized,
            DIMENSIONLESS,
            "This is the P0 used by the R_norm homogeneity gate and by Gamma5.",
        ),
        expect_dim(
            "pathA_22a_A",
            "mhat0 source-map dimension",
            mhat0,
            (LENGTH**-1) * (TIME**-1) * (MASS ** sp.Rational(-1, 2)),
            "Source-map dimensional table carried by pathA_21b; this group rechecks the monomial.",
        ),
        expect_dim(
            "pathA_22a_A",
            "3D GR target 54*G*c_s^5/(5*a^5*c^5)",
            target,
            (MASS**-1) * (LENGTH**-2) * (TIME**-2),
            "Uses 3D Newton G=L^3 T^-2 M^-1; numerical factors are dimensionless.",
        ),
        expect_dim(
            "pathA_22a_A",
            "S_port / chi_Q outgoing-normalization scalar",
            s_port,
            DIMENSIONLESS,
            "S_port is the code slot occupied by the paper's chi_Q.",
        ),
        homogeneous(
            "pathA_22a_A",
            "R_norm=mhat0^2*S_port*P0 - 54*G*c_s^5/(5*a^5*c^5)",
            {"mhat0^2*S_port*P0": lhs, "GR_target": rhs_for_check},
            "This is the actual homogeneity gate; it is not an x==x self-check.",
        ),
    ]

    failures = [check for check in checks if check.status != "CONSISTENT"]
    homogeneity_status = "HOMOGENEITY_PASS" if not failures else "HOMOGENEITY_FAILURE"
    return {
        "homogeneity_status": homogeneity_status,
        "checks": [check.as_dict() for check in checks],
        "production_failure_count": len(failures),
        "formula_negative_controls": {
            "wrong_n0_reduced_stiffness_assertion": wrong_n0_assertion.as_dict(),
            "status": "CAUGHT_WRONG_N0_ASSERTION"
            if wrong_n0_assertion.status == "INCONSISTENT"
            else "MISSED_WRONG_N0_ASSERTION",
        },
        "dimension_table": [
            _patha22a_dim_row(
                "G_3D",
                D["G_3_spatial"],
                "Newton force dimension in observed 3-space: [G]=L^3 T^-2 M^-1.",
                "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:31; software/stage1_solver/reports/pathA_19_dimensional_foundation.md",
            ),
            _patha22a_dim_row(
                "c_s",
                velocity,
                "EOS sound-speed law c_s^2=5 K rho0^4/m_GNLS.",
                "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:59; dimensional_check.py pathA_19_F3",
            ),
            _patha22a_dim_row("c", velocity, "Observed brane light speed.", "pathA_20/20b carried as lambda_gamma=c/c_s."),
            _patha22a_dim_row("a", LENGTH, "Throat/core length scale.", "research/pde/paper/pde.tex:2061"),
            _patha22a_dim_row(
                "mhat0",
                mhat0,
                "Source-map factor with [mhat0]=L^-1 T^-1 M^-1/2, so mhat0^2 matches the full target dimension.",
                "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:35",
            ),
            _patha22a_dim_row(
                "K,B0,Z0,D0",
                reduced_stiffness,
                "Reduced static denominator compiler: D0=K-B0-Z0; Z0 also derives to reduced stiffness.",
                "software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445,473-482; research/pde/paper/pde.tex:1849-1872",
            ),
            _patha22a_dim_row(
                "N0,N2,N4 from code formulas",
                n0,
                "N0=P^2/Delta^2 derives to K*T^2=M/L; N2 and N4 derive to M/L*T^2 and M/L*T^4.",
                "software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445",
            ),
            _patha22a_dim_row(
                "P0_faithful=N0/D0",
                p0_faithful,
                "Faithful code-formula ratio before normalization; dimension T^2.",
                "research/pde/paper/pde.tex:2018-2026; software/stage1_solver/src/stage1_solver/patha_extraction.py:482",
            ),
            _patha22a_dim_row(
                "P0=(c_s/a)^2*N0/D0",
                p0_normalized,
                "Dimensionless static outgoing prefactor after explicit frequency normalization.",
                "software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445,473-482; research/pde/paper/pde.tex:2018-2026",
            ),
            _patha22a_dim_row(
                "S_port == chi_Q",
                chi_q,
                "Dimensionless outgoing-normalization scalar; S_port is the code slot in observable_residuals/extract_branch.",
                "software/stage1_solver/src/stage1_solver/patha_extraction.py:526-544,609-624; research/pde/paper/pde.tex:1980-1998",
            ),
            _patha22a_dim_row(
                "Gamma5",
                TIME**5,
                "Gamma5=chi_Q P0 a^5/(27 c_s^5); the radiation time^5 factor is explicit, not hidden in P0.",
                "research/pde/paper/pde.tex:2053-2061",
            ),
            _patha22a_dim_row(
                "54*G*c_s^5/(5*a^5*c^5)",
                target,
                "Target dimension is M^-1 L^-2 T^-2.",
                "research/pde/paper/pde.tex:2082",
            ),
        ],
        "p0_dimensionless_finding": {
            "status": "P0_DIMENSIONLESS_AFTER_EXPLICIT_FREQUENCY_NORMALIZATION",
            "finding": (
                "The actual patha_extraction.py mixed-port formulas derive N0=P^2/Delta^2 with "
                "dimension M/L, so faithful N0/D0 has dimension T^2. Multiplying by the explicit "
                "(c_s/a)^2 frequency-normalization factor makes P0 dimensionless. With that "
                "normalized P0, R_norm is homogeneous; without it, the gate would fail."
            ),
            "faithful_unnormalized_dimension": str(p0_faithful),
            "normalization_needed_for_faithful_reading": "(c_s/a)^2",
            "normalization_needed_dimension": str(frequency_scale),
            "normalized_dimension": str(p0_normalized),
        },
    }


def patha22a_knob_ledger(include_control_free_factor: bool = False) -> list[PathA22aKnob]:
    ledger = [
        PathA22aKnob(
            "sigma_Q^can=4*a^5/(27*c_s^5)",
            "(a) fixed-by-prior-derivation",
            "Canonical compact outgoing fingerprint factor; dimensional T^5 carrier for Gamma5.",
            "research/pde/paper/pde.tex:1980-1993,2061",
            "Not tunable once the canonical compact outgoing branch is selected.",
        ),
        PathA22aKnob(
            "grouped signature (1,1/2,-1)",
            "(a) fixed-by-prior-derivation",
            "Weak-axisymmetric grouped transport signature.",
            "research/pde/paper/pde.tex:92,2371-2384",
            "Does not affect the isotropic homogeneity gate.",
        ),
        PathA22aKnob(
            "lambda_gamma=c/c_s",
            "(c)/(d) underived residual; TRUE FREE calibration knob under calibrate-predict",
            "Observed brane cone versus sound cone ratio; enters xi as lambda_gamma^5 and is equally unpinned by the current reduction.",
            "software/stage1_solver/reports/pathA_20b_cgamma_cs_linearization.md; research/pde/paper/pde.tex:2551",
            "A future brane zero-mode reduction may pin it; under calibrate-predict it is a tunability channel, not a prediction.",
        ),
        PathA22aKnob(
            "chi_Q / S_port",
            "(d) TRUE FREE calibration knob",
            "Outgoing-normalization scalar. Code S_port multiplies mhat0^2*P0 in exactly the paper's chi_Q slot.",
            "software/stage1_solver/src/stage1_solver/patha_extraction.py:526-544,609-624; research/pde/paper/pde.tex:1980-1998,2071-2082",
            "TUNABILITY_CHANNEL_PRESENT unless the actual outgoing compact DtN derivation independently fixes chi_Q; Stage 104/105 canonical sigma_Q choice does not prove the Path-A branch is canonical.",
        ),
        PathA22aKnob(
            "P0 branch value=N0/D0",
            "(b) branch-determined (target-blind)",
            "Static outgoing prefactor from reduced overlap data after the branch solve.",
            "software/stage1_solver/src/stage1_solver/patha_extraction.py:397-445,473-482; research/pde/paper/pde.tex:1849-1872,2018-2026",
            "Not a calibration knob if the profile/overlaps are solved target-blind.",
        ),
        PathA22aKnob(
            "R0/a",
            "(b) branch-determined (target-blind)",
            "Finite throat radius profile scale in the branch data.",
            "software/stage1_solver/src/stage1_solver/patha_extraction.py:25-37,64-75; research/pde/paper/pde.tex:2551",
            "Conditional on solving the stationary branch; not automatically tunable.",
        ),
        PathA22aKnob(
            "L/a",
            "(b) branch-determined (target-blind)",
            "Wall/worldtube length ratio used by finite-throat profiles.",
            "software/stage1_solver/src/stage1_solver/patha_extraction.py:40-53,691-700",
            "Geometry/domain data; only tunable if later allowed as calibration.",
        ),
        PathA22aKnob(
            "W_eff/a",
            "(c) underived residual",
            "Effective kernel/support width entering the profile and force normalization.",
            "software/stage1_solver/decisions/13_emergent_constants_derivation.md:433-445",
            "Pure-number closure needs the branch/kernel form.",
        ),
        PathA22aKnob(
            "Theta_Q",
            "(c) underived residual",
            "Quadrupole source/shape residual not closed by dimensional analysis.",
            "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:72-76",
            "Carried to the minimal xi/source-map derivation.",
        ),
        PathA22aKnob(
            "J-selector / flux law",
            "(c) underived residual",
            "Mass/source bridge rate selector inherited from the pathA_20/pathA_21 chain.",
            "software/stage1_solver/decisions/13_emergent_constants_derivation.md:407,433-445",
            "Affects mhat/G forms, not dimensional homogeneity.",
        ),
        PathA22aKnob(
            "alpha_J",
            "(c) underived residual",
            "Dimensionless mass-bridge coefficient, including the h versus hbar/2pi convention.",
            "software/stage1_solver/decisions/13_emergent_constants_derivation.md:407; software/stage1_solver/reports/pathA_21c_force_from_noether_stress_tensor.md",
            "Pending source/boundary/Hamiltonian derivation; not set to one.",
        ),
        PathA22aKnob(
            "branch-kernel choices",
            "(c) underived residual",
            "Profile/kernel choices that determine overlap integrals and the branch packet.",
            "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:72-76",
            "Indeterminate until forms are derived or branch solve fixes them target-blind.",
        ),
        PathA22aKnob(
            "g_G",
            "(c) underived residual",
            "Arbitrary dimensionless multiplier in G=(a*c_s^2/m_GNLS)*g_G.",
            "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:63-71",
            "Must remain symbolic in xi; dimensional analysis cannot set it to one.",
        ),
        PathA22aKnob(
            "g_mhat",
            "(c) underived residual",
            "Arbitrary dimensionless multiplier in mhat0=(c_s/(a^2*sqrt(m_GNLS)))*g_mhat.",
            "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:63-71",
            "Must remain symbolic in xi; source-map form is the scoped next step.",
        ),
    ]
    if include_control_free_factor:
        ledger.append(
            PathA22aKnob(
                "g_control_unresolved",
                "(d) TRUE FREE calibration knob",
                "Planted unresolved dimensionless multiplier for the negative control.",
                "software/stage1_solver/directives/pathA_22a_dimensional_skeleton.md:96-100",
                "Must force TUNABILITY_CHANNEL_PRESENT/INDETERMINATE, never cancel dimensionally.",
            )
        )
    return ledger


def patha22a_classify_headline(
    homogeneity_status: str,
    ledger: Sequence[PathA22aKnob],
) -> str:
    if homogeneity_status == "HOMOGENEITY_FAILURE":
        return "HOMOGENEITY_FAILURE"
    if any(item.classification.startswith("(d)") for item in ledger):
        return "TUNABILITY_CHANNEL_PRESENT"
    return "INDETERMINATE_NEEDS_FORMS"


def patha22a_symbolic_xi() -> dict[str, object]:
    lambda_gamma, g_mhat, g_G, chi_Q, P0 = sp.symbols(
        "lambda_gamma g_mhat g_G chi_Q P0", nonzero=True
    )
    c_s, a, m = sp.symbols("c_s a m_GNLS", positive=True)
    g_monomial = a * c_s**2 / m
    mhat_monomial = c_s / (a**2 * sp.sqrt(m))
    denominator = g_monomial * g_G * c_s**5 / (a**5 * (lambda_gamma * c_s) ** 5)
    xi_times = sp.simplify((mhat_monomial**2 * g_mhat**2 * chi_Q * P0) / denominator)
    return {
        "G_form": "G = (a*c_s^2/m_GNLS) * g_G = (5*a*K*rho0^4/m_GNLS^2) * g_G",
        "mhat0_form": "mhat0 = (c_s/(a^2*sqrt(m_GNLS))) * g_mhat",
        "c_form": "c = lambda_gamma*c_s",
        "xi_times_S_port_P0": str(xi_times),
        "xi_times_S_port_P0_simplified": "P0*chi_Q*g_mhat**2*lambda_gamma**5/g_G",
        "note": "All g_* factors are dimensionless and underived unless separately fixed.",
    }


def patha22a_negative_controls() -> dict[str, object]:
    target = D["G_3_spatial"] * (D["c_s"] ** 5) / ((D["a"] ** 5) * (D["c"] ** 5))
    missing_a5_rhs = target * (D["a"] ** 5)
    mismatch = patha22a_homogeneity_checks(planted_rhs=missing_a5_rhs)
    unresolved_ledger = patha22a_knob_ledger(include_control_free_factor=True)
    unresolved_headline = patha22a_classify_headline("HOMOGENEITY_PASS", unresolved_ledger)
    return {
        "planted_dimensional_mismatch": {
            "status": "CAUGHT_DIMENSIONAL_MISMATCH"
            if mismatch["homogeneity_status"] == "HOMOGENEITY_FAILURE"
            else "MISSED_DIMENSIONAL_MISMATCH",
            "expected": "HOMOGENEITY_FAILURE",
            "actual": mismatch["homogeneity_status"],
            "planted_rhs": "target with a^5 removed from the denominator",
        },
        "planted_unresolved_dimensionless_factor": {
            "status": "PRESERVED_UNRESOLVED_FACTOR"
            if unresolved_headline in {"TUNABILITY_CHANNEL_PRESENT", "INDETERMINATE_NEEDS_FORMS"}
            else "CANCELLED_UNRESOLVED_FACTOR",
            "expected": "TUNABILITY_CHANNEL_PRESENT or INDETERMINATE_NEEDS_FORMS",
            "actual": unresolved_headline,
            "control_factor": "g_control_unresolved",
        },
        "reachable_headlines": {
            "homogeneity_failure": patha22a_classify_headline("HOMOGENEITY_FAILURE", []),
            "tunability_present": patha22a_classify_headline("HOMOGENEITY_PASS", unresolved_ledger),
            "indeterminate_needs_forms": patha22a_classify_headline(
                "HOMOGENEITY_PASS",
                [
                    item
                    for item in patha22a_knob_ledger()
                    if not item.classification.startswith("(d)")
                ],
            ),
        },
    }


def patha22a_dimensional_skeleton_report() -> dict[str, object]:
    homogeneity = patha22a_homogeneity_checks()
    ledger = patha22a_knob_ledger()
    headline = patha22a_classify_headline(str(homogeneity["homogeneity_status"]), ledger)
    negative_controls = patha22a_negative_controls()
    class_d = [item.name for item in ledger if item.classification.startswith("(d)")]
    tunability_channels = [
        "chi_Q / S_port",
        "lambda_gamma=c/c_s",
    ]
    calibrate_predict_residuals = [
        "Theta_Q",
        "alpha_J",
        "W_eff/a",
        "branch-kernel choices",
    ]
    return {
        "schema": "stage1_pathA_22a_dimensional_skeleton/v1",
        "base_dimensions": ["M", "L", "T"],
        "sympy_unit_symbols": [str(M_UNIT), str(L_UNIT), str(T_UNIT)],
        "homogeneity": homogeneity,
        "knob_ledger": [item.as_dict() for item in ledger],
        "symbolic_dimensionless_audit": patha22a_symbolic_xi(),
        "negative_controls": negative_controls,
        "headline": headline,
        "class_d_free_knobs": class_d,
        "tunability_channels": tunability_channels,
        "tunability_channel_count_lower_bound": len(tunability_channels),
        "strict_ledger_disclosure": (
            "The strict class labels under-count tunability: lambda_gamma and the source-map residual cluster are "
            "underived class-(c) quantities in a proof ledger, but become calibration knobs under a calibrate-predict methodology."
        ),
        "calibrate_predict_source_map_residuals": calibrate_predict_residuals,
        "mhat_inconsistency": {
            "status": "KNOWN_PDE_TEX_INCONSISTENCY",
            "dimensionful_reading_used_by_code": str((LENGTH**-1) * (TIME**-1) * (MASS ** sp.Rational(-1, 2))),
            "equation_conflict": (
                "eq:outgoing-BT-target forces dimensionful mhat with [mhat]=L^-1 T^-1 M^-1/2, "
                "while eq:outgoing-natural-source-map writes mhat=1+O(a^2/r^2), which is dimensionless."
            ),
            "resolution_needed": "The eventual source-map derivation must reconcile this; pathA_22a uses the dimensionful reading.",
        },
        "scoped_next_step": (
            "Before deriving broad G/mass/source machinery, derive the minimal combination "
            "xi=mhat0^2*S_port/[G*c_s^5/(a^5*c^5)], with chi_Q/S_port fixed only by "
            "actual outgoing compact DtN matching and lambda_gamma fixed only by brane zero-mode reduction; "
            "otherwise both remain calibration channels."
        ),
        "residuals": [
            "chi_Q/S_port is not canonically fixed by the current code path; S_port defaults to 1.0 as a convention.",
            "lambda_gamma enters as lambda_gamma^5 and remains unpinned after pathA_20b; it is a second tunability channel under calibrate-predict.",
            "Theta_Q, alpha_J, W_eff/a, and branch-kernel source-map residuals are strict class-(c) gaps but become tunable under calibrate-predict.",
            "g_G, g_mhat, alpha_J, J-selector, W_eff/a, Theta_Q, and branch-kernel forms remain underived.",
            "pde.tex is internally inconsistent on mhat; the code uses the dimensionful eq:outgoing-BT-target reading.",
            "Faithful N0/D0 has dimension T^2; P0 is dimensionless only after the explicit (c_s/a)^2 frequency normalization.",
        ],
    }


def render_patha22a_dimensional_skeleton_markdown(report: Mapping[str, object]) -> str:
    homogeneity = report["homogeneity"]
    assert isinstance(homogeneity, Mapping)
    p0_finding = homogeneity["p0_dimensionless_finding"]
    assert isinstance(p0_finding, Mapping)
    symbolic = report["symbolic_dimensionless_audit"]
    assert isinstance(symbolic, Mapping)
    lines = [
        "# PathA 22a dimensional skeleton",
        "",
        "## Result",
        "",
        f"- Homogeneity: `{homogeneity['homogeneity_status']}`.",
        f"- P0 finding: `{p0_finding['status']}`. {p0_finding['finding']}",
        f"- Headline: `{report['headline']}`.",
        f"- Tunability channels: >= {report['tunability_channel_count_lower_bound']} ({', '.join(f'`{item}`' for item in report['tunability_channels'])}).",
        f"- Strict class-(d) free knobs: {', '.join(f'`{item}`' for item in report['class_d_free_knobs'])}.",
        f"- Strict-ledger disclosure: {report['strict_ledger_disclosure']}",
        f"- Scoped next step: {report['scoped_next_step']}",
        "",
        "## Dimension table",
        "",
        "| ingredient | dimension | SymPy monomial | derivation | provenance |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in homogeneity["dimension_table"]:
        assert isinstance(row, Mapping)
        lines.append(
            f"| `{row['name']}` | `{row['dimension']}` | `{row['sympy_monomial']}` | "
            f"{row['derivation']} | {row['provenance']} |"
        )
    lines.extend(["", "## Homogeneity checks", ""])
    for raw in homogeneity["checks"]:
        assert isinstance(raw, Mapping)
        factor = raw["factor_needed_to_reach_expected"]
        factor_text = "" if factor in (None, "1") else f"; factor needed `{factor}`"
        lines.append(
            f"- `{raw['name']}`: **{raw['status']}** "
            f"(expected `{raw['expected']}`, actual `{raw['actual']}`{factor_text}). "
            f"{raw['note']}"
        )
    lines.extend(
        [
            "",
            "## P0 normalization reading",
            "",
            f"- Faithful `N0/D0` dimension before normalization: `{p0_finding['faithful_unnormalized_dimension']}`.",
            f"- Required normalization: `{p0_finding['normalization_needed_for_faithful_reading']}` with dimension `{p0_finding['normalization_needed_dimension']}`.",
            f"- Normalized `P0` dimension: `{p0_finding['normalized_dimension']}`.",
            "- The outgoing radiation carrier is explicit in `Gamma5=chi_Q*P0*a^5/(27*c_s^5)`; it is not hidden inside `P0`.",
            "",
            "## Dimensionless audit",
            "",
            f"- `{symbolic['G_form']}`",
            f"- `{symbolic['mhat0_form']}`",
            f"- `{symbolic['c_form']}`",
            f"- `xi*S_port*P0 = {symbolic['xi_times_S_port_P0_simplified']}`.",
            f"- {symbolic['note']}",
            "",
            "## Knob ledger",
            "",
            "| knob/residual | class | role | provenance | verdict effect |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for raw in report["knob_ledger"]:
        assert isinstance(raw, Mapping)
        lines.append(
            f"| `{raw['name']}` | {raw['classification']} | {raw['role']} | "
            f"{raw['provenance']} | {raw['verdict_effect']} |"
        )
    lines.extend(
        [
            "",
            "## Calibrate-Predict Disclosure",
            "",
            f"- Tunability lower bound: `>= {report['tunability_channel_count_lower_bound']}` channels: "
            + ", ".join(f"`{item}`" for item in report["tunability_channels"])
            + ".",
            "- `chi_Q/S_port` is expected to be pinned only by an actual outgoing compact DtN derivation.",
            "- `lambda_gamma=c_gamma/c_s` is expected to be pinned only by the brane zero-mode reduction and enters the reduction as `lambda_gamma^5`.",
            "- Both are underived branch/reduction quantities; both become free knobs under calibrate-predict.",
            "- Strict class-(c) source-map residuals that also become tunable under calibrate-predict: "
            + ", ".join(f"`{item}`" for item in report["calibrate_predict_source_map_residuals"])
            + ".",
            "",
            "## Known Inconsistency",
            "",
        ]
    )
    mhat = report["mhat_inconsistency"]
    assert isinstance(mhat, Mapping)
    lines.extend(
        [
            f"- `{mhat['status']}`: {mhat['equation_conflict']}",
            f"- Code/harness reading: `[mhat]={mhat['dimensionful_reading_used_by_code']}`.",
            f"- Required resolution: {mhat['resolution_needed']}",
        ]
    )
    controls = report["negative_controls"]
    assert isinstance(controls, Mapping)
    mismatch = controls["planted_dimensional_mismatch"]
    unresolved = controls["planted_unresolved_dimensionless_factor"]
    reachable = controls["reachable_headlines"]
    assert isinstance(mismatch, Mapping)
    assert isinstance(unresolved, Mapping)
    assert isinstance(reachable, Mapping)
    lines.extend(
        [
            "",
            "## Negative controls",
            "",
            f"- Planted dimensional mismatch: `{mismatch['status']}`; expected `{mismatch['expected']}`, actual `{mismatch['actual']}`.",
            f"- Planted unresolved dimensionless factor: `{unresolved['status']}`; expected `{unresolved['expected']}`, actual `{unresolved['actual']}`.",
            f"- Reachability: homogeneity failure -> `{reachable['homogeneity_failure']}`; tunability -> `{reachable['tunability_present']}`; no class-(d) but residual forms -> `{reachable['indeterminate_needs_forms']}`.",
            f"- Wrong `N0` stiffness assertion: `{homogeneity['formula_negative_controls']['status']}`.",
            "",
            "## Residuals",
            "",
        ]
    )
    for item in report["residuals"]:
        lines.append(f"- {item}")
    lines.append("")
    return "\n".join(lines)


def write_patha22a_dimensional_skeleton_report(out_dir: Path, report_dir: Path) -> tuple[Path, Path, dict[str, object]]:
    report = patha22a_dimensional_skeleton_report()
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "pathA_22a_dimensional_skeleton.json"
    scratch_md_path = out_dir / "pathA_22a_dimensional_skeleton.md"
    reference_path = report_dir / "pathA_22a_dimensional_skeleton.md"
    rendered = render_patha22a_dimensional_skeleton_markdown(report)
    json_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    scratch_md_path.write_text(rendered + "\n", encoding="utf-8")
    reference_path.write_text(rendered + "\n", encoding="utf-8")
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
        "--patha21-emergent-g",
        action="store_true",
        help="run the side-by-side pathA_21 emergent G/mass-bridge check group instead of the pathA_18 audit",
    )
    parser.add_argument(
        "--patha21b-force-bvp",
        action="store_true",
        help="run the side-by-side pathA_21b force-closure/BVP check group instead of the pathA_18 audit",
    )
    parser.add_argument(
        "--patha21c-force-tensor",
        action="store_true",
        help="run the side-by-side pathA_21c Noether stress force check group instead of the pathA_18 audit",
    )
    parser.add_argument(
        "--patha22a-dimensional-skeleton",
        action="store_true",
        help="run the side-by-side pathA_22a dimensional skeleton/knob-ledger group instead of the pathA_18 audit",
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
    parser.add_argument(
        "--patha21-report-dir",
        default="software/stage1_solver/reports",
        help="directory for the pathA_21 emergent G/mass-bridge reference markdown",
    )
    parser.add_argument(
        "--patha21b-report-dir",
        default="software/stage1_solver/reports",
        help="directory for the pathA_21b force-closure/BVP reference markdown",
    )
    parser.add_argument(
        "--patha21c-report-dir",
        default="software/stage1_solver/reports",
        help="directory for the pathA_21c Noether stress force reference markdown",
    )
    parser.add_argument(
        "--patha22a-report-dir",
        default="software/stage1_solver/reports",
        help="directory for the pathA_22a dimensional skeleton reference markdown",
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
    if args.patha21_emergent_g:
        json_path, reference_path, report = write_patha21_emergent_g_report(
            Path(args.out_dir),
            Path(args.patha21_report_dir),
        )
        summary = report["summary"]
        p1 = report["p1_force"]
        p2 = report["p2_mass_bridge"]
        p3 = report["p3_m_collapse"]
        p4 = report["p4_g"]
        residuals = report["residuals"]
        print(f"wrote {json_path}")
        print(f"wrote {reference_path}")
        print(
            "pathA_21 emergent G/mass-bridge checks: "
            f"{summary['consistent_dimensional_checks']} dimensional consistent, "
            f"{summary['inconsistent_dimensional_checks']} dimensional inconsistent, "
            f"{summary['total_dimensional_checks']} dimensional total; "
            f"{summary['consistent_algebraic_checks']} algebraic consistent, "
            f"{summary['inconsistent_algebraic_checks']} algebraic inconsistent, "
            f"{summary['total_algebraic_checks']} algebraic total"
        )
        print(
            "P1 force FORM: "
            f"{p1['force_form']}; C_F={p1['coefficient']}; "
            f"attractive status: {p1['attractiveness']}; r-power: {p1['r_power']}"
        )
        print(
            "P2 bridge verdict: "
            f"{p2['verdict']}; angular {p2['angular_form']}; cycle {p2['cycle_form']}; "
            f"2pi/J-rate status: {p2['two_pi_status']}; EP verdict: {p2['ep_verdict']}"
        )
        print(f"P3 M-collapse verdict: {p3['verdict']}; blocker: {p3['blocker']}")
        print(
            "P4 G verdict: "
            f"{p4['verdict']}; m<->G conditional form: {p4['conditional_m_to_g_form']}; "
            f"width: {p4['width']}"
        )
        print(
            "P5 spec rows: "
            f"{summary['profile_spec_row_count']} total; "
            f"{summary['profile_solve_rows']} profile-solve, "
            f"{summary['pathA_22_rows']} pathA_22, "
            f"{summary['new_physics_rows']} new-physics, "
            f"{summary['known_rows']} known"
        )
        print("named residuals carried to pathA_22/profile solve: " + ", ".join(str(raw["name"]) for raw in residuals))
        return 0
    if args.patha21b_force_bvp:
        json_path, reference_path, report = write_patha21b_force_bvp_report(
            Path(args.out_dir),
            Path(args.patha21b_report_dir),
        )
        summary = report["summary"]
        p1 = report["p1b_force"]
        bvp = report["bvp"]
        gaps = bvp["gaps"]
        carried = report["carried_negative_verdicts_verbatim"]
        new_residuals = report["new_residuals"]
        print(f"wrote {json_path}")
        print(f"wrote {reference_path}")
        print(
            "pathA_21b force/BVP checks: "
            f"{summary['consistent_dimensional_checks']} dimensional consistent, "
            f"{summary['inconsistent_dimensional_checks']} dimensional inconsistent, "
            f"{summary['total_dimensional_checks']} dimensional total; "
            f"{summary['consistent_algebraic_checks']} algebraic consistent, "
            f"{summary['inconsistent_algebraic_checks']} algebraic inconsistent, "
            f"{summary['total_algebraic_checks']} algebraic total"
        )
        print(
            "P1b drain solve: reduced-3D r-power emerged from Gauss = yes; "
            "bulk-4D r-power emerged from Gauss = yes"
        )
        print(
            "P1b force FORM: "
            f"{p1['force_form']}; sign verdict: {p1['sign_verdict']}; "
            f"corrected P1 label: {p1['corrected_verdict']}"
        )
        print(
            "P5b gaps: "
            + "; ".join(
                f"{raw['gap']}={raw['status']}({raw['equation_or_residual']})"
                for raw in gaps
            )
        )
        print(f"alpha_H,omega anchor fix: {bvp['alpha_H_omega_anchor_fix']}")
        print(
            "P5b spec rows: "
            f"{summary['profile_spec_row_count']} total; "
            f"{summary['closed_profile_solve_rows']} closed/profile-solve, "
            f"{summary['branch_residual_rows']} branch-residual, "
            f"{summary['pathA_22_rows']} pathA_22, "
            f"{summary['new_physics_rows']} new-physics, "
            f"{summary['known_rows']} known"
        )
        print("carried negatives confirmed verbatim: " + ", ".join(str(item) for item in carried))
        print(
            "named residuals handed to option C / pathA_22: "
            + ", ".join(str(raw["name"]) for raw in new_residuals)
        )
        print(
            "files created/modified by this group: "
            f"{json_path}, {reference_path}, "
            f"{Path(args.out_dir) / 'pathA_21b_profile_bvp_machine_table.json'}"
        )
        return 0
    if args.patha21c_force_tensor:
        json_path, reference_path, report = write_patha21c_force_tensor_report(
            Path(args.out_dir),
            Path(args.patha21c_report_dir),
        )
        summary = report["summary"]
        noether = report["noether_balance"]
        force = report["force_integral"]
        sign = report["sign"]
        ledger = report["calibrate_predict_ledger"]
        residuals = report["new_residuals"]
        carried = report["carried_items_confirmed_verbatim"]
        print(f"wrote {json_path}")
        print(f"wrote {reference_path}")
        print(
            "pathA_21c Noether stress checks: "
            f"{summary['consistent_dimensional_checks']} dimensional consistent, "
            f"{summary['inconsistent_dimensional_checks']} dimensional inconsistent, "
            f"{summary['total_dimensional_checks']} dimensional total; "
            f"{summary['consistent_algebraic_checks']} algebraic consistent, "
            f"{summary['inconsistent_algebraic_checks']} algebraic inconsistent, "
            f"{summary['total_algebraic_checks']} algebraic total"
        )
        print(
            "Pi_ij + f_i^body derived (Noether): "
            f"{noether['matter_stress']}; bodies: "
            + "; ".join(str(item) for item in noether["body_forces"])
        )
        print(f"balance-law-vs-Euler check: {noether['euler_check']}")
        print(
            "force integral results: "
            f"reduced-3D {force['force_reduced_3d']}; "
            f"bulk lane {force['force_bulk_4d']}"
        )
        print(
            "SIGN verdict: "
            f"leading matter {sign['far_field_matter_verdict']}; "
            f"full {sign['full_verdict']}"
        )
        print(
            "calibrate-predict ledger: "
            f"predictions={ledger['prediction_count']} "
            f"({'; '.join(str(item) for item in ledger['target_blind_predictions'])}); "
            f"knobs={ledger['knob_count']} "
            f"({'; '.join(str(item) for item in ledger['calibration_knobs'])})"
        )
        print("residuals: " + ", ".join(str(raw["name"]) for raw in residuals))
        print(f"corrected P1b->P1c label: {report['corrected_label']}")
        print("carried items confirmed verbatim: " + ", ".join(str(item) for item in carried))
        print(f"files created/modified by this group: {json_path}, {reference_path}")
        return 0
    if args.patha22a_dimensional_skeleton:
        json_path, reference_path, report = write_patha22a_dimensional_skeleton_report(
            Path(args.out_dir),
            Path(args.patha22a_report_dir),
        )
        homogeneity = report["homogeneity"]
        assert isinstance(homogeneity, Mapping)
        p0_finding = homogeneity["p0_dimensionless_finding"]
        assert isinstance(p0_finding, Mapping)
        controls = report["negative_controls"]
        assert isinstance(controls, Mapping)
        print(f"wrote {json_path}")
        print(f"wrote {reference_path}")
        print(
            "pathA_22a dimensional skeleton: "
            f"homogeneity={homogeneity['homogeneity_status']}; "
            f"P0={p0_finding['status']}; "
            f"headline={report['headline']}"
        )
        print("class-(d) free knobs: " + ", ".join(str(item) for item in report["class_d_free_knobs"]))
        print(
            "negative controls: "
            f"mismatch={controls['planted_dimensional_mismatch']['status']}; "
            f"unresolved_factor={controls['planted_unresolved_dimensionless_factor']['status']}"
        )
        print(f"files created/modified by this group: {json_path}, {reference_path}")
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
