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
    args = parser.parse_args(list(argv) if argv is not None else None)
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
