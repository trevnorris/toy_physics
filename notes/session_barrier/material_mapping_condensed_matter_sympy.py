#!/usr/bin/env python3
"""
Material-science mappings for the damped moving-throat cold-crossing branch.

This script uses SymPy to derive three requested mappings:
1) lattice shedding rate -> electron-phonon threshold (lambda_ep * omega_D),
2) chi_lambda trigger -> harmonic-trap geometry/stiffness map,
3) Korringa spin-lattice limit -> Tmax.

Important modeling note
-----------------------
For a pure harmonic trap V_lattice(r) = (1/2) k_eff r^2, the logarithmic gradient
used in chi_lambda is
    d ln(V_lattice)/dr = 2/r,
so chi_lambda alone constrains the geometry ratio r/lambda but does NOT fix k_eff.
To obtain a stiffness in eV/Å^2, this script adds the minimal extra condition
that the *absolute* lattice force gradient matches the reduced-model force at the
turning point. That is the only extra assumption introduced here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import sympy as sp


@dataclass(frozen=True)
class PriorRunInputs:
    """Numerical anchors imported from the earlier generated reports."""
    gamma_lattice_req: float = 4.79562976   # lattice share at gamma_safe
    gamma_vac_req: float = 1.59854325
    t_cross_red: float = 1.82169718         # characteristic cold-crossing time
    lambda_th_red: float = 0.42826825       # lambda_th(r_turn)
    r_turn_red: float = 0.39096144          # lowered-barrier outer turning point
    V_turn_red: float = 2.5                 # V_eff(r_turn) = E_sub on the turning point
    v0_cold: float = 2.6


def derive_symbolic_map(inputs: PriorRunInputs) -> Dict[str, sp.Expr]:
    # ------------------------------------------------------------------
    # Symbols
    # ------------------------------------------------------------------
    zeta_ep = sp.symbols("zeta_ep", positive=True)      # O(1) transfer efficiency
    t_star = sp.symbols("t_star", positive=True)        # seconds per reduced time unit
    t_cross_phys = sp.symbols("t_cross_phys", positive=True)  # physical crossing time
    lambda_ep, omega_D = sp.symbols("lambda_ep omega_D", positive=True)

    E_star_eV = sp.symbols("E_star_eV", positive=True)      # eV per reduced energy unit
    lambda_phys_A = sp.symbols("lambda_phys_A", positive=True)  # Å, physical localization width
    a_int_A = sp.symbols("a_int_A", positive=True)          # Å, representative interstitial spacing
    r_phys = sp.symbols("r_phys", positive=True)
    k_eff = sp.symbols("k_eff", positive=True)

    Kcorr = sp.symbols("Kcorr", positive=True)  # Korringa constant [K*s]

    # Fixed previous-run numbers as exact SymPy floats/rationals
    gamma_lat = sp.Float(inputs.gamma_lattice_req)
    t_cross_red = sp.Float(inputs.t_cross_red)
    lambda_th = sp.Float(inputs.lambda_th_red)
    r_turn = sp.Float(inputs.r_turn_red)
    V_turn = sp.Float(inputs.V_turn_red)

    # ------------------------------------------------------------------
    # Task 1: gamma_lattice -> lambda_ep * omega_D
    # gamma_lattice^phys = gamma_lattice^red / t_star
    # Model: gamma_lattice^phys = zeta_ep * lambda_ep * omega_D
    # ------------------------------------------------------------------
    prod_threshold = sp.simplify(gamma_lat / (zeta_ep * t_star))
    prod_threshold_cross = sp.simplify(prod_threshold.subs(t_star, t_cross_phys / t_cross_red))
    dimless_cross_threshold = sp.simplify(prod_threshold_cross * t_cross_phys)
    # This last expression is the clean reduced-to-physical matching invariant:
    # lambda_ep * omega_D * t_cross_phys >= gamma_lattice_red * t_cross_red / zeta_ep

    # ------------------------------------------------------------------
    # Task 2: chi_lambda geometry test + force-matched stiffness
    # For harmonic trap, d ln(V)/dr = 2/r and chi_lambda = 2 lambda / r.
    # ------------------------------------------------------------------
    V_lattice = sp.Rational(1, 2) * k_eff * r_phys**2
    dlnV_lattice = sp.simplify(sp.diff(sp.log(V_lattice), r_phys))
    chi_lattice = sp.simplify(lambda_phys_A * dlnV_lattice)
    # Threshold chi_lambda >= 1  ->  r_phys <= 2 lambda_phys
    r_threshold = sp.solve_univariate_inequality(chi_lattice >= 1, r_phys)

    # Use the previous reduced branch to evaluate chi on the turning point
    ell_star_A = sp.simplify(lambda_phys_A / lambda_th)      # Å per reduced length unit
    r_turn_phys_A = sp.simplify(r_turn * ell_star_A)
    chi_on_turn = sp.simplify((2 * lambda_phys_A / r_turn_phys_A))

    # Minimal extra condition to obtain k_eff:
    # match the absolute harmonic-trap force to the reduced-model force at r_turn.
    # lambda_th = |V / V'| on the earlier branch, so |V'|_turn = V_turn / lambda_th in reduced units.
    dV_red_turn = sp.simplify(V_turn / lambda_th)
    dV_phys_turn_per_A = sp.simplify(E_star_eV * dV_red_turn / ell_star_A)  # eV/Å
    k_eff_required = sp.simplify(dV_phys_turn_per_A / r_turn_phys_A)        # eV/Å^2

    # Same stiffness written against an interstitial spacing a_int = 2 lambda_phys
    k_eff_interstitial = sp.simplify(k_eff_required.subs(lambda_phys_A, a_int_A / 2))

    # ------------------------------------------------------------------
    # Task 3: Korringa Tmax
    # Safe if T1 >= t_cross_phys with T1 * T = Kcorr.
    # ------------------------------------------------------------------
    Tmax = sp.simplify(Kcorr / t_cross_phys)
    Tmax_reduced_map = sp.simplify(Tmax.subs(t_cross_phys, t_cross_red * t_star))

    return {
        "prod_threshold": prod_threshold,
        "prod_threshold_cross": prod_threshold_cross,
        "dimless_cross_threshold": dimless_cross_threshold,
        "V_lattice": V_lattice,
        "dlnV_lattice": dlnV_lattice,
        "chi_lattice": chi_lattice,
        "r_threshold": r_threshold,
        "r_turn_phys_A": r_turn_phys_A,
        "chi_on_turn": chi_on_turn,
        "k_eff_required": k_eff_required,
        "k_eff_interstitial": k_eff_interstitial,
        "Tmax": Tmax,
        "Tmax_reduced_map": Tmax_reduced_map,
        "symbols": {
            "zeta_ep": zeta_ep,
            "t_star": t_star,
            "t_cross_phys": t_cross_phys,
            "E_star_eV": E_star_eV,
            "lambda_phys_A": lambda_phys_A,
            "a_int_A": a_int_A,
            "Kcorr": Kcorr,
        },
    }


def numeric_summary(inputs: PriorRunInputs, exprs: Dict[str, sp.Expr]) -> Dict[str, float]:
    lambda_phys_A = exprs["symbols"]["lambda_phys_A"]
    E_star_eV = exprs["symbols"]["E_star_eV"]
    a_int_A = exprs["symbols"]["a_int_A"]

    # Reference numeric evaluations
    keff_ref_1A_1eV = float(exprs["k_eff_required"].subs({lambda_phys_A: 1.0, E_star_eV: 1.0}))
    keff_ref_05A_1eV = float(exprs["k_eff_required"].subs({lambda_phys_A: 0.5, E_star_eV: 1.0}))
    keff_ref_15A_1eV = float(exprs["k_eff_required"].subs({lambda_phys_A: 1.5, E_star_eV: 1.0}))
    keff_interstitial_2A_1eV = float(exprs["k_eff_interstitial"].subs({a_int_A: 2.0, E_star_eV: 1.0}))
    dimless_cross_threshold = float(sp.N(exprs["dimless_cross_threshold"].subs({exprs["symbols"]["zeta_ep"]: 1.0})))
    chi_on_turn = float(sp.N(exprs["chi_on_turn"]))

    return {
        "gamma_lattice_req": inputs.gamma_lattice_req,
        "t_cross_red": inputs.t_cross_red,
        "lambda_th_red": inputs.lambda_th_red,
        "r_turn_red": inputs.r_turn_red,
        "V_turn_red": inputs.V_turn_red,
        "dimless_cross_threshold": dimless_cross_threshold,
        "chi_on_turn": chi_on_turn,
        "keff_ref_1A_1eV": keff_ref_1A_1eV,
        "keff_ref_05A_1eV": keff_ref_05A_1eV,
        "keff_ref_15A_1eV": keff_ref_15A_1eV,
        "keff_interstitial_2A_1eV": keff_interstitial_2A_1eV,
    }


def write_report(outpath: Path, inputs: PriorRunInputs, exprs: Dict[str, sp.Expr], nums: Dict[str, float]) -> None:
    zeta_ep = exprs["symbols"]["zeta_ep"]
    t_star = exprs["symbols"]["t_star"]
    t_cross_phys = exprs["symbols"]["t_cross_phys"]
    E_star_eV = exprs["symbols"]["E_star_eV"]
    lambda_phys_A = exprs["symbols"]["lambda_phys_A"]
    a_int_A = exprs["symbols"]["a_int_A"]
    Kcorr = exprs["symbols"]["Kcorr"]

    lines = []
    lines.append("Material mapping for the damped moving-throat cold-crossing branch")
    lines.append("=" * 92)
    lines.append("")
    lines.append("Inputs carried from earlier reduced runs")
    lines.append("-" * 92)
    lines.append(f"gamma_lattice,req  = {inputs.gamma_lattice_req:.8f}")
    lines.append(f"t_cross,red        = {inputs.t_cross_red:.8f}")
    lines.append(f"lambda_th,red      = {inputs.lambda_th_red:.8f}")
    lines.append(f"r_turn,red         = {inputs.r_turn_red:.8f}")
    lines.append(f"V_eff(r_turn),red  = {inputs.V_turn_red:.8f}")
    lines.append("")

    lines.append("Task 1 — electron-phonon threshold")
    lines.append("-" * 92)
    lines.append("Model:")
    lines.append("  gamma_lattice^phys = zeta_ep * lambda_ep * omega_D")
    lines.append("  gamma_lattice^phys = gamma_lattice^red / t_star")
    lines.append("")
    lines.append("Exact threshold:")
    lines.append(f"  (lambda_ep * omega_D)_min = {sp.simplify(exprs['prod_threshold'])}")
    lines.append("")
    lines.append("Using t_cross^phys = t_cross,red * t_star, this is equivalently:")
    lines.append(f"  (lambda_ep * omega_D)_min = {sp.simplify(exprs['prod_threshold_cross'])}")
    lines.append(f"  (lambda_ep * omega_D * t_cross^phys)_min = {sp.simplify(exprs['dimless_cross_threshold'])}")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("  For zeta_ep = 1, the lattice drain must turn over about")
    lines.append(f"  {nums['dimless_cross_threshold']:.8f} times during one physical crossing event.")
    lines.append("")

    lines.append("Task 2 — harmonic-trap geometry and stiffness")
    lines.append("-" * 92)
    lines.append("Pure harmonic trap:")
    lines.append(f"  V_lattice(r) = {sp.simplify(exprs['V_lattice'])}")
    lines.append(f"  d ln(V_lattice)/dr = {sp.simplify(exprs['dlnV_lattice'])}")
    lines.append(f"  chi_lambda,lattice = {sp.simplify(exprs['chi_lattice'])}")
    lines.append(f"  Threshold chi_lambda >= 1 implies: {exprs['r_threshold']}")
    lines.append("")
    lines.append("So chi_lambda alone constrains geometry (r <= 2 lambda_phys), not stiffness.")
    lines.append("Evaluated on the earlier reduced turning point:")
    lines.append(f"  r_turn^phys = {sp.simplify(exprs['r_turn_phys_A'])}")
    lines.append(f"  chi_lambda,lattice(r_turn) = {sp.simplify(exprs['chi_on_turn'])}")
    lines.append(f"                              ≈ {nums['chi_on_turn']:.8f}")
    lines.append("")
    lines.append("Minimal extra assumption to obtain k_eff:")
    lines.append("  match the absolute harmonic-trap force to the reduced-model force at r_turn.")
    lines.append("That gives")
    lines.append(f"  k_eff,req = {sp.simplify(exprs['k_eff_required'])}  [eV/Å^2]")
    lines.append("")
    lines.append("If lambda_phys is identified with one-half of an interstitial spacing a_int,")
    lines.append("  lambda_phys = a_int / 2,")
    lines.append("then")
    lines.append(f"  k_eff,req = {sp.simplify(exprs['k_eff_interstitial'])}  [eV/Å^2]")
    lines.append("")
    lines.append("Reference numbers for E_star = 1 eV:")
    lines.append(f"  lambda_phys = 0.5 Å  -> k_eff ≈ {nums['keff_ref_05A_1eV']:.8f} eV/Å^2")
    lines.append(f"  lambda_phys = 1.0 Å  -> k_eff ≈ {nums['keff_ref_1A_1eV']:.8f} eV/Å^2")
    lines.append(f"  lambda_phys = 1.5 Å  -> k_eff ≈ {nums['keff_ref_15A_1eV']:.8f} eV/Å^2")
    lines.append(f"  a_int = 2.0 Å        -> k_eff ≈ {nums['keff_interstitial_2A_1eV']:.8f} eV/Å^2")
    lines.append("")

    lines.append("Task 3 — Korringa depolarization limit")
    lines.append("-" * 92)
    lines.append("Korringa relation and spin-survival condition:")
    lines.append("  T1 * T = Kcorr,    safe if T1 >= t_cross^phys")
    lines.append("Hence")
    lines.append(f"  T_max = {sp.simplify(exprs['Tmax'])}")
    lines.append("")
    lines.append("Using the earlier reduced crossing time:")
    lines.append(f"  T_max = {sp.simplify(exprs['Tmax_reduced_map'])}")
    lines.append("")
    lines.append("So the spin-aligned branch survives only if the operating temperature stays")
    lines.append("below Kcorr divided by the physical crossing time.")
    lines.append("")
    lines.append("Bottom line")
    lines.append("-" * 92)
    lines.append("1) The clean database-comparison target for the lattice drain is")
    lines.append(f"     lambda_ep * omega_D >= {sp.simplify(exprs['prod_threshold_cross'])}")
    lines.append("   or, equivalently,")
    lines.append(f"     lambda_ep * omega_D * t_cross^phys >= {sp.simplify(exprs['dimless_cross_threshold'])}.")
    lines.append("2) The harmonic-trap chi_lambda test by itself fixes only the geometry ratio")
    lines.append("   r_turn / lambda_phys. A stiffness emerges only after force matching, which")
    lines.append("   gives the explicit eV/Å^2 formula above.")
    lines.append("3) The thermal depolarization ceiling is simply")
    lines.append(f"     T_max = {sp.simplify(exprs['Tmax'])}.")
    lines.append("")

    outpath.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    outdir = Path("/mnt/data")
    outdir.mkdir(parents=True, exist_ok=True)

    inputs = PriorRunInputs()
    exprs = derive_symbolic_map(inputs)
    nums = numeric_summary(inputs, exprs)

    script_out = outdir / "material_mapping_condensed_matter_sympy.py"
    # The script is already being written externally; leave a breadcrumb report instead.
    report_out = outdir / "material_mapping_condensed_matter_report.txt"

    write_report(report_out, inputs, exprs, nums)
    print(f"Wrote report to {report_out}")


if __name__ == "__main__":
    main()
