#!/usr/bin/env python3
"""
throat_kappaPV_closure.py

Path-A closure for κ_PV in the reduced throat model (v1 assumptions).

This script does *not* do parameter sweeps. Instead it:

  1) Assumes the brane-healing term is negligible in the ρ-response (Eb << Ew,Ef,Epv).
  2) Uses the power-law energy sectors used in throat_kappa_v1.py:

        Ew(a,ρ)  = Aw * c_s(ρ) / a
        Ef(a,ρ)  = Af / (ρ a^2)
        Epv(a,ρ) = Apv * P(ρ) * V(a)   with V = π Λ a^3
        (Eb ignored)

     with polytrope P = K ρ^n and c_s^2 = dP/dρ = n K ρ^(n-1).
     K is chosen so c_s(rho0) = cs0.

  3) Uses equilibrium (∂F/∂a=0) to eliminate Epv in favor of Ew and Ef (virial relation).
  4) Solves directly for the unique energy partition (x = Ef/Ew) that yields a target κ_PV,
     given κ_rho (from cavitation-mass renormalization) and the EOS exponent n.

Outputs:
  - x = Ef/Ew needed to hit κ_PV target
  - energy fractions (wave, flow, PV)
  - predicted slopes: d ln F / d ln ρ, d ln a / d ln ρ
  - implied β = κ_rho + κ_add + κ_PV
  - if you provide (Aw, Apv, rho0, cs0, Lambda), it also prints the corresponding
    equilibrium a(rho0) and the Af needed (so the system is fully specified without scanning).

This is intended to "close the knobs" for κ_PV work: once (n,Λ,κ_rho,κ_add) are frozen
by earlier papers, κ_PV=3/2 fixes a unique internal partition under Path A.

Example:
  python throat_kappaPV_closure.py --n 5 --Lambda 1.85 --kappa_rho 1 --kappa_add 0.5 --target_kappaPV 1.5 --Aw 1 --Apv 1
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class Inputs:
    # EOS and background
    n: float = 5.0
    rho0: float = 1.0
    cs0: float = 1.0

    # Geometry
    Lambda: float = 1.85

    # 1PN bookkeeping (from the orbital paper decomposition)
    kappa_rho: float = 1.0
    kappa_add: float = 0.5
    target_kappaPV: float = 1.5

    # Coefficients for the reduced energy model (only ratios matter for κ_PV)
    Aw: float = 1.0
    Apv: float = 1.0  # if Epv is literally "P*V", Apv=1 is the natural choice

    # Numerics
    x_min: float = 1e-6
    x_max: float = 1e6
    max_iter: int = 200
    tol: float = 1e-12


def fractions_from_x(x: float) -> Tuple[float, float, float]:
    """
    With Eb neglected, equilibrium implies 3Epv = Ew + 2Ef.
    Let x = Ef/Ew. Then:
        Ew : Ef : Epv = 1 : x : (1+2x)/3
    Total = (4+5x)/3, so fractions are:
        fw = 3/(4+5x)
        ff = 3x/(4+5x)
        fpv = (1+2x)/(4+5x)
    """
    denom = 4.0 + 5.0 * x
    fw = 3.0 / denom
    ff = 3.0 * x / denom
    fpv = (1.0 + 2.0 * x) / denom
    return fw, ff, fpv


def exponents(n: float) -> Tuple[float, float, float]:
    """
    ρ-exponents for each energy sector under v1 assumptions:
      Ew ∝ c_s(ρ) ∝ ρ^((n-1)/2)
      Ef ∝ ρ^{-1}
      Epv ∝ P(ρ) ∝ ρ^n
    """
    a_w = 0.5 * (n - 1.0)
    a_f = -1.0
    a_pv = n
    return a_w, a_f, a_pv


def dlnF_dlnrho_from_x(x: float, n: float) -> float:
    fw, ff, fpv = fractions_from_x(x)
    a_w, a_f, a_pv = exponents(n)
    return a_w * fw + a_f * ff + a_pv * fpv


def kappaPV_from_x(x: float, n: float, kappa_rho: float) -> float:
    # by convention used in the search scripts: κ_PV ≈ dlnF/dlnρ - κ_rho
    return dlnF_dlnrho_from_x(x, n) - kappa_rho


def dlnA_dlnrho_from_x(x: float, n: float) -> float:
    """
    Log-slope of equilibrium radius a(ρ) under Eb neglected.

    From differentiating the virial relation:
      3Epv = Ew + 2Ef

    with Ew ∝ ρ^a_w a^{-1}, Ef ∝ ρ^{-1} a^{-2}, Epv ∝ ρ^n a^{3}.

    Result (after eliminating Epv using the virial relation) becomes:
      d ln a / d ln ρ = - ( -a_w + 2x + n(1+2x) ) / (4 + 10x)
    """
    a_w, _, _ = exponents(n)
    num = (-a_w) + 2.0 * x + n * (1.0 + 2.0 * x)
    den = 4.0 + 10.0 * x
    return - num / den


def solve_x_for_target(inp: Inputs) -> float:
    """
    Solve κ_PV(x) = target_kappaPV on x>0 by bisection in log-space.
    """
    def f(logx: float) -> float:
        x = 10.0 ** logx
        return kappaPV_from_x(x, inp.n, inp.kappa_rho) - inp.target_kappaPV

    lo = math.log10(inp.x_min)
    hi = math.log10(inp.x_max)
    flo = f(lo)
    fhi = f(hi)

    # Expand bracket if needed (should be rare with generous defaults)
    expand = 0
    while flo * fhi > 0 and expand < 30:
        lo -= 1.0
        hi += 1.0
        flo = f(lo)
        fhi = f(hi)
        expand += 1

    if flo * fhi > 0:
        raise RuntimeError(
            "Failed to bracket root for x. "
            f"Try widening --x_min/--x_max. f(lo)={flo:+.3e}, f(hi)={fhi:+.3e}"
        )

    for _ in range(inp.max_iter):
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if abs(fmid) < inp.tol:
            return 10.0 ** mid
        if flo * fmid <= 0:
            hi = mid
            fhi = fmid
        else:
            lo = mid
            flo = fmid

    return 10.0 ** (0.5 * (lo + hi))


def K_from_cs0(n: float, rho0: float, cs0: float) -> float:
    # cs0^2 = n K rho0^(n-1)
    return (cs0 ** 2) / (n * (rho0 ** (n - 1.0)))


def pressure(n: float, rho: float, rho0: float, cs0: float) -> float:
    K = K_from_cs0(n, rho0, cs0)
    return K * (rho ** n)


def equilibrium_scale_a_and_Af(inp: Inputs, x: float) -> Tuple[float, float]:
    """
    Given x=Ef/Ew and chosen (Aw,Apv,rho0,cs0,Lambda,n), compute:

      a0  such that Epv = (1+2x)/3 * Ew at rho0,
      Af  such that Ef/Ew = x at rho0.

    Under v1 definitions:
      Ew0  = Aw * cs0 / a0
      Ef0  = Af / (rho0 a0^2)
      Epv0 = Apv * P(rho0) * π Λ a0^3

    Solve:
      Epv0 = (1+2x)/3 * Ew0
      Ef0/Ew0 = x
    """
    n = inp.n
    rho0 = inp.rho0
    cs0 = inp.cs0
    Aw = inp.Aw
    Apv = inp.Apv
    Lam = inp.Lambda

    # Use K such that cs(rho0)=cs0
    K = K_from_cs0(n, rho0, cs0)
    P0 = K * (rho0 ** n)

    # Epv0 = Apv * P0 * π Λ a^3
    # Ew0  = Aw * cs0 / a
    # Condition Epv0 = ((1+2x)/3) Ew0  =>  Apv P0 π Λ a^4 = ((1+2x)/3) Aw cs0
    a4 = ((1.0 + 2.0 * x) / 3.0) * (Aw * cs0) / (Apv * P0 * math.pi * Lam)
    if a4 <= 0:
        raise RuntimeError("Computed a^4 <= 0 (check inputs).")
    a0 = a4 ** 0.25

    # Ef/Ew = x => (Af/(rho0 a^2)) / (Aw cs0 / a) = x
    # => Af = x * Aw * cs0 * rho0 * a
    Af = x * Aw * cs0 * rho0 * a0

    return a0, Af


def pretty_ratio(r: float) -> str:
    # try to recognize simple rationals
    from fractions import Fraction
    frac = Fraction(r).limit_denominator(512)
    if abs(float(frac) - r) < 5e-6 and frac.denominator <= 128:
        return f"{frac.numerator}/{frac.denominator}"
    return f"{r:.6g}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, default=5.0)
    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--cs0", type=float, default=1.0)
    ap.add_argument("--Lambda", type=float, default=1.85)

    ap.add_argument("--kappa_rho", type=float, default=1.0)
    ap.add_argument("--kappa_add", type=float, default=0.5)
    ap.add_argument("--target_kappaPV", type=float, default=1.5)

    ap.add_argument("--Aw", type=float, default=1.0)
    ap.add_argument("--Apv", type=float, default=1.0)

    ap.add_argument("--x_min", type=float, default=1e-6)
    ap.add_argument("--x_max", type=float, default=1e6)

    args = ap.parse_args()
    inp = Inputs(
        n=args.n,
        rho0=args.rho0,
        cs0=args.cs0,
        Lambda=args.Lambda,
        kappa_rho=args.kappa_rho,
        kappa_add=args.kappa_add,
        target_kappaPV=args.target_kappaPV,
        Aw=args.Aw,
        Apv=args.Apv,
        x_min=args.x_min,
        x_max=args.x_max,
    )

    x = solve_x_for_target(inp)

    fw, ff, fpv = fractions_from_x(x)
    dlnF = dlnF_dlnrho_from_x(x, inp.n)
    kPV = dlnF - inp.kappa_rho
    beta = inp.kappa_rho + inp.kappa_add + kPV
    dlnA = dlnA_dlnrho_from_x(x, inp.n)

    a0, Af = equilibrium_scale_a_and_Af(inp, x)

    # energy breakdown at rho0
    K = K_from_cs0(inp.n, inp.rho0, inp.cs0)
    P0 = K * (inp.rho0 ** inp.n)
    V0 = math.pi * inp.Lambda * (a0 ** 3)

    Ew0 = inp.Aw * inp.cs0 / a0
    Ef0 = Af / (inp.rho0 * a0 ** 2)
    Epv0 = inp.Apv * P0 * V0
    Etot = Ew0 + Ef0 + Epv0

    print("=== Path-A κ_PV closure (Eb neglected) ===")
    print(f"inputs: n={inp.n:g}, Lambda={inp.Lambda:g}, rho0={inp.rho0:g}, cs0={inp.cs0:g}")
    print(f"1PN ledger: kappa_rho={inp.kappa_rho:g}, kappa_add={inp.kappa_add:g}, target_kappaPV={inp.target_kappaPV:g}")
    print()
    print("---- Derived universal partition ----")
    print(f"x = Ef/Ew = {x:.10g}   (~ {pretty_ratio(x)})")
    print(f"fractions: wave={fw:.6f}, flow={ff:.6f}, PV={fpv:.6f} (sum={fw+ff+fpv:.6f})")
    print()
    print("---- Predicted response slopes at rho0 ----")
    print(f"d ln F / d ln rho = {dlnF:.10g}   (~ {pretty_ratio(dlnF)})")
    print(f"κ_PV_est          = {kPV:.10g}")
    print(f"β_est             = {beta:.10g}   (β = κ_rho + κ_add + κ_PV)")
    print(f"d ln a / d ln rho = {dlnA:.10g}   (~ {pretty_ratio(dlnA)})")
    print()
    print("---- One concrete (a0, Af) closure given (Aw,Apv) ----")
    print(f"chosen: Aw={inp.Aw:g}, Apv={inp.Apv:g}")
    print(f"implied: a0={a0:.10g},  Af={Af:.10g}")
    print()
    print("check at rho0:")
    print(f"  Ew0={Ew0:.10g}, Ef0={Ef0:.10g}, Epv0={Epv0:.10g}, Etot={Etot:.10g}")
    print(f"  Ef/Ew={Ef0/Ew0:.10g}, Epv/Ew={Epv0/Ew0:.10g}")
    print(f"  fractions (computed) = ({Ew0/Etot:.6f}, {Ef0/Etot:.6f}, {Epv0/Etot:.6f})")

    # helpful: show the "nice integer ratio" if close
    # scale Ew:Ef:Epv
    rE = np.array([Ew0, Ef0, Epv0], dtype=float)
    rE = rE / np.min(rE)
    # attempt to integerize
    from fractions import Fraction
    ints = [Fraction(v).limit_denominator(64) for v in rE]
    if all(abs(float(f) - v) < 5e-5 for f, v in zip(ints, rE)):
        # common denom
        den = np.lcm.reduce([f.denominator for f in ints])
        num = [int(f.numerator * (den // f.denominator)) for f in ints]
        g = math.gcd(math.gcd(num[0], num[1]), num[2])
        num = [k // g for k in num]
        print()
        print(f"  Ew:Ef:Epv ≈ {num[0]}:{num[1]}:{num[2]} (integer ratio)")

    print()
    print("NOTE:")
    print("  If your parameter search is producing many 'solutions', that's expected: once Eb is negligible,")
    print("  κ_PV matching fixes only the *dimensionless partition* (x and the fractions). Apv (and Aw) set")
    print("  the overall scale (a0, total energy) but do not change κ_PV. To remove degeneracy, decide what")
    print("  Apv physically *means* (PV vs enthalpy vs internal-energy), and fix it to an O(1) value.")


if __name__ == "__main__":
    main()

"""
Run with --n 5 --Lambda 1.85 --kappa_rho 1 --kappa_add 0.5 --target_kappaPV 1.5 --Aw 1 --Apv 1

Output:
=== Path-A κ_PV closure (Eb neglected) ===
inputs: n=5, Lambda=1.85, rho0=1, cs0=1
1PN ledger: kappa_rho=1, kappa_add=0.5, target_kappaPV=1.5

---- Derived universal partition ----
x = Ef/Ew = 0.1818181818   (~ 2/11)
fractions: wave=0.611111, flow=0.111111, PV=0.277778 (sum=1.000000)

---- Predicted response slopes at rho0 ----
d ln F / d ln rho = 2.5   (~ 5/2)
κ_PV_est          = 1.5
β_est             = 3   (β = κ_rho + κ_add + κ_PV)
d ln a / d ln rho = -0.890625   (~ -57/64)

---- One concrete (a0, Af) closure given (Aw,Apv) ----
chosen: Aw=1, Apv=1
implied: a0=0.7907813725,  Af=0.1437784314

check at rho0:
  Ew0=1.264572023, Ef0=0.229922186, Epv0=0.574805465, Etot=2.069299674
  Ef/Ew=0.1818181818, Epv/Ew=0.4545454545
  fractions (computed) = (0.611111, 0.111111, 0.277778)

  Ew:Ef:Epv ≈ 11:2:5 (integer ratio)
"""
