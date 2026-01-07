#!/usr/bin/env python3
"""
step4_outer_match.py — Paper 7 Step 4 (Python port)

Inner: cylinder throat DtN (mode-resolved) with scale-aware damping.
Outer: 3D spherical Helmholtz DtN impedance (radiating or conservative-static).
Selection: fixed-q resonance-family sampling with controlled detunes.
Sanity: outer static limits, passivity sign proxy, PN locality check (monopole only).

Examples:
  python step4_outer_match.py --cli
  python step4_outer_match.py --cli --scan
  python step4_outer_match.py --cli --scan --robust
  python step4_outer_match.py --cli --scan --robust --export step4_scan.tsv

Notes:
- Uses SciPy special functions for Bessel and spherical Hankel functions.
- Avoids Mathematica's symbolic-series blowups by keeping PN check numeric + compact.
"""

from __future__ import annotations

import argparse
import math
import cmath
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple, Optional

import numpy as np
import scipy.special as sp
import sympy as sy


# ----------------------------
# Formatting helpers
# ----------------------------
def num(x, prec=12) -> str:
    if x is None:
        return "None"
    if isinstance(x, str):
        return x
    try:
        if isinstance(x, complex):
            return f"{x.real:.{prec}g}{'+' if x.imag >= 0 else '-'}{abs(x.imag):.{prec}g}j"
        return f"{float(x):.{prec}g}"
    except Exception:
        return str(x)


def tsv_row(items: Sequence[object]) -> str:
    return "\t".join(str(it) for it in items)


# ----------------------------
# Model parameters
# ----------------------------
@dataclass
class Params:
    # geometry/wave
    c: float = 1.0
    a: float = 1.0
    Lovera: float = 1.85

    # inner BCs
    wallBC: str = "Neumann"     # Neumann or Dirichlet
    bottomBC: str = "Neumann"   # Neumann or Dirichlet

    # inner mode
    mPick: int = 0
    nPick: int = 0

    # damping (scale-aware)
    gammaMode: str = "RelativeOmega"  # RelativeOmega or Absolute
    deltaRel: float = 0.02
    gammaFloor: float = 1e-6
    gammaAbs: float = 0.01

    # sampling
    qFixedList: Tuple[int, ...] = (0, 1, 2)
    detuneRel: float = 0.02
    detuneSigns: Tuple[int, ...] = (-1, +1)

    # outer model
    outerMode: str = "Radiating"  # Radiating or ConservativeStatic
    lList: Tuple[int, ...] = (0,)
    wOuter: Dict[int, float] = None
    w2Model: str = "Constant"     # Constant or KaSquared
    w2Base: float = 0.01

    # mismatch functional
    mismatchMode: str = "RelL2"   # RelL2 or AbsL2
    epsNorm: float = 1e-12

    # scan
    LoveraMin: float = 1.20
    LoveraMax: float = 4.0
    LoveraStep: float = 0.05

    # robustness
    scaleList: Tuple[float, ...] = (0.95, 1.00, 1.05)

    # PN locality check (monopole only)
    doPN: bool = True
    pnOrder: int = 8
    pnBandFrac: float = 0.25
    pnTol: float = 1e-4
    pnPrintCoeffs: bool = True


# ----------------------------
# Core math
# ----------------------------
class Step4Model:
    def __init__(self, p: Params):
        self.p = p
        if p.wOuter is None:
            p.wOuter = {0: 1.0, 2: 0.0}
        self._root_cache: Dict[Tuple[str, int, int], float] = {}

    # ---- Bessel roots ----
    def radial_root(self, wallBC: str, m: int, n: int) -> float:
        key = (wallBC, m, n)
        if key in self._root_cache:
            return self._root_cache[key]

        bc = wallBC.lower()
        if bc.startswith("dir"):
            if n < 1:
                raise ValueError("Dirichlet wall: n must be >= 1")
            x = float(sp.jn_zeros(m, n)[-1])
        else:
            # Neumann: zeros of derivative J_m'
            if m == 0 and n == 0:
                x = 0.0
            else:
                if m == 0:
                    if n < 1:
                        raise ValueError("Neumann wall with m=0: n must be 0 (special) or >=1")
                    x = float(sp.jnp_zeros(0, n)[-1])      # nth positive root
                else:
                    if n < 0:
                        raise ValueError("Neumann wall: n must be >=0")
                    x = float(sp.jnp_zeros(m, n + 1)[-1])  # (n+1)th positive root

        self._root_cache[key] = x
        return x

    def lam(self, m: int, n: int) -> float:
        return self.radial_root(self.p.wallBC, m, n) / self.p.a

    # ---- Complex omega / branch control ----
    def omega_eff(self, omega: float, deltaRel: Optional[float] = None) -> complex:
        d = self.p.deltaRel if deltaRel is None else deltaRel
        if self.p.gammaMode == "RelativeOmega":
            return omega * (1.0 + 1j * d) + 1j * self.p.gammaFloor
        if self.p.gammaMode == "Absolute":
            return omega + 1j * self.p.gammaAbs
        return omega * (1.0 + 1j * d) + 1j * self.p.gammaFloor

    @staticmethod
    def sqrt_branch(z: complex) -> complex:
        s = cmath.sqrt(z)
        return -s if s.imag < 0 else s

    def kappa(self, omega: float, m: int, n: int, Lovera: float, deltaRel: Optional[float] = None) -> complex:
        w = self.omega_eff(omega, deltaRel=deltaRel)
        k = w / self.p.c
        lam = self.lam(m, n)
        return self.sqrt_branch(k * k - lam * lam)

    def Z_in(self, omega: float, m: int, n: int, Lovera: float, deltaRel: Optional[float] = None) -> complex:
        L = Lovera * self.p.a
        kap = self.kappa(omega, m, n, Lovera=Lovera, deltaRel=deltaRel)
        if self.p.bottomBC == "Neumann":
            return -kap * cmath.tan(kap * L)
        if self.p.bottomBC == "Dirichlet":
            return kap / cmath.tan(kap * L)  # kap*cot
        raise ValueError("bottomBC must be Neumann or Dirichlet")

    # ---- Resonance predictor ----
    def resonance_omega(self, m: int, n: int, q: int, Lovera: float) -> float:
        L = Lovera * self.p.a
        lam = self.lam(m, n)
        if self.p.bottomBC == "Neumann":
            kap = (q + 0.5) * math.pi / L
        else:
            kap = q * math.pi / L
        return self.p.c * math.sqrt(lam * lam + kap * kap)

    def sample_freqs(self, m: int, n: int, Lovera: float, detuneRel: Optional[float] = None) -> List[float]:
        det = self.p.detuneRel if detuneRel is None else detuneRel
        freqs: List[float] = []
        for q in self.p.qFixedList:
            wres = self.resonance_omega(m, n, q, Lovera)
            for sgn in self.p.detuneSigns:
                w = wres * (1.0 + det * float(sgn))
                if w > 0:
                    freqs.append(w)
        return freqs

    # ---- Outer DtN ----
    @staticmethod
    def hankel1(l: int, x: complex) -> complex:
        return sp.spherical_jn(l, x) + 1j * sp.spherical_yn(l, x)

    @staticmethod
    def hankel1_prime(l: int, x: complex) -> complex:
        return sp.spherical_jn(l, x, derivative=True) + 1j * sp.spherical_yn(l, x, derivative=True)

    def Z_out_l(self, omega: float, l: int, deltaRel: Optional[float] = None) -> complex:
        if self.p.outerMode == "ConservativeStatic":
            return -(l + 1) / self.p.a

        w = self.omega_eff(omega, deltaRel=deltaRel)
        k = w / self.p.c
        x = k * self.p.a
        h = self.hankel1(l, x)
        hp = self.hankel1_prime(l, x)
        return k * hp / h

    def w_l(self, l: int, omega: float) -> float:
        if l == 0:
            return float(self.p.wOuter.get(0, 0.0))
        if l == 2 and self.p.w2Model == "KaSquared":
            ka = (omega / self.p.c) * self.p.a
            return float(self.p.w2Base * (ka * ka))
        return float(self.p.wOuter.get(l, 0.0))

    def Z_out_eff(self, omega: float, deltaRel: Optional[float] = None) -> complex:
        tot = 0.0 + 0.0j
        for l in self.p.lList:
            w = self.w_l(l, omega)
            if w == 0.0:
                continue
            tot += w * self.Z_out_l(omega, l, deltaRel=deltaRel)
        return tot

    # ---- Mismatch ----
    def mismatch(self, Lovera: float, detuneRel: Optional[float] = None, deltaRel: Optional[float] = None) -> float:
        m, n = self.p.mPick, self.p.nPick
        freqs = self.sample_freqs(m, n, Lovera, detuneRel=detuneRel)
        if not freqs:
            return float("nan")

        Zin = np.array([self.Z_in(w, m, n, Lovera, deltaRel=deltaRel) for w in freqs], dtype=complex)
        Zout = np.array([self.Z_out_eff(w, deltaRel=deltaRel) for w in freqs], dtype=complex)
        diffs = np.abs(Zin - Zout) ** 2

        if self.p.mismatchMode == "AbsL2":
            return float(np.mean(diffs))

        denom = float(np.mean(np.abs(Zout) ** 2) + self.p.epsNorm)
        return float(np.mean(diffs) / denom)

    # ---- Sanity ----
    def passivity_proxy(self, omegas: Sequence[float], Lovera: float) -> List[Tuple[float, float, float]]:
        m, n = self.p.mPick, self.p.nPick
        out = []
        for w in omegas:
            Zin = self.Z_in(w, m, n, Lovera)
            out.append((w, w * Zin.imag, -w * Zin.imag))
        return out

    def pn_locality_check(self) -> Dict[str, object]:
        """Conservative monopole check only (lambda=0)."""
        if not self.p.doPN:
            return {"enabled": False, "reason": "disabled"}

        lam = self.lam(self.p.mPick, self.p.nPick)
        if abs(lam) > 1e-14:
            return {"enabled": False, "reason": "lambda != 0 (not monopole)"}

        L = self.p.Lovera * self.p.a
        c = self.p.c

        # omega0 definition consistent with Mathematica script
        if self.p.bottomBC == "Neumann":
            omega0 = (math.pi / 2.0) * c / L
        else:
            omega0 = math.pi * c / L

        omega_max = self.p.pnBandFrac * omega0

        # exact conservative function
        def Z_exact(w: float) -> float:
            k = w / c
            if self.p.bottomBC == "Neumann":
                return float(-k * math.tan(k * L))
            return float(k / math.tan(k * L))

        # series via sympy (compact)
        w = sy.Symbol("w", real=True)
        if self.p.bottomBC == "Neumann":
            Zs = -(w / c) * sy.tan((w / c) * L)
        else:
            Zs = (w / c) / sy.tan((w / c) * L)  # cot

        ser_expr = sy.series(Zs, w, 0, self.p.pnOrder + 1).removeO().expand()

        # coefficients (only nonzero ones)
        coeffs = []
        for pwr in range(0, self.p.pnOrder + 1):
            coeff = sy.N(sy.expand(ser_expr).coeff(w, pwr), 30)
            if abs(complex(coeff)) > 1e-16:
                coeffs.append((pwr, float(coeff)))

        # numeric error check
        grid = np.linspace(1e-6 * omega0, omega_max, 120)
        errs = []
        f_ser = sy.lambdify(w, ser_expr, "numpy")
        for wi in grid:
            ze = Z_exact(float(wi))
            zs = float(f_ser(float(wi)))
            if abs(ze) < 1e-14:
                errs.append(0.0)
            else:
                errs.append(abs(ze - zs) / abs(ze))
        max_err = float(np.max(errs))
        return {
            "enabled": True,
            "pass": max_err <= self.p.pnTol,
            "maxRelErr": max_err,
            "omega0": omega0,
            "omegaMaxTest": omega_max,
            "coeffs": coeffs,
        }

    # ---- Scans ----
    def run_scan(self) -> Tuple[List[Tuple[float, float, int, float, float]], Tuple[float, float]]:
        """Return rows + (bestL, bestJ)."""
        grid = np.arange(self.p.LoveraMin, self.p.LoveraMax + 1e-12, self.p.LoveraStep)
        rows = []
        best = (float("nan"), float("inf"))

        m, n = self.p.mPick, self.p.nPick
        for Lr in grid:
            freqs = self.sample_freqs(m, n, Lr)
            if not freqs:
                J = float("nan")
                rows.append((float(Lr), J, 0, float("nan"), float("nan")))
                continue

            J = self.mismatch(float(Lr))
            rows.append((float(Lr), float(J), len(freqs), min(freqs), max(freqs)))
            if np.isfinite(J) and J < best[1]:
                best = (float(Lr), float(J))

        return rows, best

    def run_robust(self) -> List[Tuple[float, float, float]]:
        """Scale detuneRel and deltaRel together and rescan."""
        base_det = self.p.detuneRel
        base_del = self.p.deltaRel
        results = []

        for s in self.p.scaleList:
            det = base_det * s
            dlt = base_del * s

            # compute best over grid using scaled params but without mutating global state
            grid = np.arange(self.p.LoveraMin, self.p.LoveraMax + 1e-12, self.p.LoveraStep)
            bestL, bestJ = float("nan"), float("inf")
            for Lr in grid:
                try:
                    J = self.mismatch(float(Lr), detuneRel=det, deltaRel=dlt)
                except Exception:
                    J = float("nan")
                if np.isfinite(J) and J < bestJ:
                    bestL, bestJ = float(Lr), float(J)

            results.append((float(s), float(bestL), float(bestJ)))

        return results


# ----------------------------
# CLI
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", action="store_true", help="print CLI summary")
    ap.add_argument("--scan", action="store_true", help="run L/a scan")
    ap.add_argument("--robust", action="store_true", help="run robustness scan (requires --scan for context)")
    ap.add_argument("--export", type=str, default="", help="export TSV path (requires --scan)")
    ap.add_argument("--fast", action="store_true", help="coarser scan step (0.1)")

    # quick knobs (so you don't have to edit file for common changes)
    ap.add_argument("--m", type=int, default=None)
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--wallBC", type=str, default=None, choices=["Neumann", "Dirichlet"])
    ap.add_argument("--bottomBC", type=str, default=None, choices=["Neumann", "Dirichlet"])
    ap.add_argument("--outerMode", type=str, default=None, choices=["Radiating", "ConservativeStatic"])
    ap.add_argument("--q", type=str, default=None, help="q list like '0,1,2'")
    ap.add_argument("--detune", type=float, default=None)
    ap.add_argument("--delta", type=float, default=None)

    ap.add_argument("--l", type=str, default=None, help="outer multipoles list like '0' or '0,2' or '2'")
    ap.add_argument("--wouter", type=str, default=None, help="outer weights like '0:1,2:0.01' (keys are l)")
    ap.add_argument("--w2Model", type=str, default=None, choices=["Constant", "KaSquared"])
    ap.add_argument("--w2Base", type=float, default=None)
    ap.add_argument("--mismatch", type=str, default=None, choices=["RelL2", "AbsL2"])

    args = ap.parse_args()

    p = Params()
    if args.fast:
        p.LoveraStep = 0.10
    if args.m is not None:
        p.mPick = args.m
    if args.n is not None:
        p.nPick = args.n
    if args.wallBC is not None:
        p.wallBC = args.wallBC
    if args.bottomBC is not None:
        p.bottomBC = args.bottomBC
    if args.outerMode is not None:
        p.outerMode = args.outerMode
    if args.q is not None:
        p.qFixedList = tuple(int(x.strip()) for x in args.q.split(",") if x.strip() != "")
    if args.detune is not None:
        p.detuneRel = float(args.detune)
    if args.delta is not None:
        p.deltaRel = float(args.delta)


    # outer multipoles / weights (optional)
    if args.l is not None:
        p.lList = tuple(int(x.strip()) for x in args.l.split(",") if x.strip() != "")
    if args.wouter is not None:
        d: Dict[int, float] = {}
        for part in args.wouter.split(","):
            part = part.strip()
            if not part:
                continue
            if ":" not in part:
                raise ValueError(f"--wouter entries must be like 'l:weight' (got '{part}')")
            k, v = part.split(":", 1)
            d[int(k.strip())] = float(v.strip())
        p.wOuter = d
    if args.w2Model is not None:
        p.w2Model = args.w2Model
    if args.w2Base is not None:
        p.w2Base = float(args.w2Base)
    if args.mismatch is not None:
        p.mismatchMode = args.mismatch
    # auto-default nPick similar to Mathematica behavior if not provided explicitly
    # (you can still override with --n)
    if args.n is None:
        p.nPick = 1 if p.wallBC == "Dirichlet" else 0

    model = Step4Model(p)

    # Guard: if all selected multipoles have zero weight, Z_out_eff == 0 and RelL2 will blow up.
    eff_w = []
    for ell in p.lList:
        w = float(p.wOuter.get(ell, 0.0))
        if ell == 2 and p.w2Model == "KaSquared":
            # weight becomes omega-dependent; treat base as nonzero if w2Base>0
            w = float(p.w2Base)
        eff_w.append(w)
    if all(abs(w) == 0.0 for w in eff_w):
        raise ValueError(f"All outer weights are zero for lList={p.lList}. "
                         f"Set --wouter (e.g. --wouter '2:1,0:0') or change lList.")


    if args.cli:
        L = p.Lovera * p.a
        lam = model.lam(p.mPick, p.nPick)
        print("[CLI mode enabled]")
        print("=== Step 4: Inner/Outer DtN Matching Summary (Python) ===")
        print(f"Inner: wallBC={p.wallBC} bottomBC={p.bottomBC}  a={num(p.a,6)}  L={num(L,6)}  (L/a)={num(L/p.a,6)}")
        print(f"Picked inner mode: (m,n)=({p.mPick},{p.nPick})  lambda={num(lam,12)}")
        print(f"Damping: gammaMode={p.gammaMode}  deltaRel={num(p.deltaRel,6)}  gammaFloor={num(p.gammaFloor,6)}  gammaAbs={num(p.gammaAbs,6)}")
        print(f"Sampling: qFixedList={p.qFixedList}  detuneRel={num(p.detuneRel,6)}  detuneSigns={p.detuneSigns}")
        print(f"Outer: mode={p.outerMode}  lList={p.lList}  wOuter={p.wOuter}  w2Model={p.w2Model}  w2Base={num(p.w2Base,6)}")
        print(f"Mismatch: mode={p.mismatchMode}")

        print("--- Outer sanity (static limits) ---")
        print(f"Zout0_static = {num(-1.0/p.a,12)}   (expected -1/a)")
        print(f"Zout2_static = {num(-3.0/p.a,12)}   (expected -3/a)")

        print("--- Passivity sign proxy (sample ω points; shows convention) ---")
        samples = [0.2, 0.5, 1.0, 2.0, 3.0]
        print(tsv_row(["ω", "ω Im(Zin)", "-ω Im(Zin)"]))
        for w, p1, p2 in model.passivity_proxy(samples, p.Lovera):
            print(tsv_row([num(w,6), num(p1,12), num(p2,12)]))

        pn = model.pn_locality_check()
        print("--- PN locality check (conservative, monopole only if λ=0) ---")
        if pn.get("enabled"):
            print(f"Result: {'PASS' if pn['pass'] else 'FAIL'}  maxRelErr={num(pn['maxRelErr'],12)}  ωmaxTest={num(pn['omegaMaxTest'],12)}  ω0={num(pn['omega0'],12)}")
            if p.pnPrintCoeffs:
                print("Series coefficients (power p -> coeff):")
                print(tsv_row(["p", "coeff"]))
                for pwr, coeff in pn["coeffs"]:
                    if pwr == 0:
                        continue
                    print(tsv_row([pwr, num(coeff,18)]))
        else:
            print(f"(skipped) reason: {pn.get('reason')}")

    if args.scan:
        rows, best = model.run_scan()
        if args.cli:
            print("--- Scan table: L/a, J(L/a), Nω, ωmin, ωmax ---")
            print(tsv_row(["L/a", "J", "Nω", "ωmin", "ωmax"]))
        for Lr, J, nW, wmin, wmax in rows:
            if args.cli:
                print(tsv_row([num(Lr,4), "NaN" if not np.isfinite(J) else num(J,12), nW, num(wmin,10), num(wmax,10)]))
        if args.cli:
            print("--- Best point (scan) ---")
            print(f"best L/a = {num(best[0],6)}   J = {num(best[1],12)}   qFixedList={p.qFixedList}")

        if args.export:
            path = args.export
            with open(path, "w", encoding="utf-8") as f:
                f.write(tsv_row(["Lovera", "J", "Nω", "ωmin", "ωmax"]) + "\n")
                for Lr, J, nW, wmin, wmax in rows:
                    f.write(tsv_row([Lr, J, nW, wmin, wmax]) + "\n")
            if args.cli:
                print(f"[exported] {path}")

    if args.robust:
        res = model.run_robust()
        if args.cli:
            print("--- Robustness (scale detuneRel and deltaRel together) ---")
            print(tsv_row(["scale s", "best L/a", "best J"]))
            for s, bestL, bestJ in res:
                print(tsv_row([num(s,4), num(bestL,8), "NaN" if not np.isfinite(bestJ) else num(bestJ,12)]))
            finiteLs = [r[1] for r in res if np.isfinite(r[1])]
            if finiteLs:
                print(f"Drift Δ(L/a) = {num(max(finiteLs) - min(finiteLs),8)}")


if __name__ == "__main__":
    main()
