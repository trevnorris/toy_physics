(* ::Package:: *)

(*
========================================================
Paper VIII (referee add-on v3): moving defect branches + covariant propagator (KK tower)
========================================================

Purpose (referee-walkable):
  (1) Make explicit that the brane-to-brane response depends only on the Lorentz
      scalar k^2 (in momentum space), hence has *exact* 3+1 Lorentz invariance
      (no preferred-frame artifacts) even when KK corrections are included.
  (2) Provide a compact expression for the brane effective propagator as a KK sum:
         D_eff(k^2) = mu0 * Sum_n c_n / (k^2 + m_n^2 + i eps)
      with the Gaussian-localization spectrum:
         m_n^2 = 2 n / lambdaConf^2,
         c_n = f_n(0)^2 / ∫dw Z(w) f_n(w)^2,
      (odd n have c_n = 0).
  (3) Connect to configuration-space intuition:
      - static limit gives Coulomb + Yukawa tower
      - retarded solutions exist mode-by-mode (massless + massive Lorentz-covariant Green functions)

What this establishes for the paper:
  - The EM sector is exactly Lorentz invariant in the *brane* (t,x,y,z) directions:
    Z(w) breaks boosts mixing w with the brane, but preserves the full SO(1,3) on-brane.
  - Therefore, moving-defect-branch fields computed from this sector cannot contain a preferred
    direction in the brane (no anisotropic v/c artifacts from Z(w) itself).
  - Deviations from standard EM appear only as Lorentz-covariant Yukawa/KK corrections.

Notes:
  - We keep everything symbolic; this is a "structure" proof for referees.
  - No symbol names contain underscores ( _ ).

Conventions:
  - Brane metric signature: (-,+,+,+)
  - k2 := k_mu k^mu = -k0^2 + kx^2 + ky^2 + kz^2
========================================================
*)

Print["\n========================================"];
Print["Paper VIII add-on v3: moving defect branches + covariant propagator (KK tower)"];
Print["========================================\n"];

ClearAll["Global`*"];

(* ---------- Assumptions + helpers ---------- *)
$Assumptions = {
  Element[{mu0, lambdaConf, betaB, k0, kx, ky, kz, epsReg, qStar, eStar, etaQ}, Reals],
  mu0 > 0, lambdaConf > 0, eStar > 0, etaQ^2 == 1,
  0 <= betaB < 1,
  epsReg > 0
};

Simp[expr_] := FullSimplify[Together[expr], Assumptions -> $Assumptions];
ZeroQAssume[expr_] := TrueQ[FullSimplify[expr == 0, Assumptions -> $Assumptions]];

(* ---------- (A) Brane Lorentz invariance: k^2 is invariant ---------- *)

Print["--- (A) Brane Lorentz invariance: k^2 is invariant under boosts ---\n"];

eta4 = DiagonalMatrix[{-1, 1, 1, 1}];

gammaB[bb_] := 1/Sqrt[1 - bb^2];

(* Boost along x: contravariant vector transform k'^mu = Lambda^mu_nu k^nu *)
LambdaUp = {
  {gammaB[betaB], -betaB*gammaB[betaB], 0, 0},
  {-betaB*gammaB[betaB], gammaB[betaB], 0, 0},
  {0, 0, 1, 0},
  {0, 0, 0, 1}
};

kVec = {k0, kx, ky, kz};
kPrime = Simp[LambdaUp . kVec];

k2 = Simp[kVec . (eta4 . kVec)];
k2p = Simp[kPrime . (eta4 . kPrime)];

Print["k2 = -k0^2 + kx^2 + ky^2 + kz^2 -> ", k2];
Print["k2' after boost -> ", k2p];
Print["Check k2' - k2 (should be 0): ", Simp[k2p - k2]];

If[ZeroQAssume[k2p - k2],
  Print["OK: k^2 is invariant under the brane Lorentz boost."],
  Print["WARNING: k^2 invariance check failed (unexpected)."]
];

(* ---------- (B) Gaussian-localized KK data (spectrum + brane couplings) ---------- *)

Print["\n--- (B) Gaussian-localized KK data: m_n^2 and c_n ---\n"];

Zint = Simp[Integrate[Exp[-w^2/lambdaConf^2], {w, -Infinity, Infinity}]];
mu0eff = Simp[mu0/Zint];

Print["Zint -> ", Zint];
Print["mu0eff -> ", mu0eff];
Print["Fixed defect branch relation: qStar = etaQ * eStar"];
Print["Canonical brane charge relation: qEff = qStar / Sqrt[Zint]"];

(* KK spectrum for the Gaussian Z(w) = exp(-w^2/lambda^2) (as validated in add-on v2):
     f_n(w) = HermiteH[n, w/lambdaConf],    m_n^2 = 2n/lambdaConf^2
   and (standard orthogonality) normClosed[n] = ∫dw Z f_n^2 = 2^n lambdaConf Sqrt[Pi] n!
*)
mn2[n_Integer?NonNegative] := 2 n/lambdaConf^2;
normClosed[n_Integer?NonNegative] := 2^n lambdaConf Sqrt[Pi] Factorial[n];
f0[n_Integer?NonNegative] := HermiteH[n, 0];

cCoef[n_Integer?NonNegative] := Simp[(f0[n]^2)/normClosed[n]];

Print["m_n^2 = 2 n / lambdaConf^2"];
Print["c_n = f_n(0)^2 / norm_n, with norm_n = 2^n lambdaConf Sqrt[Pi] n!"];

(* Spot-check: odd n should vanish *)
oddCheck = Table[{n, Simp[cCoef[n]]}, {n, 1, 9, 2}];
Print["Odd-n c_n should be 0. Spot-check (n,c_n): ", oddCheck];
Print["Interpretation: odd modes decouple from centered brane sources, but odd-w core structure need not be absent microscopically."];

(* ---------- (C) Brane-to-brane propagator depends only on k^2 ---------- *)

Print["\n--- (C) Brane-to-brane propagator D_eff(k^2) depends only on k^2 ---\n"];

(* Effective (Feynman) propagator kernel in momentum space:
     D_eff(k^2) = mu0 * Sum_n c_n / (k^2 + m_n^2 + i eps)
   For classical retarded solutions, replace the i eps prescription accordingly;
   the key point here is Lorentz covariance: dependence only on k^2.
*)
Deff[k2in_, nmax_Integer?NonNegative] := mu0 * Sum[cCoef[n]/(k2in + mn2[n] + I epsReg), {n, 0, nmax}];

Print["Example truncations (nmax=0..6), showing dependence only on k2:"];
Do[
  Print["  nmax=", nmax, " : Deff(k2) = ", Simp[Deff[k2sym, nmax]]],
  {nmax, 0, 6}
];

(* Explicit invariance check under boost: Deff(k'^2) - Deff(k^2) == 0 *)
invCheck = Simp[Deff[k2p, 6] - Deff[k2, 6]];
Print["\nCheck Deff(k2') - Deff(k2) at nmax=6 (should be 0): ", invCheck];

If[ZeroQAssume[invCheck],
  Print["OK: D_eff is a function of the Lorentz scalar k^2 only (manifest brane Lorentz invariance)."],
  Print["WARNING: propagator invariance check failed (unexpected)."]
];

(* ---------- (D) Configuration-space intuition: static potential is isotropic ---------- *)

Print["\n--- (D) Static limit (k0=0): Coulomb + Yukawa tower, isotropic in r ---\n"];

(* In the static limit, each massive mode gives Yukawa exp(-m r)/(4 pi r).
   The brane potential for a fixed defect branch qStar becomes:
     A0(r) = (mu0 qStar)/(4 pi r) * Sum_n c_n exp(-m_n r)
   with c_0 = 1/Zint, so the leading term is Coulomb with mu0eff.
*)
mn[n_Integer?NonNegative] := Sqrt[mn2[n]];

A0tower[r_, nmax_Integer?NonNegative] := (mu0*qStar)/(4*Pi*r) * Sum[cCoef[n] Exp[-mn[n] r], {n, 0, nmax}];

A0zero[r_] := (mu0eff*qStar)/(4*Pi*r);

Print["A0_zero(r) = mu0eff qStar/(4 pi r) -> ", A0zero[r]];
Print["First correction (keep n=0,2): A0(r) = ", Simp[A0tower[r, 2]]];
Print["Relative correction Delta(r) = (A0 - A0zero)/A0zero (nmax=2): ", Simp[(A0tower[r, 2] - A0zero[r])/A0zero[r]]];

Print["Note: all static corrections depend only on r (no angular dependence), hence no brane anisotropy."];

(* ---------- (E) Retarded solutions exist mode-by-mode (structure) ---------- *)

Print["\n--- (E) Retarded solutions exist mode-by-mode (structure statement) ---\n"];

Print["Each KK mode on the brane is Lorentz covariant:"];
Print["  - n=0 behaves like standard Maxwell (massless) with retarded Green G0(x-x')."];
Print["  - n>0 modes behave like massive (Proca-type) vectors with Klein-Gordon operator (Box + m_n^2) on each component."];
Print["  The 4D retarded Green function for (Box + m^2) depends only on the invariant interval (t^2 - r^2) and theta(t),"];
Print["  so the sum over modes preserves exact brane Lorentz covariance."];

Print["\nCanonical (3+1) retarded Green functions (c=1) for reference (printed as strings):"];
Print["  Massless wave:   G0_ret(t,r) = theta(t) * delta(t-r) / (4 pi r)"];
Print["  Massive KG:      Gm_ret(t,r) = theta(t) * [ delta(t^2-r^2)/(2 pi) - (m/(4 pi)) theta(t^2-r^2) J1(m sqrt(t^2-r^2))/sqrt(t^2-r^2) ]"];

Print["\nMoving-defect-branch potentials can be written covariantly as:"];
Print["  A_mu(x) = ∫ d^4x' D_eff_ret(x-x') J_mu(x')"];
Print["and in momentum space:"];
Print["  A_mu(k) = D_eff(k^2) J_mu(k),  with D_eff depending only on k^2."];
Print["This is the cleanest way to see there is no preferred-frame artifact in the brane EM sector for a fixed charge branch."];

Print["\n========== End add-on v3 ==========\n"];

(*"
Output:

========================================
Paper VIII add-on v3: moving charges + covariant propagator (KK tower)
========================================

--- (A) Brane Lorentz invariance: k^2 is invariant under boosts ---

k2 = -k0^2 + kx^2 + ky^2 + kz^2 -> -k0^2 + kx^2 + ky^2 + kz^2
k2' after boost -> -k0^2 + kx^2 + ky^2 + kz^2
Check k2' - k2 (should be 0): 0
OK: k^2 is invariant under the brane Lorentz boost.

--- (B) Gaussian-localized KK data: m_n^2 and c_n ---

Zint -> lambdaConf*Sqrt[Pi]
mu0eff -> mu0/(lambdaConf*Sqrt[Pi])
m_n^2 = 2 n / lambdaConf^2
c_n = f_n(0)^2 / norm_n, with norm_n = 2^n lambdaConf Sqrt[Pi] n!
Odd-n c_n should be 0. Spot-check (n,c_n): {{1, 0}, {3, 0}, {5, 0}, {7, 0}, {9, 0}}

--- (C) Brane-to-brane propagator D_eff(k^2) depends only on k^2 ---

Example truncations (nmax=0..6), showing dependence only on k2:
  nmax=0 : Deff(k2) = mu0/((I*epsReg + k2sym)*lambdaConf*Sqrt[Pi])
  nmax=1 : Deff(k2) = mu0/((I*epsReg + k2sym)*lambdaConf*Sqrt[Pi])
  nmax=2 : Deff(k2) = ((8 - (3*I)*epsReg*lambdaConf^2 - 3*k2sym*lambdaConf^2)*mu0)/(2*(4*(I*epsReg + k2sym)*lambdaConf + (epsReg - I*k2sym)^2*lambdaConf^3)*Sqrt[Pi])
  nmax=3 : Deff(k2) = ((8 - (3*I)*epsReg*lambdaConf^2 - 3*k2sym*lambdaConf^2)*mu0)/(2*(4*(I*epsReg + k2sym)*lambdaConf + (epsReg - I*k2sym)^2*lambdaConf^3)*Sqrt[Pi])
  nmax=4 : Deff(k2) = ((256*I + 5*(epsReg - I*k2sym)*lambdaConf^2*(28 - (3*I)*epsReg*lambdaConf^2 - 3*k2sym*lambdaConf^2))*mu0)/(8*(epsReg - I*k2sym)*lambdaConf*(4*I + (epsReg - I*k2sym)*lambdaConf^2)*(8*I + (epsReg - I*k2sym)*lambdaConf^2)*Sqrt[Pi])
  nmax=5 : Deff(k2) = ((256*I + 5*(epsReg - I*k2sym)*lambdaConf^2*(28 - (3*I)*epsReg*lambdaConf^2 - 3*k2sym*lambdaConf^2))*mu0)/(8*(epsReg - I*k2sym)*lambdaConf*(4*I + (epsReg - I*k2sym)*lambdaConf^2)*(8*I + (epsReg - I*k2sym)*lambdaConf^2)*Sqrt[Pi])
  nmax=6 : Deff(k2) = ((-6144 + 4032*(I*epsReg + k2sym)*lambdaConf^2 + 700*(epsReg - I*k2sym)^2*lambdaConf^4 - (35*I)*(epsReg - I*k2sym)^3*lambdaConf^6)*mu0)/(16*(epsReg - I*k2sym)*lambdaConf*(4*I + (epsReg - I*k2sym)*lambdaConf^2)*(8*I + (epsReg - I*k2sym)*lambdaConf^2)*(12*I + (epsReg - I*k2sym)*lambdaConf^2)*Sqrt[Pi])

Check Deff(k2') - Deff(k2) at nmax=6 (should be 0): 0
OK: D_eff is a function of the Lorentz scalar k^2 only (manifest brane Lorentz invariance).

--- (D) Static limit (k0=0): Coulomb + Yukawa tower, isotropic in r ---

A0_zero(r) = mu0eff qStar/(4 pi r) -> (mu0*qStar)/(4*lambdaConf*Pi^(3/2)*r)
First correction (keep n=0,2): A0(r) = ((2 + E^((-2*r)/lambdaConf))*mu0*qStar)/(8*lambdaConf*Pi^(3/2)*r)
Relative correction Delta(r) = (A0 - A0zero)/A0zero (nmax=2): 1/(2*E^((2*r)/lambdaConf))
Note: all static corrections depend only on r (no angular dependence), hence no brane anisotropy.

--- (E) Retarded solutions exist mode-by-mode (structure statement) ---

Each KK mode on the brane is Lorentz covariant:
  - n=0 behaves like standard Maxwell (massless) with retarded Green G0(x-x').
  - n>0 modes behave like massive (Proca-type) vectors with Klein-Gordon operator (Box + m_n^2) on each component.
  The 4D retarded Green function for (Box + m^2) depends only on the invariant interval (t^2 - r^2) and theta(t),
  so the sum over modes preserves exact brane Lorentz covariance.

Canonical (3+1) retarded Green functions (c=1) for reference (printed as strings):
  Massless wave:   G0_ret(t,r) = theta(t) * delta(t-r) / (4 pi r)
  Massive KG:      Gm_ret(t,r) = theta(t) * [ delta(t^2-r^2)/(2 pi) - (m/(4 pi)) theta(t^2-r^2) J1(m sqrt(t^2-r^2))/sqrt(t^2-r^2) ]

Moving-defect-branch potentials can be written covariantly as:
  A_mu(x) = ∫ d^4x' D_eff_ret(x-x') J_mu(x')
and in momentum space:
  A_mu(k) = D_eff(k^2) J_mu(k),  with D_eff depending only on k^2.
This is the cleanest way to see there is no preferred-frame artifact in the brane EM sector.

========== End add-on v3 ==========
"*)
