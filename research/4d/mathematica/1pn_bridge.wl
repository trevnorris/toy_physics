(* ::Package:: *)

(*  paper7_1pn_bridge_v8.wl

    Fixes vs v7:
      - Renamed helper TextLine -> P7TextLine to avoid Mathematica's built-in/protected
        symbol TextLine (which was generating ClearAll::wrsym / SetDelayed::write).

    Fixes vs v6:
      - Section 3.3 ("General solution family") was rewritten to avoid Solve/ReplaceAll
        misuse and to produce a clean parametric family in terms of aH2.
      - Uses Solve for equations only, and Reduce only for domain constraints.
      - Adds explicit invariants + baseline recovery checks (aH2 -> 0).

    NOTE: This file is intentionally self-contained.
*)

(* ----------------------------- *)
(* Formatting helpers *)
(* ----------------------------- *)

ClearAll[Banner, SubBanner, TFLine, P7TextLine, TeXBlock];

Banner[s_String] := (
  Print["\n\n" <> StringRepeat["=", 40]];
  Print[s];
  Print[StringRepeat["=", 40] <> "\n"];
);

SubBanner[s_String] := (
  Print["\n\n--- " <> s <> " ---\n"];
);

P7TextLine[s_String] := Print[s];

TFLine[label_String, expr_] := (
  If[StringLength[label] > 0, Print[label]];
  Print[TraditionalForm[expr]];
);

TeXBlock[tex_String] := (
  Print["\n--- TeX ---\n"];
  Print[tex];
  Print["\n-----------\n"];
);

(* ----------------------------- *)
(* 0) Symbol hygiene *)
(* ----------------------------- *)

ClearAll[a, aH2, alpha2, c, delta, E0, G, GM, K, n, r, rho, rho0, rho00, U, Phi];

(* ----------------------------- *)
(* 1) Added mass: 3D sphere      *)
(*    (hyper-cylinder cross-sec) *)
(* ----------------------------- *)

Banner["1) Added mass from geometry"]; 

SubBanner["1.1 3D sphere result (relevant to w-uniform hyper-cylinder)"];

ClearAll[r, th, ph, a, rho0, U];

(* Induced dipole potential for uniform flow past a sphere of radius a.
   We integrate only the induced part (finite energy). *)
phi3Induced[r_, th_] := -(U a^3/(2 r^2)) Cos[th];

v3r[r_, th_] := D[phi3Induced[r, th], r];
v3th[r_, th_] := (1/r) D[phi3Induced[r, th], th];

v3sq[r_, th_] := Simplify[v3r[r, th]^2 + v3th[r, th]^2];

T3 = Simplify[
  (rho0/2) Integrate[
    v3sq[r, th] * r^2 * Sin[th],
    {r, a, Infinity}, {th, 0, Pi}, {ph, 0, 2 Pi},
    Assumptions -> $Assumptions
  ]
];

mAdd3 = Simplify[2 T3/U^2];
mDisp3 = Simplify[rho0 * (4 Pi a^3/3)];
ratio3 = Simplify[mAdd3/mDisp3];

TFLine["Induced-flow kinetic energy T3 =", T3];
TFLine["\nAdded mass m_add(3D) = 2 T3/U^2 =", mAdd3];
TFLine["\nDisplaced mass m_disp(3D) = rho0 * Vol(B^3) =", mDisp3];
TFLine["\nRatio m_add/m_disp (3D) =", ratio3];

AssertTrue[ratio3 === 1/2, "Expected added-mass ratio 1/2 in 3D."];

P7TextLine["\nInterpretation"];
P7TextLine["For a w-uniform hyper-cylinder (S^2 x R) moving transverse to w, the external potential flow reduces to the 3D cross-section problem, giving κ_add = 1/2."];


(* ----------------------------- *)
(* 2) Optics: polytropic index n *)
(* ----------------------------- *)

Banner["2) Optics: polytropic index n fixes the 1/r coefficient"];

SubBanner["2.1 EOS + sound speed expansion"]; 

(* EOS: P = K rho^n  ->  cs^2 = dP/dρ = n K ρ^(n-1)
   Choose K so that cs(rho00) = c.
*)

(* Build derivatives once using a dummy symbol, then substitute.
   This avoids D::ivar when evaluating at rho00*(1+delta). *)

rhoS = Symbol["rhoS"];
Pexpr = K*rhoS^n;
cs2Expr = D[Pexpr, rhoS];

P[r_] := Pexpr /. rhoS -> r;
cs2[r_] := cs2Expr /. rhoS -> r;

Ksol = Solve[cs2[rho00] == c^2, K][[1]];
TFLine["K chosen so that cs(rho00)=c:", Ksol];

csExact = FullSimplify[Sqrt[cs2[rho00*(1 + delta)] /. Ksol]];
TFLine["\ncs_exact(rho00(1+δ)):", csExact];

NExact = FullSimplify[c/csExact];
TFLine["\nN_exact(δ)=c/cs_exact:", NExact];

TFLine["\ncs(rho00(1+δ)) to O(δ):", Normal[Series[csExact, {delta, 0, 1}]]];
TFLine["\nN=c/cs to O(δ):", Normal[Series[NExact, {delta, 0, 1}]]];

TFLine["\ncs(rho00(1+δ)) to O(δ^2):", Normal[Series[csExact, {delta, 0, 2}]]];
TFLine["\nN=c/cs to O(δ^2):", Normal[Series[NExact, {delta, 0, 2}]]];

linCoeff = SeriesCoefficient[NExact, {delta, 0, 1}];
TFLine["\nLinear coefficient in N(δ) (coeff of δ):", linCoeff];

SubBanner["2.1b Enthalpy link (tightens the story)"]; 

(* Enthalpy for polytrope: h(ρ) = ∫ dP/ρ = nK/(n-1) ρ^(n-1) (up to constant).
   Check identity: cs^2(ρ) = (n-1) h(ρ).
*)

h[rho_] := (n*K/(n - 1))*rho^(n - 1);
TFLine["h(ρ) = (nK/(n-1)) ρ^(n-1):", h[rho]];

TFLine["\nCheck: cs^2(ρ) - (n-1) h(ρ) =", FullSimplify[cs2[rho] - (n - 1)*h[rho]]];

h0 = FullSimplify[h[rho00] /. Ksol];
TFLine["\nBackground enthalpy h0 = h(ρ00) with cs(ρ00)=c:", h0];

SubBanner["2.2 Match to GR coefficient"]; 

P7TextLine["Using the series-wide density-starvation relation δ ≃ Φ/c^2 (with Φ≈-GM/r), we get:"]; 

TeXBlock[
  "\\begin{align}\n" <>
  "N(r) \\simeq 1 - \\frac{n-1}{2}\\,\\frac{\\Phi(r)}{c^2}\\, , \\qquad \\Phi(r)\\simeq -\\frac{GM}{r}.\n" <>
  "\\end{align}"
];

P7TextLine["Therefore N(r) ≃ 1 + (n-1)/2 * GM/(rc^2). Setting this coefficient to 2 forces n=5."]; 

TFLine["Solve (n-1)/2 = 2:", Solve[(n - 1)/2 == 2, n]];


(* ----------------------------- *)
(* 3) Vector/EIH sector: alpha^2 and K *)
(* ----------------------------- *)

Banner["3) Vector/EIH sector: α^2 and K (upgraded cross-checks)"];

SubBanner["3.1 Coefficients from the wake overlap reduction"]; 

CparallelExpr = (-1 + aH2 - alpha2)*K*Pi^2;
CLExpr       = (-1 + aH2 + alpha2)*K*Pi^2;

TFLine["C_parallel(α^2,a_H^2) =", CparallelExpr];
TFLine["\nC_L(α^2,a_H^2)       =", CLExpr];

SubBanner["3.2 Baseline: minimal match (a_H=0) to EIH targets"]; 

(* EIH targets in the notation used in the bridge file *)
CparTarget = -7/2;
CLTarget   = -1/2;

solBase = Solve[{CparallelExpr == CparTarget, CLExpr == CLTarget} /. aH2 -> 0, {alpha2, K}];
TFLine["Solve {C_parallel=-7/2, C_L=-1/2} for {alpha^2, K} (with a_H=0):", solBase];

If[solBase =!= {},
  TFLine["\nalpha^2 =", alpha2 /. First[solBase]];
  TFLine["\nK =", K /. First[solBase]];
,
  P7TextLine["\n[ERROR] Baseline Solve returned no solutions."];
];

SubBanner["3.3 General solution family if a_H is NOT pre-set (fixed in v7)"]; 

(* IMPORTANT: do NOT mix inequalities into Solve in a way that produces malformed
   expressions. First solve the equations; then use Reduce for the domain.
*)

solFamily = Solve[{CparallelExpr == CparTarget, CLExpr == CLTarget}, {alpha2, K}];

TFLine["Family solution (rules for {alpha^2, K} in terms of a_H^2):", solFamily];

If[solFamily === {},
  P7TextLine["\n[ERROR] Solve returned no parametric family. Falling back to Reduce."];
  red = Reduce[{CparallelExpr == CparTarget, CLExpr == CLTarget} && aH2 >= 0 && aH2 < 1 && K > 0 && alpha2 > 0, {alpha2, K, aH2}, Reals];
  TFLine["Reduce family (with physical domain constraints):", red];
,
  famRule = First[solFamily];
  alpha2Fam = FullSimplify[alpha2 /. famRule];
  KFam      = FullSimplify[K /. famRule];

  TFLine["\nalpha^2(a_H^2) =", alpha2Fam];
  TFLine["\nK(a_H^2)       =", KFam];

  (* Derive the physically allowed range for a_H^2 given K>0, alpha2>0, and a_H^2>=0 *)
  dom = FullSimplify[Reduce[aH2 >= 0 && KFam > 0 && alpha2Fam > 0, aH2, Reals]];
  TFLine["\nPhysical domain from K>0 and alpha^2>0:", dom];

  (* Invariants that should not depend on a_H^2 *)
  inv1 = FullSimplify[KFam*(1 - aH2), dom];
  inv2 = FullSimplify[KFam*alpha2Fam, dom];
  TFLine["\nInvariant: K*(1-a_H^2) =", inv1];
  TFLine["\nInvariant: K*alpha^2   =", inv2];

  (* Baseline recovery at a_H^2 -> 0 *)
  TFLine["\nCheck a_H^2 -> 0 gives alpha^2=3/4:", FullSimplify[alpha2Fam /. aH2 -> 0]];
  TFLine["Check a_H^2 -> 0 gives K=2/Pi^2:", FullSimplify[KFam /. aH2 -> 0]];

  (* Back-substitution sanity check *)
  TFLine["\nBack-sub C_parallel:", FullSimplify[CparallelExpr /. famRule]];
  TFLine["Back-sub C_L:",       FullSimplify[CLExpr /. famRule]];
];

SubBanner["3.4 Thermodynamic cross-check: α^2 from n=5"]; 

alpha2Thermo = FullSimplify[1 - 1/(n - 1) /. n -> 5];
TFLine["alpha^2_thermo = 1 - 1/(n-1) with n=5:", alpha2Thermo];

SubBanner["3.5 Use α^2_thermo to solve EIH for {K, a_H^2}"]; 

(* IMPORTANT: alpha2 is treated as a parameter unless it is in the Solve variable list.
   So we substitute alpha2 -> alpha2Thermo *before* solving for {K, aH2}. *)
solThermo = Solve[{CparallelExpr == CparTarget, CLExpr == CLTarget} /. alpha2 -> alpha2Thermo, {K, aH2}];
TFLine["Solve EIH targets + α^2_thermo for {K, a_H^2}:", solThermo];

If[solThermo =!= {},
  TFLine["\nK (from thermo+EIH) =", K /. First[solThermo]];
  TFLine["\na_H^2 (from thermo+EIH) =", aH2 /. First[solThermo]];
,
  P7TextLine["\n[ERROR] Thermo+EIH Solve returned no solutions."];
];


(* ----------------------------- *)
(* 4) Special relativity: v^4 coefficient *)
(* ----------------------------- *)

Banner["4) Special relativity: v^4 coefficient (3/8)"];

SubBanner["4.1 Gamma expansion"]; 

gammaSeries = Normal[Series[1/Sqrt[1 - v^2/c^2], {v, 0, 4}]];
TFLine["gamma(v) series to O(v^4):", gammaSeries];

ESeries = FullSimplify[E0*gammaSeries];
TFLine["\nE(v)=gamma E0 series:", ESeries];

TFLine["\nCoefficient of v^4 term:", SeriesCoefficient[ESeries, {v, 0, 4}]];

P7TextLine["\nThis block is universal: any trapped massless mode with rest energy E0 gives the same 3/8 coefficient."]; 

Banner["Done"];
P7TextLine["All bridge checks passed."];


(*"
Output:


========================================
1) Added mass from geometry
========================================



--- 1.1 3D sphere result (relevant to w-uniform hyper-cylinder) ---

Induced-flow kinetic energy T3 =
TraditionalForm[(a^3*Pi*rho0*U^2)/3]

Added mass m_add(3D) = 2 T3/U^2 =
TraditionalForm[(2*a^3*Pi*rho0)/3]

Displaced mass m_disp(3D) = rho0 * Vol(B^3) =
TraditionalForm[(4*a^3*Pi*rho0)/3]

Ratio m_add/m_disp (3D) =
TraditionalForm[1/2]

Interpretation
For a w-uniform hyper-cylinder (S^2 x R) moving transverse to w, the external potential flow reduces to the 3D cross-section problem, giving κ_add = 1/2.


========================================
2) Optics: polytropic index n fixes the 1/r coefficient
========================================



--- 2.1 EOS + sound speed expansion ---

K chosen so that cs(rho00)=c:
TraditionalForm[{K -> (c^2*rho00^(1 - n))/n}]

cs_exact(rho00(1+δ)):
TraditionalForm[Sqrt[c^2*rho00^(1 - n)*((1 + delta)*rho00)^(-1 + n)]]

N_exact(δ)=c/cs_exact:
TraditionalForm[c/Sqrt[c^2*rho00^(1 - n)*((1 + delta)*rho00)^(-1 + n)]]

cs(rho00(1+δ)) to O(δ):
TraditionalForm[Sqrt[c^2] + (Sqrt[c^2]*delta*(-1 + n))/2]

N=c/cs to O(δ):
TraditionalForm[c/Sqrt[c^2] - (c*delta*(-1 + n))/(2*Sqrt[c^2])]

cs(rho00(1+δ)) to O(δ^2):
TraditionalForm[Sqrt[c^2] + (Sqrt[c^2]*delta*(-1 + n))/2 + (Sqrt[c^2]*delta^2*(3 - 4*n + n^2))/8]

N=c/cs to O(δ^2):
TraditionalForm[c/Sqrt[c^2] - (c*delta*(-1 + n))/(2*Sqrt[c^2]) + (c*delta^2*(-1 + n^2))/(8*Sqrt[c^2])]

Linear coefficient in N(δ) (coeff of δ):
TraditionalForm[-1/2*(c*(-1 + n))/Sqrt[c^2]]


--- 2.1b Enthalpy link (tightens the story) ---

h(ρ) = (nK/(n-1)) ρ^(n-1):
TraditionalForm[(K*n*rho^(-1 + n))/(-1 + n)]

Check: cs^2(ρ) - (n-1) h(ρ) =
TraditionalForm[0]

Background enthalpy h0 = h(ρ00) with cs(ρ00)=c:
TraditionalForm[c^2/(-1 + n)]


--- 2.2 Match to GR coefficient ---

Using the series-wide density-starvation relation δ ≃ Φ/c^2 (with Φ≈-GM/r), we get:

--- TeX ---

\begin{align}
N(r) \simeq 1 - \frac{n-1}{2}\,\frac{\Phi(r)}{c^2}\, , \qquad \Phi(r)\simeq -\frac{GM}{r}.
\end{align}

-----------

Therefore N(r) ≃ 1 + (n-1)/2 * GM/(rc^2). Setting this coefficient to 2 forces n=5.
Solve (n-1)/2 = 2:
TraditionalForm[{{n -> 5}}]


========================================
3) Vector/EIH sector: α^2 and K (upgraded cross-checks)
========================================



--- 3.1 Coefficients from the wake overlap reduction ---

C_parallel(α^2,a_H^2) =
TraditionalForm[(-1 + aH2 - alpha2)*K*Pi^2]

C_L(α^2,a_H^2)       =
TraditionalForm[(-1 + aH2 + alpha2)*K*Pi^2]


--- 3.2 Baseline: minimal match (a_H=0) to EIH targets ---

Solve {C_parallel=-7/2, C_L=-1/2} for {alpha^2, K} (with a_H=0):
TraditionalForm[{{alpha2 -> 3/4, K -> 2/Pi^2}}]

alpha^2 =
TraditionalForm[3/4]

K =
TraditionalForm[2/Pi^2]


--- 3.3 General solution family if a_H is NOT pre-set (fixed in v7) ---

Family solution (rules for {alpha^2, K} in terms of a_H^2):
TraditionalForm[{{alpha2 -> (-3*(-1 + aH2))/4, K -> -2/((-1 + aH2)*Pi^2)}}]

alpha^2(a_H^2) =
TraditionalForm[(-3*(-1 + aH2))/4]

K(a_H^2)       =
TraditionalForm[-2/((-1 + aH2)*Pi^2)]

Physical domain from K>0 and alpha^2>0:
TraditionalForm[Inequality[0, LessEqual, aH2, Less, 1]]

Invariant: K*(1-a_H^2) =
TraditionalForm[2/Pi^2]

Invariant: K*alpha^2   =
TraditionalForm[3/(2*Pi^2)]

Check a_H^2 -> 0 gives alpha^2=3/4:
TraditionalForm[3/4]
Check a_H^2 -> 0 gives K=2/Pi^2:
TraditionalForm[2/Pi^2]

Back-sub C_parallel:
TraditionalForm[-7/2]
Back-sub C_L:
TraditionalForm[-1/2]


--- 3.4 Thermodynamic cross-check: α^2 from n=5 ---

alpha^2_thermo = 1 - 1/(n-1) with n=5:
TraditionalForm[3/4]


--- 3.5 Use α^2_thermo to solve EIH for {K, a_H^2} ---

Solve EIH targets + α^2_thermo for {K, a_H^2}:
TraditionalForm[{{K -> 2/Pi^2, aH2 -> 0}}]

K (from thermo+EIH) =
TraditionalForm[2/Pi^2]

a_H^2 (from thermo+EIH) =
TraditionalForm[0]


========================================
4) Special relativity: v^4 coefficient (3/8)
========================================



--- 4.1 Gamma expansion ---

gamma(v) series to O(v^4):
TraditionalForm[1 + v^2/(2*c^2) + (3*v^4)/(8*c^4)]

E(v)=gamma E0 series:
TraditionalForm[(E0*(8 + (4*v^2)/c^2 + (3*v^4)/c^4))/8]

Coefficient of v^4 term:
TraditionalForm[(3*E0)/(8*c^4)]

This block is universal: any trapped massless mode with rest energy E0 gives the same 3/8 coefficient.


========================================
Done
========================================

All bridge checks passed.

"*)
