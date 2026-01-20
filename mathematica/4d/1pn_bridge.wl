(* ::Package:: *)

(*******************************************************
 Paper 7 (Hard-Mode 4D) -> Papers 1–6 bridge checks

 This Wolfram Language file is meant to be a small, self-contained
 “bridge appendix harness”: it reproduces the series’ key 1PN
 coefficients as explicit algebra/integrals, and can export clean
 TeX blocks.

 Scope:
   A) Added-mass coefficient for (i) 3D sphere (hyper-cylinder cross-section)
      and (ii) 4D ball (hypersphere/bubble) to show the topology lever.
   B) Optics: N(r) expansion from polytropic EOS, with n=5 giving coefficient 2.
   C) EIH wake-sector: solve {C_parallel, C_L} for {alpha^2, K} (minimal a_H=0).
   D) Relativistic kinetic-energy expansion: gamma series gives 3/8 v^4 term.

 v4 fixes (relative to v3)
 ------------------------
 FIX #3 (critical): Mathematica cannot Solve[] for an expression like alpha^2.
   - Solve variables must be symbols.
   - We introduce alpha2 as a stand-in for alpha^2, solve for {alpha2, K},
     then report alpha^2 = alpha2.

 This file is designed to be “referee proof”: every claimed number is
 either computed from an integral or from a symbolic solve.
*******************************************************)

(* ----------------------------- *)
(* 0) Pretty printing + asserts  *)
(* ----------------------------- *)

ClearAll[Section, Subsection, EqPrint, AssertTrue, TeXBlock];

Section[title_String] := (Print["\n\n========================================"]; Print[title]; Print["========================================\n"]);
Subsection[title_String] := (Print["\n--- ", title, " ---\n"]);

EqPrint[label_String, expr_] := (Print[label]; Print[TraditionalForm[expr]]; Print[""]);

AssertTrue[cond_, msg_String : "Assertion failed"] := If[TrueQ[cond], Null, Print[Style[msg, Red, Bold]]; Abort[]];

TeXBlock[tex_String] := (Print["\n--- TeX ---\n"]; Print[tex]; Print["\n-----------\n"]);

(*)
  Global assumptions used by Integrate/Simplify.

  IMPORTANT: keep these *minimal* and section-agnostic.
  In particular, do NOT include unrelated inequalities like v<c here,
  because Solve[] will sometimes return ConditionalExpression[..., v<c]
  even when v does not appear.
*)
$Assumptions = {a > 0, rho0 > 0, U > 0, c > 0, n > 1, K > 0, rho00 > 0};

(* ----------------------------- *)
(* 1) Added mass: 3D sphere      *)
(*    (hyper-cylinder cross-sec) *)
(* ----------------------------- *)

Section["1) Added mass from geometry"]; 

Subsection["1.1 3D sphere result (relevant to w-uniform hyper-cylinder)"];

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

EqPrint["Induced-flow kinetic energy T3 =", T3];
EqPrint["Added mass m_add(3D) = 2 T3/U^2 =", mAdd3];
EqPrint["Displaced mass m_disp(3D) = rho0 * Vol(B^3) =", mDisp3];
EqPrint["Ratio m_add/m_disp (3D) =", ratio3];

AssertTrue[ratio3 === 1/2, "Expected added-mass ratio 1/2 in 3D."];

Subsection["Interpretation"];
Print["For a w-uniform hyper-cylinder (S^2 x R) moving transverse to w, the external potential flow reduces to the 3D cross-section problem, giving κ_add = 1/2."];

(* ----------------------------- *)
(* 2) Added mass: 4D sphere      *)
(*    (bubble / hypersphere)     *)
(* ----------------------------- *)

Subsection["1.2 4D ball result (bubble/hypersphere counterfactual)"];

ClearAll[r, chi, th, ph];

(* 4D spherical coordinates:
   ds^2 = dr^2 + r^2 ( dχ^2 + sin^2χ ( dθ^2 + sin^2θ dφ^2 ) )
   Volume element: r^3 sin^2χ sinθ dr dχ dθ dφ.

   Induced dipole potential for uniform flow past a 4-ball:
     φ_induced = -(U a^4/(3 r^3)) cosχ.
*)

phi4Induced[r_, chi_] := -(U a^4/(3 r^3)) Cos[chi];

v4r[r_, chi_] := D[phi4Induced[r, chi], r];
v4chi[r_, chi_] := (1/r) D[phi4Induced[r, chi], chi];

v4sq[r_, chi_] := Simplify[v4r[r, chi]^2 + v4chi[r, chi]^2];

T4 = Simplify[
  (rho0/2) Integrate[
    v4sq[r, chi] * r^3 * Sin[chi]^2 * Sin[th],
    {r, a, Infinity}, {chi, 0, Pi}, {th, 0, Pi}, {ph, 0, 2 Pi},
    Assumptions -> $Assumptions
  ]
];

mAdd4 = Simplify[2 T4/U^2];
(* Volume of 4-ball: (π^2/2) a^4 *)
mDisp4 = Simplify[rho0 * (Pi^2/2) a^4];
ratio4 = Simplify[mAdd4/mDisp4];

EqPrint["Induced-flow kinetic energy T4 =", T4];
EqPrint["Added mass m_add(4D) = 2 T4/U^2 =", mAdd4];
EqPrint["Displaced mass m_disp(4D) = rho0 * Vol(B^4) =", mDisp4];
EqPrint["Ratio m_add/m_disp (4D) =", ratio4];

AssertTrue[ratio4 === 1/3, "Expected added-mass ratio 1/3 in 4D."];

Print["This is the clean topology discriminator: if the defect were a 4D bubble (B^4), κ_add=1/3 rather than 1/2."];

(* ----------------------------- *)
(* 3) Optics: stiffness n=5      *)
(* ----------------------------- *)

Section["2) Optics: polytropic index n fixes the 1/r coefficient"]; 

Subsection["2.1 EOS + sound speed expansion"]; 

ClearAll[rho, rho00, delta, n, K, c, s, cs, cs2, assumeOptics, csExact, NExact, csSeries1, NSeries1, csSeries2, NSeries2];

(* Local assumptions for this block only (keeps Solve[] from emitting ConditionalExpression noise). *)
assumeOptics = {rho00 > 0, n > 1, c > 0, K > 0, Element[delta, Reals], -1 < delta < 1};

P[rho_] := K rho^n;

(* FIX #1 (carried from v3):
   Do NOT define cs2[rho_] as D[P[rho], rho], because then calling cs2[rho00 (1+δ)]
   asks Mathematica to differentiate with respect to (rho00(1+δ)), which is not a valid
   differentiation variable.
   Instead, differentiate with respect to a dummy symbol and substitute.
*)
cs2[rho_] := Module[{s}, D[P[s], s] /. s -> rho];
cs[rho_] := Sqrt[cs2[rho]];

(* Fix K so that cs(rho00)=c (background sets the wave speed). *)
Ksol = First @ Assuming[assumeOptics, Solve[cs2[rho00] == c^2, K]] // FullSimplify;

(* FIX #2 (carried from v3):
   Compute N from the *exact* cs expression, then Series[] expand.
   (Do not expand cs first and then take 1/cs, because truncation can corrupt coefficients.)
*)
csExact = Assuming[assumeOptics, FullSimplify[cs[rho00 (1 + delta)] /. Ksol]];
NExact  = Assuming[assumeOptics, FullSimplify[c/csExact]];

csSeries1 = Assuming[assumeOptics, Normal @ Series[csExact, {delta, 0, 1}] // FullSimplify];
NSeries1  = Assuming[assumeOptics, Normal @ Series[NExact,  {delta, 0, 1}] // FullSimplify];

(* Upgrade #1 (carried from v3): also print O(δ^2) so the referee can see the next correction. *)
csSeries2 = Assuming[assumeOptics, Normal @ Series[csExact, {delta, 0, 2}] // FullSimplify];
NSeries2  = Assuming[assumeOptics, Normal @ Series[NExact,  {delta, 0, 2}] // FullSimplify];

EqPrint["K chosen so that cs(rho00)=c:", Ksol];
EqPrint["cs_exact(rho00(1+δ)):", csExact];
EqPrint["N_exact(δ)=c/cs_exact:", NExact];

EqPrint["cs(rho00(1+δ)) to O(δ):", csSeries1];
EqPrint["N=c/cs to O(δ):", NSeries1];

EqPrint["cs(rho00(1+δ)) to O(δ^2):", csSeries2];
EqPrint["N=c/cs to O(δ^2):", NSeries2];

(* Upgrade #2 (carried from v3): sanity checks that pin the normalization unambiguously. *)
AssertTrue[Assuming[assumeOptics, FullSimplify[csExact /. delta -> 0]] === c, "Expected cs(δ=0)=c"]; 
AssertTrue[Assuming[assumeOptics, FullSimplify[NExact  /. delta -> 0]] === 1, "Expected N(δ=0)=1"]; 

(* N = 1 - (n-1)/2 δ + ... *)
AssertTrue[
  Assuming[assumeOptics, FullSimplify[NSeries1 - (1 - (n - 1) delta/2)]] === 0,
  "Expected N = 1 - (n-1)δ/2 + ..."
];

EqPrint["Linear coefficient in N(δ) (coeff of δ):", Assuming[assumeOptics, SeriesCoefficient[NExact, {delta, 0, 1}] // FullSimplify]];

Subsection["2.1b Enthalpy link (tightens the story)"];
(* Enthalpy per unit mass for a barotrope: h(ρ)=∫ (dP/ρ). For P=K ρ^n (n>1):
     h(ρ) = (n K/(n-1)) ρ^(n-1)    (constant absorbed).
   This gives the clean identity:  c_s^2 = dP/dρ = (n-1) h(ρ).
*)
ClearAll[h];
h[rho_] := (n K/(n - 1)) rho^(n - 1);
EqPrint["h(ρ) = (nK/(n-1)) ρ^(n-1):", h[rho]];
EqPrint["Check: cs^2(ρ) - (n-1) h(ρ) =", Assuming[assumeOptics, FullSimplify[cs2[rho] - (n - 1) h[rho]]]];
EqPrint["Background enthalpy h0 = h(ρ00) with cs(ρ00)=c:", Assuming[assumeOptics, FullSimplify[h[rho00] /. Ksol]]];
AssertTrue[
  Assuming[assumeOptics, FullSimplify[(h[rho00] /. Ksol) - c^2/(n - 1)]] === 0,
  "Expected h0 = c^2/(n-1) under cs(ρ00)=c"
];

Subsection["2.2 Match to GR coefficient"];
Print["Using the series-wide density-starvation relation δ ≃ Φ/c^2 (with Φ≈-GM/r), we get:"];
TeXBlock[
  "\\begin{align}\n" <>
  "N(r) \\simeq 1 - \\frac{n-1}{2}\\,\\frac{\\Phi(r)}{c^2}\\, , \\qquad \\Phi(r)\\simeq -\\frac{GM}{r}.\n" <>
  "\\end{align}"
];
Print["Therefore N(r) ≃ 1 + (n-1)/2 * GM/(rc^2). Setting this coefficient to 2 forces n=5."];

(* Quick solve: (n-1)/2 == 2 *)
ns = Solve[(n - 1)/2 == 2, n] // Simplify;
EqPrint["Solve (n-1)/2 = 2:", ns];
AssertTrue[ns === {{n -> 5}}, "Expected n=5"]; 

(* ----------------------------- *)
(* 4) EIH wake-sector solve      *)
(* ----------------------------- *)

Section["3) Vector/EIH sector: solve for α^2 and K"]; 

Subsection["3.1 Coefficients from the wake overlap reduction"];

ClearAll[alpha, alpha2, aH];

(* These closed-form coefficients are the ones carried forward in the series.
   Minimal EIH match uses a_H=0. *)
Cpar[alpha_, aH_] := K Pi^2 (-1 + aH^2 - alpha^2);
CL  [alpha_, aH_] := K Pi^2 (-1 + aH^2 + alpha^2);

EqPrint["C_parallel(alpha,a_H) =", Cpar[alpha, aH]];
EqPrint["C_L(alpha,a_H)       =", CL[alpha, aH]];

Subsection["3.2 Minimal match (a_H=0) to EIH targets"];

(* FIX #3:
   Solve variables must be symbols. Mathematica rejects Solve[..., {alpha^2, K}]
   because alpha^2 is an expression, not a symbol.

   Workaround: introduce alpha2 ≡ alpha^2 as an independent symbol, solve for {alpha2, K},
   and then *report* alpha^2 = alpha2.
*)

eqsAlpha2 = {
    Cpar[alpha, 0] == -7/2,
    CL[alpha, 0]   == -1/2
  } /. alpha^2 -> alpha2;

solAlpha2K = Solve[eqsAlpha2, {alpha2, K}] // FullSimplify;

EqPrint["Solve {C_parallel=-7/2, C_L=-1/2} for {alpha^2, K}:", solAlpha2K];

(* Extract principal solution *)
alpha2Val = (alpha2 /. First[solAlpha2K]) // Simplify;
KVal      = (K      /. First[solAlpha2K]) // Simplify;

EqPrint["alpha^2 =", alpha2Val];
EqPrint["K =", KVal];

AssertTrue[alpha2Val === 3/4, "Expected alpha^2=3/4"]; 
AssertTrue[KVal === 2/Pi^2, "Expected K=2/pi^2"]; 

(* ----------------------------- *)
(* 5) Relativistic v^4 term      *)
(* ----------------------------- *)

Section["4) Special relativity: v^4 coefficient (3/8)"]; 

Subsection["4.1 Gamma expansion"];

ClearAll[v, c];

gamma = 1/Sqrt[1 - v^2/c^2];
seriesGamma = Normal@Series[gamma, {v, 0, 4}] // Simplify;
EqPrint["gamma(v) series to O(v^4):", seriesGamma];

(* Energy of trapped massless mode: E(v)=gamma E0. *)
ClearAll[E0];
Eseries = Expand[seriesGamma * E0];
EqPrint["E(v)=gamma E0 series:", Eseries];

(* Coefficient check *)
coeffV4 = Coefficient[Eseries, v, 4] // Simplify;
EqPrint["Coefficient of v^4 term:", coeffV4];
AssertTrue[coeffV4 === (3 E0)/(8 c^4), "Expected (3/8) E0 v^4/c^4 term."];

Print["This block is universal: any trapped massless mode with rest energy E0 gives the same 3/8 coefficient."];

(* ----------------------------- *)
(* End *)
(* ----------------------------- *)

Section["Done"]; 
Print["All bridge checks passed."];

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


--- Interpretation ---

For a w-uniform hyper-cylinder (S^2 x R) moving transverse to w, the external potential flow reduces to the 3D cross-section problem, giving κ_add = 1/2.

--- 1.2 4D ball result (bubble/hypersphere counterfactual) ---

Induced-flow kinetic energy T4 =
TraditionalForm[(a^4*Pi^2*rho0*U^2)/12]

Added mass m_add(4D) = 2 T4/U^2 =
TraditionalForm[(a^4*Pi^2*rho0)/6]

Displaced mass m_disp(4D) = rho0 * Vol(B^4) =
TraditionalForm[(a^4*Pi^2*rho0)/2]

Ratio m_add/m_disp (4D) =
TraditionalForm[1/3]

This is the clean topology discriminator: if the defect were a 4D bubble (B^4), κ_add=1/3 rather than 1/2.


========================================
2) Optics: polytropic index n fixes the 1/r coefficient
========================================


--- 2.1 EOS + sound speed expansion ---

K chosen so that cs(rho00)=c:
TraditionalForm[{K -> (c^2*rho00^(1 - n))/n}]

cs_exact(rho00(1+δ)):
TraditionalForm[c*(1 + delta)^((-1 + n)/2)]

N_exact(δ)=c/cs_exact:
TraditionalForm[(1 + delta)^(1/2 - n/2)]

cs(rho00(1+δ)) to O(δ):
TraditionalForm[c + (c*delta*(-1 + n))/2]

N=c/cs to O(δ):
TraditionalForm[(2 + delta - delta*n)/2]

cs(rho00(1+δ)) to O(δ^2):
TraditionalForm[c + (c*delta*(4 + delta*(-3 + n))*(-1 + n))/8]

N=c/cs to O(δ^2):
TraditionalForm[1 + (delta*(-1 + n)*(-4 + delta + delta*n))/8]

Linear coefficient in N(δ) (coeff of δ):
TraditionalForm[(1 - n)/2]


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
3) Vector/EIH sector: solve for α^2 and K
========================================


--- 3.1 Coefficients from the wake overlap reduction ---

C_parallel(alpha,a_H) =
TraditionalForm[(-1 + aH^2 - alpha^2)*K*Pi^2]

C_L(alpha,a_H)       =
TraditionalForm[(-1 + aH^2 + alpha^2)*K*Pi^2]


--- 3.2 Minimal match (a_H=0) to EIH targets ---

Solve {C_parallel=-7/2, C_L=-1/2} for {alpha^2, K}:
TraditionalForm[{{alpha2 -> 3/4, K -> 2/Pi^2}}]

alpha^2 =
TraditionalForm[3/4]

K =
TraditionalForm[2/Pi^2]



========================================
4) Special relativity: v^4 coefficient (3/8)
========================================


--- 4.1 Gamma expansion ---

gamma(v) series to O(v^4):
TraditionalForm[1 + v^2/(2*c^2) + (3*v^4)/(8*c^4)]

E(v)=gamma E0 series:
TraditionalForm[E0 + (E0*v^2)/(2*c^2) + (3*E0*v^4)/(8*c^4)]

Coefficient of v^4 term:
TraditionalForm[(3*E0)/(8*c^4)]

This block is universal: any trapped massless mode with rest energy E0 gives the same 3/8 coefficient.


========================================
Done
========================================

All bridge checks passed.
"*)
