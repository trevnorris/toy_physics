(*
========================================================
Paper VIII (referee add-on v2): symmetry checks
  (A) Gauge invariance (5D, Z(w)-localized Maxwell sector)
  (B) 3+1 Lorentz invariance of the reduced Maxwell sector

Changes vs v1:
  - Make gammaB an explicit function of betaB, rather than a separate symbol
    constrained by an equality in $Assumptions. This avoids "not fully simplified"
    residuals of the form (1 + (-1 + beta^2) gamma^2).
  - Add an explicit check that the covector transform used is the inverse-transpose
    of the vector transform.
========================================================
*)

Print["\n========================================"];
Print["Paper VIII add-on v2: symmetry checks (gauge + brane Lorentz)"];
Print["========================================\n"];

ClearAll["Global`*"];

(* ---------- Assumptions ---------- *)
$Assumptions = {
  Element[{mu0, xi, lambdaConf, qStar, hbar}, Reals],
  mu0 > 0, xi > 0, lambdaConf > 0, hbar > 0
};

(* Robust "is zero" helpers (avoid brittle === 0 tests). *)
zeroQ[expr_] := TrueQ[FullSimplify[expr == 0, Assumptions -> $Assumptions]];
zeroMatQ[m_] := And @@ Flatten[Map[zeroQ, m, {2}]];

(* ---------- Coordinates and metrics ---------- *)
coords5 = {t, x, y, z, w};
eta5 = DiagonalMatrix[{-1, 1, 1, 1, 1}];

coords4 = {t, x, y, z};
eta4 = DiagonalMatrix[{-1, 1, 1, 1}];

Zprof[ww_] := Exp[-ww^2/lambdaConf^2];

Print["--- Objects ---"];
Print[HoldForm[coords5], " -> ", coords5];
Print[HoldForm[eta5], " -> ", eta5];
Print[HoldForm[Zprof[w]], " -> ", Zprof[w]];

(* ---------- Helper: build F_{MN} from a covector A_M ---------- *)
MakeFLower[A_List, coords_List] := Module[{n = Length[coords]},
  Table[
    D[A[[j]], coords[[i]]] - D[A[[i]], coords[[j]]]
    , {i, 1, n}, {j, 1, n}
  ]
];

RaiseBothIndices[Flow_, eta_] := eta . Flow . eta;

FFScalar[Flow_, eta_] := Module[{Fup},
  Fup = RaiseBothIndices[Flow, eta];
  (* The 1/2 avoids double-counting antisymmetric pairs. *)
  FullSimplify[
    Sum[Flow[[i, j]]*Fup[[i, j]], {i, 1, Length[eta]}, {j, 1, Length[eta]}]/2,
    Assumptions -> $Assumptions
  ]
];

(* ===================================================================== *)
(* (A) Gauge invariance in 5D                                             *)
(* ===================================================================== *)

Print["\n--- (A) Gauge invariance check (5D) ---\n"];

(* Fields and source components (treated as independent functions) *)
A5 = {
  A0[t, x, y, z, w],
  Ax[t, x, y, z, w],
  Ay[t, x, y, z, w],
  Az[t, x, y, z, w],
  Aw[t, x, y, z, w]
};

J5 = {
  J0[t, x, y, z, w],
  Jx[t, x, y, z, w],
  Jy[t, x, y, z, w],
  Jz[t, x, y, z, w],
  Jw[t, x, y, z, w]
};

chi5 = chi[t, x, y, z, w];

(* Gauge transform: A_M -> A_M + d_M chi *)
A5g = A5 + {
  D[chi5, t],
  D[chi5, x],
  D[chi5, y],
  D[chi5, z],
  D[chi5, w]
};

(* Gauge-invariant Lagrangian density: kinetic + source, no gauge-fixing *)
Flow5 = MakeFLower[A5, coords5];
Flow5g = MakeFLower[A5g, coords5];

FF5 = FFScalar[Flow5, eta5];
FF5g = FFScalar[Flow5g, eta5];

Lint5 = -(Zprof[w]/(4*mu0))*FF5 - (A5 . J5);
Lint5g = -(Zprof[w]/(4*mu0))*FF5g - (A5g . J5);

deltaLint5 = FullSimplify[Expand[Lint5g - Lint5], Assumptions -> $Assumptions];

Print["Delta L (gauge-invariant sector) should be a total derivative plus chi*(d·J):"];
Print["deltaLint5 = ", deltaLint5];

(* Target identity: deltaL = -d_M(chi J^M) + chi d_M J^M *)
DivJ5 = D[J5[[1]], t] + D[J5[[2]], x] + D[J5[[3]], y] + D[J5[[4]], z] + D[J5[[5]], w];
DivChiJ5 = D[chi5*J5[[1]], t] + D[chi5*J5[[2]], x] + D[chi5*J5[[3]], y] + D[chi5*J5[[4]], z] + D[chi5*J5[[5]], w];

targetDelta = FullSimplify[-DivChiJ5 + chi5*DivJ5, Assumptions -> $Assumptions];

resGauge = FullSimplify[deltaLint5 - targetDelta, Assumptions -> $Assumptions];

Print["Residual (deltaLint5 - (-d(chi J) + chi dJ)) should be 0:"];
Print[resGauge];

If[zeroQ[resGauge],
  Print["OK: Gauge invariance holds up to a boundary term when d_M J^M = 0."],
  Print["WARNING: Gauge invariance residual is nonzero (unexpected)."]
];

(* Show how DivA changes under gauge transform (relevant for gauge-fixing discussion). *)
DivA5[A_List] := Module[{A0c, Axc, Ayc, Azc, Awc},
  {A0c, Axc, Ayc, Azc, Awc} = A;
  (* DivA := d_M A^M, with A^0 = -A0 and A^i = Ai *)
  FullSimplify[
    D[Axc, x] + D[Ayc, y] + D[Azc, z] + D[Awc, w] - D[A0c, t],
    Assumptions -> $Assumptions
  ]
];

Box5[phi_] := FullSimplify[
  -D[phi, {t, 2}] + D[phi, {x, 2}] + D[phi, {y, 2}] + D[phi, {z, 2}] + D[phi, {w, 2}],
  Assumptions -> $Assumptions
];

divAorig = DivA5[A5];
divAg = DivA5[A5g];
boxChi = Box5[chi5];

resDiv = FullSimplify[divAg - divAorig - boxChi, Assumptions -> $Assumptions];

Print["\nCheck: DivA(A + d chi) - DivA(A) - Box(chi) should be 0:"];
Print[resDiv];

If[zeroQ[resDiv],
  Print["OK: DivA shifts by Box(chi), as expected."],
  Print["WARNING: DivA shift check failed (unexpected)."]
];

Print["\nNote: The gauge-fixing term is not gauge invariant; this is expected and harmless."];

(* ===================================================================== *)
(* (A2) Paper VII brane matter gauge convention                            *)
(* ===================================================================== *)

Print["\n--- (A2) Paper VII brane matter gauge convention ---\n"];

psiB = psiBf[t, x, y, z];
A0m = A0mF[t, x, y, z];
Axm = AxmF[t, x, y, z];
Aym = AymF[t, x, y, z];
Azm = AzmF[t, x, y, z];
chiB = chiBf[t, x, y, z];
phaseB = Exp[I*qStar*chiB/hbar];

DtMatter[field_, a0_] := D[field, t] + I*(qStar/hbar)*a0*field;
DxMatter[field_, ax_] := D[field, x] - I*(qStar/hbar)*ax*field;
DyMatter[field_, ay_] := D[field, y] - I*(qStar/hbar)*ay*field;
DzMatter[field_, az_] := D[field, z] - I*(qStar/hbar)*az*field;

psiBg = phaseB*psiB;
A0mg = A0m - D[chiB, t];
Axmg = Axm + D[chiB, x];
Aymg = Aym + D[chiB, y];
Azmg = Azm + D[chiB, z];

resDtMatter = FullSimplify[DtMatter[psiBg, A0mg] - phaseB*DtMatter[psiB, A0m], Assumptions -> $Assumptions];
resDxMatter = FullSimplify[DxMatter[psiBg, Axmg] - phaseB*DxMatter[psiB, Axm], Assumptions -> $Assumptions];
resDyMatter = FullSimplify[DyMatter[psiBg, Aymg] - phaseB*DyMatter[psiB, Aym], Assumptions -> $Assumptions];
resDzMatter = FullSimplify[DzMatter[psiBg, Azmg] - phaseB*DzMatter[psiB, Azm], Assumptions -> $Assumptions];

Print["Residuals for Dt/Dx/Dy/Dz matter covariant derivatives (each should be 0):"];
Print["  Dt -> ", resDtMatter];
Print["  Dx -> ", resDxMatter];
Print["  Dy -> ", resDyMatter];
Print["  Dz -> ", resDzMatter];

If[And @@ (zeroQ /@ {resDtMatter, resDxMatter, resDyMatter, resDzMatter}),
  Print["OK: Paper VII matter-gauge convention is internally consistent."],
  Print["WARNING: Matter-gauge convention check failed (unexpected)."]
];

(* ===================================================================== *)
(* (B) 3+1 Lorentz invariance (effective brane Maxwell sector)            *)
(* ===================================================================== *)

Print["\n--- (B) 3+1 Lorentz invariance check (boost along x) ---\n"];

ClearAll[betaB];
$Assumptions = Join[
  $Assumptions,
  {
    Element[betaB, Reals],
    0 < betaB < 1
  }
];

(* Update helpers with updated assumptions *)
zeroQ[expr_] := TrueQ[FullSimplify[expr == 0, Assumptions -> $Assumptions]];
zeroMatQ[m_] := And @@ Flatten[Map[zeroQ, m, {2}]];

(* Explicit gamma as function of beta (avoids weak simplification). *)
gammaB := 1/Sqrt[1 - betaB^2];

(* Contravariant Lorentz boost matrix (x-direction), acting on {t,x,y,z} *)
LambdaUp = {
  {gammaB, -betaB*gammaB, 0, 0},
  {-betaB*gammaB, gammaB, 0, 0},
  {0, 0, 1, 0},
  {0, 0, 0, 1}
};

(* Covector transform matrix should be Inverse[LambdaUp]^T *)
LambdaDownFromInv = FullSimplify[Transpose[Inverse[LambdaUp]], Assumptions -> $Assumptions];

(* For a Lorentz transform: Inverse[Lambda]^T = eta . Lambda . eta *)
LambdaDownFromEta = FullSimplify[eta4 . LambdaUp . eta4, Assumptions -> $Assumptions];

checkDown = FullSimplify[LambdaDownFromInv - LambdaDownFromEta, Assumptions -> $Assumptions];

Print["Covector transform check: (Inverse[LambdaUp]^T) - (eta LambdaUp eta) should be 0:"];
Print[checkDown];
If[zeroMatQ[checkDown],
  Print["OK: Covector transform matrix matches the Lorentz form eta Lambda eta."],
  Print["WARNING: Covector transform matrices differ (unexpected)."]
];

LambdaDown = LambdaDownFromEta;

(* Sanity: Lorentz condition Lambda^T eta Lambda = eta *)
lorentzResidual = FullSimplify[Transpose[LambdaUp].eta4.LambdaUp - eta4, Assumptions -> $Assumptions];

Print["\nLorentz condition residual (should be zero matrix):"];
Print[lorentzResidual];

If[zeroMatQ[lorentzResidual],
  Print["OK: LambdaUp is Lorentz (Lambda^T eta Lambda = eta)."],
  Print["WARNING: Lorentz condition failed (unexpected)."]
];

(* Build a generic field tensor in terms of E and B (c = 1 units). *)
ClearAll[Ex, Ey, Ez, Bx, By, Bz];
Flow4 = {
  {0,  Ex,  Ey,  Ez},
  {-Ex, 0,  -Bz, By},
  {-Ey, Bz, 0,  -Bx},
  {-Ez, -By, Bx, 0}
};

(* Invariant scalar: F_{mu nu} F^{mu nu} *)
inv0 = FFScalar[Flow4, eta4];

Print["\nInvariant scalar F_{mu nu} F^{mu nu} (should reduce to 2(B^2 - E^2) in c=1):"];
Print[inv0];

(* Transform F_{mu nu} as a rank-2 covariant tensor *)
Flow4p = FullSimplify[LambdaDown . Flow4 . Transpose[LambdaDown], Assumptions -> $Assumptions];
invP = FFScalar[Flow4p, eta4];

Print["Invariant after boost (should match):"];
Print[invP];

diffInv = FullSimplify[invP - inv0, Assumptions -> $Assumptions];
Print["Difference invP - inv0 (should be 0):"];
Print[diffInv];

If[zeroQ[diffInv],
  Print["OK: F_{mu nu} F^{mu nu} is Lorentz invariant under the boost."],
  Print["WARNING: Lorentz invariance check failed (unexpected)."]
];

(* Scalar coupling invariance: A_mu J^mu *)
ClearAll[A0s, Axs, Ays, Azs, J0s, Jxs, Jys, Jzs];
Aco = {A0s, Axs, Ays, Azs}; (* covector components *)
Jcon = {J0s, Jxs, Jys, Jzs}; (* contravariant components *)

AcoP = FullSimplify[LambdaDown . Aco, Assumptions -> $Assumptions];
JconP = FullSimplify[LambdaUp . Jcon, Assumptions -> $Assumptions];

scal0 = FullSimplify[Aco . Jcon, Assumptions -> $Assumptions];
scalP = FullSimplify[AcoP . JconP, Assumptions -> $Assumptions];

diffScal = FullSimplify[scalP - scal0, Assumptions -> $Assumptions];

Print["\nScalar coupling check (A_mu J^mu):"];
Print["Difference (A'·J') - (A·J) should be 0:"];
Print[diffScal];

If[zeroQ[diffScal],
  Print["OK: A_mu J^mu is Lorentz invariant when A transforms as a covector and J as a vector."],
  Print["WARNING: Scalar coupling invariance check failed (unexpected)."]
];

Print["\n========== End symmetry checks v2 ==========\n"];

(*"
Output:


========================================
Paper VIII add-on v2: symmetry checks (gauge + brane Lorentz)
========================================

--- Objects ---
HoldForm[coords5] -> {t, x, y, z, w}
HoldForm[eta5] -> {{-1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}}
HoldForm[Zprof[w]] -> E^(-(w^2/lambdaConf^2))

--- (A) Gauge invariance check (5D) ---

Delta L (gauge-invariant sector) should be a total derivative plus chi*(d·J):
deltaLint5 = -(Jw[t, x, y, z, w]*Derivative[0, 0, 0, 0, 1][chi][t, x, y, z, w]) - Jz[t, x, y, z, w]*Derivative[0, 0, 0, 1, 0][chi][t, x, y, z, w] - Jy[t, x, y, z, w]*Derivative[0, 0, 1, 0, 0][chi][t, x, y, z, w] - Jx[t, x, y, z, w]*Derivative[0, 1, 0, 0, 0][chi][t, x, y, z, w] - J0[t, x, y, z, w]*Derivative[1, 0, 0, 0, 0][chi][t, x, y, z, w]
Residual (deltaLint5 - (-d(chi J) + chi dJ)) should be 0:
0
OK: Gauge invariance holds up to a boundary term when d_M J^M = 0.

Check: DivA(A + d chi) - DivA(A) - Box(chi) should be 0:
0
OK: DivA shifts by Box(chi), as expected.

Note: The gauge-fixing term is not gauge invariant; this is expected and harmless.

--- (A2) Paper VII brane matter gauge convention ---

Residuals for Dt/Dx/Dy/Dz matter covariant derivatives (each should be 0):
  Dt -> 0
  Dx -> 0
  Dy -> 0
  Dz -> 0
OK: Paper VII matter-gauge convention is internally consistent.

--- (B) 3+1 Lorentz invariance check (boost along x) ---

Covector transform check: (Inverse[LambdaUp]^T) - (eta LambdaUp eta) should be 0:
{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
OK: Covector transform matrix matches the Lorentz form eta Lambda eta.

Lorentz condition residual (should be zero matrix):
{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
OK: LambdaUp is Lorentz (Lambda^T eta Lambda = eta).

Invariant scalar F_{mu nu} F^{mu nu} (should reduce to 2(B^2 - E^2) in c=1):
Bx^2 + By^2 + Bz^2 - Ex^2 - Ey^2 - Ez^2
Invariant after boost (should match):
Bx^2 + By^2 + Bz^2 - Ex^2 - Ey^2 - Ez^2
Difference invP - inv0 (should be 0):
0
OK: F_{mu nu} F^{mu nu} is Lorentz invariant under the boost.

Scalar coupling check (A_mu J^mu):
Difference (A'·J') - (A·J) should be 0:
0
OK: A_mu J^mu is Lorentz invariant when A transforms as a covector and J as a vector.

========== End symmetry checks v2 ==========
"*)
