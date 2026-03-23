(* ::Package:: *)

(**************************************************)
(*
========================================
Paper VIII add-on v6: brane matter current conservation (fixed local-phase Noether extraction)
========================================

Goal:
  Provide a referee-grade, *symbolic* derivation that the brane matter model used to
  define Maxwell sources has a conserved U(1) current.

Key improvement vs v5:
  - Use Inactive[D] in the Lagrangian so a local phase probe eps(t,x,y,z) correctly
    generates terms proportional to derivatives of eps when we Activate[...].
  - Extract j^mu from the coefficient of \!\(\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)] eps\)
    and verify it matches the standard minimally-coupled Schrödinger number current.

Conventions:
  - coords4 = {t,x,y,z}
  - psi and psiC are treated as independent fields (psiC plays the role of psi* ).
  - A_\[Mu] is an external background (no need to use Maxwell EOM here).
  - Ufun[rho] is an arbitrary real function of rho = psiC psi (gauge/phase invariant).

Outputs:
  (A) E-L equations (symbolic)
  (B) Noether current from local phase probe (correctly accounting for derivative terms)
  (C) Off-shell Noether identity:  \[PartialD]_\[Mu] j^\[Mu] == (I/\[HBar])(psiC ELpsiC - psi ELpsi)
      hence on-shell continuity holds.

*)
(**************************************************)

Print["\n========================================"]; 
Print["Paper VIII add-on v6: brane matter current conservation (fixed local-phase Noether extraction)"];
Print["========================================\n"];

Needs["VariationalMethods`"];

(* --- Coordinates --- *)
coords4 = {t, x, y, z};

(* --- Fields (as functions) --- *)
psi  = psiF[t, x, y, z];
psiC = psiCF[t, x, y, z];
A0b = A0bF[t, x, y, z];
Axb = AxbF[t, x, y, z];
Ayb = AybF[t, x, y, z];
Azb = AzbF[t, x, y, z];

Print["--- Brane fields ---"]; 
Print["psi(t,x,y,z)  = ", psi];
Print["psiC(t,x,y,z) = ", psiC];
Print["A_\[Mu](t,x,y,z) = {A0,A\_x,A\_y,A\_z} = ", {A0b, Axb, Ayb, Azb}];
Print[""]; 

(* --- Parameters & assumptions --- *)
$Assumptions = {
  Element[{m, qStar, hbar, eStar, etaQ}, Reals], m > 0, hbar > 0, eStar > 0, etaQ^2 == 1
};

Print["Charge/source bookkeeping: qStar = etaQ * eStar and Jpsi^mu = qStar * j^mu."]; 
Print["Interpretation: this brane matter model is a representative localized-defect / background-neutralization closure.\n"]; 

(* Potential: generic function of rho; remains symbolic. *)
Ufun[r_] := UfunSym[r];

(* --- Inactive covariant derivatives (key trick for local phase probe) --- *)
DtI[e_]  := Inactive[D][e, t];
DxI[e_]  := Inactive[D][e, x];
DyI[e_]  := Inactive[D][e, y];
DzI[e_]  := Inactive[D][e, z];

DtPsiI[e_]  := DtI[e] + I*(qStar/hbar)*A0b*e;
DtPsiCI[e_] := DtI[e] - I*(qStar/hbar)*A0b*e;

DxPsiI[e_]  := DxI[e] - I*(qStar/hbar)*Axb*e;
DxPsiCI[e_] := DxI[e] + I*(qStar/hbar)*Axb*e;

DyPsiI[e_]  := DyI[e] - I*(qStar/hbar)*Ayb*e;
DyPsiCI[e_] := DyI[e] + I*(qStar/hbar)*Ayb*e;

DzPsiI[e_]  := DzI[e] - I*(qStar/hbar)*Azb*e;
DzPsiCI[e_] := DzI[e] + I*(qStar/hbar)*Azb*e;

rho = psiC*psi;

(* Decompose L for readability *)
LtimeI  = (I*hbar/2) ( psiC*DtPsiI[psi] - psi*DtPsiCI[psiC] );
LspaceI = -(hbar^2/(2*m)) ( DxPsiCI[psiC]*DxPsiI[psi] + DyPsiCI[psiC]*DyPsiI[psi] + DzPsiCI[psiC]*DzPsiI[psi] );
LpotI   = -Ufun[rho];

LpsiI = LtimeI + LspaceI + LpotI;
Lpsi  = Activate[LpsiI];

Print["----------------------------------------\n"]; 
Print["--- Lagrangian density pieces (Activated for display) ---"]; 
Print["L_time = "]; Print[Activate[LtimeI] // TraditionalForm];
Print["\nL_space = "]; Print[Activate[LspaceI] // TraditionalForm];
Print["\nL_pot = "]; Print[Activate[LpotI] // TraditionalForm];
Print["\nL_total = "]; Print[Lpsi // TraditionalForm];
Print["\n----------------------------------------\n"]; 

(* --- Euler–Lagrange equations (robust extraction) --- *)
EqToExpr[eq_] := If[Head[eq] === Equal, eq[[1]] - eq[[2]], eq];

Print["--- Euler–Lagrange equations from L_total ---"]; 

eomEqnPsi  = EulerEquations[Lpsi, psiF[t, x, y, z], coords4][[1]];
eomEqnPsiC = EulerEquations[Lpsi, psiCF[t, x, y, z], coords4][[1]];

Print["Head[eomEqnPsi]  -> ", Head[eomEqnPsi]];
Print["Head[eomEqnPsiC] -> ", Head[eomEqnPsiC]];

ELpsi  = FullSimplify[EqToExpr[eomEqnPsi],  Assumptions -> $Assumptions];
ELpsiC = FullSimplify[EqToExpr[eomEqnPsiC], Assumptions -> $Assumptions];

Print["\nELpsi  (variation wrt psiF)  == 0 with ELpsi  ="]; Print[ELpsi // TraditionalForm];
Print["\nELpsiC (variation wrt psiCF) == 0 with ELpsiC ="]; Print[ELpsiC // TraditionalForm];

Print["\n----------------------------------------\n"]; 

(* --- Noether current from local phase probe eps(t,x,y,z) --- *)
Print["--- Noether current extracted from a local phase probe eps(t,x,y,z) ---"]; 
Print["We probe: psi -> exp(i s eps/\[HBar]) psi,   psiC -> exp(-i s eps/\[HBar]) psiC, expand to O(s).\n"]; 

s = Unique["s"]; (* expansion parameter *)
epsF = epsN[t, x, y, z];

repPhase = {
  psiF[t, x, y, z]  -> Exp[I*s*epsF/hbar]*psiF[t, x, y, z],
  psiCF[t, x, y, z] -> Exp[-I*s*epsF/hbar]*psiCF[t, x, y, z]
};

(* IMPORTANT: apply phase probe to the Inactive Lagrangian, then Activate so product rule generates d eps terms *)
deltaLraw = Coefficient[Normal@Series[Activate[LpsiI /. repPhase], {s, 0, 1}], s, 1];

depsT = Derivative[1, 0, 0, 0][epsN][t, x, y, z];
depsX = Derivative[0, 1, 0, 0][epsN][t, x, y, z];
depsY = Derivative[0, 0, 1, 0][epsN][t, x, y, z];
depsZ = Derivative[0, 0, 0, 1][epsN][t, x, y, z];

deltaLcol = FullSimplify[
  Collect[deltaLraw, {epsF, depsT, depsX, depsY, depsZ}, Simplify],
  Assumptions -> $Assumptions
];

Print["deltaL (linearized, collected in eps and its derivatives) ="]; 
Print[deltaLcol // TraditionalForm];

(* For a true global U(1) symmetry, deltaL must vanish when eps is constant (deps=0). *)
deltaLconst = FullSimplify[deltaLcol /. {depsT -> 0, depsX -> 0, depsY -> 0, depsZ -> 0}, Assumptions -> $Assumptions];
Print["\nCheck: deltaL for constant eps (set d eps = 0) (should be 0):"]; 
Print[deltaLconst // TraditionalForm];

(* Extract current: deltaL should be - j^mu \[PartialD]_mu eps  (with our eps/\[HBar] normalization). *)
j0N = FullSimplify[-Coefficient[deltaLcol, depsT], Assumptions -> $Assumptions];
jxN = FullSimplify[-Coefficient[deltaLcol, depsX], Assumptions -> $Assumptions];
jyN = FullSimplify[-Coefficient[deltaLcol, depsY], Assumptions -> $Assumptions];
jzN = FullSimplify[-Coefficient[deltaLcol, depsZ], Assumptions -> $Assumptions];

Print["\nNoether density j0 ="]; Print[j0N // TraditionalForm];
Print["\nNoether spatial current j = {jx,jy,jz} ="]; Print[{jxN, jyN, jzN} // TraditionalForm];

(* Expected (minimally coupled) Schr number current *)
rhoExpected = rho;

jExpected = {
  (hbar/(2*I*m)) (psiC*D[psi, x] - psi*D[psiC, x]) - (qStar/m)*Axb*rho,
  (hbar/(2*I*m)) (psiC*D[psi, y] - psi*D[psiC, y]) - (qStar/m)*Ayb*rho,
  (hbar/(2*I*m)) (psiC*D[psi, z] - psi*D[psiC, z]) - (qStar/m)*Azb*rho
};

Print["\nCheck j0 - rhoExpected (should be 0):"]; 
Print[FullSimplify[j0N - rhoExpected, Assumptions -> $Assumptions] // TraditionalForm];

Print["\nCheck j - jExpected (should be {0,0,0}):"]; 
Print[FullSimplify[{jxN, jyN, jzN} - jExpected, Assumptions -> $Assumptions] // TraditionalForm];

Print["\n----------------------------------------\n"]; 

(* --- Off-shell Noether identity (continuity on-shell) --- *)
Print["--- Off-shell Noether identity and on-shell continuity ---"]; 

contResidual = FullSimplify[
  D[j0N, t] + D[jxN, x] + D[jyN, y] + D[jzN, z],
  Assumptions -> $Assumptions
];

noetherCombo = FullSimplify[(I/hbar) (psiC*ELpsiC - psi*ELpsi), Assumptions -> $Assumptions];

idResidual = FullSimplify[contResidual + noetherCombo, Assumptions -> $Assumptions];

Print["Continuity divergence  \[PartialD]_t j0 + \[Del]\[CenterDot]j  ="]; 
Print[contResidual // TraditionalForm];

Print["\nNoether combo  (I/\[HBar]) (psiC ELpsiC - psi ELpsi) ="]; 
Print[noetherCombo // TraditionalForm];

Print["\nIdentity residual (should be 0):"]; 
Print[idResidual // TraditionalForm];

zeroBool = TrueQ[FullSimplify[idResidual == 0, Assumptions -> $Assumptions]];
Print["\nZero-check boolean under $Assumptions: ", zeroBool];
If[!zeroBool,
  Print["WARNING: Identity did not simplify to 0. If needed, try FunctionExpand / stronger assumptions."],
  Print["OK: Off-shell Noether identity verified. On-shell (ELpsi=ELpsiC=0), continuity holds." ]
];

Print["\n----------------------------------------\n"]; 

(* --- Charge current used as Maxwell source (project convention) --- *)
Print["--- Charge current used as Maxwell source ---"]; 

J0 = qStar*j0N;
Jvec = qStar*{jxN, jyN, jzN};

Print["Jpsi^0 = qStar * j0 ="]; Print[J0 // TraditionalForm];
Print["\nJpsi^vec = qStar * j ="]; Print[Jvec // TraditionalForm];

Print["\nCharge continuity:  \[PartialD]_t Jpsi^0 + \[Del]\[CenterDot]Jpsi = qStar (\[PartialD]_t j0 + \[Del]\[CenterDot]j).\n"]; 
Print["Background-neutralization, when needed, belongs in Jext^0 rather than in a fluctuating qStar.\n"]; 

Print["========== End add-on v6 (matter current) =========="];

(*"
Output:

========================================
Paper VIII add-on v6: brane matter current conservation (fixed local-phase Noether extraction)
========================================

--- Brane fields ---
psi(t,x,y,z)  = psiF[t, x, y, z]
psiC(t,x,y,z) = psiCF[t, x, y, z]
A_μ(t,x,y,z) = {A0,A\_x,A\_y,A\_z} = {A0bF[t, x, y, z], AxbF[t, x, y, z], AybF[t, x, y, z], AzbF[t, x, y, z]}

----------------------------------------

--- Lagrangian density pieces (Activated for display) ---
L_time =
TraditionalForm[I/2*hbar*(-(psiF[t, x, y, z]*((-I*q*A0bF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[1, 0, 0, 0][psiCF][t, x, y, z])) + psiCF[t, x, y, z]*((I*q*A0bF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[1, 0, 0, 0][psiF][t, x, y, z]))]

L_space =
TraditionalForm[-1/2*(hbar^2*(((I*q*AzbF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[0, 0, 0, 1][psiCF][t, x, y, z])*((-I*q*AzbF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[0, 0, 0, 1][psiF][t, x, y, z]) + ((I*q*AybF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[0, 0, 1, 0][psiCF][t, x, y, z])*((-I*q*AybF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[0, 0, 1, 0][psiF][t, x, y, z]) + ((I*q*AxbF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[0, 1, 0, 0][psiCF][t, x, y, z])*((-I*q*AxbF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[0, 1, 0, 0][psiF][t, x, y, z])))/m]

L_pot =
TraditionalForm[-UfunSym[psiCF[t, x, y, z]*psiF[t, x, y, z]]]

L_total =
TraditionalForm[-UfunSym[psiCF[t, x, y, z]*psiF[t, x, y, z]] - (hbar^2*(((I*q*AzbF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[0, 0, 0, 1][psiCF][t, x, y, z])*((-I*q*AzbF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[0, 0, 0, 1][psiF][t, x, y, z]) + ((I*q*AybF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[0, 0, 1, 0][psiCF][t, x, y, z])*((-I*q*AybF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[0, 0, 1, 0][psiF][t, x, y, z]) + ((I*q*AxbF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[0, 1, 0, 0][psiCF][t, x, y, z])*((-I*q*AxbF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[0, 1, 0, 0][psiF][t, x, y, z])))/(2*m) + I/2*hbar*(-(psiF[t, x, y, z]*((-I*q*A0bF[t, x, y, z]*psiCF[t, x, y, z])/hbar + Derivative[1, 0, 0, 0][psiCF][t, x, y, z])) + psiCF[t, x, y, z]*((I*q*A0bF[t, x, y, z]*psiF[t, x, y, z])/hbar + Derivative[1, 0, 0, 0][psiF][t, x, y, z]))]

----------------------------------------

--- Euler–Lagrange equations from L_total ---
Head[eomEqnPsi]  -> Times
Head[eomEqnPsiC] -> Times

ELpsi  (variation wrt psiF)  == 0 with ELpsi  =
TraditionalForm[(-2*m*q*A0bF[t, x, y, z]*psiCF[t, x, y, z] - q^2*AxbF[t, x, y, z]^2*psiCF[t, x, y, z] - q^2*AybF[t, x, y, z]^2*psiCF[t, x, y, z] - q^2*AzbF[t, x, y, z]^2*psiCF[t, x, y, z] - 2*m*psiCF[t, x, y, z]*Derivative[1][UfunSym][psiCF[t, x, y, z]*psiF[t, x, y, z]] + I*hbar*q*psiCF[t, x, y, z]*Derivative[0, 0, 0, 1][AzbF][t, x, y, z] + (2*I)*hbar*q*AzbF[t, x, y, z]*Derivative[0, 0, 0, 1][psiCF][t, x, y, z] + hbar^2*Derivative[0, 0, 0, 2][psiCF][t, x, y, z] + I*hbar*q*psiCF[t, x, y, z]*Derivative[0, 0, 1, 0][AybF][t, x, y, z] + (2*I)*hbar*q*AybF[t, x, y, z]*Derivative[0, 0, 1, 0][psiCF][t, x, y, z] + hbar^2*Derivative[0, 0, 2, 0][psiCF][t, x, y, z] + I*hbar*q*psiCF[t, x, y, z]*Derivative[0, 1, 0, 0][AxbF][t, x, y, z] + (2*I)*hbar*q*AxbF[t, x, y, z]*Derivative[0, 1, 0, 0][psiCF][t, x, y, z] + hbar^2*Derivative[0, 2, 0, 0][psiCF][t, x, y, z] - (2*I)*hbar*m*Derivative[1, 0, 0, 0][psiCF][t, x, y, z])/(2*m)]

ELpsiC (variation wrt psiCF) == 0 with ELpsiC =
TraditionalForm[-(q*A0bF[t, x, y, z]*psiF[t, x, y, z]) - (q^2*AxbF[t, x, y, z]^2*psiF[t, x, y, z] + q^2*AybF[t, x, y, z]^2*psiF[t, x, y, z] + q^2*AzbF[t, x, y, z]^2*psiF[t, x, y, z] + 2*m*psiF[t, x, y, z]*Derivative[1][UfunSym][psiCF[t, x, y, z]*psiF[t, x, y, z]] + I*hbar*q*psiF[t, x, y, z]*Derivative[0, 0, 0, 1][AzbF][t, x, y, z] + (2*I)*hbar*q*AzbF[t, x, y, z]*Derivative[0, 0, 0, 1][psiF][t, x, y, z] - hbar^2*Derivative[0, 0, 0, 2][psiF][t, x, y, z] + I*hbar*q*psiF[t, x, y, z]*Derivative[0, 0, 1, 0][AybF][t, x, y, z] + (2*I)*hbar*q*AybF[t, x, y, z]*Derivative[0, 0, 1, 0][psiF][t, x, y, z] - hbar^2*Derivative[0, 0, 2, 0][psiF][t, x, y, z] + I*hbar*q*psiF[t, x, y, z]*Derivative[0, 1, 0, 0][AxbF][t, x, y, z] + (2*I)*hbar*q*AxbF[t, x, y, z]*Derivative[0, 1, 0, 0][psiF][t, x, y, z] - hbar^2*Derivative[0, 2, 0, 0][psiF][t, x, y, z])/(2*m) + I*hbar*Derivative[1, 0, 0, 0][psiF][t, x, y, z]]

----------------------------------------

--- Noether current extracted from a local phase probe eps(t,x,y,z) ---
We probe: psi -> exp(i s eps/ℏ) psi,   psiC -> exp(-i s eps/ℏ) psiC, expand to O(s).

deltaL (linearized, collected in eps and its derivatives) =
TraditionalForm[(Derivative[0, 0, 0, 1][epsN][t, x, y, z]*(psiF[t, x, y, z]*(2*q*AzbF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 0, 0, 1][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 0, 0, 1][psiF][t, x, y, z]) + Derivative[0, 0, 1, 0][epsN][t, x, y, z]*(psiF[t, x, y, z]*(2*q*AybF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 0, 1, 0][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 0, 1, 0][psiF][t, x, y, z]) + Derivative[0, 1, 0, 0][epsN][t, x, y, z]*(psiF[t, x, y, z]*(2*q*AxbF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 1, 0, 0][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 1, 0, 0][psiF][t, x, y, z]) - 2*m*psiCF[t, x, y, z]*psiF[t, x, y, z]*Derivative[1, 0, 0, 0][epsN][t, x, y, z])/(2*m)]

Check: deltaL for constant eps (set d eps = 0) (should be 0):
TraditionalForm[0]

Noether density j0 =
TraditionalForm[psiCF[t, x, y, z]*psiF[t, x, y, z]]

Noether spatial current j = {jx,jy,jz} =
TraditionalForm[{-1/2*(psiF[t, x, y, z]*(2*q*AxbF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 1, 0, 0][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 1, 0, 0][psiF][t, x, y, z])/m, -1/2*(psiF[t, x, y, z]*(2*q*AybF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 0, 1, 0][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 0, 1, 0][psiF][t, x, y, z])/m, -1/2*(psiF[t, x, y, z]*(2*q*AzbF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 0, 0, 1][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 0, 0, 1][psiF][t, x, y, z])/m}]

Check j0 - rhoExpected (should be 0):
TraditionalForm[0]

Check j - jExpected (should be {0,0,0}):
TraditionalForm[{0, 0, 0}]

----------------------------------------

--- Off-shell Noether identity and on-shell continuity ---
Continuity divergence  ∂_t j0 + ∇·j  =
TraditionalForm[(psiF[t, x, y, z]*(-2*q*(AzbF[t, x, y, z]*Derivative[0, 0, 0, 1][psiCF][t, x, y, z] + AybF[t, x, y, z]*Derivative[0, 0, 1, 0][psiCF][t, x, y, z] + AxbF[t, x, y, z]*Derivative[0, 1, 0, 0][psiCF][t, x, y, z]) + I*hbar*(Derivative[0, 0, 0, 2][psiCF][t, x, y, z] + Derivative[0, 0, 2, 0][psiCF][t, x, y, z] + Derivative[0, 2, 0, 0][psiCF][t, x, y, z]) + 2*m*Derivative[1, 0, 0, 0][psiCF][t, x, y, z]) - psiCF[t, x, y, z]*(2*q*(AzbF[t, x, y, z]*Derivative[0, 0, 0, 1][psiF][t, x, y, z] + AybF[t, x, y, z]*Derivative[0, 0, 1, 0][psiF][t, x, y, z] + psiF[t, x, y, z]*(Derivative[0, 0, 0, 1][AzbF][t, x, y, z] + Derivative[0, 0, 1, 0][AybF][t, x, y, z] + Derivative[0, 1, 0, 0][AxbF][t, x, y, z]) + AxbF[t, x, y, z]*Derivative[0, 1, 0, 0][psiF][t, x, y, z]) + I*hbar*(Derivative[0, 0, 0, 2][psiF][t, x, y, z] + Derivative[0, 0, 2, 0][psiF][t, x, y, z] + Derivative[0, 2, 0, 0][psiF][t, x, y, z]) - 2*m*Derivative[1, 0, 0, 0][psiF][t, x, y, z]))/(2*m)]

Noether combo  (I/ℏ) (psiC ELpsiC - psi ELpsi) =
TraditionalForm[(psiF[t, x, y, z]*(2*q*(AzbF[t, x, y, z]*Derivative[0, 0, 0, 1][psiCF][t, x, y, z] + AybF[t, x, y, z]*Derivative[0, 0, 1, 0][psiCF][t, x, y, z] + AxbF[t, x, y, z]*Derivative[0, 1, 0, 0][psiCF][t, x, y, z]) - I*hbar*(Derivative[0, 0, 0, 2][psiCF][t, x, y, z] + Derivative[0, 0, 2, 0][psiCF][t, x, y, z] + Derivative[0, 2, 0, 0][psiCF][t, x, y, z]) - 2*m*Derivative[1, 0, 0, 0][psiCF][t, x, y, z]) + psiCF[t, x, y, z]*(2*q*(AzbF[t, x, y, z]*Derivative[0, 0, 0, 1][psiF][t, x, y, z] + AybF[t, x, y, z]*Derivative[0, 0, 1, 0][psiF][t, x, y, z] + psiF[t, x, y, z]*(Derivative[0, 0, 0, 1][AzbF][t, x, y, z] + Derivative[0, 0, 1, 0][AybF][t, x, y, z] + Derivative[0, 1, 0, 0][AxbF][t, x, y, z]) + AxbF[t, x, y, z]*Derivative[0, 1, 0, 0][psiF][t, x, y, z]) + I*hbar*(Derivative[0, 0, 0, 2][psiF][t, x, y, z] + Derivative[0, 0, 2, 0][psiF][t, x, y, z] + Derivative[0, 2, 0, 0][psiF][t, x, y, z]) - 2*m*Derivative[1, 0, 0, 0][psiF][t, x, y, z]))/(2*m)]

Identity residual (should be 0):
TraditionalForm[0]

Zero-check boolean under $Assumptions: True
OK: Off-shell Noether identity verified. On-shell (ELpsi=ELpsiC=0), continuity holds.

----------------------------------------

--- Charge current used as Maxwell source ---
J0 = q * j0 =
TraditionalForm[q*psiCF[t, x, y, z]*psiF[t, x, y, z]]

Jvec = q * j =
TraditionalForm[{-1/2*(q*(psiF[t, x, y, z]*(2*q*AxbF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 1, 0, 0][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 1, 0, 0][psiF][t, x, y, z]))/m, -1/2*(q*(psiF[t, x, y, z]*(2*q*AybF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 0, 1, 0][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 0, 1, 0][psiF][t, x, y, z]))/m, -1/2*(q*(psiF[t, x, y, z]*(2*q*AzbF[t, x, y, z]*psiCF[t, x, y, z] - I*hbar*Derivative[0, 0, 0, 1][psiCF][t, x, y, z]) + I*hbar*psiCF[t, x, y, z]*Derivative[0, 0, 0, 1][psiF][t, x, y, z]))/m}]

Charge continuity:  ∂_t J0 + ∇·J = q (∂_t j0 + ∇·j).

========== End add-on v6 (matter current) ==========
"*)
