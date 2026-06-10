ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 040 — GENERALIZED SELECTED-BRANCH NORMALIZATION"];

Clear[a0, delta, z0, q, eta, xi, rU, eps, alpha];
$Assumptions =
  Element[{a0, delta, z0, q, eta, xi, rU, eps, alpha}, Reals] &&
  a0 > 0 && delta > 0 && z0 > 0 && xi > 0 && rU > 0 && eps > 0 && q != 0;

lambda0 = 2/9;

subbanner["1. Exact 2x2 selected-branch solve with generic loading ratio"];

(* Independent derivation: build the perturbed matrix, solve the         *)
(* characteristic equation det(M - alpha z z^T - lam I) = 0 for alpha    *)
(* at lam = a0 (1 - xi), then derive the eigenvector via NullSpace.     *)
mBase = DiagonalMatrix[{a0, a0 (1 + delta)}];
zVec = {z0, q z0};
mPert[alphaVal_] := mBase - alphaVal Outer[Times, zVec, zVec];
lamMinus = a0 (1 - xi);

charEq = Det[mPert[alpha] - lamMinus IdentityMatrix[2]] == 0;
alphaSol = Solve[charEq, alpha];
alphaReq = FullSimplify[alpha /. alphaSol[[1]], Assumptions -> $Assumptions];
Print["alpha_req (from Det = 0) = ", fmt[alphaReq]];

(* Eigenvector via NullSpace of (M - alpha_req z z^T - lam I).          *)
nsVec = NullSpace[mPert[alphaReq] - lamMinus IdentityMatrix[2]];
eMinusRaw = FullSimplify[nsVec[[1]], Assumptions -> $Assumptions];
(* Normalize so the first component is 1. *)
eMinus = FullSimplify[eMinusRaw/eMinusRaw[[1]], Assumptions -> $Assumptions];
eMinusFromNullSpace = eMinus;
r = eMinus[[2]];
Print["e1/e0 (from NullSpace) = ", fmt[r]];
expectZero["e1/e0 closed form", r - xi q/(delta + xi)];

(* Explicit eigenvalue/eigenvector residual check. *)
(* Convention: perturbed matrix is M - alpha z z^T, with             *)
(* M = DiagonalMatrix[{a0, a0 (1+delta)}] and z = {z0, q z0}.         *)
mBase = DiagonalMatrix[{a0, a0 (1 + delta)}];
zVec = {z0, q z0};
mPerturbed = mBase - alphaReq Outer[Times, zVec, zVec];
eMinus = {1, q xi/(delta + xi)};
eigResidual = FullSimplify[mPerturbed.eMinus - lamMinus eMinus,
  Assumptions -> $Assumptions];
expectZero["eigenvector residual row 0", eigResidual[[1]]];
expectZero["eigenvector residual row 1", eigResidual[[2]]];
eMinus = eMinusFromNullSpace;

subbanner["2. Exact overlap formulas and generalized F,G functions"];

(* Build z and s overlaps from the derived eigenvector, NOT from a      *)
(* posited r ratio. Use s = (1, eta/q) so that s1/s0 * q = eta.         *)
sVec = {1, eta/q};
zDotE = zVec.eMinus;
sDotE = sVec.eMinus;
normESq = eMinus.eMinus;
zOverlapSq = FullSimplify[zDotE^2/(z0^2 normESq), Assumptions -> $Assumptions];
sOverlapSq = FullSimplify[sDotE^2/(1^2 normESq), Assumptions -> $Assumptions];
fGeneral = FullSimplify[(a0/lamMinus) zOverlapSq sOverlapSq, Assumptions -> $Assumptions];
gGeneral = FullSimplify[(z0^2/a0) alphaReq, Assumptions -> $Assumptions];

Print["(z.e_-)^2 / z0^2 = ", fmt[zOverlapSq]];
Print["(s.e_-)^2 / s0^2 = ", fmt[sOverlapSq]];
Print["F_(q,eta) = ", fmt[fGeneral]];
Print["G_q = ", fmt[gGeneral]];

(* Cross-check against the SymPy script's claimed closed forms. *)
fExpected = FullSimplify[
  (delta + (1 + q^2) xi)^2 (delta + (1 + eta) xi)^2/
    ((1 - xi) ((delta + xi)^2 + q^2 xi^2)^2),
  Assumptions -> $Assumptions
];
gExpected = FullSimplify[xi (delta + xi)/(delta + (1 + q^2) xi), Assumptions -> $Assumptions];
expectZero["F_general - expected", fGeneral - fExpected];
expectZero["G_general - expected", gGeneral - gExpected];

subbanner["3. Split-U specialization"];

qU = FullSimplify[-Sqrt[lambda0] rU, Assumptions -> $Assumptions];
etaU = FullSimplify[lambda0 rU, Assumptions -> $Assumptions];
fU = FullSimplify[fExpected /. {q -> qU, eta -> etaU}, Assumptions -> $Assumptions];
gU = FullSimplify[gExpected /. q -> qU, Assumptions -> $Assumptions];

(* fStage18 reproduces the Stage-035 closed-form F(xi, delta) verified  *)
(* in mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl lines 50-61. *)
(* gStage19 reproduces the Stage-036 closed-form G(xi, delta) verified  *)
(* in mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl lines 41-58. *)
(* Keep these literals in sync with the upstream source of truth.      *)
fStage18 = FullSimplify[
  (9 delta + 11 xi)^4/(81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2),
  Assumptions -> $Assumptions
];
gStage19 = FullSimplify[9 xi (delta + xi)/(9 delta + 11 xi), Assumptions -> $Assumptions];

Print["F_U(xi,delta;R_U) = ", fmt[fU]];
Print["G_U(xi,delta;R_U) = ", fmt[gU]];
expectZero["F_U(R_U=1) - Stage035 F", (fU /. rU -> 1) - fStage18];
expectZero["G_U(R_U=1) - Stage036 G", (gU /. rU -> 1) - gStage19];

subbanner["4. Independent cross-check of first-order deformation about flat-U limit"];

qUEps = FullSimplify[-Sqrt[lambda0] (1 + eps), Assumptions -> $Assumptions];
etaUEps = FullSimplify[lambda0 (1 + eps), Assumptions -> $Assumptions];

hF = FullSimplify[(D[fU /. rU -> 1 + eps, eps] /. eps -> 0)/fStage18, Assumptions -> $Assumptions];
hG = FullSimplify[(D[gU /. rU -> 1 + eps, eps] /. eps -> 0)/gStage19, Assumptions -> $Assumptions];

fGeneralEps = fGeneral /. {q -> qUEps, eta -> etaUEps};
gGeneralEps = gGeneral /. q -> qUEps;
hFDirect = FullSimplify[(D[fGeneralEps, eps] /. eps -> 0)/fStage18, Assumptions -> $Assumptions];
hGDirect = FullSimplify[(D[gGeneralEps, eps] /. eps -> 0)/gStage19, Assumptions -> $Assumptions];

Print["H_F (via F_U)        = ", fmt[hF]];
Print["H_F (via F_general)  = ", fmt[hFDirect]];
Print["H_G (via G_U)        = ", fmt[hG]];
Print["H_G (via G_general)  = ", fmt[hGDirect]];

expectZero["H_F cross-check (F_U vs F_general)", hF - hFDirect];
expectZero["H_G cross-check (G_U vs G_general)", hG - hGDirect];

Print[""];
Print["Stage 040 Mathematica audit passed."];

Exit[0];
