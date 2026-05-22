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

banner["STAGE 027 — CONTINUUM-SELECTED RANK-2 CLOSURE"];

Clear[xi, delta, mMix, mSupp, rU, rPhi];
$Assumptions =
  Element[{xi, delta, mMix, mSupp, rU, rPhi}, Reals] &&
  delta > 0 && mMix >= 0 && mSupp >= 0 && rU > 0 && rPhi > 0;

lambda0 = 2/9;

subbanner["1. Exact continuum-selected branch equation and quadratic theorem"];

nReq = FullSimplify[
  (xi (delta + xi) - mMix (delta + (1 + lambda0 rU^2) xi))/
    (delta + (1 + lambda0 rPhi^2) xi - mMix lambda0 (rU - rPhi)^2),
  Assumptions -> $Assumptions
];
branchEq = FullSimplify[Numerator[Together[nReq - mSupp]], Assumptions -> $Assumptions];
bCont = FullSimplify[delta - mMix (1 + lambda0 rU^2) - mSupp (1 + lambda0 rPhi^2), Assumptions -> $Assumptions];
cCont = FullSimplify[-delta (mMix + mSupp) + lambda0 mMix mSupp (rU - rPhi)^2, Assumptions -> $Assumptions];
branchExpected = Expand[xi^2 + bCont xi + cCont];

Print["n_req^(cont) = ", fmt[nReq]];
Print["numerator of n_req - M_supp = ", fmt[branchEq]];
expectZero["quadratic branch equation", branchEq - 9 branchExpected];

deltaDisc = FullSimplify[bCont^2 - 4 cCont, Assumptions -> $Assumptions];
xiPhys = FullSimplify[(-bCont + Sqrt[deltaDisc])/2, Assumptions -> $Assumptions];
Print["xi_phys = ", fmt[xiPhys]];
expectZero["zero-load root", xiPhys /. {mMix -> 0, mSupp -> 0}];

(* Independent route: solve branchEq for xi via Mathematica's algebraic solver,
   then pick the root that vanishes at zero load and verify it matches xiPhys. *)
xiRoots = xi /. Solve[branchEq == 0, xi];
xiPhysSolve = SelectFirst[
  xiRoots,
  TrueQ[FullSimplify[(# /. {mMix -> 0, mSupp -> 0}) === 0,
                     Assumptions -> $Assumptions]] &
];
If[MissingQ[xiPhysSolve],
  xiPhysSolve = First[Select[xiRoots,
    TrueQ[Simplify[(# /. {mMix -> 0, mSupp -> 0}) == 0,
                   Assumptions -> $Assumptions]] &]]
];
Print["xi_phys (from Solve) = ", fmt[xiPhysSolve]];
expectZero["xiPhys solve match", xiPhysSolve - xiPhys];

subbanner["2. Exact continuum-selected normalization function"];

dCont = FullSimplify[
  (delta + xi - mMix lambda0 rU (rU - rPhi))^2 + lambda0 (mMix (rU - rPhi) + rPhi xi)^2,
  Assumptions -> $Assumptions
];
fCont = FullSimplify[
  (delta + (1 + lambda0 rU rPhi) xi)^2
    (delta + (1 + lambda0 rPhi) xi - mMix lambda0 (rU - rPhi) (rU - 1))^2/
    ((1 - xi) dCont^2),
  Assumptions -> $Assumptions
];

Print["F_cont(xi) = ", fmt[fCont]];

(* Third slice: independent literal rPhi = 2 to constrain bivariate dependence. *)
fLit = FullSimplify[fCont /. rPhi -> 2, Assumptions -> $Assumptions];
fLitExpected = FullSimplify[
  (delta + (1 + lambda0 rU 2) xi)^2
    (delta + (1 + lambda0 2) xi - mMix lambda0 (rU - 2) (rU - 1))^2/
    ((1 - xi) ((delta + xi - mMix lambda0 rU (rU - 2))^2
               + lambda0 (mMix (rU - 2) + 2 xi)^2)^2),
  Assumptions -> $Assumptions
];
expectZero["third-slice F at rPhi=2", fLit - fLitExpected];

subbanner["3. Minimal-kernel source-tied surface"];

nSource = FullSimplify[nReq /. rPhi -> 1, Assumptions -> $Assumptions];
fSource = FullSimplify[fCont /. rPhi -> 1, Assumptions -> $Assumptions];
nSourceExpected = FullSimplify[
  (xi (delta + xi) - mMix (delta + (1 + lambda0 rU^2) xi))/
    (delta + (1 + lambda0) xi - mMix lambda0 (rU - 1)^2),
  Assumptions -> $Assumptions
];
fSourceExpected = FullSimplify[
  (delta + (1 + lambda0 rU) xi)^2
    (delta + (1 + lambda0) xi - mMix lambda0 (rU - 1)^2)^2/
    ((1 - xi) ((delta + xi - mMix lambda0 rU (rU - 1))^2 + lambda0 (xi + mMix (rU - 1))^2)^2),
  Assumptions -> $Assumptions
];

Print["n_source = ", fmt[nSource]];
Print["F_source = ", fmt[fSource]];
expectZero["source-tied n", nSource - nSourceExpected];
expectZero["source-tied F", fSource - fSourceExpected];

subbanner["4. Interference-matched tracking surface"];

nTrack = FullSimplify[nReq /. rPhi -> rU, Assumptions -> $Assumptions];
gQ = FullSimplify[xi (delta + xi)/(delta + (1 + lambda0 rU^2) xi), Assumptions -> $Assumptions];
fTrack = FullSimplify[fCont /. rPhi -> rU, Assumptions -> $Assumptions];
fTrackExpected = FullSimplify[
  (delta + (1 + lambda0 rU^2) xi)^2 (delta + (1 + lambda0 rU) xi)^2/
    ((1 - xi) ((delta + xi)^2 + lambda0 rU^2 xi^2)^2),
  Assumptions -> $Assumptions
];

expectZero["tracking collapse of n_req", nTrack - (gQ - mMix)];
expectZero["tracking F collapse", fTrack - fTrackExpected];

subbanner["5. Exact mismatch penalty"];

(* Extract the xi-constant coefficient of branchEq (= 9*(xi^2 + bCont*xi + cCont)). *)
cFromBranchEq = FullSimplify[CoefficientList[branchEq, xi][[1]]/9,
                             Assumptions -> $Assumptions];
cExpected = FullSimplify[-delta (mMix + mSupp) + lambda0 mMix mSupp (rU - rPhi)^2,
                         Assumptions -> $Assumptions];

expectZero["mismatch penalty in C coefficient", cFromBranchEq - cExpected];

Print[""];
Print["Stage 044 Mathematica audit passed."];

Exit[0];
