ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
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

banner["STAGE 091 — GROUPED-P2 + STATIC-GEOMETRY DERIVATION"];

Clear[omega, kGeom, kPole, omegaQ];
$Assumptions = Element[{omega, kGeom, kPole, omegaQ}, Reals] && kGeom > 0 && kPole > 0 && omegaQ > 0;

kCons = FullSimplify[kGeom + kPole/(1 - omega^2/omegaQ^2), Assumptions -> $Assumptions];
series = Expand[Normal[Series[kCons, {omega, 0, 4}]]];
k0 = FullSimplify[Coefficient[series, omega, 0], Assumptions -> $Assumptions];
k2 = FullSimplify[Coefficient[series, omega, 2], Assumptions -> $Assumptions];
k4 = FullSimplify[Coefficient[series, omega, 4], Assumptions -> $Assumptions];

Print["K_Q^cons(omega) = ", fmt[kCons]];
Print["Series = ", fmt[series]];
Print["K0 = ", fmt[k0]];
Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];

expectZero["K0 - (Kgeom + Kpole)", k0 - (kGeom + kPole)];
expectZero["K2 - Kpole/OmegaQ^2", k2 - kPole/omegaQ^2];
expectZero["K4 - Kpole/OmegaQ^4", k4 - kPole/omegaQ^4];

branchIdentity = FullSimplify[k0*k4 - 4*k2^2, Assumptions -> $Assumptions];
Print["Branch identity K0*K4 - 4*K2^2 = ", fmt[branchIdentity]];

kGeomSol = FullSimplify[First[Solve[branchIdentity == 0, kGeom]], Assumptions -> $Assumptions];
kGeomForced = FullSimplify[kGeom /. kGeomSol, Assumptions -> $Assumptions];
Print["K_geom forced by branch identity = ", fmt[kGeomForced]];
expectZero["K_geom - 3*K_pole", kGeomForced - 3*kPole];

k0OnBranch = FullSimplify[k0 /. kGeomSol, Assumptions -> $Assumptions];
expectZero["K0 - 4*K_pole on branch", k0OnBranch - 4*kPole];

yHat = FullSimplify[(kCons /. kGeomSol)/k0OnBranch, Assumptions -> $Assumptions];
yHatExpected = FullSimplify[3/4 + (1/4)/(1 - omega^2/omegaQ^2), Assumptions -> $Assumptions];
Print["Normalized module on branch = ", fmt[yHat]];
expectZero["Yhat - [3/4 + 1/4/(1-omega^2/OmegaQ^2)]", yHat - yHatExpected];

rhoAlpha = FullSimplify[k0OnBranch/kGeomForced, Assumptions -> $Assumptions];
zetaReq = FullSimplify[(k0OnBranch - kGeomForced)/kGeomForced, Assumptions -> $Assumptions];
Print["rho_alpha = ", fmt[rhoAlpha]];
Print["zeta_req = ", fmt[zetaReq]];
expectZero["rho_alpha - 4/3", rhoAlpha - 4/3];
expectZero["zeta_req - 1/3", zetaReq - 1/3];

(* Independent derivation: bypass Series + Solve and recombine directly. Start
   from the branch-identity result K_geom = 3 K_pole and verify Yhat via
   Together, not via series expansion. Algebraic path is independent of the
   SymPy Series+Solve route, exercising the same bottom-line identity. *)
kConsBranchDirect = 3*kPole + kPole/(1 - omega^2/omegaQ^2);
k0BranchDirect = 4*kPole;
yHatRecomb = Together[kConsBranchDirect/k0BranchDirect];
yHatTargetRecomb = Together[3/4 + (1/4)/(1 - omega^2/omegaQ^2)];
expectZero["Yhat partial-fraction recombination", yHatRecomb - yHatTargetRecomb];

Print[""];
Print["Stage 091 Mathematica audit passed."];

Exit[0];
