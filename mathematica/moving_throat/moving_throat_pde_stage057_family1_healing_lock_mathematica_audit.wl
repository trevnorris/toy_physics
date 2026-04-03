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

banner["STAGE 057 — HEALING-LENGTH LOCK AND SUPPORT SCALE"];

Clear[lambdaEll, chiS, kappa, len, ell, mpsi, cSw, hbar];
$Assumptions =
  Element[{lambdaEll, chiS, kappa, len, ell, mpsi, cSw, hbar}, Reals] &&
  lambdaEll > 0 && chiS > 0 && kappa > 0 && len > 0 && ell > 0 && mpsi > 0 && cSw > 0 && hbar > 0;

chiFromHealing = FullSimplify[
  (mpsi*cSw*len/hbar) /. cSw -> hbar/(2*mpsi*ell),
  Assumptions -> $Assumptions
];
chiLock = FullSimplify[chiFromHealing /. (len/ell) -> lambdaEll, Assumptions -> $Assumptions];
kappaLock = FullSimplify[4*chiLock^2 + (4/5)*lambdaEll^2, Assumptions -> $Assumptions];

Print["chi_s (locked) = ", fmt[chiLock]];
Print["kappa(Lambda_ell) = ", fmt[kappaLock]];
expectZero["chi_s - Lambda_ell/2", chiLock - lambdaEll/2];
expectZero["kappa - (9/5) Lambda_ell^2", kappaLock - (9/5)*lambdaEll^2];

lambdaRef = 37;
chiRef = FullSimplify[chiLock /. lambdaEll -> lambdaRef];
kappaRef = FullSimplify[kappaLock /. lambdaEll -> lambdaRef];
alphaRef = FullSimplify[Sqrt[kappaRef]];

Print["Reference branch:"];
Print["Lambda_ell = ", fmt[lambdaRef]];
Print["chi_s = ", fmt[chiRef]];
Print["kappa = ", fmt[kappaRef]];
Print["alpha = ", fmt[alphaRef]];
Print["alpha (numeric) = ", fmt[N[alphaRef, 20]]];

expectZero["chi_ref - 37/2", chiRef - 37/2];
expectZero["kappa_ref - 12321/5", kappaRef - 12321/5];

Print[""];
Print["Stage 057 Mathematica audit passed."];

Exit[0];
