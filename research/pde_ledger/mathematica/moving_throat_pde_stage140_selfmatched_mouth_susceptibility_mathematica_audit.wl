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

banner["STAGE 140 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE"];

Clear[lM, a, ell, rhoW, tM, hbar, cSound, mPsi, tHat];
$Assumptions =
  Element[{lM, a, ell, rhoW, tM, hbar, cSound, mPsi, tHat}, Reals] &&
  lM > 0 && a > 0 && ell > 0 && rhoW > 0 && tM > 0 && hbar > 0 && cSound > 0 && mPsi > 0 && tHat > 0;

jS = 4*Pi*a^2*ell/3;
hWall = mPsi*cSound^2/rhoW;
thetaSigma = hWall*jS;
kS = 3*Pi*a^2*hbar^2/(5*mPsi*rhoW*ell);
gS = tM*jS;
sigma0 = FullSimplify[lM*gS^2/(kS*thetaSigma), Assumptions -> $Assumptions];

Print["Sigma_0 = ", fmt[sigma0]];

sigma0Hat = FullSimplify[sigma0 /. tM -> tHat*hbar*cSound/(rhoW*ell*Sqrt[lM]), Assumptions -> $Assumptions];
Print["Sigma_0 in terms of That = ", fmt[sigma0Hat]];
expectZero["Sigma_0_hat - 20 That^2/9", sigma0Hat - (20/9)*tHat^2];

mSNat = SetPrecision[1.6685425296562397, 30];
mSComp = SetPrecision[1.80594111095636, 30];
tHatNat = N[Sqrt[9*mSNat/20], 30];
tHatComp = N[Sqrt[9*mSComp/20], 30];

Print["That_nat = ", fmt[tHatNat]];
Print["That_comp = ", fmt[tHatComp]];
Print["fractional traction enhancement = ", fmt[N[tHatComp/tHatNat - 1, 20]]];

Module[{diff1, diff2, diff3, tol},
  tol = 10^-12;
  diff1 = N[tHatNat - SetPrecision[0.866512630228382, 30], 30];
  diff2 = N[tHatComp - SetPrecision[0.901484054174206, 30], 30];
  diff3 = N[(tHatComp/tHatNat - 1) - SetPrecision[0.0403588161624, 30], 30];
  If[Abs[diff1] < tol, pass["That_nat matches notes"], fail["That_nat matches notes", diff1]];
  If[Abs[diff2] < tol, pass["That_comp matches notes"], fail["That_comp matches notes", diff2]];
  If[Abs[diff3] < tol, pass["fractional enhancement matches notes"], fail["fractional enhancement matches notes", diff3]];
];

Print[""];
Print["Stage 140 Mathematica audit passed."];

Exit[0];
