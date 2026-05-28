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

banner["STAGE 160 — BARE MIXED-PORT SLIPPAGE THEOREM"];

Clear[rStar, deltaR, dKappa0, dGamma0];
$Assumptions = Element[{rStar, deltaR, dKappa0, dGamma0}, Reals];

kappa0Canon = (1 + rStar)/3;
gamma0Canon = (1 + rStar)/9;
(* Total differential of kappa_W = kappa_0 / (1 + r_c) at the canonical point. *)
deltaKappaW = Together[
  dKappa0/(1 + rStar) - (kappa0Canon/(1 + rStar)^2) * deltaR
];
deltaGammaW = Together[
  dGamma0/(1 + rStar) - (gamma0Canon/(1 + rStar)^2) * deltaR
];

Print["dκ_W = ", fmt[FullSimplify[deltaKappaW]]];
Print["dγ_W = ", fmt[FullSimplify[deltaGammaW]]];

identity = deltaGammaW - (1/3)*deltaKappaW - (dGamma0 - (1/3)*dKappa0)/(1 + rStar);
expectZero["exact compensated-branch slippage identity", identity];

deltaGammaWGate = FullSimplify[(dGamma0 - (1/3)*dKappa0)/(1 + rStar)];
Print["dγ_W under dκ_W = 0 = ", fmt[deltaGammaWGate]];
expectZero["pure-scale harmlessness", deltaGammaWGate /. dGamma0 -> dKappa0/3];

banner["Tangential DtN susceptibility and final defect law"];

Clear[upsilonPi, dSigma0, dS, sigmaStar];
$Assumptions = Element[{upsilonPi, dSigma0, dS, sigmaStar, rc}, Reals];

dPiTan = 0.832409471081635*dSigma0 - 1.16275838754222*dS;
dgWTan = FullSimplify[upsilonPi*dPiTan/(1 + rc)];
deltaQ = FullSimplify[-9*sigmaStar*dgWTan/(1 - sigmaStar)];
nQm1 = FullSimplify[9*sigmaStar*dgWTan/(1 - sigmaStar)];

Print["dPi_tan = ", fmt[dPiTan]];
Print["dγ_W = ", fmt[dgWTan]];
Print["Δ_Q = ", fmt[deltaQ]];
Print["N_Q - 1 = ", fmt[nQm1]];

banner["Carry-forward formulas"];
Print["1) dγ_W - (1/3)dκ_W = (dγ_0 - (1/3)dκ_0)/(1+r_c)"];
Print["2) with dκ_W = 0: dγ_W = (dγ_0 - (1/3)dκ_0)/(1+r_c)"];
Print["3) if dγ_0 = dκ_0/3, then dγ_W = 0"];
Print["4) if dγ_0 - dκ_0/3 = Upsilon_Pi dPi_tan, then"];
Print["   Δ_Q = -9 σ_* Upsilon_Pi dPi_tan / [(1-σ_*)(1+r_c)]"];
Print["   N_Q-1 = +9 σ_* Upsilon_Pi dPi_tan / [(1-σ_*)(1+r_c)]"];

Print[""];
Print["Stage 160 Mathematica audit passed."];

Exit[0];
