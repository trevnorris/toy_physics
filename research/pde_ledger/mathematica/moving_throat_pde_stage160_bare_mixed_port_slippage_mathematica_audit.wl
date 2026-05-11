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

banner["STAGE 143 — BARE MIXED-PORT SLIPPAGE THEOREM"];

Clear[eps, rc, drc, dk0, dg0];
$Assumptions = Element[{eps, rc, drc, dk0, dg0}, Reals];

k0Star = (1 + rc)/3;
g0Star = (1 + rc)/9;
kW = Normal[Series[(k0Star + eps*dk0)/(1 + rc + eps*drc), {eps, 0, 1}]];
gW = Normal[Series[(g0Star + eps*dg0)/(1 + rc + eps*drc), {eps, 0, 1}]];
dkW = Expand[Coefficient[kW, eps, 1]];
dgW = Expand[Coefficient[gW, eps, 1]];

Print["dκ_W = ", fmt[FullSimplify[dkW]]];
Print["dγ_W = ", fmt[FullSimplify[dgW]]];

identity = dgW - 1/3*dkW - (dg0 - 1/3*dk0)/(1 + rc);
expectZero["exact compensated-branch slippage identity", identity];

dgWGate = FullSimplify[(dg0 - 1/3*dk0)/(1 + rc)];
Print["dγ_W under dκ_W = 0 = ", fmt[dgWGate]];
expectZero["pure-scale harmlessness", dgWGate /. dg0 -> dk0/3];

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
