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

banner["STAGE 035 — NON-TWIN ASYMMETRY THRESHOLD"];

Clear[piTr, cMix, eps, kW, kPhi0, omega0];
$Assumptions =
  Element[{piTr, cMix, eps, kW, kPhi0, omega0}, Reals] &&
  piTr > 0 && cMix > 0 && 0 < eps < 1 && kW > 0 && kPhi0 > 0 && omega0 > 0;

sReq = FullSimplify[piTr/cMix, Assumptions -> $Assumptions];
zetaReq = FullSimplify[(sReq - 1)/(1 + eps (sReq - 2)), Assumptions -> $Assumptions];

Print["S_req = ", fmt[sReq]];
Print["zeta_req = ", fmt[zetaReq]];
expectZero["zeta_req at Pi_tr = C_mix", zetaReq /. piTr -> cMix];
expectZero["zeta_req at Pi_tr = 2 C_mix minus 1", (zetaReq /. piTr -> 2 cMix) - 1];

dZdPi = FullSimplify[D[zetaReq, piTr], Assumptions -> $Assumptions];
dZExpected = FullSimplify[cMix (1 - eps)/(cMix - eps (2 cMix - piTr))^2, Assumptions -> $Assumptions];
Print["d zeta_req / d Pi_tr = ", fmt[dZdPi]];
expectZero["dzeta_req/dPi - expected", dZdPi - dZExpected];

deltaZeta = FullSimplify[zetaReq - 1, Assumptions -> $Assumptions];
deltaExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(cMix - eps (2 cMix - piTr)), Assumptions -> $Assumptions];
Print["Delta_zeta = ", fmt[deltaZeta]];
expectZero["Delta_zeta - expected", deltaZeta - deltaExpected];

zetaPhys = FullSimplify[kW omega0^2/kPhi0, Assumptions -> $Assumptions];
omega0ReqSq = FullSimplify[zetaReq kPhi0/kW, Assumptions -> $Assumptions];
kPhi0Req = FullSimplify[kW omega0^2/zetaReq, Assumptions -> $Assumptions];

Print["zeta_0^(phys) = ", fmt[zetaPhys]];
Print["Omega_(0,req)^2 = ", fmt[omega0ReqSq]];
Print["K_(phi,0)^(req) = ", fmt[kPhi0Req]];
expectZero["threshold equality at fixed stiffness", (zetaPhys /. omega0^2 -> omega0ReqSq) - zetaReq];
expectZero["threshold equality at fixed overlap", (zetaPhys /. kPhi0 -> kPhi0Req) - zetaReq];

zetaTwin = FullSimplify[zetaPhys /. {omega0 -> 1, kPhi0 -> kW}, Assumptions -> $Assumptions];
omegaReqEqualStiff = FullSimplify[Sqrt[zetaReq], Assumptions -> $Assumptions];
kPhiReqEqualOv = FullSimplify[kW/zetaReq, Assumptions -> $Assumptions];
softFrac = FullSimplify[1 - kPhiReqEqualOv/kW, Assumptions -> $Assumptions];
softExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(piTr - cMix), Assumptions -> $Assumptions];

expectZero["zeta_0^(twin) - 1", zetaTwin - 1];
Print["Required overlap boost at equal stiffness = ", fmt[omegaReqEqualStiff]];
Print["Required softened stiffness at equal overlap = ", fmt[kPhiReqEqualOv]];
Print["Exact softening fraction = ", fmt[softFrac]];
expectZero["softening fraction - expected", softFrac - softExpected];

Print[""];
Print["Stage 035 Mathematica audit passed."];

Exit[0];
