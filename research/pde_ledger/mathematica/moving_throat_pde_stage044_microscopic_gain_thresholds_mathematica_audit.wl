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

banner["STAGE 044 — MICROSCOPIC GAIN THRESHOLDS"];

Clear[chiSigma, lambdaPhi, kX, tX, ell, kappa, eta, peReq, alpha];
$Assumptions =
  Element[{chiSigma, lambdaPhi, kX, tX, ell, kappa, eta, peReq, alpha}, Reals] &&
  chiSigma > 0 && lambdaPhi > 0 && kX > 0 && tX > 0 && ell > 0 &&
  kappa > 0 && eta > 0 && peReq > 0;

alpha = Sqrt[kappa];
delta0 = FullSimplify[eta*(Cosh[alpha] - 1)/(kappa*(alpha*Sinh[alpha] + eta*Cosh[alpha])), Assumptions -> $Assumptions];
deltaInf = FullSimplify[(Cosh[alpha] + (eta/alpha)*Sinh[alpha] - 1)/(alpha*Sinh[alpha] + eta*Cosh[alpha]), Assumptions -> $Assumptions];
xiMicro = FullSimplify[chiSigma*lambdaPhi^2*ell^2/tX, Assumptions -> $Assumptions];
gMicro = FullSimplify[chiSigma*lambdaPhi^2/kX, Assumptions -> $Assumptions];
gFail = FullSimplify[peReq/(kappa*deltaInf), Assumptions -> $Assumptions];
gSuff = FullSimplify[peReq/(kappa*delta0), Assumptions -> $Assumptions];

Print["Delta_0 = ", fmt[delta0]];
Print["Delta_inf = ", fmt[deltaInf]];
Print["Xi_micro = ", fmt[xiMicro]];
Print["G_micro = ", fmt[gMicro]];
Print["G_fail = ", fmt[gFail]];
Print["G_suff = ", fmt[gSuff]];
expectZero["Xi_micro - kappa G_micro", FullSimplify[xiMicro /. tX -> kX*ell^2/kappa, Assumptions -> $Assumptions] - kappa*gMicro];

chiFail = FullSimplify[tX*peReq/(lambdaPhi^2*ell^2*deltaInf), Assumptions -> $Assumptions];
chiSuff = FullSimplify[tX*peReq/(lambdaPhi^2*ell^2*delta0), Assumptions -> $Assumptions];
lambda2Fail = FullSimplify[tX*peReq/(chiSigma*ell^2*deltaInf), Assumptions -> $Assumptions];
lambda2Suff = FullSimplify[tX*peReq/(chiSigma*ell^2*delta0), Assumptions -> $Assumptions];
expectZero["chi_fail from G_fail", FullSimplify[chiFail /. tX -> kX*ell^2/kappa, Assumptions -> $Assumptions] - (kX/lambdaPhi^2)*gFail];
expectZero["chi_suff from G_suff", FullSimplify[chiSuff /. tX -> kX*ell^2/kappa, Assumptions -> $Assumptions] - (kX/lambdaPhi^2)*gSuff];
expectZero["Lambda^2_fail from G_fail", FullSimplify[lambda2Fail /. tX -> kX*ell^2/kappa, Assumptions -> $Assumptions] - (kX/chiSigma)*gFail];
expectZero["Lambda^2_suff from G_suff", FullSimplify[lambda2Suff /. tX -> kX*ell^2/kappa, Assumptions -> $Assumptions] - (kX/chiSigma)*gSuff];

delta0K0 = FullSimplify[Limit[delta0, kappa -> 0, Direction -> "FromAbove"], Assumptions -> eta > 0];
deltaInfK0 = FullSimplify[Limit[deltaInf, kappa -> 0, Direction -> "FromAbove"], Assumptions -> eta > 0];
Print["lim_{kappa->0+} Delta_0 = ", fmt[delta0K0]];
Print["lim_{kappa->0+} Delta_inf = ", fmt[deltaInfK0]];
expectZero["Delta0 soft-support limit - 1/2", delta0K0 - 1/2];
expectZero["Delta_inf soft-support limit - 1", deltaInfK0 - 1];
expectZero["kappa G_fail soft-support limit - Pe_req", FullSimplify[Limit[kappa*gFail, kappa -> 0, Direction -> "FromAbove"], Assumptions -> peReq > 0 && eta > 0] - peReq];
expectZero["kappa G_suff soft-support limit - 2 Pe_req", FullSimplify[Limit[kappa*gSuff, kappa -> 0, Direction -> "FromAbove"], Assumptions -> peReq > 0 && eta > 0] - 2*peReq];

delta0EtaInf = FullSimplify[Limit[delta0, eta -> Infinity], Assumptions -> kappa > 0];
deltaInfEtaInf = FullSimplify[Limit[deltaInf, eta -> Infinity], Assumptions -> kappa > 0];
gFailInf = FullSimplify[Limit[gFail, eta -> Infinity], Assumptions -> kappa > 0 && peReq > 0];
gSuffInf = FullSimplify[Limit[gSuff, eta -> Infinity], Assumptions -> kappa > 0 && peReq > 0];

Print["lim_{eta->Infinity} Delta_0 = ", fmt[delta0EtaInf]];
Print["lim_{eta->Infinity} Delta_inf = ", fmt[deltaInfEtaInf]];
Print["G_fail^(inf) = ", fmt[gFailInf]];
Print["G_suff^(inf) = ", fmt[gSuffInf]];
expectZero["Delta0 eta->Infinity formula", delta0EtaInf - (1 - Sech[alpha])/kappa];
expectZero["Delta_inf eta->Infinity formula", deltaInfEtaInf - Tanh[alpha]/alpha];
expectZero["G_fail^(inf) formula", gFailInf - peReq/(alpha*Tanh[alpha])];
expectZero["G_suff^(inf) formula", gSuffInf - peReq/(1 - Sech[alpha])];

z = Symbol["z"];
gFailInfZ = FullSimplify[gFailInf /. kappa -> z^2, Assumptions -> z > 0 && peReq > 0];
gSuffInfZ = FullSimplify[gSuffInf /. kappa -> z^2, Assumptions -> z > 0 && peReq > 0];
expectZero["stiff-support compliant-mouth limit: sqrt(kappa) G_fail -> Pe_req", FullSimplify[Limit[z*gFailInfZ, z -> Infinity], Assumptions -> peReq > 0] - peReq];
expectZero["stiff-support compliant-mouth limit: G_suff -> Pe_req", FullSimplify[Limit[gSuffInfZ, z -> Infinity], Assumptions -> peReq > 0] - peReq];

Print[""];
Print["Stage 044 Mathematica audit passed."];

Exit[0];
