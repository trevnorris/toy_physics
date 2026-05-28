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

banner["STAGE 163 — OFF-FAMILY NORMAL COORDINATE"];

Clear[r, g, dg, dr];
$Assumptions = Element[{r, g, dg, dr}, Reals] && r > 0 && g > 0;

s = Sqrt[1 + r^2];
gMinus = r - s/2;
gPrime = D[gMinus, r];
fComp = 1 + r^2 - 4*(g - r)^2;
rComp = (g - r)^2/(1 + r^2);
deltaPerp = dg - gPrime*dr;

(* Independent route: implicit-function derivative of g_-(r) from F(g,r)=0 *)
gPrimeImplicit = -((D[fComp, r])/(D[fComp, g])) /. g -> gMinus;
expectZero["gPrime matches implicit-function route", gPrime - gPrimeImplicit];

dF = (D[fComp, g] /. g -> gMinus)*dg + (D[fComp, r] /. g -> gMinus)*dr;
dR = (D[rComp, g] /. g -> gMinus)*dg + (D[rComp, r] /. g -> gMinus)*dr;

expectZero["delta F - 4 s delta_perp", dF - 4*s*deltaPerp];
expectZero["delta R + delta_perp/s", dR + deltaPerp/s];
Print["g_-'(r) = ", fmt[FullSimplify[gPrime]]];

Clear[dlnKs, dlnKq, dlnLam, dlnGs, dlnGq];
$Assumptions = Element[{r, dlnKs, dlnKq, dlnLam, dlnGs, dlnGq}, Reals] && r > 0;

deltaR = r*(dlnLam - 1/2*dlnKs - 1/2*dlnKq);
deltaG = gMinus*(dlnGq - dlnGs + 1/2*dlnKs - 1/2*dlnKq);
deltaPerpMicro = FullSimplify[deltaG - gPrime*deltaR];
deltaPerpExpected = FullSimplify[
  gMinus*(dlnGq - dlnGs - dlnLam + dlnKs) + (dlnKs + dlnKq - 2*dlnLam)/(4*s)
];
expectZero["microscopic delta_perp identity", deltaPerpMicro - deltaPerpExpected];
Print["delta_perp microscopic form = ", fmt[deltaPerpExpected]];

(* Independent route: derive delta r, delta g via chain-rule on Log[...] *)
Clear[Ks, Kq, lam, gs, gq, eps];
rExpr  = lam/Sqrt[Ks*Kq];
gExpr  = gq*Sqrt[Ks]/(gs*Sqrt[Kq]);
pertRule = {
  Ks -> Ks (1 + eps dlnKs),
  Kq -> Kq (1 + eps dlnKq),
  lam -> lam (1 + eps dlnLam),
  gs -> gs (1 + eps dlnGs),
  gq -> gq (1 + eps dlnGq)
};
deltaRSeries = Coefficient[Series[rExpr /. pertRule, {eps, 0, 1}] // Normal, eps];
deltaGSeries = Coefficient[Series[gExpr /. pertRule, {eps, 0, 1}] // Normal, eps];
deltaRSubst = FullSimplify[deltaRSeries] /. {lam -> r*Sqrt[Ks*Kq]};
deltaGSubst = FullSimplify[deltaGSeries] /. {gq -> gMinus*gs*Sqrt[Kq]/Sqrt[Ks]};
expectZero["delta r series matches hand form", deltaRSubst - r*(dlnLam - dlnKs/2 - dlnKq/2)];
expectZero["delta g series matches hand form", deltaGSubst - gMinus*(dlnGq - dlnGs + dlnKs/2 - dlnKq/2)];
deltaPerpSeries = FullSimplify[deltaGSubst - gPrime*deltaRSubst];
expectZero["microscopic delta_perp via series route",
  deltaPerpSeries -
    (gMinus*(dlnGq - dlnGs - dlnLam + dlnKs) + (dlnKs + dlnKq - 2*dlnLam)/(4*s))];

Clear[Ks, Kq, lam, gs, gq, eps, pertRule, rExpr, gExpr, deltaRSeries, deltaGSeries, deltaRSubst, deltaGSubst, deltaPerpSeries];

Clear[sigmaStar, dkapW, dgamW];
$Assumptions = Element[{r, sigmaStar, dkapW, dgamW, deltaPerp}, Reals] && r > 0;

deltaC = 16*sigmaStar*deltaPerp/s;
dE2 = (deltaC - 9*sigmaStar*dkapW)/(27*(1 - sigmaStar));
dE4 = (5*deltaC - 72*sigmaStar*dkapW)/(243*(1 - sigmaStar));
deltaQ = (deltaC - 27*sigmaStar*dgamW)/(3*(1 - sigmaStar));

Print["delta C = ", fmt[FullSimplify[deltaC]]];
Print["delta E2 = ", fmt[FullSimplify[dE2]]];
Print["delta E4 = ", fmt[FullSimplify[dE4]]];
Print["Delta_Q = ", fmt[FullSimplify[deltaQ]]];
expectZero["delta C - 4 sigma_star deltaF/(1+r^2)", deltaC - sigmaStar*dF*4/(1 + r^2)];

Clear[sigma0, sStar, dSigma0, dS];
$Assumptions = Element[{r, sigma0, sStar, dSigma0, dS, deltaPerp}, Reals] && r > 0;

rStar = 1/4;
dRFromPerp = -deltaPerp/s;
dPi = Expand[(1 - rStar*sStar)*dSigma0 - sigma0*(rStar*dS + sStar*dRFromPerp)];
dPiExpected = Expand[(1 - sStar/4)*dSigma0 - sigma0*dS/4 + sigma0*sStar*deltaPerp/s];
expectZero["delta Pi tangent/normal split", dPi - dPiExpected];
Print["delta Pi = ", fmt[FullSimplify[dPiExpected]]];

rf1 = SetPrecision[1.77799353547498, 30];
gf1 = SetPrecision[0.758035078944663, 30];
sigma0Can = SetPrecision[4.651033550168876, 30];
sCan = SetPrecision[0.6703621156734617, 30];
sf1 = Sqrt[1 + rf1^2];

banner["Family-1 numerical coefficients"];
Print["4 sqrt(1+r_*^2) = ", fmt[N[4*sf1, 20]]];
Print["-1/sqrt(1+r_*^2) = ", fmt[N[-1/sf1, 20]]];
Print["g_* = ", fmt[gf1]];
Print["1/(4 sqrt(1+r_*^2)) = ", fmt[N[1/(4*sf1), 20]]];
Print["Sigma0_can S_can / sqrt(1+r_*^2) = ", fmt[N[sigma0Can*sCan/sf1, 20]]];
Print["16 / sqrt(1+r_*^2) = ", fmt[N[16/sf1, 20]]];

Print[""];
Print["Carry-forward formulas:"];
Print["  delta_perp = delta g - g_-'(r_*) delta r"];
Print["  delta F = 4 sqrt(1+r_*^2) delta_perp"];
Print["  delta R_q = -delta_perp/sqrt(1+r_*^2)"];
Print["  delta_perp = g_* dln[(g_q K_s)/(g_s lam)] + [4 sqrt(1+r_*^2)]^{-1} dln[(K_s K_q)/lam^2]"];
Print["  delta C, delta E2, delta E4, Delta_Q are all linear in delta_perp"];

Print[""];
Print["Stage 163 Mathematica audit passed."];

Exit[0];
