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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 170 — LINEAR GROUPED-P2 DIRECT OUTLET MAP"];

Clear[
  D0, dD0, dD2, dD4, N0, dN0, sigma, P0, dkappa, dgamma, eps,
  aD0, aD2, aN0, bD0, bD2, bN0
];

$Assumptions =
  Element[
    {
      D0, dD0, dD2, dD4, N0, dN0, sigma, P0, dkappa, dgamma, eps,
      aD0, aD2, aN0, bD0, bD2, bN0
    },
    Reals
  ] && D0 != 0 && N0 != 0 && sigma != 0 && sigma != 1;

u2 = 1/9;
u4 = 4/81;
D2 = -u2*D0;
D4 = -D0/27;

u2Full = -(D2 + eps*dD2)/(D0 + eps*dD0);
u4Full = ((D2 + eps*dD2)^2 - (D0 + eps*dD0)*(D4 + eps*dD4))/(D0 + eps*dD0)^2;
P0Full = (N0 + eps*dN0)/(D0 + eps*dD0);

(* Independent route: first-order coefficient via D[..., eps] /. eps -> 0, a    *)
(* different mechanism than the SymPy series().coeff(eps,1) this once mirrored.  *)
du2 = FullSimplify[(D[u2Full, eps] /. eps -> 0), Assumptions -> $Assumptions];
du4 = FullSimplify[(D[u4Full, eps] /. eps -> 0), Assumptions -> $Assumptions];
dP0 = FullSimplify[(D[P0Full /. N0 -> P0*D0, eps] /. eps -> 0), Assumptions -> $Assumptions];

banner["Linear grouped conservative/output transport"];
expectZero["delta u2 + (dD2 + dD0/9)/D0", du2 + (dD2 + dD0/9)/D0];
expectZero["delta u4 + (dD4 + 2 dD2/9 + 5 dD0/81)/D0", du4 + (dD4 + 2*dD2/9 + 5*dD0/81)/D0];
expectZero["delta P0 - (dN0 - P0 dD0)/D0", dP0 - (dN0 - P0*dD0)/D0];

du2Hyb = -sigma*dkappa/(3*(1 - sigma));
dP0OverP0Hyb = -9*sigma*dgamma/(1 - sigma);

(* Invert directly (no du2sym/dP0sym placeholder idiom — that was a SymPy tell). *)
dkappaFromdu2 = FullSimplify[
  dkappa /. First[Solve[du2 == du2Hyb, dkappa]],
  Assumptions -> $Assumptions
];
dgammaFromdP0 = FullSimplify[
  dgamma /. First[Solve[dP0/P0 == dP0OverP0Hyb, dgamma]],
  Assumptions -> $Assumptions
];

banner["Direct outlet coefficients"];
expectZero[
  "delta kappa_W - 3(1-sigma)(dD2 + dD0/9)/(sigma D0)",
  dkappaFromdu2 - 3*(1 - sigma)*(dD2 + dD0/9)/(sigma*D0)
];
expectZero[
  "delta gamma_W + (1-sigma)(dN0 - P0 dD0)/(9 sigma N0)",
  (dgammaFromdP0 /. P0 -> N0/D0) + (1 - sigma)*(dN0 - (N0/D0)*dD0)/(9*sigma*N0)
];

du4FromHyb = -8*sigma*dkappa/(27*(1 - sigma));
du4FromKappa = FullSimplify[du4FromHyb /. dkappa -> dkappaFromdu2, Assumptions -> $Assumptions];

banner["Even one-parameter consistency"];
expectZero["delta u4 - (8/9) delta u2", du4FromKappa - (8/9)*du2];
relation = FullSimplify[
  dD4 /. First[Solve[du4 == (8/9)*du2, dD4, Reals]],
  Assumptions -> $Assumptions
];
expectZero["delta D4 - (2/3) delta D2 - dD0/27", relation - ((2/3)*dD2 + dD0/27)];

aKappa = FullSimplify[3*(1 - sigma)*(aD2 + aD0/9)/(sigma*D0), Assumptions -> $Assumptions];
bKappa = FullSimplify[3*(1 - sigma)*(bD2 + bD0/9)/(sigma*D0), Assumptions -> $Assumptions];
aGamma = FullSimplify[-(1 - sigma)*(aN0 - P0*aD0)/(9*sigma*N0), Assumptions -> $Assumptions];
bGamma = FullSimplify[-(1 - sigma)*(bN0 - P0*bD0)/(9*sigma*N0), Assumptions -> $Assumptions];
P0Ref = FullSimplify[N0/D0, Assumptions -> $Assumptions];

aKappaFromMap = FullSimplify[Expand[dkappaFromdu2 /. {dD2 -> aD2, dD0 -> aD0}], Assumptions -> $Assumptions];
bKappaFromMap = FullSimplify[Expand[dkappaFromdu2 /. {dD2 -> bD2, dD0 -> bD0}], Assumptions -> $Assumptions];
aGammaFromMap = FullSimplify[Expand[dgammaFromdP0 /. {dN0 -> aN0, dD0 -> aD0}], Assumptions -> $Assumptions];
bGammaFromMap = FullSimplify[Expand[dgammaFromdP0 /. {dN0 -> bN0, dD0 -> bD0}], Assumptions -> $Assumptions];

banner["Grouped trace/anomaly transport"];
Print["a_kappa = ", fmt[aKappa]];
Print["b_kappa = ", fmt[bKappa]];
Print["a_gamma = ", fmt[aGamma]];
Print["b_gamma = ", fmt[bGamma]];
Print["a_kappa from map = ", fmt[aKappaFromMap]];
Print["b_kappa from map = ", fmt[bKappaFromMap]];
Print["a_gamma from map = ", fmt[aGammaFromMap]];
Print["b_gamma from map = ", fmt[bGammaFromMap]];

expectZero["trace kappa coefficient", aKappaFromMap - aKappa];
expectZero["anomaly kappa coefficient", bKappaFromMap - bKappa];
expectZero[
  "trace gamma coefficient",
  FullSimplify[(aGammaFromMap /. P0 -> P0Ref) - (aGamma /. P0 -> P0Ref), Assumptions -> $Assumptions]
];
expectZero[
  "anomaly gamma coefficient",
  FullSimplify[(bGammaFromMap /. P0 -> P0Ref) - (bGamma /. P0 -> P0Ref), Assumptions -> $Assumptions]
];

(* ----------------------------------------------------------------------- *)
(* 5. Weak-axisymmetric branch: signature (1,1/2,-1) and scalar amplitudes  *)
(*    (paper Sec. 5 / card Checks item 2)                                   *)
(* ----------------------------------------------------------------------- *)
(* Lane-scaled grouped bundle defects delta D_(A,n)=eps*lambda_A*D_n^(1),    *)
(* delta N_(A,0)=eps*lambda_A*N_0^(1) with lambda=(1,1/2,-1) feed the SAME    *)
(* linear outlet maps verified in Sec. 2; output must inherit the signature  *)
(* and collapse to kappa1, gamma1 with the closed forms boxed in notes Sec.5. *)
$Assumptions = $Assumptions && Element[{D01, D21, N01, epsL}, Reals];
kappaMap[dD2x_, dD0x_] := 3*(1 - sigma)*(dD2x + dD0x/9)/(sigma*D0);
gammaMap[dN0x_, dD0x_] := -(1 - sigma)*(dN0x - P0*dD0x)/(9*sigma*N0);
kappa1 = 3*(1 - sigma)*(D21 + D01/9)/(sigma*D0);
gamma1 = -(1 - sigma)*(N01 - P0*D01)/(9*sigma*N0);
dkW20 = kappaMap[epsL*1*D21, epsL*1*D01];
dkW21 = kappaMap[epsL*(1/2)*D21, epsL*(1/2)*D01];
dkW22 = kappaMap[epsL*(-1)*D21, epsL*(-1)*D01];
dgW20 = gammaMap[epsL*1*N01, epsL*1*D01];
dgW21 = gammaMap[epsL*(1/2)*N01, epsL*(1/2)*D01];
dgW22 = gammaMap[epsL*(-1)*N01, epsL*(-1)*D01];

banner["Weak-axisymmetric signature (1, 1/2, -1) and scalar amplitudes"];
expectZero["delta kappa_W^(20) - eps kappa1", dkW20 - epsL*kappa1];
expectZero["delta kappa_W^(21) - (eps/2) kappa1", dkW21 - epsL*(1/2)*kappa1];
expectZero["delta kappa_W^(22) + eps kappa1", dkW22 + epsL*kappa1];
expectZero["delta gamma_W^(20) - eps gamma1", dgW20 - epsL*gamma1];
expectZero["delta gamma_W^(21) - (eps/2) gamma1", dgW21 - epsL*(1/2)*gamma1];
expectZero["delta gamma_W^(22) + eps gamma1", dgW22 + epsL*gamma1];
expectZero["kappa signature: 21 = (1/2) 20", dkW21 - (1/2)*dkW20];
expectZero["kappa signature: 22 = -20", dkW22 + dkW20];
expectZero["gamma signature: 21 = (1/2) 20", dgW21 - (1/2)*dgW20];
expectZero["gamma signature: 22 = -20", dgW22 + dgW20];

Print[""];
Print["Carry-forward formulas:"];
Print["  delta kappa_W^(A) = 3(1-sigma_*) [delta D_(A,2) + delta D_(A,0)/9] / (sigma_* D0)"];
Print["  delta gamma_W^(A) = -(1-sigma_*) [delta N_(A,0) - P0 delta D_(A,0)] / (9 sigma_* N0)"];
Print["  direct even-pole consistency: delta D_(A,4) = (2/3) delta D_(A,2) + delta D_(A,0)/27"];
Print["  so the full linear grouped-anisotropy outlet problem collapses to the pair"];
Print["      K_A := delta D_(A,2) + delta D_(A,0)/9"];
Print["      G_A := delta N_(A,0) - P0 delta D_(A,0)"];

Exit[0];
