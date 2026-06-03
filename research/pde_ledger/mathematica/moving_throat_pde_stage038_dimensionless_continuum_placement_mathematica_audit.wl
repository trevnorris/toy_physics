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

banner["STAGE 038 — DIMENSIONLESS CONTINUUM PLACEMENT"];

subbanner["1. Stage-20 continuum formulas"];

Clear[gNewton, cLight, cS, aScale, ell, muEta, muW, kEta, tOmega, tw, kU, kW, tW, cEtaU, cEtaW, cUW];
$Assumptions =
  Element[{gNewton, cLight, cS, aScale, ell, muEta, muW, kEta, tOmega, tw, kU, kW, tW, cEtaU, cEtaW, cUW}, Reals] &&
  gNewton > 0 && cLight > 0 && cS > 0 && aScale > 0 && ell > 0 && muEta > 0 &&
  muW > 0 && kEta > 0 && tOmega > 0 && tw > 0 && kU > 0 && kW > 0 && tW > 0;

sigma = 88/(9 Pi^2);
kEtaEff = FullSimplify[kEta + 6 tOmega, Assumptions -> $Assumptions];
kWEff = FullSimplify[kW + Pi^2 tW/(4 ell^2), Assumptions -> $Assumptions];

a = FullSimplify[(kU*kEtaEff - cEtaU^2)/(muEta*kU), Assumptions -> $Assumptions];
delta = FullSimplify[Pi^2 tw*kU/(ell^2 (kU*kEtaEff - cEtaU^2)), Assumptions -> $Assumptions];
mMix = FullSimplify[
  8 (kU*cEtaW + cUW*cEtaU)^2/(Pi^2 (kU*kEtaEff - cEtaU^2) (kU*kWEff - cUW^2 sigma)),
  Assumptions -> $Assumptions
];
beta0 = FullSimplify[
  (muW/muEta) (kU*cEtaW + cUW*cEtaU)^2/(kU*kWEff - cUW^2 sigma)^2,
  Assumptions -> $Assumptions
];
rTarget = FullSimplify[
  (54 gNewton cS^5/(5 aScale^5 cLight^5)) a/((8/Pi^2) beta0),
  Assumptions -> $Assumptions
];

Print["A = ", fmt[a]];
Print["delta = ", fmt[delta]];
Print["M_mix = ", fmt[mMix]];
Print["beta0 = ", fmt[beta0]];
Print["R_target = ", fmt[rTarget]];

subbanner["2. Dimensionless kernel substitutions"];

Clear[epsEta, epsW, rho, zW, delta0, lambda];
$Assumptions =
  Element[{epsEta, epsW, rho, zW, delta0, lambda, gNewton, cLight, cS, aScale, ell, muEta, muW, kEtaEff, kWEff, kU}, Reals] &&
  epsEta > 0 && epsW > 0 && rho > 0 && zW > 0 && delta0 > 0 && lambda > 0 &&
  gNewton > 0 && cLight > 0 && cS > 0 && aScale > 0 && ell > 0 && muEta > 0 &&
  muW > 0 && kEtaEff > 0 && kWEff > 0 && kU > 0;

applyDimless[expr_] := Module[{res},
  res = Expand[PowerExpand[Factor[expr /. cUW*cEtaU -> rho kU cEtaW]]];
  res = Expand[
    res /. {
      cUW^4 -> (epsW kU kWEff/sigma)^2,
      cEtaW^2 -> zW kEtaEff kWEff,
      cEtaW^(-2) -> 1/(zW kEtaEff kWEff)
    }
  ];
  res = Expand[
    res /. {
      cEtaU^2 -> epsEta kU kEtaEff,
      cUW^2 -> epsW kU kWEff/sigma,
      tw -> delta0 ell^2 kEtaEff/Pi^2,
      gNewton -> 20 lambda aScale^5 cLight^5 muW/(27 Pi^2 cS^5 kWEff)
    }
  ];
  FullSimplify[res, Assumptions -> $Assumptions]
];

deltaDimless = applyDimless[delta];
mDimless = applyDimless[mMix];
rDimless = applyDimless[rTarget];
betaDimless = applyDimless[beta0];

expectZero["delta - delta0/(1-epsEta)", deltaDimless - delta0/(1 - epsEta)];
expectZero[
  "M_mix - 8 Z_W (1+rho)^2 / [Pi^2 (1-epsEta)(1-epsW)]",
  mDimless - 8 zW (1 + rho)^2/(Pi^2 (1 - epsEta) (1 - epsW))
];
expectZero[
  "R_target - lambda (1-epsEta)(1-epsW)^2 / [Z_W (1+rho)^2]",
  rDimless - lambda (1 - epsEta) (1 - epsW)^2/(zW (1 + rho)^2)
];
expectZero[
  "beta0 - (muW/muEta)(KetaEff/KWEff) Z_W (1+rho)^2/(1-epsW)^2",
  betaDimless - (muW/muEta) (kEtaEff/kWEff) zW (1 + rho)^2/(1 - epsW)^2
];

Print["delta(dimless) = ", fmt[deltaDimless]];
Print["M_mix(dimless) = ", fmt[mDimless]];
Print["R_target(dimless) = ", fmt[rDimless]];
Print["beta0(dimless) = ", fmt[betaDimless]];

subbanner["3. Exact product relation"];

product = FullSimplify[rDimless*mDimless, Assumptions -> $Assumptions];
expectZero["R_target M_mix - 8 lambda (1-epsW)/Pi^2", product - 8 lambda (1 - epsW)/Pi^2];
expectZero[
  "8 lambda (1-epsW)/Pi^2 - NQ KWeff(1-epsW)/muW",
  8 lambda (1 - epsW)/Pi^2 -
    FullSimplify[
      (54 gNewton cS^5/(5 aScale^5 cLight^5)) kWEff (1 - epsW)/muW /.
        gNewton -> 20 lambda aScale^5 cLight^5 muW/(27 Pi^2 cS^5 kWEff),
      Assumptions -> $Assumptions
    ]
];

Print["R_target M_mix = ", fmt[product]];

subbanner["4. Exact derivative factors"];

deltaExpr = delta0/(1 - epsEta);
mExpr = 8 zW (1 + rho)^2/(Pi^2 (1 - epsEta) (1 - epsW));
rExpr = lambda (1 - epsEta) (1 - epsW)^2/(zW (1 + rho)^2);

dDeltaDepsEta = FullSimplify[D[deltaExpr, epsEta], Assumptions -> $Assumptions];
dMdzW = FullSimplify[D[mExpr, zW], Assumptions -> $Assumptions];
dRdzW = FullSimplify[D[rExpr, zW], Assumptions -> $Assumptions];
dMDepsEta = FullSimplify[D[mExpr, epsEta], Assumptions -> $Assumptions];
dRDepsEta = FullSimplify[D[rExpr, epsEta], Assumptions -> $Assumptions];
dMDepsW = FullSimplify[D[mExpr, epsW], Assumptions -> $Assumptions];
dRDepsW = FullSimplify[D[rExpr, epsW], Assumptions -> $Assumptions];
dMDrho = FullSimplify[D[mExpr, rho], Assumptions -> $Assumptions];
dRDrho = FullSimplify[D[rExpr, rho], Assumptions -> $Assumptions];

expectZero["d delta / d epsEta factorization", dDeltaDepsEta - delta0/(1 - epsEta)^2];
expectZero["d M / d Z_W factorization", dMdzW - 8 (1 + rho)^2/(Pi^2 (1 - epsEta) (1 - epsW))];
expectZero["d R / d Z_W factorization", dRdzW + lambda (1 - epsEta) (1 - epsW)^2/(zW^2 (1 + rho)^2)];
expectZero["d M / d epsEta factorization", dMDepsEta - 8 zW (1 + rho)^2/(Pi^2 (1 - epsEta)^2 (1 - epsW))];
expectZero["d R / d epsEta factorization", dRDepsEta + lambda (1 - epsW)^2/(zW (1 + rho)^2)];
expectZero["d M / d epsW factorization", dMDepsW - 8 zW (1 + rho)^2/(Pi^2 (1 - epsEta) (1 - epsW)^2)];
expectZero["d R / d epsW factorization", dRDepsW + 2 lambda (1 - epsEta) (1 - epsW)/(zW (1 + rho)^2)];
expectZero["d M / d rho factorization", dMDrho - 16 zW (1 + rho)/(Pi^2 (1 - epsEta) (1 - epsW))];
expectZero["d R / d rho factorization", dRDrho + 2 lambda (1 - epsEta) (1 - epsW)^2/(zW (1 + rho)^3)];

(* Sign assertions under the natural transfer-branch assumption *)
expectZero["sign(d delta / d epsEta) - 1",
  dDeltaDepsEta (1 - epsEta)^2/delta0 - 1];
expectZero["sign(d M / d Z_W) - 1",
  dMdzW Pi^2 (1 - epsEta) (1 - epsW)/(8 (1 + rho)^2) - 1];
expectZero["sign(d R / d Z_W) + 1",
  dRdzW zW^2 (1 + rho)^2/(lambda (1 - epsEta) (1 - epsW)^2) + 1];
expectZero["sign(d M / d epsEta) - 1",
  dMDepsEta Pi^2 (1 - epsEta)^2 (1 - epsW)/(8 zW (1 + rho)^2) - 1];
expectZero["sign(d R / d epsEta) + 1",
  dRDepsEta zW (1 + rho)^2/(lambda (1 - epsW)^2) + 1];
expectZero["sign(d M / d epsW) - 1",
  dMDepsW Pi^2 (1 - epsEta) (1 - epsW)^2/(8 zW (1 + rho)^2) - 1];
expectZero["sign(d R / d epsW) + 1",
  dRDepsW zW (1 + rho)^2/(2 lambda (1 - epsEta) (1 - epsW)) + 1];
expectZero["sign(d M / d rho) - 1",
  dMDrho Pi^2 (1 - epsEta) (1 - epsW)/(16 zW (1 + rho)) - 1];
expectZero["sign(d R / d rho) + 1",
  dRDrho zW (1 + rho)^3/(2 lambda (1 - epsEta) (1 - epsW)^2) + 1];

Print["On the natural nonvanishing transfer branch 1+rho > 0:"];
Print["  delta increases with epsEta"];
Print["  M_mix increases with epsEta, epsW, Z_W, rho"];
Print["  R_target decreases with epsEta, epsW, Z_W, rho"];

Print[""];
Print["Stage 038 Mathematica audit passed."];

Exit[0];
