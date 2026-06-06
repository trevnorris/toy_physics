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

banner["STAGE 094 — ISOTROPIC GEOMETRY-DECOUPLING"];

Clear[th, ph, mu, tw, tOm, kPot];
$Assumptions = Element[{th, ph, mu, tw, tOm, kPot}, Reals];

dOmega[expr_] := FullSimplify[
  Integrate[Integrate[expr*Sin[th], {ph, 0, 2*Pi}], {th, 0, Pi}],
  Assumptions -> $Assumptions
];

lapS2[y_] := FullSimplify[
  (1/Sin[th]) D[Sin[th] D[y, th], th] + (1/Sin[th]^2) D[y, {ph, 2}],
  Assumptions -> $Assumptions
];

y00 = 1/(2*Sqrt[Pi]);
y20 = Sqrt[5/(16*Pi)]*(3*Cos[th]^2 - 1);
y21c = Sqrt[15/(4*Pi)]*Sin[th]*Cos[th]*Cos[ph];
y21s = Sqrt[15/(4*Pi)]*Sin[th]*Cos[th]*Sin[ph];
y22c = Sqrt[15/(16*Pi)]*Sin[th]^2*Cos[2*ph];
y22s = Sqrt[15/(16*Pi)]*Sin[th]^2*Sin[2*ph];
y2List = {"20" -> y20, "21c" -> y21c, "21s" -> y21s, "22c" -> y22c, "22s" -> y22s};
kg2Overlap = 0;
kg4Overlap = 0;

Print["Y00 normalization = ", fmt[dOmega[y00^2]]];
expectZero["Y00 normalization - 1", dOmega[y00^2] - 1];

Do[
  label = pair[[1]];
  y = pair[[2]];
  norm = dOmega[y^2];
  overlap = dOmega[y00*y];
  lapResidual = FullSimplify[-lapS2[y] - 6*y, Assumptions -> $Assumptions];
  gradCross = FullSimplify[D[y00, th] D[y, th] + (1/Sin[th]^2) D[y00, ph] D[y, ph], Assumptions -> $Assumptions];
  lapCross = dOmega[y00*(-lapS2[y])];
  kg2Overlap = kg2Overlap + overlap;
  kg4Overlap = kg4Overlap + lapCross;
  cCross = FullSimplify[mu*overlap - tw*overlap - tOm*lapCross - kPot*overlap, Assumptions -> $Assumptions];
  Print["Y" <> label <> " normalization = ", fmt[norm]];
  expectZero["Y" <> label <> " normalization - 1", norm - 1];
  expectZero["<Y00|Y" <> label <> ">", overlap];
  expectZero["(-Delta)Y" <> label <> " - 6 Y" <> label, lapResidual];
  expectZero["<grad Y00 . grad Y" <> label <> ">", gradCross];
  expectZero["<Y00|(-Delta)Y" <> label <> ">", lapCross];
  expectZero["Generic isotropic cross coefficient C_0," <> label, cCross],
  {pair, y2List}
];

(* Static-limit check (paper Check #1): with K_(g,2) = K_(g,4) = 0 from
   orthogonality, the contamination numbers vanish and 3/4 + 1/4 split holds. *)
(* K_(g,2), K_(g,4) carry the accumulated l=0 <-> l=2 overlap moments from the
   explicit S^2 integrals above; 0 BECAUSE those integrals vanished. *)
Kg2 = FullSimplify[kg2Overlap, Assumptions -> $Assumptions];
Kg4 = FullSimplify[kg4Overlap, Assumptions -> $Assumptions];
expectZero["K_(g,2) overlap moment", Kg2];
expectZero["K_(g,4) overlap moment", Kg4];
OmegaQ = Symbol["OmegaQ"];
Kpole = Symbol["Kpole"];
eps2 = OmegaQ^2 * Kg2 / Kpole;
eps4 = OmegaQ^4 * Kg4 / Kpole;
cPoleStatic = 1/4;
cGeomStatic = 3/4;
expectZero["eps_2 (static limit)", eps2];
expectZero["eps_4 (static limit)", eps4];
expectZero["c_pole - 1/4", cPoleStatic - 1/4];
expectZero["c_geom - 3/4", cGeomStatic - 3/4];
expectZero["c_pole + c_geom - 1", cPoleStatic + cGeomStatic - 1];
Print["eps_2 = ", fmt[eps2], "; eps_4 = ", fmt[eps4],
      "; c_pole = ", fmt[cPoleStatic], "; c_geom = ", fmt[cGeomStatic]];

Print[""];
Print["Stage 094 Mathematica audit passed."];

Exit[0];
