(*
  Stage 006 Mathematica audit.
  This script independently checks the projected Maxwell vector rearrangements,
  Gaussian projection identities, a concrete moving-throat bulk potential, and
  two sign-mutation guards for the unit 006 projected vector system.
*)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

zeroLikeQ[expr_] := If[ListQ[expr], And @@ (zeroLikeQ /@ expr), TrueQ[expr === 0]];

expectZero[label_String, expr_] := Module[{res},
  res = FullSimplify[expr, Assumptions -> $Assumptions];
  Print[label, " residual = ", fmt[res]];
  If[!zeroLikeQ[res],
    Print["FAIL: ", label]; Exit[1],
    Print["PASS: ", label]
  ];
];

expectNonzero[label_String, expr_] := Module[{res},
  res = FullSimplify[expr, Assumptions -> $Assumptions];
  Print[label, " residual = ", fmt[res]];
  If[zeroLikeQ[res],
    Print["FAIL: ", label]; Exit[1],
    Print["PASS: ", label]
  ];
];

Clear[t, x, y, z, w, ww, mu0];
$Assumptions = Element[{t, x, y, z, mu0}, Reals] && mu0 != 0;

W[w_] := Exp[-w^2]/Sqrt[Pi];
Wp[w_] := D[W[ww], ww] /. ww -> w;
ZZ[w_] := Exp[-w^2];

Pg[f_] := FullSimplify[
  Integrate[W[w] f, {w, -Infinity, Infinity}, Assumptions -> $Assumptions],
  Assumptions -> $Assumptions
];

Pgp[f_] := FullSimplify[
  Integrate[Wp[w] f, {w, -Infinity, Infinity}, Assumptions -> $Assumptions],
  Assumptions -> $Assumptions
];

boundary[f_] := FullSimplify[
  Limit[W[w] f, w -> Infinity] - Limit[W[w] f, w -> -Infinity],
  Assumptions -> $Assumptions
];

spaceCoords = {x, y, z};
braneCoords = {t, x, y, z};
bulkCoords = {t, x, y, z, w};
eps3 = LeviCivitaTensor[3];

divergence3[v_] := Sum[D[v[[i]], spaceCoords[[i]]], {i, {1, 2, 3}}];
curl3[v_] := Table[
  Sum[eps3[[i, j, k]] D[v[[k]], spaceCoords[[j]]], {j, {1, 2, 3}}, {k, {1, 2, 3}}],
  {i, 1, 3}
];
ampereCurl3[v_] := Table[
  Sum[eps3[[k, j, i]] D[v[[k]], spaceCoords[[j]]], {j, {1, 2, 3}}, {k, {1, 2, 3}}],
  {i, 1, 3}
];

Print["STAGE 006 -- projected Maxwell vector audit"];

(* === M1: homogeneous Bianchi rearrangement === *)
Evec = Table[Symbol["E" <> ToString[i]][t, x, y, z], {i, 1, 3}];
Bvec = Table[Symbol["B" <> ToString[i]][t, x, y, z], {i, 1, 3}];

fieldF[0, 0] := 0;
fieldF[i_Integer, 0] /; 1 <= i <= 3 := Evec[[i]];
fieldF[0, i_Integer] /; 1 <= i <= 3 := -Evec[[i]];
fieldF[i_Integer, j_Integer] /; 1 <= i <= 3 && 1 <= j <= 3 :=
  Sum[eps3[[k, i, j]] Bvec[[k]], {k, {1, 2, 3}}];

timeCycleResiduals = Table[
  1/2 Sum[
    eps3[[i, j, k]] (
      D[fieldF[j, k], t] +
      D[fieldF[k, 0], spaceCoords[[j]]] +
      D[fieldF[0, j], spaceCoords[[k]]]
    ),
    {j, {1, 2, 3}}, {k, {1, 2, 3}}
  ] - (D[Bvec[[i]], t] + curl3[Evec][[i]]),
  {i, 1, 3}
];
divCycleResidual = 1/2 Sum[
  eps3[[i, j, k]] D[fieldF[j, k], spaceCoords[[i]]],
  {i, {1, 2, 3}}, {j, {1, 2, 3}}, {k, {1, 2, 3}}
] - divergence3[Bvec];
expectZero["M1 Faraday rearrangement", timeCycleResiduals];
expectZero["M1 divB rearrangement", divCycleResidual];

(* === M2: inhomogeneous projected rearrangement === *)
Dflux = Table[Symbol["Dflux" <> ToString[i]][t, x, y, z], {i, 1, 3}];
Hflux = Table[Symbol["Hflux" <> ToString[i]][t, x, y, z], {i, 1, 3}];
leak = Table[Symbol["leak" <> ToString[i]][t, x, y, z], {i, 0, 3}];
gauge = Table[Symbol["gauge" <> ToString[i]][t, x, y, z], {i, 0, 3}];
rho = rhoProj[t, x, y, z];
J = Table[Symbol["J" <> ToString[i]][t, x, y, z], {i, 1, 3}];

fluxG[0, 0] := 0;
fluxG[i_Integer, 0] /; 1 <= i <= 3 := Dflux[[i]];
fluxG[0, i_Integer] /; 1 <= i <= 3 := -Dflux[[i]];
fluxG[i_Integer, j_Integer] /; 1 <= i <= 3 && 1 <= j <= 3 :=
  Sum[eps3[[k, i, j]] Hflux[[k]], {k, {1, 2, 3}}];

projectedInhom[nu_Integer] := Sum[
  D[fluxG[mu, nu], braneCoords[[mu + 1]]],
  {mu, {0, 1, 2, 3}}
] + leak[[nu + 1]] + gauge[[nu + 1]];

gaussRearranged = projectedInhom[0] - (divergence3[Dflux] + leak[[1]] + gauge[[1]]);
ampereRearranged = Table[
  projectedInhom[i] - (ampereCurl3[Hflux][[i]] - D[Dflux[[i]], t] + leak[[i + 1]] + gauge[[i + 1]]),
  {i, 1, 3}
];
expectZero["M2 Gauss rearrangement", gaussRearranged];
expectZero["M2 Ampere rearrangement", ampereRearranged];

(* === M3: Gaussian projection identities === *)
gaussianSample = ZZ[w] w;
gaussianBoundary = boundary[gaussianSample];
ibpResidual = Pg[D[gaussianSample, w]] - (gaussianBoundary - Pgp[gaussianSample]);
vectorLeakMoment = -Pgp[gaussianSample];
expectZero["M3 boundary discharge", gaussianBoundary];
expectZero["M3 projected IBP relation", ibpResidual];
expectZero["M3 leakage normalization", vectorLeakMoment - 1/(2 Sqrt[2])];

(* === M4: concrete bulk potential === *)
potential = {
  x z (1 + w^2),
  t y (1 + w^2),
  t z (1 + w^2),
  x y (1 + w^2),
  0
};

twoForm[a_Integer, b_Integer] := D[potential[[b + 1]], bulkCoords[[a + 1]]] -
  D[potential[[a + 1]], bulkCoords[[b + 1]]];

spatialTwoForm[v_Integer] := 1/2 Sum[
  eps3[[v, j, k]] twoForm[j, k],
  {j, {1, 2, 3}}, {k, {1, 2, 3}}
];

Eprojected = Table[Pg[twoForm[i, 0]], {i, 1, 3}];
Bprojected = Table[Pg[spatialTwoForm[i]], {i, 1, 3}];
Dprojected = Table[Pg[ZZ[w] twoForm[i, 0]], {i, 1, 3}];
Hprojected = Table[Pg[ZZ[w] spatialTwoForm[i]], {i, 1, 3}];
leakProjected = Table[-Pgp[ZZ[w] twoForm[4, nu]], {nu, 0, 3}];

bulkCurrent[nu_Integer] := FullSimplify[
  (1/mu0) Sum[D[ZZ[w] twoForm[mu, nu], bulkCoords[[mu + 1]]], {mu, {0, 1, 2, 3, 4}}],
  Assumptions -> $Assumptions
];

expectZero["M4 projected Bianchi divB", divergence3[Bprojected]];
expectZero["M4 projected Bianchi Faraday", D[Bprojected, t] + curl3[Eprojected]];
expectZero[
  "M4 projected Gauss law",
  divergence3[Dprojected] + leakProjected[[1]] - mu0 Pg[bulkCurrent[0]]
];
expectZero[
  "M4 projected Ampere law",
  ampereCurl3[Hprojected] - D[Dprojected, t] + Rest[leakProjected] -
    mu0 Table[Pg[bulkCurrent[i]], {i, 1, 3}]
];

(* === M5: adversarial sign mutations === *)
ibpSignMutation = Pg[D[gaussianSample, w]] - (gaussianBoundary + Pgp[gaussianSample]);
faradaySignMutation = D[Bprojected[[1]], t] - D[Eprojected[[3]], y] + D[Eprojected[[2]], z];
expectNonzero["M5 IBP sign mutation", ibpSignMutation];
expectNonzero["M5 concrete Faraday sign mutation", faradaySignMutation];

Print["STATUS: PASS"];
Exit[0];
