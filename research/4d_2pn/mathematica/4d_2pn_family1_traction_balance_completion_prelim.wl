(* ::Package:: *)

ClearAll[section, pass, fail, mu, P2, densityLM, gradThetaDensity, zBase, tProf,
  b0, b2, t0, t2, t4, c0, c2, c4, targets, eqns, sol, zBaseSol, tProfSol, deltaT,
  pProfSol, localChannel, operatorChannel, keys12, residualsLocal, residualsOp,
  baseLegCoeffs, tLegCoeffs, pLegCoeffs, normalizeExpr];

section[s_String] := Print["\n" <> StringRepeat["=", 78] <> "\n" <> s <> "\n" <> StringRepeat["=", 78]];
pass[s_String] := Print["PASS: " <> s];
fail[s_String] := Print["FAIL: " <> s];
normalizeExpr[expr_] := FullSimplify[Together[Expand[expr]]];

$Assumptions = -1 < mu < 1;

section["1) Carried-forward support targets"];

targets = <|
  {0, 0} -> 4/45,
  {1, 1} -> 2/7,
  {1, 0} -> 1/4,
  {2, 0} -> 4/9,
  {2, 1} -> 2/3,
  {2, 2} -> 8/3
|>;

KeyValueMap[Print["K_", #1, " = ", #2] &, targets];
Print["\nK_(0,0)=4/45 is kept as the separate monopole wall channel."];
Print["We solve the exact local traction-balance completion on the l=1,2 support sector."];

section["2) Minimal local traction-balance ansatz"];

P2 = LegendreP[2, mu];
zBase = b0 + b2 mu^2;
tProf = t0 + t2 mu^2 + t4 mu^4;

Print["zBase(mu) = b0 + b2 mu^2"];
Print["t(mu)     = t0 + t2 mu^2 + t4 mu^4"];
Print["\nChannel formula on the l=1,2 support sector:"];
Print["K_(l,m) = <zBase>_(l,m) + (l(l+1)-2) <t>_(l,m)"];


densityLM[1, 1] := 3 (1 - mu^2)/4;
densityLM[1, 0] := 3 mu^2/2;
densityLM[2, 0] := 5 (3 mu^2 - 1)^2/8;
densityLM[2, 1] := 15 mu^2 (1 - mu^2)/4;
densityLM[2, 2] := 15 (1 - mu^2)^2/16;

gradThetaDensity[1, 1] := 3 mu^2/4;
gradThetaDensity[1, 0] := 3 (1 - mu^2)/2;
gradThetaDensity[2, 0] := 45 mu^2 (1 - mu^2)/2;
gradThetaDensity[2, 1] := 15 mu^4 - 15 mu^2 + 15/4;
gradThetaDensity[2, 2] := 15 mu^2 (1 - mu^2)/4;

keys12 = {{1, 1}, {1, 0}, {2, 0}, {2, 1}, {2, 2}};

eqns = Table[
  With[{ell = key[[1]], m = key[[2]], w = densityLM[key[[1]], key[[2]]]},
    normalizeExpr[
      Integrate[Expand[zBase w], {mu, -1, 1}] + (ell (ell + 1) - 2) Integrate[Expand[tProf w], {mu, -1, 1}]
    ] == targets[key]
  ],
  {key, keys12}
];

sol = First@Solve[eqns, {b0, b2, t0, t2, t4}, Rationals];
zBaseSol = normalizeExpr[zBase /. sol];
tProfSol = normalizeExpr[tProf /. sol];

Print["Exact solution:"];
Do[Print["  ", sym, " = ", normalizeExpr[sym /. sol]], {sym, {b0, b2, t0, t2, t4}}];

Print["\nSolved profiles:"];
Print["  zBase(mu) = ", zBaseSol];
Print["  t(mu)     = ", tProfSol];

baseLegCoeffs = Association @ First @ Solve[
   Thread[CoefficientList[Expand[c0 + c2 P2 - zBaseSol], mu] == 0],
   {c0, c2}
];

tLegCoeffs = Association @ First @ Solve[
   Thread[CoefficientList[Expand[c0 + c2 P2 + c4 P2^2 - tProfSol], mu] == 0],
   {c0, c2, c4}
];

Print["\nLegendre-basis form:"];
Print["  zBase(mu) = ", baseLegCoeffs[c0], " + (", baseLegCoeffs[c2], ") P2(mu)"];
Print["  t(mu)     = ", tLegCoeffs[c0], " + (", tLegCoeffs[c2], ") P2(mu) + (", tLegCoeffs[c4], ") P2(mu)^2"];

section["3) Equivalent local wall-energy form"];

Print["For the operator zBase + 1/2{L-2,t}, with L=-Delta_S, the local quadratic form is"];
Print["  E = 1/2 ∫ dΩ [ p(mu) Psi^2 + t(mu) |∇_S Psi|^2 ]"];
Print["with p(mu) = zBase(mu) - 2 t(mu) - 1/2 Delta_S t(mu)."];

deltaT = normalizeExpr[D[(1 - mu^2) D[tProfSol, mu], mu]];
pProfSol = normalizeExpr[zBaseSol - 2 tProfSol - (1/2) deltaT];

Print["Delta_S t(mu) = ", deltaT];
Print["p(mu)         = ", pProfSol];

pLegCoeffs = Association @ First @ Solve[
   Thread[CoefficientList[Expand[c0 + c2 P2 + c4 P2^2 - pProfSol], mu] == 0],
   {c0, c2, c4}
];

Print["p(mu)         = ", pLegCoeffs[c0], " + (", pLegCoeffs[c2], ") P2(mu) + (", pLegCoeffs[c4], ") P2(mu)^2"];

section["4) Verification"];

localChannel[ell_Integer, m_Integer] := Module[{w, gt, gp},
  w = densityLM[ell, m];
  gt = gradThetaDensity[ell, m];
  gp = normalizeExpr[Expand[(m^2/(1 - mu^2)) w]];
  normalizeExpr[Integrate[Expand[pProfSol w + tProfSol (gt + gp)], {mu, -1, 1}]]
];

operatorChannel[ell_Integer, m_Integer] := Module[{w},
  w = densityLM[ell, m];
  normalizeExpr[Integrate[Expand[zBaseSol w], {mu, -1, 1}] + (ell (ell + 1) - 2) Integrate[Expand[tProfSol w], {mu, -1, 1}]]
];

residualsLocal = Table[key -> normalizeExpr[localChannel[key[[1]], key[[2]]] - targets[key]], {key, keys12}];
residualsOp = Table[key -> normalizeExpr[operatorChannel[key[[1]], key[[2]]] - targets[key]], {key, keys12}];

Scan[(Print["Local residual ", First[#], " = ", Last[#]]) &, residualsLocal];
Scan[(Print["Operator residual ", First[#], " = ", Last[#]]) &, residualsOp];

If[AllTrue[residualsLocal, Last[#] === 0 &] &&
   AllTrue[residualsOp, Last[#] === 0 &],
  pass["Exact l=1,2 support-sector reconstruction holds in both operator and local wall-energy forms."],
  fail["Residual check failed."]
];

Print["\nMonopole channel:"];
Print["  K_(0,0) = ", targets[{0, 0}], " is kept separate."];
Print["  The local traction-balance completion derived here supplements that monopole wall mode."];

section["5) Numerical sign structure"];

Print["t(0) = ", N[tProfSol /. mu -> 0, 16]];
Print["t(1) = ", N[tProfSol /. mu -> 1, 16]];
Print["p(0) = ", N[pProfSol /. mu -> 0, 16]];
Print["p(1) = ", N[pProfSol /. mu -> 1, 16]];

section["6) Conclusion"];

Print["Main result:"];
Print["- The passed l=1,2 support-sector data admit an exact minimal local traction-balance completion with one base profile zBase(mu) and one tangential wall-stress profile t(mu)."];
Print["- Both profiles sit in the Family-1 low-order basis {1, mu^2, mu^4}."];
Print["- The resulting local wall-energy form reproduces all dipole and quadrupole support channels exactly, while K_(0,0)=4/45 remains the separate monopole wall channel from the earlier PDE steps."];
