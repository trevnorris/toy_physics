(* ::Package:: *)

ClearAll[section, sWall, gLocal, gSharp, xTurn, uTurn, nuK, fFull, c0Full, c2Full,
  poleData, maxRelErr];

section[title_String] := Module[{},
  Print["\n" <> StringRepeat["=", 78]];
  Print[title];
  Print[StringRepeat["=", 78]];
];

pass[name_String, cond_] := Print[If[TrueQ[cond], "PASS: ", "FAIL: "], name];

sWall[x_] := (1 + Tanh[x])/2;

gLocal[x_?NumericQ, alpha_?NumericQ, p_Integer] := Module[{arg},
  arg = -x - alpha*sWall[x]^p;
  If[arg > 0, arg^(1/4), 0]
];

gSharp[x_?NumericQ] := If[x < 0, (-x)^(1/4), 0];

xTurn[alpha_?NumericQ, p_Integer] :=
  x /. FindRoot[x + alpha*sWall[x]^p == 0,
    {x, SetPrecision[-alpha - 1, 50], SetPrecision[0.1, 50]},
    WorkingPrecision -> 50, AccuracyGoal -> 20, PrecisionGoal -> 20];

nuK[alpha_?NumericQ, p_Integer, k_Integer] := Module[{xt, integrand, res1, res2},
  xt = xTurn[alpha, p];
  integrand[x_?NumericQ] := x^k (gLocal[x, alpha, p] - gSharp[x]);
  res1 = NIntegrate[integrand[x], {x, -Infinity, xt},
    WorkingPrecision -> 50, AccuracyGoal -> 16, PrecisionGoal -> 16,
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
    MaxRecursion -> 25
  ];
  res2 = NIntegrate[integrand[x], {x, xt, 0},
    WorkingPrecision -> 50, AccuracyGoal -> 16, PrecisionGoal -> 16,
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
    MaxRecursion -> 25
  ];
  N[res1 + res2, 30]
];

uTurn[eps_?NumericQ, alpha_?NumericQ, p_Integer] := N[1 + eps*xTurn[alpha, p], 50];

fFull[u_?NumericQ, eps_?NumericQ, alpha_?NumericQ, p_Integer] := Module[{x, arg},
  x = (Abs[u] - 1)/eps;
  arg = 1 - u^2 - 2 eps alpha sWall[x]^p;
  If[arg > 0, arg^(1/4), 0]
];

c0Full[eps_?NumericQ, alpha_?NumericQ, p_Integer] := Module[{uMax},
  uMax = uTurn[eps, alpha, p];
  Quiet[
    NIntegrate[fFull[u, eps, alpha, p], {u, 0, uMax},
      WorkingPrecision -> 50, AccuracyGoal -> 16, PrecisionGoal -> 16,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
      MaxRecursion -> 25
    ],
    {NIntegrate::slwcon, NIntegrate::ncvb}
  ]
];

c2Full[eps_?NumericQ, alpha_?NumericQ, p_Integer] := Module[{uMax},
  uMax = uTurn[eps, alpha, p];
  Quiet[
    NIntegrate[u^2 fFull[u, eps, alpha, p], {u, 0, uMax},
      WorkingPrecision -> 50, AccuracyGoal -> 16, PrecisionGoal -> 16,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
      MaxRecursion -> 25
    ],
    {NIntegrate::slwcon, NIntegrate::ncvb}
  ]
];

poleData[hNum_?MatrixQ, gNum_List, sigmaStar_?NumericQ, mHat_?MatrixQ] := Module[
  {vals, vecs, order, valsS, vecsS, vecsN, residues, lamEff},
  {vals, vecs} = Eigensystem[LinearSolve[mHat, hNum]];
  order = Ordering[vals];
  valsS = N[vals[[order]], 30];
  vecsS = N[vecs[[order]], 30];
  vecsN = Map[#/Sqrt[#.mHat.#] &, vecsS];
  residues = MapThread[((gNum.#1)^2)/(sigmaStar #2) &, {vecsN, valsS}];
  lamEff = Total[residues]/Total[residues/valsS];
  <|
    "Values" -> valsS,
    "Vectors" -> vecsN,
    "Residues" -> residues,
    "LamEff" -> lamEff
  |>
];

maxRelErr[vals_List, residues_List, lamEff_?NumericQ] := Module[{sGrid, exact, pade},
  sGrid = N@Subdivide[0, 0.1 First[vals], 400];
  exact = Table[Total[MapThread[#1/(1 - s/#2) &, {residues, vals}]], {s, sGrid}];
  pade = Table[Total[residues]/(1 - s/lamEff), {s, sGrid}];
  Max[Abs[(pade - exact)/exact]]
];

Module[
  {
    c0, c2, aCap, bCap,
    epsRep = SetPrecision[0.05, 50], alphaRep = SetPrecision[1, 50], pRep = 2,
    nu0Rep, xStarRep, dcRep, c0As, c2As, mllAs, c0Ex, c2Ex, mllEx,
    a0, lam, rho, beta, sigma, a, l, vGeom, aGeom, pVac, kappaB, eGeom,
    h0, g0, hBar, gBar, delta0, x01, lamEM, rhoEx, betaEx, target, sigmaStar,
    hNum, gNum, baseMHat, baseData, baseErr, rCap, capMHat, capData, capErr,
    omRatioCap, v0Num, rSide, maaSide, rFull, fullMHat, fullData, fullErr,
    omRatioFull
  },

  c0 = FullSimplify[Sqrt[Pi] Gamma[5/4]/(2 Gamma[7/4])];
  c2 = FullSimplify[(2/7) c0];
  aCap = FullSimplify[2^(1/4)/c0];
  bCap = FullSimplify[5 2^(1/4)/(28 c0)];

  section["1) Filled-to-endcap soft-cap scaling"];
  Print["Baseline n = 5 filled-to-endcap TF weight:"];
  Print["  rho_0(u) ∝ (1 - u^2)^(1/4),   u = 2 w / L"];
  Print[""];
  Print["Near the cap, u = 1 + eps_z x, so"];
  Print["  1 - u^2 = eps_z (-2 x - eps_z x^2)."];
  Print[""];
  Print["Therefore a genuine thin-endcap layer must scale as"];
  Print["  V_cap / mu_c = 2 eps_z alpha_z S(x)^p,"];
  Print["not O(1) in mu_c.  This is the key difference from the sidewall."];
  Print[""];
  Print["Reduced local profile:"];
  Print["  g_{alpha,p}(x) = (-x - alpha S(x)^p)_+^(1/4)"];
  Print[""];
  Print["So the endcap correction is weaker than the sidewall by one extra"];
  Print["power of the vanishing TF profile at the cap, namely eps_z^(5/4)."];

  section["2) Exact asymptotic form of the endcap correction"];
  Print["Define the defect moments"];
  Print["  nu_k(alpha,p) = Integral x^k [g_{alpha,p}(x) - (-x)_+^(1/4)] dx"];
  Print[""];
  Print["Then for n = 5:"];
  Print["  c0^cap = c0 + 2^(1/4) eps_z^(5/4) nu_0 + O(eps_z^(9/4))"];
  Print["  c2^cap = c2 + 2^(1/4) eps_z^(5/4) nu_0 + O(eps_z^(9/4))"];
  Print[""];
  Print["Therefore"];
  Print["  rho_eff^(TF+cap) / rho_eff^TF = 1 + A_cap nu_0 eps_z^(5/4) + ..."];
  Print["  M_LL^(TF+cap) = 1/14 + B_cap nu_0 eps_z^(5/4) + ..."];
  Print[""];
  Print["A_cap = ", N[aCap, 30]];
  Print["B_cap = ", N[bCap, 30]];

  nu0Rep = nuK[alphaRep, pRep, 0];
  xStarRep = xTurn[alphaRep, pRep];
  dcRep = 2^(1/4) epsRep^(5/4) nu0Rep;
  c0As = N[c0 + dcRep, 30];
  c2As = N[c2 + dcRep, 30];
  mllAs = N[c2As/(4 c0As), 30];
  c0Ex = c0Full[epsRep, alphaRep, pRep];
  c2Ex = c2Full[epsRep, alphaRep, pRep];
  mllEx = N[c2Ex/(4 c0Ex), 30];

  section["3) Representative steep-cap branch and direct full-profile check"];
  Print["Representative endcap branch:"];
  Print["  eps_z = ", epsRep];
  Print["  alpha_z = ", alphaRep];
  Print["  p_z = ", pRep];
  Print[""];
  Print["Local turning point x_* solving x + alpha S(x)^p = 0:"];
  Print["  x_* = ", N[xStarRep, 30]];
  Print[""];
  Print["Defect moment nu_0 = ", N[nu0Rep, 30]];
  Print[""];
  Print["Asymptotic vs direct full-profile values:"];
  Print["  c0 asymptotic = ", c0As];
  Print["  c0 direct     = ", N[c0Ex, 30]];
  Print["  relative error = ", N[(c0As - c0Ex)/c0Ex, 20]];
  Print[""];
  Print["  c2 asymptotic = ", c2As];
  Print["  c2 direct     = ", N[c2Ex, 30]];
  Print["  relative error = ", N[(c2As - c2Ex)/c2Ex, 20]];
  Print[""];
  Print["  M_LL asymptotic = ", mllAs];
  Print["  M_LL direct     = ", mllEx];
  Print["  relative error  = ", N[(mllAs - mllEx)/mllEx, 20]];
  Print[""];
  Print["So the leading eps_z^(5/4) asymptotic is already very accurate at eps_z = 0.05."];

  (* Carried-forward geometry response *)
  a0 = Unique["a0"]; lam = Unique["Lam"]; rho = Unique["rho"]; beta = Unique["beta"]; sigma = Unique["Sigma"];
  a = Unique["a"]; l = Unique["L"];
  vGeom = (4 Pi/3) a^3 l;
  aGeom = 4 Pi a^2 l + (8 Pi/3) a^3;
  pVac = rho sigma/a0^4;
  kappaB = beta sigma/a0;
  eGeom = Expand[pVac vGeom + sigma/a0^3 aGeom + kappaB a^2/l];
  h0 = D[eGeom, {{a, l}, 2}] /. {a -> a0, l -> lam a0};
  g0 = {D[vGeom, a], D[vGeom, l]} /. {a -> a0, l -> lam a0};
  hBar = FullSimplify[a0^2 h0/sigma];
  gBar = FullSimplify[a0 g0/((vGeom /. {a -> a0, l -> lam a0}))];
  delta0 = FullSimplify[(gBar . Inverse[hBar] . gBar)/sigma];

  x01 = SetPrecision[2.40482555769577276862163187933, 50];
  lamEM = N[Sqrt[2] Pi/x01, 40];
  rhoEx = 1/10; betaEx = 12; target = N[109/280, 40];
  sigmaStar = N[(delta0 /. {a0 -> 1, lam -> lamEM, rho -> rhoEx, beta -> betaEx, sigma -> 1})/target, 40];
  hNum = N[hBar /. {a0 -> 1, lam -> lamEM, rho -> rhoEx, beta -> betaEx}, 30];
  gNum = Flatten @ N[gBar /. {a0 -> 1, lam -> lamEM}, 30];
  v0Num = N[(vGeom /. {a -> 1, l -> lamEM}), 30];

  baseMHat = N[{{3/5, 0}, {0, 1/14}}, 30];
  baseData = poleData[hNum, gNum, sigmaStar, baseMHat];
  baseErr = maxRelErr[baseData["Values"], baseData["Residues"], baseData["LamEff"]];

  rCap = N[c0Ex/c0, 30];
  capMHat = N[{{3/5, 0}, {0, mllEx}}, 30];
  capData = poleData[hNum, gNum, sigmaStar, capMHat];
  capErr = maxRelErr[capData["Values"], capData["Residues"], capData["LamEff"]];
  omRatioCap = N[(capData["Values"]/rCap)/baseData["Values"], 20];

  section["4) Dynamic monopole response with the endcap correction"];
  Print["Baseline TF branch:"];
  Print["  lambda_- = ", baseData["Values"][[1]]];
  Print["  lambda_+ = ", baseData["Values"][[2]]];
  Print["  residues = ", baseData["Residues"]];
  Print["  lambda_eff = ", baseData["LamEff"]];
  Print["  max relative Pade error = ", baseErr];
  Print[""];
  Print["Endcap-corrected branch:"];
  Print["  rho_eff factor R_cap = c0^cap / c0 = ", rCap];
  Print["  mHat_cap = ", capMHat];
  Print[""];
  Print["  lambda_- = ", capData["Values"][[1]]];
  Print["  lambda_+ = ", capData["Values"][[2]]];
  Print["  residues = ", capData["Residues"]];
  Print["  residue fractions = ", N[capData["Residues"]/Total[capData["Residues"]], 20]];
  Print["  lambda_eff = ", capData["LamEff"]];
  Print["  max relative Pade error = ", capErr];
  Print[""];
  Print["Physical pole-squared ratios relative to the sharp-wall TF baseline:"];
  Print["  Omega_-^2 / Omega_-^2(sharp) = ", omRatioCap[[1]]];
  Print["  Omega_+^2 / Omega_+^2(sharp) = ", omRatioCap[[2]]];
  Print[""];
  Print["Static sum of residues = ", Total[capData["Residues"]], " (target 109/280 = ", target, ")"];

  (* Leading separated-order full-wall composite using the carried-forward sidewall branch *)
  rSide = 0.9060975247692787;
  maaSide = 0.5623811549096673;
  rFull = N[rSide*rCap, 30];
  fullMHat = N[{{maaSide, 0}, {0, mllEx}}, 30];
  fullData = poleData[hNum, gNum, sigmaStar, fullMHat];
  fullErr = maxRelErr[fullData["Values"], fullData["Residues"], fullData["LamEff"]];
  omRatioFull = N[(fullData["Values"]/rFull)/baseData["Values"], 20];

  section["5) Leading separated-order full-wall composite branch"];
  Print["Carried-forward sidewall branch (from the previous step):"];
  Print["  R_side  = ", rSide];
  Print["  M_aa    = ", maaSide];
  Print[""];
  Print["New endcap branch:"];
  Print["  R_cap   = ", rCap];
  Print["  M_LL    = ", mllEx];
  Print[""];
  Print["Leading separated-order composite:"];
  Print["  R_full = R_side * R_cap = ", rFull];
  Print["  mHat_full = ", fullMHat];
  Print[""];
  Print["  lambda_- = ", fullData["Values"][[1]]];
  Print["  lambda_+ = ", fullData["Values"][[2]]];
  Print["  residues = ", fullData["Residues"]];
  Print["  residue fractions = ", N[fullData["Residues"]/Total[fullData["Residues"]], 20]];
  Print["  lambda_eff = ", fullData["LamEff"]];
  Print["  max relative Pade error = ", fullErr];
  Print[""];
  Print["Physical pole-squared ratios relative to the sharp-wall TF baseline:"];
  Print["  Omega_-^2 / Omega_-^2(sharp) = ", omRatioFull[[1]]];
  Print["  Omega_+^2 / Omega_+^2(sharp) = ", omRatioFull[[2]]];
  Print[""];
  Print["This is the first full wall-completed monopole breathing branch in the"];
  Print["current Family-1 program, at least to separated leading order in the"];
  Print["sidewall and endcap thickness parameters."];

  section["6) Interpretation"];
  Print["1) The endcap layer is parametrically weaker than the sidewall because the"];
  Print["   filled-to-endcap TF profile already vanishes at the cap."];
  Print[""];
  Print["2) The correct cap scaling is W_cap / mu_c = O(eps_z), and the first"];
  Print["   nontrivial correction is O(eps_z^(5/4)) on the n = 5 branch."];
  Print[""];
  Print["3) Even after adding the cap correction, the monopole wall channel remains"];
  Print["   an excellent positive two-pole Stieltjes response with a one-pole"];
  Print["   reduction that is accurate well below the lower pole."];
  Print[""];
  Print["4) Combining the carried-forward sidewall branch with the new cap branch"];
  Print["   gives a near-final full-wall dynamic monopole response.  The remaining"];
  Print["   gap is the fully coupled sidewall-cap derivation beyond separated order."];

  section["7) Verification"];
  pass[
    "Representative endcap asymptotic c0 correction stays within 1% of the direct full-profile value",
    Abs[N[(c0As - c0Ex)/c0Ex, 20]] < 10^-2
  ];
  pass[
    "Representative endcap asymptotic c2 correction stays within 1% of the direct full-profile value",
    Abs[N[(c2As - c2Ex)/c2Ex, 20]] < 10^-2
  ];
  pass[
    "Representative endcap asymptotic M_LL stays within 1% of the direct full-profile value",
    Abs[N[(mllAs - mllEx)/mllEx, 20]] < 10^-2
  ];
  pass[
    "Endcap-corrected monopole residues stay positive",
    Min[capData["Residues"]] > 0
  ];
  pass[
    "Endcap-corrected monopole residues still sum to 109/280",
    Abs[N[Total[capData["Residues"]] - target, 20]] < 10^-12
  ];
  pass[
    "Endcap-corrected one-pole Pade error stays below 1e-3 on the low-frequency band",
    capErr < 10^-3
  ];
  pass[
    "Separated-order full-wall monopole residues stay positive",
    Min[fullData["Residues"]] > 0
  ];
  pass[
    "Separated-order full-wall monopole residues still sum to 109/280",
    Abs[N[Total[fullData["Residues"]] - target, 20]] < 10^-12
  ];
  pass[
    "Separated-order full-wall one-pole Pade error stays below 1e-3 on the low-frequency band",
    fullErr < 10^-3
  ];

  Family1EndcapSoftwallResults = <|
    "RepresentativeBranch" -> <|
      "epsZ" -> epsRep,
      "alphaZ" -> alphaRep,
      "pZ" -> pRep,
      "nu0" -> nu0Rep,
      "xTurn" -> xStarRep,
      "c0Asymptotic" -> c0As,
      "c0Direct" -> c0Ex,
      "c2Asymptotic" -> c2As,
      "c2Direct" -> c2Ex,
      "MLLAsymptotic" -> mllAs,
      "MLLDirect" -> mllEx
    |>,
    "BaselineData" -> baseData,
    "BaselineError" -> baseErr,
    "CapFactor" -> rCap,
    "CapData" -> capData,
    "CapError" -> capErr,
    "CapPoleRatios" -> omRatioCap,
    "FullFactor" -> rFull,
    "FullData" -> fullData,
    "FullError" -> fullErr,
    "FullPoleRatios" -> omRatioFull
  |>;

  Print["Key exported symbol: Family1EndcapSoftwallResults."];
];

(*"
Output:


==============================================================================
1) Filled-to-endcap soft-cap scaling
==============================================================================
Baseline n = 5 filled-to-endcap TF weight:
  rho_0(u) ∝ (1 - u^2)^(1/4),   u = 2 w / L

Near the cap, u = 1 + eps_z x, so
  1 - u^2 = eps_z (-2 x - eps_z x^2).

Therefore a genuine thin-endcap layer must scale as
  V_cap / mu_c = 2 eps_z alpha_z S(x)^p,
not O(1) in mu_c.  This is the key difference from the sidewall.

Reduced local profile:
  g_{alpha,p}(x) = (-x - alpha S(x)^p)_+^(1/4)

So the endcap correction is weaker than the sidewall by one extra
power of the vanishing TF profile at the cap, namely eps_z^(5/4).

==============================================================================
2) Exact asymptotic form of the endcap correction
==============================================================================
Define the defect moments
  nu_k(alpha,p) = Integral x^k [g_{alpha,p}(x) - (-x)_+^(1/4)] dx

Then for n = 5:
  c0^cap = c0 + 2^(1/4) eps_z^(5/4) nu_0 + O(eps_z^(9/4))
  c2^cap = c2 + 2^(1/4) eps_z^(5/4) nu_0 + O(eps_z^(9/4))

Therefore
  rho_eff^(TF+cap) / rho_eff^TF = 1 + A_cap nu_0 eps_z^(5/4) + ...
  M_LL^(TF+cap) = 1/14 + B_cap nu_0 eps_z^(5/4) + ...

A_cap = 1.3606190066912236182756960790094114888931090988606056803949`30.
B_cap = 0.2429676797662899318349457283945377658737694819393938714991`30.

==============================================================================
3) Representative steep-cap branch and direct full-profile check
==============================================================================
Representative endcap branch:
  eps_z = 0.05000000000000000277555756156289135105907917022705078125`50.
  alpha_z = 1.`50.
  p_z = 2

Local turning point x_* solving x + alpha S(x)^p = 0:
  x_* = -0.1720646550263600131524020623998347021416427990314126189252`30.

Defect moment nu_0 = -0.129717121094554977624409185396406052639451568605557971578`30.

Asymptotic vs direct full-profile values:
  c0 asymptotic = 0.8703719198752304850128091226758594147985878012473424646553`30.
  c0 direct     = 0.8703582632533208735917762478962476635195924316532794385425`30.
  relative error = 0.00001569080513875308904577782490903477`20.

  c2 asymptotic = 0.2460725021866305301402282679398786020169794794779083782703`30.
  c2 direct     = 0.2461235845668365322952753581469143125434274532380978150235`30.
  relative error = -0.00020754768502134498814903878106410022`20.

  M_LL asymptotic = 0.0706802737334131551577823344204524518913281962776625658844`29.69897000433602
  M_LL direct     = 0.0706960555664884249145380213260314649250027532901543262391`30.
  relative error  = -0.00022323498742341026446783374922772638`20.

So the leading eps_z^(5/4) asymptotic is already very accurate at eps_z = 0.05.

==============================================================================
4) Dynamic monopole response with the endcap correction
==============================================================================
Baseline TF branch:
  lambda_- = 5.9255625769268585838222385594374512585639450952200111942518`30.
  lambda_+ = 237.91117494303324065318581486681708803894`30.
  residues = {0.0026280028657021827179081311195504674743576298665594278846`27.357675499791714, 0.3866577114200121029963775831661638182397013620131685287482`29.09691001262318}
  lambda_eff = 188.17695898017126349792842186993058185251`27.949704504763947
  max relative Pade error = 0.00007100969970122965

Endcap-corrected branch:
  rho_eff factor R_cap = c0^cap / c0 = 0.9958113945614278400150315119581863189245647077529426983699`30.
  mHat_cap = {{0.6000000000000000000000000000000000000000000000000000000001`30., 0}, {0, 0.0706960555664884249145380213260314649250027532901543262391`30.}}

  lambda_- = 5.9743204945455909174564671786352624767651340024390619490172`30.
  lambda_+ = 238.41451639640104966595731700374101049039`30.
  residues = {0.0025851672043531042688485264913050708348536525990567167356`27.354082151651184, 0.3867005470813611814454371877944092148795721158225923448356`29.09691001262318}
  residue fractions = {0.00664079648824650637869346254647174159`20., 0.99335920351175349362130653745352825841`20.}
  lambda_eff = 189.46289723924532568038193028445062914108`27.953147544026734
  max relative Pade error = 0.00006983937917456873

Physical pole-squared ratios relative to the sharp-wall TF baseline:
  Omega_-^2 / Omega_-^2(sharp) = 1.01246923727752781965038674301492711931`20.
  Omega_+^2 / Omega_+^2(sharp) = 1.00633079228437766651961604440700012432`20.

Static sum of residues = 0.3892857142857142857142857142857142857144257684216490615712`28.96315333239801 (target 109/280 = 0.3892857142857142857142857142857142857142857142857142857144`40.)

==============================================================================
5) Leading separated-order full-wall composite branch
==============================================================================
Carried-forward sidewall branch (from the previous step):
  R_side  = 0.9060975247692787
  M_aa    = 0.5623811549096673

New endcap branch:
  R_cap   = 0.9958113945614278400150315119581863189245647077529426983699`30.
  M_LL    = 0.0706960555664884249145380213260314649250027532901543262391`30.

Leading separated-order composite:
  R_full = R_side * R_cap = 0.9023022397491534
  mHat_full = {{0.5623811549096673, 0}, {0, 0.0706960555664884249145380213260314649250027532901543262391`30.}}

  lambda_- = 6.052355928868017
  lambda_+ = 251.0829615183736
  residues = {0.0028552853638614327, 0.3864304289218528}
  residue fractions = {0.007334677998910103, 0.9926653220010898}
  lambda_eff = 193.59559675191792
  max relative Pade error = 0.00007722575556500872

Physical pole-squared ratios relative to the sharp-wall TF baseline:
  Omega_-^2 / Omega_-^2(sharp) = 1.1319906403274251
  Omega_+^2 / Omega_+^2(sharp) = 1.1696350262065394

This is the first full wall-completed monopole breathing branch in the
current Family-1 program, at least to separated leading order in the
sidewall and endcap thickness parameters.

==============================================================================
6) Interpretation
==============================================================================
1) The endcap layer is parametrically weaker than the sidewall because the
   filled-to-endcap TF profile already vanishes at the cap.

2) The correct cap scaling is W_cap / mu_c = O(eps_z), and the first
   nontrivial correction is O(eps_z^(5/4)) on the n = 5 branch.

3) Even after adding the cap correction, the monopole wall channel remains
   an excellent positive two-pole Stieltjes response with a one-pole
   reduction that is accurate well below the lower pole.

4) Combining the carried-forward sidewall branch with the new cap branch
   gives a near-final full-wall dynamic monopole response.  The remaining
   gap is the fully coupled sidewall-cap derivation beyond separated order.

==============================================================================
7) Verification
==============================================================================
PASS: Representative endcap asymptotic c0 correction stays within 1% of the direct full-profile value
PASS: Representative endcap asymptotic c2 correction stays within 1% of the direct full-profile value
PASS: Representative endcap asymptotic M_LL stays within 1% of the direct full-profile value
PASS: Endcap-corrected monopole residues stay positive
PASS: Endcap-corrected monopole residues still sum to 109/280
PASS: Endcap-corrected one-pole Pade error stays below 1e-3 on the low-frequency band
PASS: Separated-order full-wall monopole residues stay positive
PASS: Separated-order full-wall monopole residues still sum to 109/280
PASS: Separated-order full-wall one-pole Pade error stays below 1e-3 on the low-frequency band
Key exported symbol: Family1EndcapSoftwallResults.
"*)
