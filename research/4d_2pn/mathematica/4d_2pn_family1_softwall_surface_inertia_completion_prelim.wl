
ClearAll["Global`*"];

section[title_String] := Module[{},
  Print["\n" <> StringRepeat["=", 78]];
  Print[title];
  Print[StringRepeat["=", 78]];
];

pass[name_String, cond_] := Print[If[TrueQ[cond], "PASS: ", "FAIL: "], name];

(* ---------------------------------------------------------------------- *)
(* 2PN Family-1 soft-wall surface-inertia completion                      *)
(* ---------------------------------------------------------------------- *)

SmoothStep[x_?NumericQ] := (1 + Tanh[x])/2;

XiStar[alpha_?NumericQ, p_?NumericQ] := ArcTanh[2 alpha^(-1/p) - 1];

FProfile[x_?NumericQ, alpha_?NumericQ, p_?NumericQ] := Module[{val},
  val = 1 - alpha SmoothStep[x]^p;
  If[val > 0, N[val^(1/4), 30], 0.]
];

DefectIntegrand[x_?NumericQ, alpha_?NumericQ, p_?NumericQ, k_Integer?NonNegative] :=
  x^k (FProfile[x, alpha, p] - If[x < 0, 1., 0.]);

WallMoment[alpha_?NumericQ, p_?NumericQ, k_Integer?NonNegative] :=
  WallMoment[alpha, p, k] = Module[{xs, f, res1, res2},
    xs = XiStar[alpha, p];
    f[x_?NumericQ] := DefectIntegrand[x, alpha, p, k];
    If[xs >= 0,
      res1 = NIntegrate[f[x], {x, -Infinity, 0},
        WorkingPrecision -> 30, AccuracyGoal -> 10, PrecisionGoal -> 10,
        MaxRecursion -> 12, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
      ];
      res2 = NIntegrate[f[x], {x, 0, xs},
        WorkingPrecision -> 30, AccuracyGoal -> 10, PrecisionGoal -> 10,
        MaxRecursion -> 12, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
      ];
      N[res1 + res2, 20],
      res1 = NIntegrate[f[x], {x, -Infinity, xs},
        WorkingPrecision -> 30, AccuracyGoal -> 10, PrecisionGoal -> 10,
        MaxRecursion -> 12, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
      ];
      res2 = NIntegrate[f[x], {x, xs, 0},
        WorkingPrecision -> 30, AccuracyGoal -> 10, PrecisionGoal -> 10,
        MaxRecursion -> 12, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
      ];
      N[res1 + res2, 20]
    ]
];

GaussLegendreRule[n_Integer?Positive] := GaussLegendreRule[n] = Module[
  {x, roots, weights, dPn},
  roots = x /. NSolve[LegendreP[n, x] == 0, x, Reals, WorkingPrecision -> 40];
  roots = Sort[N[roots, 30]];
  dPn[z_?NumericQ] := N[D[LegendreP[n, x], x] /. x -> z, 30];
  weights = N[2/((1 - #^2) dPn[#]^2) & /@ roots, 30];
  {roots, weights}
];

SharpC0N5[] := N[Sqrt[Pi] Gamma[5/4]/(2 Gamma[7/4]), 30];

AveragedWallMoment[alpha0_?NumericQ, p_?NumericQ, k_Integer?NonNegative, nGL_Integer:20] :=
  AveragedWallMoment[alpha0, p, k, nGL] = Module[
    {nodes, weights, c0, vals, wgt},
    {nodes, weights} = GaussLegendreRule[nGL];
    wgt = (1 - nodes^2)^(1/4);
    c0 = SharpC0N5[];
    vals = WallMoment[alpha0/(1 - #^2), p, k] & /@ nodes;
    N[(1/2 Total[weights * wgt * vals])/c0, 20]
];

(* Universal thin-wall series *)
eps = Symbol["eps"];
mb0 = Symbol["mb0"];
mb1 = Symbol["mb1"];
mb2 = Symbol["mb2"];

J2Series = 1/3 + eps mb0 + 2 eps^2 mb1 + eps^3 mb2;
J4Series = 1/5 + eps mb0 + 4 eps^2 mb1 + 6 eps^3 mb2;
RMassSeries = Expand[3 J2Series];
MaaSeries = Normal@Series[Expand[J4Series/J2Series], {eps, 0, 3}];

section["1) Universal Family-1 thin-wall formulas"];
Print["J2 = "];
Print[TraditionalForm[J2Series]];
Print[""];
Print["J4 = "];
Print[TraditionalForm[J4Series]];
Print[""];
Print["R_mass = 3 J2 = "];
Print[TraditionalForm[RMassSeries]];
Print[""];
Print["M_aa = J4/J2 through O(eps^3) = "];
Print[TraditionalForm[MaaSeries]];

section["2) Worked geometry point carried forward from the monopole branch"];

x01 = N[2.40482555769577276862163187933, 30];
LamEM = N[Sqrt[2] Pi/x01, 30];

a0 = Symbol["a0"];
Lam = Symbol["Lam"];
Sigma = Symbol["Sigma"];
rhoGeom = Symbol["rhoGeom"];
betaGeom = Symbol["betaGeom"];
a = Symbol["a"];
L = Symbol["L"];

V = (4 Pi/3) a^3 L;
A = 4 Pi a^2 L + (8 Pi/3) a^3;

sigma = Sigma/a0^3;
Pvac = rhoGeom Sigma/a0^4;
kappaB = betaGeom Sigma/a0;

EGeom = Expand[Pvac V + sigma A + kappaB a^2/L];
H = D[EGeom, {{a, L}, 2}];
g = {D[V, a], D[V, L]};

subs0 = {a -> a0, L -> Lam a0};
H0 = Simplify[H /. subs0];
V0 = Simplify[V /. subs0];
g0 = Simplify[g /. subs0];

hBar = Simplify[H0 a0^2/Sigma];
gBar = Simplify[g0/V0 a0];

hNum = N[hBar /. {a0 -> 1, Lam -> LamEM, rhoGeom -> 1/10, betaGeom -> 12}, 30];
gNum = N[gBar /. {a0 -> 1, Lam -> LamEM}, 30];

deltaUnit = N[(gBar . Inverse[hBar] . gBar) /. {a0 -> 1, Lam -> LamEM, rhoGeom -> 1/10, betaGeom -> 12}, 30];
SigmaStar = N[deltaUnit/(109/280), 30];

Print["Lam_EM = ", LamEM];
Print["Sigma*  = ", SigmaStar];
Print[""];
Print["hBar = "];
Print[MatrixForm[hNum]];
Print[""];
Print["gBar = ", gNum];

EigResiduesForMetric[Maa_?NumericQ, Mll_?NumericQ] := Module[
  {M, vals, vecs, ord, vecsN, valsN, residues},
  M = {{Maa, 0.}, {0., Mll}};
  {vals, vecs} = Eigensystem[LinearSolve[M, hNum]];
  ord = Ordering[Re[vals]];
  valsN = N[vals[[ord]], 20];
  vecsN = N[vecs[[ord]], 20];
  vecsN = (#/Sqrt[# . M . #] &) /@ vecsN;
  residues = N[((gNum . #)^2/(SigmaStar #2)) & @@@ Transpose[{vecsN, valsN}], 20];
  {valsN, residues}
];

PadePole[vals_List, residues_List] := N[Total[residues]/Total[residues/vals], 20];

MaxPadeRelativeError[vals_List, residues_List] := Module[
  {sgrid, exact, pade, lamEff},
  lamEff = PadePole[vals, residues];
  sgrid = N[Subdivide[0, 0.1 vals[[1]], 400], 20];
  exact = Table[Total[MapThread[#1/(1 - s/#2) &, {residues, vals}]], {s, sgrid}];
  pade = Table[Total[residues]/(1 - s/lamEff), {s, sgrid}];
  Max[Abs[(pade - exact)/exact]]
];

{baseVals, baseResidues} = EigResiduesForMetric[3/5, 1/14];
baseLamEff = PadePole[baseVals, baseResidues];
baseErr = MaxPadeRelativeError[baseVals, baseResidues];

Print[""];
Print["Baseline sharp-wall poles = ", baseVals];
Print["Baseline residues        = ", baseResidues];
Print["Baseline Pade pole       = ", baseLamEff];
Print["Baseline max rel. err.   = ", baseErr];

PoleShiftCoefficients[] := Module[
  {h = 10.^-7, valsPlus, valsMinus, dVals},
  valsPlus = First@EigResiduesForMetric[3/5 + h, 1/14];
  valsMinus = First@EigResiduesForMetric[3/5 - h, 1/14];
  dVals = (valsPlus - valsMinus)/(2 h);
  N[(dVals/baseVals) (6/5) - 3, 20]
];

coeffs = PoleShiftCoefficients[];
Print[""];
Print["Leading physical pole-shift coefficients c_{-,+} = ", coeffs];

EvaluateCase[p_?NumericQ, alpha0_?NumericQ, epsR_?NumericQ, nGL_Integer:20] := Module[
  {bm0Loc, bm1Loc, bm2Loc, RMass, Maa, vals, residues, lamEff, maxErr, omegaRatios, xstar0},
  bm0Loc = AveragedWallMoment[alpha0, p, 0, nGL];
  bm1Loc = AveragedWallMoment[alpha0, p, 1, nGL];
  bm2Loc = AveragedWallMoment[alpha0, p, 2, nGL];

  RMass = N[1 + 3 epsR bm0Loc + 6 epsR^2 bm1Loc + 3 epsR^3 bm2Loc, 20];
  Maa = N[
    3/5 + (6/5) epsR bm0Loc +
    epsR^2 ((42/5) bm1Loc - (18/5) bm0Loc^2) +
    epsR^3 ((54/5) bm0Loc^3 - (162/5) bm0Loc bm1Loc + (81/5) bm2Loc),
    20
  ];

  {vals, residues} = EigResiduesForMetric[Maa, 1/14];
  lamEff = PadePole[vals, residues];
  maxErr = MaxPadeRelativeError[vals, residues];
  omegaRatios = N[(vals/baseVals)/RMass, 20];
  xstar0 = XiStar[alpha0, p];

  <|
    "p" -> p,
    "alpha0" -> alpha0,
    "epsR" -> epsR,
    "threshold" -> N[2^p, 20],
    "xstar0" -> N[xstar0, 20],
    "mbar0" -> bm0Loc,
    "mbar1" -> bm1Loc,
    "mbar2" -> bm2Loc,
    "RMass" -> RMass,
    "Maa" -> Maa,
    "vals" -> vals,
    "residues" -> residues,
    "omegaRatios" -> omegaRatios,
    "lamEff" -> lamEff,
    "maxErr" -> maxErr
  |>
];

section["3) Representative steep-wall cases"];

cases = {
  EvaluateCase[2., 10., 0.05, 20],
  EvaluateCase[4., 10., 0.05, 20]
};

Do[
  Print["\nCase p=", case["p"], ", alpha0=", case["alpha0"], ", epsR=", case["epsR"]];
  Print["  threshold 2^p = ", case["threshold"]];
  Print["  xi*(0)        = ", case["xstar0"]];
  Print["  mbar0         = ", case["mbar0"]];
  Print["  mbar1         = ", case["mbar1"]];
  Print["  mbar2         = ", case["mbar2"]];
  Print["  RMass         = ", case["RMass"]];
  Print["  Maa           = ", case["Maa"]];
  Print["  poles         = ", case["vals"]];
  Print["  residues      = ", case["residues"]];
  Print["  Omega^2 ratios= ", case["omegaRatios"]];
  Print["  Pade pole     = ", case["lamEff"]];
  Print["  max rel. err. = ", case["maxErr"]];
, {case, cases}];

section["4) Summary"];
Print["1) The Family-1 radial soft wall gives a derived boundary-layer correction to the"];
Print["   monopole inertia metric, not a new phenomenological completion term."];
Print[""];
Print["2) At O(epsR), the same averaged wall moment mbar0 controls both"];
Print["      RMass = 1 + 3 epsR mbar0 + ..."];
Print["   and"];
Print["      Maa   = 3/5 + (6/5) epsR mbar0 + ..."];
Print[""];
Print["3) On the representative steep-wall branch (p=2, alpha0=10, epsR=0.05), the"];
Print["   physical monopole poles move upward while the static 109/280 closure remains"];
Print["   unchanged, and the one-pole Pade reduction stays at ~10^-4 relative error."];
Print[""];
Print["4) The next natural tightening is to add the endcap soft-wall layer, or tie the"];
Print["   same wall profile directly to the earlier tangential traction/support law."];

section["5) Verification"];

target109280 = N[109/280, 30];

pass[
  "Sharp-wall baseline residues sum to 109/280",
  Abs[N[Total[baseResidues] - target109280, 30]] < 10^-12
];

pass[
  "Sharp-wall baseline one-pole Pade error stays below 1e-3",
  baseErr < 10^-3
];

Do[
  pass[
    "Surface-inertia case p=" <> ToString[case["p"]] <> " keeps positive residues",
    Min[case["residues"]] > 0
  ];
  pass[
    "Surface-inertia case p=" <> ToString[case["p"]] <> " preserves the static 109/280 closure",
    Abs[N[Total[case["residues"]] - target109280, 30]] < 10^-12
  ];
  pass[
    "Surface-inertia case p=" <> ToString[case["p"]] <> " keeps one-pole Pade error below 1e-3",
    case["maxErr"] < 10^-3
  ];
  pass[
    "Surface-inertia case p=" <> ToString[case["p"]] <> " shifts both physical poles upward",
    Min[case["omegaRatios"]] > 1
  ];
,
  {case, cases}
];

Family1SurfaceInertiaResults = <|
  "RMassSeries" -> RMassSeries,
  "MaaSeries" -> MaaSeries,
  "BaseValues" -> baseVals,
  "BaseResidues" -> baseResidues,
  "BaseLamEff" -> baseLamEff,
  "BaseError" -> baseErr,
  "PoleShiftCoefficients" -> coeffs,
  "RepresentativeCases" -> cases
|>;

Print["Key exported symbol: Family1SurfaceInertiaResults."];

(*"
Output:


==============================================================================
1) Universal Family-1 thin-wall formulas
==============================================================================
J2 = 
TraditionalForm[1/3 + eps*mb0 + 2*eps^2*mb1 + eps^3*mb2]

J4 = 
TraditionalForm[1/5 + eps*mb0 + 4*eps^2*mb1 + 6*eps^3*mb2]

R_mass = 3 J2 = 
TraditionalForm[1 + 3*eps*mb0 + 6*eps^2*mb1 + 3*eps^3*mb2]

M_aa = J4/J2 through O(eps^3) = 
TraditionalForm[3/5 + (6*eps*mb0)/5 - (6*eps^2*(3*mb0^2 - 7*mb1))/5 + (27*eps^3*(2*mb0^3 - 6*mb0*mb1 + 3*mb2))/5]

==============================================================================
2) Worked geometry point carried forward from the monopole branch
==============================================================================
Lam_EM = 1.8474865771201280510433744839396122428696573796896114673325`29.3810835788064
Sigma*  = 0.20761432918354888540403250898444970518`26.476476754460833

hBar = 
MatrixForm[{{114.33174685298724550956098370044657032947`29.632621016158737, 19.35786732628000931888935542023829255417`29.51986243386294}, {19.35786732628000931888935542023829255417`29.51986243386294, 3.80598757845104928960045813753644681536`28.903962324086738}}]

gBar = {3.`30., 0.5412759217762790697885051933117400179826769616175860802256`29.3810835788064}

Baseline sharp-wall poles = {5.925562576926863, 237.91117494303325}
Baseline residues        = {0.0026280028657021397, 0.3866577114200121}
Baseline Pade pole       = 188.17695898017192
Baseline max rel. err.   = 0.00007100969970122967

Leading physical pole-shift coefficients c_{-,+} = {-3.408286213590945, -4.5917137855849415}

==============================================================================
3) Representative steep-wall cases
==============================================================================

Case p=2., alpha0=10., epsR=0.05
  threshold 2^p = 4.
  xi*(0)        = -0.3855810692154256
  mbar0         = -0.65088873122065008379221256309790053717`19.9996325479206
  mbar1         = 0.25050113075080546135567455398106278544`19.999473638705286
  mbar2         = -0.15564575968474289279578537218256416667`19.999327217351524
  RMass         = 0.9060658401182828
  Maa           = 0.5623671912556893
  poles         = {6.002317874599835, 250.58594815325637}
  residues      = {0.0028985744267034873, 0.3863871398590109}
  Omega^2 ratios= {1.117968701458279, 1.1624709663217863}
  Pade pole     = 192.25469161621535
  max rel. err. = 0.00007841090874976033

Case p=4., alpha0=10., epsR=0.05
  threshold 2^p = 16.
  xi*(0)        = 0.12533484020558097
  mbar0         = -0.06871900887756139638433062476173914115`19.8715058298112
  mbar1         = 0.02098204877157089886189677040803353417`19.999203774519213
  mbar2         = -0.00650709388379880272996267095677777998`19.999099784835735
  RMass         = 0.9900044392397329
  Maa           = 0.5962672063713681
  poles         = {5.933096894056334, 239.0965532765938}
  residues      = {0.002653972006883264, 0.386631742278831}
  Omega^2 ratios= {1.0113808123389363, 1.0151292266168788}
  Pade pole     = 188.57368458791117
  max rel. err. = 0.00007171939497249488

==============================================================================
4) Summary
==============================================================================
1) The Family-1 radial soft wall gives a derived boundary-layer correction to the
   monopole inertia metric, not a new phenomenological completion term.

2) At O(epsR), the same averaged wall moment mbar0 controls both
      RMass = 1 + 3 epsR mbar0 + ...
   and
      Maa   = 3/5 + (6/5) epsR mbar0 + ...

3) On the representative steep-wall branch (p=2, alpha0=10, epsR=0.05), the
   physical monopole poles move upward while the static 109/280 closure remains
   unchanged, and the one-pole Pade reduction stays at ~10^-4 relative error.

4) The next natural tightening is to add the endcap soft-wall layer, or tie the
   same wall profile directly to the earlier tangential traction/support law.

==============================================================================
5) Verification
==============================================================================
PASS: Sharp-wall baseline residues sum to 109/280
PASS: Sharp-wall baseline one-pole Pade error stays below 1e-3
PASS: Surface-inertia case p=2. keeps positive residues
PASS: Surface-inertia case p=2. preserves the static 109/280 closure
PASS: Surface-inertia case p=2. keeps one-pole Pade error below 1e-3
PASS: Surface-inertia case p=2. shifts both physical poles upward
PASS: Surface-inertia case p=4. keeps positive residues
PASS: Surface-inertia case p=4. preserves the static 109/280 closure
PASS: Surface-inertia case p=4. keeps one-pole Pade error below 1e-3
PASS: Surface-inertia case p=4. shifts both physical poles upward
Key exported symbol: Family1SurfaceInertiaResults.
"*)
