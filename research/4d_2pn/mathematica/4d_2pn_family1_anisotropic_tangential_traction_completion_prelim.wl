(* ::Package:: *)

ClearAll[
  section, pass, fail,
  A, Btan0, Btan2, q, r, mu,
  K00expr, K11expr, K10expr, K20expr, K21expr, K22expr,
  eqns, solFamily1, residuals, sigmaExpr, bTanExpr, sigmaMin, bTanMin,
  beta0, beta2, identity, k00FromSupport, monopoleAdd
];

section[s_String] := Print["\n" <> StringRepeat["=", 78] <> "\n" <> s <> "\n" <> StringRepeat["=", 78]];
pass[s_String] := Print["PASS: " <> s];
fail[s_String] := Print["FAIL: " <> s];

section["1) Minimal anisotropic tangential-traction completion"];

Print["Strict Family-1 boundary-layer ansatz:"];
Print["  isotropic normal penetration moment   A"];
Print["  axisymmetric tangential wall moment   B_tan(mu) = Btan0 + Btan2 mu^2"];
Print["  flared mouth profile                  sigma(mu) = 1 - q mu^2 + r mu^4"];

K00expr =
  7 A q^2/15 - 26 A q r/35 - 2 A q/3 + 23 A r^2/63 + 2 A r/5 + A;

K11expr =
  11 A q^2/35 - 2 A q r/5 - 2 A q/5 + 13 A r^2/77 + 6 A r/35 + A +
  8 Btan0 q^2/35 - 32 Btan0 q r/105 + 32 Btan0 r^2/231 + 2 Btan0 +
  8 Btan2 q^2/105 - 32 Btan2 q r/231 + 32 Btan2 r^2/429 + 4 Btan2/5;

K10expr =
  27 A q^2/35 - 10 A q r/7 - 6 A q/5 + 25 A r^2/33 + 6 A r/7 + A -
  16 Btan0 q^2/35 + 64 Btan0 q r/105 - 64 Btan0 r^2/231 + 2 Btan0 -
  16 Btan2 q^2/105 + 64 Btan2 q r/231 - 64 Btan2 r^2/429 + 2 Btan2/5;

K20expr =
  13 A q^2/21 - 94 A q r/77 - 22 A q/21 + 6185 A r^2/9009 + 6 A r/7 + A -
  16 Btan0 q^2/7 + 320 Btan0 q r/77 - 320 Btan0 r^2/143 + 6 Btan0 -
  80 Btan2 q^2/77 + 320 Btan2 q r/143 - 192 Btan2 r^2/143 + 18 Btan2/7;

K21expr =
  13 A q^2/21 - 230 A q r/231 - 6 A q/7 + 205 A r^2/429 + 10 A r/21 + A +
  8 Btan0 q^2/21 - 96 Btan0 q r/77 + 800 Btan0 r^2/1001 + 6 Btan0 +
  24 Btan2 q^2/77 - 800 Btan2 q r/1001 + 224 Btan2 r^2/429 + 16 Btan2/7;

K22expr =
  5 A q^2/21 - 58 A q r/231 - 2 A q/7 + 25 A r^2/273 + 2 A r/21 + A +
  16 Btan0 q^2/21 - 64 Btan0 q r/77 + 320 Btan0 r^2/1001 + 6 Btan0 +
  16 Btan2 q^2/77 - 320 Btan2 q r/1001 + 64 Btan2 r^2/429 + 10 Btan2/7;

Print["\nExact projected channel formulas:"];
Print["  K_(0,0) = ", K00expr];
Print["  K_(1,1) = ", K11expr];
Print["  K_(1,0) = ", K10expr];
Print["  K_(2,0) = ", K20expr];
Print["  K_(2,1) = ", K21expr];
Print["  K_(2,2) = ", K22expr];

section["2) Exact l=1,2 support solve on the Family-1-like branch"];

eqns = {
  K11expr == 2/7,
  K10expr == 1/4,
  K20expr == 4/9,
  K21expr == 2/3,
  K22expr == 8/3
};

solFamily1 = Quiet@FindRoot[
  eqns,
  {{A, -0.281923219302567},
   {Btan0, 0.648914709228733},
   {Btan2, -1.085434603876673},
   {q, 2.370915717168574},
   {r, 2.758343474052255}},
  WorkingPrecision -> 50,
  AccuracyGoal -> 25,
  PrecisionGoal -> 25,
  MaxIterations -> 200
];

Print["Selected Family-1-like physical branch:"];
Do[Print["  ", sym, " = ", N[sym /. solFamily1, 20]], {sym, {A, Btan0, Btan2, q, r}}];

residuals = Association[
  "K11-2/7" -> N[(K11expr - 2/7) /. solFamily1, 30],
  "K10-1/4" -> N[(K10expr - 1/4) /. solFamily1, 30],
  "K20-4/9" -> N[(K20expr - 4/9) /. solFamily1, 30],
  "K21-2/3" -> N[(K21expr - 2/3) /. solFamily1, 30],
  "K22-8/3" -> N[(K22expr - 8/3) /. solFamily1, 30]
];

KeyValueMap[Print["  ", #1, " = ", #2] &, residuals];

sigmaExpr = 1 - q mu^2 + r mu^4 /. solFamily1;
bTanExpr = Btan0 + Btan2 mu^2 /. solFamily1;

sigmaMin = NMinimize[{sigmaExpr, -1 <= mu <= 1}, mu, WorkingPrecision -> 40];
bTanMin = NMinimize[{bTanExpr, -1 <= mu <= 1}, mu, WorkingPrecision -> 40];

Print["\nsigma(mu) minimum on [-1,1] = ", sigmaMin];
Print["B_tan(mu) minimum on [-1,1] = ", bTanMin];

beta0 = N[(Btan0 + Btan2/3) /. solFamily1, 20];
beta2 = N[(2 Btan2/3) /. solFamily1, 20];

Print["\nLegendre form of the selected tangential wall profile:"];
Print["  B_tan(mu) = beta0 + beta2 P2(mu)"];
Print["  beta0 = ", beta0];
Print["  beta2 = ", beta2];

If[AllTrue[Values[residuals], Abs[#] < 10^-20 &],
  pass["Exact l=1,2 support-sector reconstruction achieved on the Family-1-like branch."],
  fail["Support-sector residuals are not sufficiently small."]
];

section["3) Universal monopole prediction"];

identity = FullSimplify[
  Expand[K00expr - (K11expr + 1/2 K10expr - 1/10 K20expr - 1/5 K21expr - 1/5 K22expr)]
];

Print["Exact structural identity:"];
Print["  K_(0,0) - [ K_(1,1) + 1/2 K_(1,0) - 1/10 K_(2,0) - 1/5 K_(2,1) - 1/5 K_(2,2) ] = ", identity];

k00FromSupport = FullSimplify[2/7 + 1/2 (1/4) - 1/10 (4/9) - 1/5 (2/3) - 1/5 (8/3)];
monopoleAdd = FullSimplify[4/45 - k00FromSupport];

Print["\nSo once the l=1,2 support targets are matched exactly, the raw boundary-layer monopole is fixed to"];
Print["  K_(0,0)^BL = ", k00FromSupport];
Print["and the separate monopole wall add-on needed to recover K_(0,0)=4/45 is"];
Print["  Delta K_(0,0)^mono = ", monopoleAdd];

Print["\nNumerically:"];
Print["  K_(0,0)^BL = ", N[k00FromSupport, 20]];
Print["  Delta K_(0,0)^mono = ", N[monopoleAdd, 20]];

If[identity === 0 && k00FromSupport === -757/2520 && monopoleAdd === 109/280,
  pass["Universal monopole prediction and add-on are exact."],
  fail["Monopole identity or exact rational values did not match."]
];

section["4) Conclusion"];

Print["Main result:"];
Print["- The smallest genuine strict boundary-layer / soft-wall completion that works is to keep the normal penetration moment isotropic while promoting the tangential wall moment to the axisymmetric profile B_tan(mu)=Btan0+Btan2 mu^2."];
Print["- That single extra tangential wall-stress degree of freedom reproduces the entire l=1,2 support sector exactly."];
Print["- The monopole remains separate, but now with an exact universal prediction K_(0,0)^BL = -757/2520 and exact required monopole add-on 109/280."];

(*"
Output:


==============================================================================
1) Minimal anisotropic tangential-traction completion
==============================================================================
Strict Family-1 boundary-layer ansatz:
  isotropic normal penetration moment   A
  axisymmetric tangential wall moment   B_tan(mu) = Btan0 + Btan2 mu^2
  flared mouth profile                  sigma(mu) = 1 - q mu^2 + r mu^4

Exact projected channel formulas:
  K_(0,0) = A - (2*A*q)/3 + (7*A*q^2)/15 + (2*A*r)/5 - (26*A*q*r)/35 + (23*A*r^2)/63
  K_(1,1) = A + 2*Btan0 + (4*Btan2)/5 - (2*A*q)/5 + (11*A*q^2)/35 + (8*Btan0*q^2)/35 + (8*Btan2*q^2)/105 + (6*A*r)/35 - (2*A*q*r)/5 - (32*Btan0*q*r)/105 - (32*Btan2*q*r)/231 + (13*A*r^2)/77 + (32*Btan0*r^2)/231 + (32*Btan2*r^2)/429
  K_(1,0) = A + 2*Btan0 + (2*Btan2)/5 - (6*A*q)/5 + (27*A*q^2)/35 - (16*Btan0*q^2)/35 - (16*Btan2*q^2)/105 + (6*A*r)/7 - (10*A*q*r)/7 + (64*Btan0*q*r)/105 + (64*Btan2*q*r)/231 + (25*A*r^2)/33 - (64*Btan0*r^2)/231 - (64*Btan2*r^2)/429
  K_(2,0) = A + 6*Btan0 + (18*Btan2)/7 - (22*A*q)/21 + (13*A*q^2)/21 - (16*Btan0*q^2)/7 - (80*Btan2*q^2)/77 + (6*A*r)/7 - (94*A*q*r)/77 + (320*Btan0*q*r)/77 + (320*Btan2*q*r)/143 + (6185*A*r^2)/9009 - (320*Btan0*r^2)/143 - (192*Btan2*r^2)/143
  K_(2,1) = A + 6*Btan0 + (16*Btan2)/7 - (6*A*q)/7 + (13*A*q^2)/21 + (8*Btan0*q^2)/21 + (24*Btan2*q^2)/77 + (10*A*r)/21 - (230*A*q*r)/231 - (96*Btan0*q*r)/77 - (800*Btan2*q*r)/1001 + (205*A*r^2)/429 + (800*Btan0*r^2)/1001 + (224*Btan2*r^2)/429
  K_(2,2) = A + 6*Btan0 + (10*Btan2)/7 - (2*A*q)/7 + (5*A*q^2)/21 + (16*Btan0*q^2)/21 + (16*Btan2*q^2)/77 + (2*A*r)/21 - (58*A*q*r)/231 - (64*Btan0*q*r)/77 - (320*Btan2*q*r)/1001 + (25*A*r^2)/273 + (320*Btan0*r^2)/1001 + (64*Btan2*r^2)/429

==============================================================================
2) Exact l=1,2 support solve on the Family-1-like branch
==============================================================================
Selected Family-1-like physical branch:
  A = -0.28192321930256801306217652511668432095`20.
  Btan0 = 0.64891470922873163178469476960102677809`20.
  Btan2 = -1.08543460387667108583027478499355652474`20.
  q = 2.37091571716857638961663327439521543074`20.
  r = 2.75834347405224663495172884331598732546`20.
  K11-2/7 = 0``48.644578961902404
  K10-1/4 = 0``48.295444607504045
  K20-4/9 = 0``47.62926904811592
  K21-2/3 = 0``48.05745389626773
  K22-8/3 = 0``48.32488137360163

sigma(mu) minimum on [-1,1] = {0.4905238061543064337863356228519762886551107971869820347949`40., {mu -> -0.6555697230448126001326550714083302247219097175702208906147`40.}}
B_tan(mu) minimum on [-1,1] = {-0.4365198946479394540455800153925297466515206722420924556133`40., {mu -> 1.`40.}}

Legendre form of the selected tangential wall profile:
  B_tan(mu) = beta0 + beta2 P2(mu)
  beta0 = 0.28710317460317460317460317460317460317`20.
  beta2 = -0.72362306925111405722018318999570434983`20.
PASS: Exact l=1,2 support-sector reconstruction achieved on the Family-1-like branch.

==============================================================================
3) Universal monopole prediction
==============================================================================
Exact structural identity:
  K_(0,0) - [ K_(1,1) + 1/2 K_(1,0) - 1/10 K_(2,0) - 1/5 K_(2,1) - 1/5 K_(2,2) ] = 0

So once the l=1,2 support targets are matched exactly, the raw boundary-layer monopole is fixed to
  K_(0,0)^BL = -757/2520
and the separate monopole wall add-on needed to recover K_(0,0)=4/45 is
  Delta K_(0,0)^mono = 109/280

Numerically:
  K_(0,0)^BL = -0.30039682539682539682539682539682539683`20.
  Delta K_(0,0)^mono = 0.38928571428571428571428571428571428572`20.
PASS: Universal monopole prediction and add-on are exact.

==============================================================================
4) Conclusion
==============================================================================
Main result:
- The smallest genuine strict boundary-layer / soft-wall completion that works is to keep the normal penetration moment isotropic while promoting the tangential wall moment to the axisymmetric profile B_tan(mu)=Btan0+Btan2 mu^2.
- That single extra tangential wall-stress degree of freedom reproduces the entire l=1,2 support sector exactly.
- The monopole remains separate, but now with an exact universal prediction K_(0,0)^BL = -757/2520 and exact required monopole add-on 109/280.
"*)
