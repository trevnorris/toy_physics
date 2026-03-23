(***
  4D -> candidate 2PN master harness
  ----------------------------------
  Purpose:
    A single PASS/FAIL style Wolfram Language script that carries forward the
    established 1PN derivation logic and then pushes the candidate 2PN sectors
    we discussed into one place.

  What this script does:
    1) Re-runs the key 1PN regressions that must remain frozen:
         kappa_rho = 1,
         n = 5,
         C_self^(1PN) = -3/2 in the one-body scalar/optical ledger,
         static 1PN pair-count coefficient = -1/2,
         kappa_add = 1/2,
         kappa_PV = 3/2,
         beta_1PN = 3.
    2) Exactifies the n = 5 Bernoulli / barotropic closure.
    3) Builds the minimal exact 2PN self-sector worldline and extracts its
       2PN monomials.
    4) Shows why a fully Bernoulli-exactified mass prefactor fails the 1PN
       static regression.
    5) Extends local mass scaling to quadratic order and extracts the 2PN
       static hierarchy (2-body, selected 3-body, selected 4-body terms).
    6) Derives the generic cubic adiabatic-elimination response formula for a
       single internal fast variable.
    7) Assembles a self+static-only 2PN comparable-mass candidate and reduces
       it to the test-mass limit.
    8) Compares that test-mass candidate to an exact Schwarzschild-isotropic
       coordinate target through 2PN.
    9) States the worldtube finite-size gate that any claimed universal 1/r^4
       term must pass.

  Scope note:
    This script does NOT claim to close the full comparable-mass conservative
    2PN problem. The 2PN wake / cross sector remains open here. What this file
    does is make the self/static/response/test-mass bookkeeping explicit and
    machine-checkable.
***)

ClearAll["Global`*"];
Print["--- 4D gravity + candidate 2PN master harness ---"];

(* ---------------------------------------------------------------------- *)
(* Helpers                                                                *)
(* ---------------------------------------------------------------------- *)

passCount = 0;
failCount = 0;
skipCount = 0;

section[name_String] := Print["\n=== ", name, " ==="];
pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);
skip[name_String, msg_String : ""] := (
  skipCount++;
  If[msg === "",
    Print["SKIP: ", name],
    Print["SKIP: ", name, "\n  ", msg]
  ]
);
info[msg_String] := Print["INFO: ", msg];

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet @ FullSimplify[res == 0, assum]]
];

zeroArrayQ[arr_, assum_: True] := Module[{res, flat},
  res = Quiet @ FullSimplify[arr, assum];
  flat = Flatten[{res}];
  And @@ (TrueQ[# === 0] || TrueQ[Quiet @ FullSimplify[# == 0, assum]] & /@ flat)
];

checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  If[zeroScalarQ[expr, assum], pass[name], fail[name, res]]
];

checkVector[name_String, expr_List, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  If[zeroArrayQ[expr, assum], pass[name], fail[name, res]]
];

checkEqScalar[name_String, lhs_, rhs_, assum_: True] := checkScalar[name, lhs - rhs, assum];
checkEqVector[name_String, lhs_List, rhs_List, assum_: True] := checkVector[name, lhs - rhs, assum];
checkCondition[name_String, cond_] := If[TrueQ[cond], pass[name], fail[name, cond]];

PNSeries[expr_, vars_List, order_Integer] := Module[{lam, rules},
  rules = (# -> lam #) & /@ vars;
  Expand[Normal[Series[Expand[expr /. rules], {lam, 0, order}]] /. lam -> 1]
];

dot3[a_List, b_List] := Total[a b];

positiveMassAssumptions = cLight > 0 && Gconst > 0 && rAB > 0 && rAC > 0 && rBC > 0 && rAD > 0 && rBD > 0 && rCD > 0 &&
  mA > 0 && mB > 0 && mC > 0 && mD > 0 && Mcentral > 0 && rho0 > 0 && aBody > 0 && rSep > 0;

(* ---------------------------------------------------------------------- *)
section["CARRY-FORWARD / frozen 1PN regressions"];
(* ---------------------------------------------------------------------- *)

alphaN[n_] := (n - 1)/2;
nSolution = nEOS /. First[Solve[alphaN[nEOS] == 2, nEOS]];

checkEqScalar[
  "Optical weak-field matching still fixes n = 5",
  nSolution,
  5
];

LscRaw1PN = -M0 (1 + kappaRhoVar epsPhi) cLight^2 Sqrt[1 - epsV2/(1 + (nEOS - 1) epsPhi)];
LscScaled1PN = Expand[LscRaw1PN/(M0 cLight^2)];
LscPNScaled1PN = Expand[PNSeries[LscScaled1PN, {epsPhi, epsV2}, 2]];

coeffNewtonPhi = FullSimplify[Coefficient[LscPNScaled1PN /. epsV2 -> 0, epsPhi], cLight > 0];
coeffV2 = FullSimplify[Coefficient[LscPNScaled1PN /. epsPhi -> 0, epsV2], cLight > 0];
coeffV4 = FullSimplify[Coefficient[LscPNScaled1PN /. epsPhi -> 0, epsV2, 2], cLight > 0];
coeffPhiV2 = FullSimplify[Coefficient[Coefficient[LscPNScaled1PN, epsPhi], epsV2], cLight > 0];

checkEqScalar[
  "Scalar-sector Newtonian coupling coefficient remains -kappa_rho",
  coeffNewtonPhi,
  -kappaRhoVar,
  cLight > 0
];

kappaRhoSolution = kappaRhoVar /. First[Solve[coeffNewtonPhi == -1, kappaRhoVar]];

checkEqScalar[
  "Newtonian matching still fixes kappa_rho = 1 (historically q = 1)",
  kappaRhoSolution,
  1,
  cLight > 0
];

checkEqScalar[
  "Free kinetic coefficient remains 1/2",
  coeffV2,
  1/2,
  cLight > 0
];

checkEqScalar[
  "Free-particle quartic kinetic coefficient remains 1/8 in the Lagrangian",
  coeffV4,
  1/8,
  cLight > 0
];

checkEqScalar[
  "Scalar+optical 1PN self coefficient remains kappa_rho/2 - (n-1)/2",
  coeffPhiV2,
  kappaRhoVar/2 - (nEOS - 1)/2,
  cLight > 0
];

checkEqScalar[
  "With kappa_rho = 1 and n = 5 the one-body Phi v^2 coefficient is still -3/2",
  coeffPhiV2 /. {kappaRhoVar -> kappaRhoSolution, nEOS -> nSolution},
  -3/2,
  cLight > 0
];

PhiLocA = -(Gconst mB)/rAB;
PhiLocB = -(Gconst mA)/rAB;
LpairPotential1PN = Expand[-(1/2) (mA (1 + kappaRhoSolution PhiLocA/cLight^2) PhiLocA + mB (1 + kappaRhoSolution PhiLocB/cLight^2) PhiLocB)];
invC = Unique["invC"];
LpairPotential1PNInv = Expand[LpairPotential1PN /. cLight -> 1/invC];

staticCoeff1PN = FullSimplify[
  Coefficient[LpairPotential1PNInv, invC, 2] / (Gconst^2 mA mB (mA + mB)/rAB^2),
  Gconst > 0 && rAB > 0 && mA > 0 && mB > 0
];

checkEqScalar[
  "Pair-counted local mass scaling still gives the 1PN static coefficient -1/2",
  staticCoeff1PN,
  -1/2,
  Gconst > 0 && rAB > 0 && mA > 0 && mB > 0
];

kappaAdd[d_] := 1/(d - 1);
checkEqScalar["Added-mass closure still gives kappa_add = 1/2 in 3D", kappaAdd[3], 1/2];

EwExpr = rho0^((nEOS - 1)/2)/aRad;
EfExpr = xRatio EwExpr;
EpvExpr = EwExpr (1 + 2 xRatio)/3;
FtotExpr = EwExpr + EfExpr + EpvExpr;

dlnFgeneral = FullSimplify[
  (((nEOS - 1)/2) EwExpr - EfExpr + nEOS EpvExpr)/FtotExpr,
  rho0 > 0 && aRad > 0
];

dlnFn5 = FullSimplify[dlnFgeneral /. nEOS -> 5, rho0 > 0 && aRad > 0];

checkEqScalar[
  "At n = 5 the adiabatic response slope remains (11+7x)/(4+5x)",
  dlnFn5,
  (11 + 7 xRatio)/(4 + 5 xRatio),
  rho0 > 0 && aRad > 0
];

xSolution = xRatio /. First[Solve[dlnFn5 == 5/2, xRatio]];
kappaPV = FullSimplify[(dlnFn5 /. xRatio -> xSolution) - 1, rho0 > 0 && aRad > 0];
beta1PN = FullSimplify[kappaRhoSolution + kappaAdd[3] + kappaPV, rho0 > 0 && aRad > 0];

checkEqScalar[
  "Adiabatic internal response closure still gives x = 2/11",
  xSolution,
  2/11,
  rho0 > 0 && aRad > 0
];

checkEqScalar[
  "Adiabatic internal response closure still gives kappa_PV = 3/2",
  kappaPV,
  3/2,
  rho0 > 0 && aRad > 0
];

checkEqScalar[
  "The 1PN ledger still closes to beta = 3",
  beta1PN,
  3,
  rho0 > 0 && aRad > 0
];

(* ---------------------------------------------------------------------- *)
section["EXACT / n = 5 Bernoulli continuation of the barotropic sector"];
(* ---------------------------------------------------------------------- *)

PEOS[rho_] := KEOS rho^5;
UEOS[rho_] := (KEOS/4) rho^5;
hEOS[rho_] := D[UEOS[rhoTmp], rhoTmp] /. rhoTmp -> rho;
cs2EOS[rho_] := ((1/mMatter) D[PEOS[rhoTmp], rhoTmp]) /. rhoTmp -> rho;

matterUnitRule = mMatter -> 1;
KEOSRule = First[Solve[(cs2EOS[rho0] /. matterUnitRule) == cLight^2, KEOS]];
h0 = FullSimplify[(hEOS[rho0] /. KEOSRule) /. matterUnitRule, cLight > 0 && rho0 > 0];

checkEqScalar[
  "Background sound-speed normalization gives h0 = c^2/4 in reduced matter units",
  h0,
  cLight^2/4,
  cLight > 0 && rho0 > 0
];

rhoBernoulli = rho0 (1 + 4 PhiB/cLight^2)^(1/4);
bernoulliDomain = cLight > 0 && rho0 > 0 && 1 + 4 PhiB/cLight^2 > 0;

checkEqScalar[
  "Exact n=5 Bernoulli density profile solves h(rho) = h0 + Phi",
  FullSimplify[((hEOS[rhoBernoulli] /. KEOSRule) /. matterUnitRule) - (h0 + PhiB), bernoulliDomain],
  0,
  bernoulliDomain
];

cs2Bernoulli = FullSimplify[((cs2EOS[rho] /. KEOSRule) /. matterUnitRule /. rho -> rhoBernoulli), bernoulliDomain];
NoptBernoulli = FullSimplify[cLight/Sqrt[cs2Bernoulli], bernoulliDomain];

checkEqScalar[
  "Exact Bernoulli closure gives cs^2/c^2 = 1 + 4 Phi/c^2",
  cs2Bernoulli/cLight^2,
  1 + 4 PhiB/cLight^2,
  bernoulliDomain
];

checkEqScalar[
  "Exact Bernoulli closure gives N(Phi) = (1 + 4 Phi/c^2)^(-1/2)",
  NoptBernoulli,
  (1 + 4 PhiB/cLight^2)^(-1/2),
  bernoulliDomain
];

NoptU = FullSimplify[NoptBernoulli /. PhiB -> -Upot, cLight > 0 && Upot >= 0 && 1 - 4 Upot/cLight^2 > 0];
NoptUSeries2PN = Expand[PNSeries[NoptU, {Upot}, 2]];

checkEqScalar[
  "Bernoulli optical index expands as 1 + 2 U/c^2 + 6 U^2/c^4 through 2PN",
  NoptUSeries2PN,
  1 + 2 Upot/cLight^2 + 6 Upot^2/cLight^4,
  cLight > 0 && Upot >= 0
];

(* ---------------------------------------------------------------------- *)
section["CONTROLLED / minimal Bernoulli exactification of the self sector through 2PN"];
(* ---------------------------------------------------------------------- *)

LselfExactMinimal = -M0 (1 + PhiB/cLight^2) cLight^2 Sqrt[1 - v2/(cLight^2 (1 + 4 PhiB/cLight^2))];
LselfExactMinimalScaled = Expand[(LselfExactMinimal/(M0 cLight^2)) /. {PhiB -> cLight^2 epsPhi, v2 -> cLight^2 epsV2}];
Lself2PNScaled = Expand[PNSeries[LselfExactMinimalScaled, {epsPhi, epsV2}, 3]];
Lself2PNScaledU = Expand[Lself2PNScaled /. epsPhi -> -epsU];

checkEqScalar[
  "Minimal exactification still reproduces the 1PN self coefficient -3/2",
  Coefficient[Coefficient[Lself2PNScaled, epsPhi], epsV2],
  -3/2
];

checkEqScalar[
  "Minimal exactification gives the universal free 2PN kinetic coefficient 1/16",
  Coefficient[Lself2PNScaledU, epsV2, 3],
  1/16
];

checkEqScalar[
  "Minimal exactification gives the candidate U v^4 coefficient 7/8",
  Coefficient[Coefficient[Lself2PNScaledU, epsU], epsV2, 2],
  7/8
];

checkEqScalar[
  "Minimal exactification gives the candidate U^2 v^2 coefficient 6",
  Coefficient[Coefficient[Lself2PNScaledU, epsU, 2], epsV2],
  6
];

Lself2PNPhysical = Expand[M0 cLight^2 (Lself2PNScaledU /. {epsU -> Upot/cLight^2, epsV2 -> v2/cLight^2})];
info["Lself2PNScaledU stores the dimensionless one-body candidate through 2PN in U = -Phi variables."];

(* ---------------------------------------------------------------------- *)
section["CAUTION / fully Bernoulli-exactified mass dressing fails a frozen 1PN regression"];
(* ---------------------------------------------------------------------- *)

LselfExactAggressiveScaled = -(1 + 4 epsPhi)^(1/4) Sqrt[1 - epsV2/(1 + 4 epsPhi)];
LselfAggressive1PN = Expand[PNSeries[LselfExactAggressiveScaled, {epsPhi, epsV2}, 2]];
aggressiveStatic1PN = Coefficient[LselfAggressive1PN /. epsV2 -> 0, epsPhi, 2];

checkEqScalar[
  "Fully Bernoulli-exactified mass dressing would generate +3/2 epsPhi^2 at 1PN",
  aggressiveStatic1PN,
  3/2
];

info["That +3/2 static coefficient is incompatible with the already-frozen pair-counted 1PN static ledger, so only the optical denominator is exactified here."];

(* ---------------------------------------------------------------------- *)
section["EXACT / quadratic local mass scaling and the 2PN static hierarchy"];
(* ---------------------------------------------------------------------- *)

mEffA = mA (1 + kappaRhoSolution PhiLocA/cLight^2 + lambdaRho PhiLocA^2/cLight^4);
mEffB = mB (1 + kappaRhoSolution PhiLocB/cLight^2 + lambdaRho PhiLocB^2/cLight^4);
LpairPotential2PN = Expand[-(1/2) (mEffA PhiLocA + mEffB PhiLocB)];
LpairPotential2PNInv = Expand[LpairPotential2PN /. cLight -> 1/invC];

staticCoeff2PNTwoBody = FullSimplify[
  Coefficient[LpairPotential2PNInv, invC, 4] / (Gconst^3 mA mB (mA^2 + mB^2)/rAB^3),
  Gconst > 0 && rAB > 0 && mA > 0 && mB > 0
];

checkEqScalar[
  "Quadratic local mass scaling gives lambda_rho/2 on the 2PN two-body static term",
  staticCoeff2PNTwoBody,
  lambdaRho/2,
  Gconst > 0 && rAB > 0 && mA > 0 && mB > 0
];

PhiA3 = -Gconst (mB/rAB + mC/rAC);
PhiB3 = -Gconst (mA/rAB + mC/rBC);
PhiC3 = -Gconst (mA/rAC + mB/rBC);
Lstatic3Body2PN = Expand[-(lambdaRho/(2 cLight^4)) (mA PhiA3^3 + mB PhiB3^3 + mC PhiC3^3)];
Lstatic3Body2PNSub = Expand[Lstatic3Body2PN /. {rAB -> 1/uAB, rAC -> 1/uAC, rBC -> 1/uBC}];

coeff3AB2AC = FullSimplify[
  Coefficient[Lstatic3Body2PNSub, uAB^2 uAC] / ((3 Gconst^3 lambdaRho mA mB^2 mC)/(2 cLight^4)),
  positiveMassAssumptions && lambdaRho != 0
];
coeff3ABAC2 = FullSimplify[
  Coefficient[Lstatic3Body2PNSub, uAB uAC^2] / ((3 Gconst^3 lambdaRho mA mB mC^2)/(2 cLight^4)),
  positiveMassAssumptions && lambdaRho != 0
];
coeff3ABCmixed = FullSimplify[Coefficient[Lstatic3Body2PNSub, uAB uAC uBC], positiveMassAssumptions];

checkEqVector[
  "Selected 3-body cubic-static monomials carry the expected unit-normalized coefficients",
  {coeff3AB2AC, coeff3ABAC2},
  {1, 1},
  positiveMassAssumptions
];

checkEqScalar[
  "Quadratic local mass scaling generates no pure uAB uAC uBC triple-product term at cubic order",
  coeff3ABCmixed,
  0,
  positiveMassAssumptions
];

PhiA4 = -Gconst (mB/rAB + mC/rAC + mD/rAD);
PhiB4 = -Gconst (mA/rAB + mC/rBC + mD/rBD);
PhiC4 = -Gconst (mA/rAC + mB/rBC + mD/rCD);
PhiD4 = -Gconst (mA/rAD + mB/rBD + mC/rCD);
Lstatic4Body2PN = Expand[-(lambdaRho/(2 cLight^4)) (mA PhiA4^3 + mB PhiB4^3 + mC PhiC4^3 + mD PhiD4^3)];
Lstatic4Body2PNSub = Expand[Lstatic4Body2PN /. {rAB -> 1/uAB, rAC -> 1/uAC, rAD -> 1/uAD, rBC -> 1/uBC, rBD -> 1/uBD, rCD -> 1/uCD}];

coeff4ABACAD = FullSimplify[
  Coefficient[Lstatic4Body2PNSub, uAB uAC uAD] / ((3 Gconst^3 lambdaRho mA mB mC mD)/cLight^4),
  positiveMassAssumptions && lambdaRho != 0
];

checkEqScalar[
  "The first genuine 4-body cubic-static monomial has the expected unit-normalized coefficient",
  coeff4ABACAD,
  1,
  positiveMassAssumptions
];

(* ---------------------------------------------------------------------- *)
section["CONTROLLED / generic cubic adiabatic-elimination closure for one fast internal mode"];
(* ---------------------------------------------------------------------- *)

Hgen = H0 + H1 epsResp + (H2/2) epsResp^2 + (H3/6) epsResp^3
  + (Kstiff/2) daResp^2 + fMix daResp epsResp
  + (gMix/2) daResp^2 epsResp + (hMix/2) daResp epsResp^2
  + (Ccub/6) daResp^3;

daAnsatz = Aresp epsResp + Bresp epsResp^2;
stationarySeries = Expand[Normal[Series[(D[Hgen, daResp] /. daResp -> daAnsatz), {epsResp, 0, 2}]]];
stationaryEq1 = Coefficient[stationarySeries, epsResp, 1];
stationaryEq2 = Coefficient[stationarySeries, epsResp, 2];

AResponse = FullSimplify[-fMix/Kstiff, Kstiff != 0];
BResponse = FullSimplify[-(gMix AResponse + hMix/2 + Ccub AResponse^2/2)/Kstiff, Kstiff != 0];
ABSolution = {Aresp -> AResponse, Bresp -> BResponse};

checkEqScalar[
  "First-order fast-mode response coefficient is A = -f/K",
  AResponse,
  -fMix/Kstiff,
  Kstiff != 0
];

checkEqScalar[
  "Second-order fast-mode response coefficient is B = -(g A + h/2 + C A^2/2)/K",
  BResponse,
  -((gMix (-fMix/Kstiff)) + hMix/2 + Ccub (-fMix/Kstiff)^2/2)/Kstiff,
  Kstiff != 0
];

checkEqScalar[
  "Stationary condition vanishes through O(eps^2) after substituting A and B",
  stationarySeries /. ABSolution,
  0,
  Kstiff != 0
];

HeffSeriesExplicit = Expand @ PNSeries[
  H0 + H1 epsResp + (H2/2) epsResp^2 + (H3/6) epsResp^3
    + (Kstiff/2) (AResponse epsResp + BResponse epsResp^2)^2
    + fMix (AResponse epsResp + BResponse epsResp^2) epsResp
    + (gMix/2) (AResponse epsResp + BResponse epsResp^2)^2 epsResp
    + (hMix/2) (AResponse epsResp + BResponse epsResp^2) epsResp^2
    + (Ccub/6) (AResponse epsResp + BResponse epsResp^2)^3,
  {epsResp},
  3
];
HeffExpectedSeries = Expand[
  H0 + H1 epsResp
  + (1/2) (H2 - fMix^2/Kstiff) epsResp^2
  + (1/6) (H3 - (3 hMix fMix)/Kstiff + (3 gMix fMix^2)/Kstiff^2 - (Ccub fMix^3)/Kstiff^3) epsResp^3
];
HeffQuadratic = H2 - fMix^2/Kstiff;
HeffCubic = H3 - (3 hMix fMix)/Kstiff + (3 gMix fMix^2)/Kstiff^2 - (Ccub fMix^3)/Kstiff^3;
HeffQuadraticDirect = FullSimplify[
  Expand[H2 + Kstiff AResponse^2 + 2 fMix AResponse],
  Kstiff != 0
];
HeffCubicDirect = FullSimplify[
  Expand[H3 + 6 Kstiff AResponse BResponse + 6 fMix BResponse + 3 gMix AResponse^2 + 3 hMix AResponse + Ccub AResponse^3],
  Kstiff != 0
];

checkScalar[
  "Direct adiabatic-elimination series matches the expected through O(eps^3)",
  Together[Expand[HeffSeriesExplicit - HeffExpectedSeries]],
  Kstiff != 0
];

checkEqScalar[
  "Adiabatic elimination shifts the quadratic coefficient to H2 - f^2/K",
  HeffQuadraticDirect,
  HeffQuadratic,
  Kstiff != 0
];

checkEqScalar[
  "Adiabatic elimination shifts the cubic coefficient to H3 - 3 h f/K + 3 g f^2/K^2 - C f^3/K^3",
  HeffCubicDirect,
  HeffCubic,
  Kstiff != 0
];

info["HeffCubic is the natural one-DOF symbolic template for a second-order response contribution to the 2PN ledger."];

(* ---------------------------------------------------------------------- *)
section["FULL / self+static-only comparable-mass 2PN candidate and test-mass reduction"];
(* ---------------------------------------------------------------------- *)

Clear[vA2, vB2, vAB, vAn, vBn];

L1PNTwoBodyFrozen = (
  (1/2) mA vA2 + (1/2) mB vB2
  + (mA vA2^2 + mB vB2^2)/(8 cLight^2)
  + Gconst mA mB/rAB
  + (Gconst mA mB)/(cLight^2 rAB)
      ((3/2) (vA2 + vB2) - (7/2) vAB - (1/2) vAn vBn)
  - (Gconst^2 mA mB (mA + mB))/(2 cLight^2 rAB^2)
);

L2PNSelfStaticCandidate = (
  (mA vA2^3 + mB vB2^3)/(16 cLight^4)
  + (7 Gconst mA mB)/(8 cLight^4 rAB) (vA2^2 + vB2^2)
  + (6 Gconst^2 mA mB)/(cLight^4 rAB^2) (mB vA2 + mA vB2)
  + (lambdaRho Gconst^3 mA mB (mA^2 + mB^2))/(2 cLight^4 rAB^3)
);

L2BodyCandidateThrough2PN = Expand[L1PNTwoBodyFrozen + L2PNSelfStaticCandidate];
info["L2BodyCandidateThrough2PN is the explicit self+static-only comparable-mass candidate through 2PN; no 2PN wake/cross sector has been added."];

L2BodyCandidateScaledForTestMass = Expand[(L2BodyCandidateThrough2PN/cLight^2) /. {
    mB -> Mcentral,
    rAB^(-1) -> (cLight^2 epsU)/(Gconst Mcentral),
    rAB^(-2) -> (cLight^4 epsU^2)/(Gconst^2 Mcentral^2),
    rAB^(-3) -> (cLight^6 epsU^3)/(Gconst^3 Mcentral^3),
    vA2 -> cLight^2 epsV2,
    vB2 -> 0,
    vAB -> 0,
    vAn -> 0,
    vBn -> 0
  }];

LtestMassCandidate2PNScaledFromTwoBody = Expand[
  FullSimplify[
    Coefficient[L2BodyCandidateScaledForTestMass, mA, 1],
    positiveMassAssumptions && epsU >= 0 && epsV2 >= 0
  ]
];

LstaticTestMass2PNScaled = Expand[-epsU^2/2 + (lambdaRho/2) epsU^3];
LtestMassCandidate2PNScaled = Expand[Lself2PNScaledU + LstaticTestMass2PNScaled];
LtestMassCandidate2PNScaledDyn = Expand[LtestMassCandidate2PNScaled + 1];

checkEqScalar[
  "Two-body reduction agrees with the direct self+static test-mass construction through 2PN",
  LtestMassCandidate2PNScaledFromTwoBody,
  LtestMassCandidate2PNScaledDyn
];

checkEqScalar[
  "Self+static test-mass candidate has the expected U^3 coefficient lambda_rho/2",
  SeriesCoefficient[LtestMassCandidate2PNScaled /. epsV2 -> 0, {epsU, 0, 3}],
  lambdaRho/2
];

checkEqScalar[
  "Self+static test-mass candidate has the expected U^2 v^2 coefficient 6",
  SeriesCoefficient[SeriesCoefficient[LtestMassCandidate2PNScaled, {epsU, 0, 2}], {epsV2, 0, 1}],
  6
];

(* ---------------------------------------------------------------------- *)
section["TARGET / exact Schwarzschild-isotropic test-mass target through 2PN"];
(* ---------------------------------------------------------------------- *)

uIso = Ucentral/cLight^2;
LisoExactPerUnit = -cLight^2 Sqrt[((1 - uIso/2)/(1 + uIso/2))^2 - (1 + uIso/2)^4 (v2/cLight^2)];
LisoScaled = Expand[(LisoExactPerUnit/cLight^2) /. {Ucentral -> cLight^2 epsU, v2 -> cLight^2 epsV2}];
Liso2PNScaled = Expand[PNSeries[LisoScaled, {epsU, epsV2}, 3]];

checkEqScalar[
  "Isotropic Schwarzschild target has the 1PN static coefficient -1/2",
  Coefficient[Liso2PNScaled /. epsV2 -> 0, epsU, 2],
  -1/2
];

checkEqScalar[
  "Isotropic Schwarzschild target has the 2PN static coefficient +1/4",
  Coefficient[Liso2PNScaled /. epsV2 -> 0, epsU, 3],
  1/4
];

checkEqScalar[
  "Isotropic Schwarzschild target has the 2PN U^2 v^2 coefficient 2",
  Coefficient[Coefficient[Liso2PNScaled, epsU, 2], epsV2],
  2
];

Liso2PNScaledDyn = Expand[Liso2PNScaled + 1];
ResidualCandidateToIso = Expand[LtestMassCandidate2PNScaledDyn - Liso2PNScaledDyn];
lambdaRhoStaticResidual = FullSimplify[SeriesCoefficient[ResidualCandidateToIso /. epsV2 -> 0, {epsU, 0, 3}]];
lambdaRhoIsoFitRules = Quiet @ Solve[lambdaRhoStaticResidual == 0, lambdaRho];
lambdaRhoIsoFit = If[Length[lambdaRhoIsoFitRules] >= 1,
  FullSimplify[lambdaRho /. First[lambdaRhoIsoFitRules]],
  Missing["NoSolution"]
];
ResidualCandidateToIsoStaticMatched = If[MissingQ[lambdaRhoIsoFit],
  Missing["NoSolution"],
  Expand[ResidualCandidateToIso /. lambdaRho -> lambdaRhoIsoFit]
];

If[MissingQ[lambdaRhoIsoFit],
  fail["Pure-static isotropic target matching would set lambda_rho = 1/2", "no solution from Solve"],
  checkEqScalar[
    "Pure-static isotropic target matching would set lambda_rho = 1/2",
    lambdaRhoIsoFit,
    1/2
  ]
];

If[MissingQ[ResidualCandidateToIsoStaticMatched],
  fail["After fixing lambda_rho by the pure-static isotropic match, the remaining residual is 4 U^2 v^2", "no residual because the lambda fit failed"],
  checkEqScalar[
    "After fixing lambda_rho by the pure-static isotropic match, the remaining residual is 4 U^2 v^2",
    ResidualCandidateToIsoStaticMatched,
    4 epsU^2 epsV2
  ]
];

info["Within this candidate assembly, once lambda_rho is fixed by the pure-static isotropic test-mass match, the unresolved test-mass discrepancy sits entirely in the U^2 v^2 sector."];

(* ---------------------------------------------------------------------- *)
section["WORLDTUBE / finite-size gate on any claimed universal 1/r^4 term"];
(* ---------------------------------------------------------------------- *)

FNewtonScale = Gconst Mcentral/rSep^2;
FQuadScale = FNewtonScale (aBody/rSep)^2;
F2PNUniversalScale = FNewtonScale ((Gconst Mcentral/rSep)/cLight^2)^2;
ratioFiniteToUniversal2PN = FullSimplify[
  FQuadScale/F2PNUniversalScale,
  cLight > 0 && Gconst > 0 && Mcentral > 0 && aBody > 0 && rSep > 0
];

checkEqScalar[
  "Finite-size / universal-2PN force ratio collapses to (a c^2 / (G M))^2",
  ratioFiniteToUniversal2PN,
  aBody^2 cLight^4/(Gconst^2 Mcentral^2),
  cLight > 0 && Gconst > 0 && Mcentral > 0 && aBody > 0 && rSep > 0
];

info["Universal 2PN domination therefore requires a << G M / c^2, not merely a << r."];

(* ---------------------------------------------------------------------- *)
section["SUMMARY"];
(* ---------------------------------------------------------------------- *)

TwoPNResults = <|
  "kappaRhoSolution" -> kappaRhoSolution,
  "nSolution" -> nSolution,
  "kappaAdd" -> kappaAdd[3],
  "kappaPV" -> kappaPV,
  "beta1PN" -> beta1PN,
  "BernoulliDensity" -> rhoBernoulli,
  "BernoulliSoundSpeedSq" -> cs2Bernoulli,
  "BernoulliIndex" -> NoptBernoulli,
  "Lself2PNScaledU" -> Lself2PNScaledU,
  "StaticCoeff1PN" -> staticCoeff1PN,
  "StaticCoeff2PNTwoBody" -> staticCoeff2PNTwoBody,
  "ThreeBodyCoeff_uAB2uAC" -> coeff3AB2AC,
  "ThreeBodyCoeff_uABuAC2" -> coeff3ABAC2,
  "ThreeBodyPureTripleProductCoeff" -> coeff3ABCmixed,
  "FourBodyCoeff_uABuACuAD" -> coeff4ABACAD,
  "HeffQuadratic" -> HeffQuadratic,
  "HeffCubic" -> HeffCubic,
  "L2BodyCandidateThrough2PN" -> L2BodyCandidateThrough2PN,
  "LtestMassCandidate2PNScaled" -> LtestMassCandidate2PNScaled,
  "LtestMassCandidate2PNScaledDyn" -> LtestMassCandidate2PNScaledDyn,
  "LtestMassCandidate2PNScaledFromTwoBody" -> LtestMassCandidate2PNScaledFromTwoBody,
  "Liso2PNScaled" -> Liso2PNScaled,
  "Liso2PNScaledDyn" -> Liso2PNScaledDyn,
  "lambdaRhoIsoFit" -> lambdaRhoIsoFit,
  "ResidualCandidateToIsoStaticMatched" -> ResidualCandidateToIsoStaticMatched,
  "FiniteSizeToUniversal2PNRatio" -> ratioFiniteToUniversal2PN
|>;

Print["Key exported symbols: TwoPNResults, Lself2PNScaledU, L2BodyCandidateThrough2PN, LtestMassCandidate2PNScaled, LtestMassCandidate2PNScaledDyn, LtestMassCandidate2PNScaledFromTwoBody, Liso2PNScaled, ResidualCandidateToIsoStaticMatched, lambdaRhoIsoFit."];
Print["Passes: ", passCount];
Print["Fails : ", failCount];
Print["Skips : ", skipCount];

If[failCount == 0,
  Print["\nALL CHECKS PASSED."],
  Print["\nSOME CHECKS FAILED. Inspect the residuals above."]
];

(*"
Output:

--- 4D gravity + candidate 2PN master harness ---

=== CARRY-FORWARD / frozen 1PN regressions ===
PASS: Optical weak-field matching still fixes n = 5
PASS: Scalar-sector Newtonian coupling coefficient remains -kappa_rho
PASS: Newtonian matching still fixes kappa_rho = 1 (historically q = 1)
PASS: Free kinetic coefficient remains 1/2
PASS: Free-particle quartic kinetic coefficient remains 1/8 in the Lagrangian
PASS: Scalar+optical 1PN self coefficient remains kappa_rho/2 - (n-1)/2
PASS: With kappa_rho = 1 and n = 5 the one-body Phi v^2 coefficient is still -3/2
PASS: Pair-counted local mass scaling still gives the 1PN static coefficient -1/2
PASS: Added-mass closure still gives kappa_add = 1/2 in 3D
PASS: At n = 5 the adiabatic response slope remains (11+7x)/(4+5x)
PASS: Adiabatic internal response closure still gives x = 2/11
PASS: Adiabatic internal response closure still gives kappa_PV = 3/2
PASS: The 1PN ledger still closes to beta = 3

=== EXACT / n = 5 Bernoulli continuation of the barotropic sector ===
PASS: Background sound-speed normalization gives h0 = c^2/4 in reduced matter units
PASS: Exact n=5 Bernoulli density profile solves h(rho) = h0 + Phi
PASS: Exact Bernoulli closure gives cs^2/c^2 = 1 + 4 Phi/c^2
PASS: Exact Bernoulli closure gives N(Phi) = (1 + 4 Phi/c^2)^(-1/2)
PASS: Bernoulli optical index expands as 1 + 2 U/c^2 + 6 U^2/c^4 through 2PN

=== CONTROLLED / minimal Bernoulli exactification of the self sector through 2PN ===
PASS: Minimal exactification still reproduces the 1PN self coefficient -3/2
PASS: Minimal exactification gives the universal free 2PN kinetic coefficient 1/16
PASS: Minimal exactification gives the candidate U v^4 coefficient 7/8
PASS: Minimal exactification gives the candidate U^2 v^2 coefficient 6
INFO: Lself2PNScaledU stores the dimensionless one-body candidate through 2PN in U = -Phi variables.

=== CAUTION / fully Bernoulli-exactified mass dressing fails a frozen 1PN regression ===
PASS: Fully Bernoulli-exactified mass dressing would generate +3/2 epsPhi^2 at 1PN
INFO: That +3/2 static coefficient is incompatible with the already-frozen pair-counted 1PN static ledger, so only the optical denominator is exactified here.

=== EXACT / quadratic local mass scaling and the 2PN static hierarchy ===
PASS: Quadratic local mass scaling gives lambda_rho/2 on the 2PN two-body static term
PASS: Selected 3-body cubic-static monomials carry the expected unit-normalized coefficients
PASS: Quadratic local mass scaling generates no pure uAB uAC uBC triple-product term at cubic order
PASS: The first genuine 4-body cubic-static monomial has the expected unit-normalized coefficient

=== CONTROLLED / generic cubic adiabatic-elimination closure for one fast internal mode ===
PASS: First-order fast-mode response coefficient is A = -f/K
PASS: Second-order fast-mode response coefficient is B = -(g A + h/2 + C A^2/2)/K
PASS: Stationary condition vanishes through O(eps^2) after substituting A and B
PASS: Direct adiabatic-elimination series matches the expected through O(eps^3)
PASS: Adiabatic elimination shifts the quadratic coefficient to H2 - f^2/K
PASS: Adiabatic elimination shifts the cubic coefficient to H3 - 3 h f/K + 3 g f^2/K^2 - C f^3/K^3
INFO: HeffCubic is the natural one-DOF symbolic template for a second-order response contribution to the 2PN ledger.

=== FULL / self+static-only comparable-mass 2PN candidate and test-mass reduction ===
INFO: L2BodyCandidateThrough2PN is the explicit self+static-only comparable-mass candidate through 2PN; no 2PN wake/cross sector has been added.
PASS: Two-body reduction agrees with the direct self+static test-mass construction through 2PN
PASS: Self+static test-mass candidate has the expected U^3 coefficient lambda_rho/2
PASS: Self+static test-mass candidate has the expected U^2 v^2 coefficient 6

=== TARGET / exact Schwarzschild-isotropic test-mass target through 2PN ===
PASS: Isotropic Schwarzschild target has the 1PN static coefficient -1/2
PASS: Isotropic Schwarzschild target has the 2PN static coefficient +1/4
PASS: Isotropic Schwarzschild target has the 2PN U^2 v^2 coefficient 2
PASS: Pure-static isotropic target matching would set lambda_rho = 1/2
PASS: After fixing lambda_rho by the pure-static isotropic match, the remaining residual is 4 U^2 v^2
INFO: Within this candidate assembly, once lambda_rho is fixed by the pure-static isotropic test-mass match, the unresolved test-mass discrepancy sits entirely in the U^2 v^2 sector.

=== WORLDTUBE / finite-size gate on any claimed universal 1/r^4 term ===
PASS: Finite-size / universal-2PN force ratio collapses to (a c^2 / (G M))^2
INFO: Universal 2PN domination therefore requires a << G M / c^2, not merely a << r.

=== SUMMARY ===
Key exported symbols: TwoPNResults, Lself2PNScaledU, L2BodyCandidateThrough2PN, LtestMassCandidate2PNScaled, LtestMassCandidate2PNScaledDyn, LtestMassCandidate2PNScaledFromTwoBody, Liso2PNScaled, ResidualCandidateToIsoStaticMatched, lambdaRhoIsoFit.
Passes: 42
Fails : 0
Skips : 0

ALL CHECKS PASSED.
"*)
