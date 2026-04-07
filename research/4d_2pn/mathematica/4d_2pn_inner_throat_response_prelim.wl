(**
  4D -> preliminary 2PN inner-throat response prototype
  -----------------------------------------------------
  Purpose:
    A compact Wolfram Language companion to the SymPy prototype that isolates the
    missing one-body 2PN response slot left open by 4d_gravity_2pn_master_harness.wl.

  What this script does:
    1) Rebuilds the current self+static test-mass candidate through 2PN.
    2) Shows the exact isotropic Schwarzschild residual is +4 U^2 v^2.
    3) Introduces a minimal denominator correction
           D(U)=1-4U/c^2+chi (U/c^2)^2
       and solves chi = 8.
    4) Verifies that chi = 8 reproduces the exact isotropic test-mass target through 2PN.
    5) Extends the 1DOF adiabatic throat closure to second order:
         - solves a(rho) about the operating point,
         - derives F_eq(rho),
         - isolates the pure internal response factor,
         - composes with exact n=5 Bernoulli density.
    6) Shows that a minimal quadratic geometry dressing
           D_eff = (1 - 4u)(1 + mu (delta a)^2)
       can realize chi = 8 without disturbing the frozen 1PN slot.

  Scope note:
    This file is a targeted prototype. It does not close the full comparable-mass 2PN
    wake/cross sector.
**)

ClearAll["Global`*"];
Print["--- 4D preliminary 2PN inner-throat response prototype ---"];

passCount = 0;
failCount = 0;
section[name_String] := Print["\n=== ", name, " ==="];
pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);
info[msg_String] := Print["INFO: ", msg];

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet @ FullSimplify[res == 0, assum]]
];
checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  If[zeroScalarQ[expr, assum], pass[name], fail[name, res]]
];
checkEqScalar[name_String, lhs_, rhs_, assum_: True] := checkScalar[name, lhs - rhs, assum];

PNSeries[expr_, vars_List, order_Integer] := Module[{lam, rules},
  rules = (# -> lam #) & /@ vars;
  Expand[Normal[Series[Expand[expr /. rules], {lam, 0, order}]] /. lam -> 1]
];

(* ---------------------------------------------------------------------- *)
section["Baseline one-body 2PN residual"]; 
(* ---------------------------------------------------------------------- *)

Lbaseline = -(1 - epsU) Sqrt[1 - epsV2/(1 - 4 epsU)] - epsU^2/2 + epsU^3/4;
Lbaseline2PN = Expand[PNSeries[Lbaseline, {epsU, epsV2}, 3]];

LisoExact = -Sqrt[((1 - epsU/2)/(1 + epsU/2))^2 - (1 + epsU/2)^4 epsV2];
Liso2PN = Expand[PNSeries[LisoExact, {epsU, epsV2}, 3]];

ResidualBaseline = Expand[Lbaseline2PN - Liso2PN];

checkEqScalar[
  "Current self+static one-body candidate misses the isotropic target only by +4 U^2 v^2",
  ResidualBaseline,
  4 epsU^2 epsV2
];

(* ---------------------------------------------------------------------- *)
section["Minimal denominator correction and isotropic solve"]; 
(* ---------------------------------------------------------------------- *)

Lchi = -(1 - epsU) Sqrt[1 - epsV2/(1 - 4 epsU + chiResp epsU^2)] - epsU^2/2 + epsU^3/4;
Lchi2PN = Expand[PNSeries[Lchi, {epsU, epsV2}, 3]];
chiU2V2Coeff = FullSimplify[Coefficient[Coefficient[Lchi2PN, epsU, 2], epsV2], True];
chiSolution = chiResp /. First[Solve[chiU2V2Coeff == Coefficient[Coefficient[Liso2PN, epsU, 2], epsV2], chiResp]];

checkEqScalar[
  "With D = 1 - 4u + chi u^2, the U^2 v^2 coefficient becomes 6 - chi/2",
  chiU2V2Coeff,
  6 - chiResp/2
];

checkEqScalar[
  "Exact isotropic one-body matching fixes chi = 8",
  chiSolution,
  8
];

checkEqScalar[
  "Substituting chi = 8 reproduces the isotropic one-body target through 2PN",
  Expand[Lchi2PN /. chiResp -> chiSolution],
  Liso2PN
];

info["The dynamic denominator D(u)=1-4u+8u^2 is the unique minimal one-body fix within this separation of static vs dynamic sectors."];

(* ---------------------------------------------------------------------- *)
section["Second-order 1DOF throat closure around the operating point"]; 
(* ---------------------------------------------------------------------- *)

Fclosure = 11 rhoVar^2/aVar + 2/(rhoVar aVar^2) + 5 rhoVar^5 aVar^3;
aAnsatz = 1 + A1 epsRho + A2 epsRho^2 + A3 epsRho^3;
stationarySeries = Expand[Normal[Series[(D[Fclosure, aVar] /. {rhoVar -> 1 + epsRho, aVar -> aAnsatz}), {epsRho, 0, 3}]]];
A1Solution = A1 /. First[Solve[Coefficient[stationarySeries, epsRho, 1] == 0, A1]];
A2Solution = A2 /. First[Solve[(Coefficient[stationarySeries, epsRho, 2] /. A1 -> A1Solution) == 0, A2]];
A3Solution = A3 /. First[Solve[(Coefficient[stationarySeries, epsRho, 3] /. {A1 -> A1Solution, A2 -> A2Solution}) == 0, A3]];
closureRule = {A1 -> A1Solution, A2 -> A2Solution, A3 -> A3Solution};

checkEqScalar[
  "Second-order closure reproduces the known breathing slope A1 = -57/64",
  A1Solution,
  -57/64
];

FeqSeries = Expand[Normal[Series[(Fclosure /. {rhoVar -> 1 + epsRho, aVar -> aAnsatz} /. closureRule), {epsRho, 0, 3}]]];
FeqRatio = FullSimplify[FeqSeries/(FeqSeries /. epsRho -> 0)];
RPVeps = Expand[Normal[Series[FeqRatio/(1 + epsRho), {epsRho, 0, 3}]]];

checkEqScalar[
  "Pure internal response factor has the linear coefficient +3/2 in eps = delta rho / rho0",
  Coefficient[RPVeps, epsRho, 1],
  3/2
];

info["RPVeps stores the pure internal response factor after dividing out the baseline density scaling kappa_rho = 1."];

(* ---------------------------------------------------------------------- *)
section["Compose with exact n=5 Bernoulli density"]; 
(* ---------------------------------------------------------------------- *)

rhoBernoulliU = Expand[Normal[Series[(1 - 4 epsU)^(1/4), {epsU, 0, 3}]]];
aOfU = Expand[Normal[Series[(aAnsatz /. closureRule /. epsRho -> (rhoBernoulliU - 1)), {epsU, 0, 3}]]];
RPVU = Expand[Normal[Series[(RPVeps /. epsRho -> (rhoBernoulliU - 1)), {epsU, 0, 3}]]];

deltaA = Expand[aOfU - 1];

checkEqScalar[
  "Bernoulli-composed throat breathing starts as a(u) = 1 + 57/64 u + ...",
  Coefficient[aOfU, epsU, 1],
  57/64
];

checkEqScalar[
  "Bernoulli-composed pure response factor starts as R_PV(u) = 1 - 3/2 u + ...",
  Coefficient[RPVU, epsU, 1],
  -3/2
];

(* ---------------------------------------------------------------------- *)
section["Minimal quadratic geometry dressing"]; 
(* ---------------------------------------------------------------------- *)

muSolution = FullSimplify[chiSolution/Coefficient[Expand[deltaA^2], epsU, 2]];
Deff = Expand[Normal[Series[(1 - 4 epsU) (1 + muSolution deltaA^2), {epsU, 0, 2}]]];
DrawResonance = Expand[Normal[Series[(1 - 4 epsU)/aOfU^2, {epsU, 0, 2}]]];

checkEqScalar[
  "A quadratic geometry dressing D_eff = (1-4u)(1 + mu deltaa^2) can realize chi = 8",
  Deff,
  1 - 4 epsU + 8 epsU^2
];

info[StringJoin["mu = ", ToString[TraditionalForm[muSolution]], " is the required quadratic geometry coefficient in the minimal dressing ansatz."]];
info[StringJoin["Raw resonance proxy c_s^2/L^2 gives u^2 coefficient ", ToString[TraditionalForm[Coefficient[DrawResonance, epsU, 2]]], ", so the sign is right but the magnitude is not yet enough."]];

(* ---------------------------------------------------------------------- *)
section["Summary"]; 
(* ---------------------------------------------------------------------- *)

PrelimResponseResults = <|
  "BaselineResidual" -> ResidualBaseline,
  "Lbaseline2PN" -> Lbaseline2PN,
  "Liso2PN" -> Liso2PN,
  "Lchi2PN" -> Lchi2PN,
  "chiSolution" -> chiSolution,
  "FeqRatio" -> FeqRatio,
  "RPVeps" -> RPVeps,
  "rhoBernoulliU" -> rhoBernoulliU,
  "aOfU" -> aOfU,
  "RPVU" -> RPVU,
  "deltaA" -> deltaA,
  "muSolution" -> muSolution,
  "Deff" -> Deff,
  "DrawResonance" -> DrawResonance
|>;

Print["Key exported symbol: PrelimResponseResults."];
Print["Passes: ", passCount];
Print["Fails : ", failCount];
If[failCount == 0,
  Print["\nALL CHECKS PASSED."],
  Print["\nSOME CHECKS FAILED. Inspect the residuals above."]
];

(*"
Output:

--- 4D preliminary 2PN inner-throat response prototype ---

=== Baseline one-body 2PN residual ===
PASS: Current self+static one-body candidate misses the isotropic target only by +4 U^2 v^2

=== Minimal denominator correction and isotropic solve ===
PASS: With D = 1 - 4u + chi u^2, the U^2 v^2 coefficient becomes 6 - chi/2
PASS: Exact isotropic one-body matching fixes chi = 8
PASS: Substituting chi = 8 reproduces the isotropic one-body target through 2PN
INFO: The dynamic denominator D(u)=1-4u+8u^2 is the unique minimal one-body fix within this separation of static vs dynamic sectors.

=== Second-order 1DOF throat closure around the operating point ===
PASS: Second-order closure reproduces the known breathing slope A1 = -57/64
PASS: Pure internal response factor has the linear coefficient +3/2 in eps = delta rho / rho0
INFO: RPVeps stores the pure internal response factor after dividing out the baseline density scaling kappa_rho = 1.

=== Compose with exact n=5 Bernoulli density ===
PASS: Bernoulli-composed throat breathing starts as a(u) = 1 + 57/64 u + ...
PASS: Bernoulli-composed pure response factor starts as R_PV(u) = 1 - 3/2 u + ...

=== Minimal quadratic geometry dressing ===
PASS: A quadratic geometry dressing D_eff = (1-4u)(1 + mu deltaa^2) can realize chi = 8
INFO: mu = DisplayForm[FormBox[FractionBox[32768, 3249], TraditionalForm]] is the required quadratic geometry coefficient in the minimal dressing ansatz.
INFO: Raw resonance proxy c_s^2/L^2 gives u^2 coefficient DisplayForm[FormBox[FractionBox[324075, 65536], TraditionalForm]], so the sign is right but the magnitude is not yet enough.

=== Summary ===
Key exported symbol: PrelimResponseResults.
Passes: 9
Fails : 0

ALL CHECKS PASSED.
"*)
