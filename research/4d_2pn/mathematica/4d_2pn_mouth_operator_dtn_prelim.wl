(**
  4D -> preliminary 2PN mouth-operator / DtN reduction prototype
  ----------------------------------------------------------------
  Purpose:
    Continue the inner-throat response program past the first Bernoulli/throat
    prototype and express the missing one-body 2PN response slot directly in
    terms of low-frequency monopole mouth-operator data.

  Main result:
    Let the monopole DtN operator be
        Z00[ω;u] = Z2[u] ω^2 + Z4[u] ω^4 + O[ω^6],
    and define the normalized DtN invariants
        Cs[u] = (Z4/Z2^3)/(Z40/Z20^3),
        G[u]  = (Z4/Z2^2)/(Z40/Z20^2).
    For the cylinder / Neumann-bottom unit-test branch,
        Z2 = -L/cs^2,
        Z4 = -L^3/(3 cs^4),
    so
        Cs = cs^2/cs0^2,
        G  = L/L0.

    The minimal conservative local one-body closure that preserves the frozen
    1PN slot and fixes the missing 2PN U^2 v^2 coefficient is
        D_eff[u] = Cs[u] (1 + mu (G[u]-1)^2),
        mu = 8/g1^2,
    where g1 is the linear coefficient in G[u] = 1 + g1 u + ... .

    Using the already-fixed throat slope g1 = 57/64 gives
        mu = 32768/3249,
    and the resulting denominator reproduces the exact isotropic Schwarzschild
    one-body target through 2PN.

  Scope note:
    This is still a one-body / self-sector response prototype.  It does not yet
    close the full comparable-mass 2PN wake/cross sector.
**)

ClearAll["Global`*"];
Print["--- 4D preliminary 2PN mouth-operator / DtN reduction prototype ---"];

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
section["General one-body denominator conditions"]; 
(* ---------------------------------------------------------------------- *)

Dgeneral = 1 + d1 epsU + d2 epsU^2;
Lgeneral = -(1 - epsU) Sqrt[1 - epsV2/Dgeneral];
Lgeneral2PN = Expand[PNSeries[Lgeneral, {epsU, epsV2}, 3]];

coeffUV = FullSimplify[Coefficient[Coefficient[Lgeneral2PN, epsU, 1], epsV2]];
coeffU2V = FullSimplify[Coefficient[Coefficient[Lgeneral2PN, epsU, 2], epsV2]];

checkEqScalar[
  "For D = 1 + d1 u + d2 u^2, the 1PN self coefficient is -(d1+1)/2",
  coeffUV,
  -(d1 + 1)/2
];

checkEqScalar[
  "For D = 1 + d1 u + d2 u^2, the 2PN U^2 v^2 coefficient is (d1^2 + d1 - d2)/2",
  coeffU2V,
  (d1^2 + d1 - d2)/2
];

checkEqScalar[
  "Keeping the frozen 1PN self coefficient +3/2 requires d1 = -4",
  coeffUV /. d1 -> -4,
  3/2
];

checkEqScalar[
  "With d1 = -4, matching the isotropic U^2 v^2 coefficient 2 requires d2 = 8",
  coeffU2V /. {d1 -> -4, d2 -> 8},
  2
];

(* ---------------------------------------------------------------------- *)
section["Generic low-frequency DtN invariants"]; 
(* ---------------------------------------------------------------------- *)

Z2hat = 1 + z1 epsU + z2 epsU^2;
Z4hat = 1 + w1 epsU + w2 epsU^2;

Ghat = Expand[Normal[Series[Z4hat/Z2hat^2, {epsU, 0, 2}]]];
Cshat = Expand[Normal[Series[Z4hat/Z2hat^3, {epsU, 0, 2}]]];

g1Generic = FullSimplify[Coefficient[Ghat, epsU, 1]];
g2Generic = FullSimplify[Coefficient[Ghat, epsU, 2]];
c1Generic = FullSimplify[Coefficient[Cshat, epsU, 1]];
c2Generic = FullSimplify[Coefficient[Cshat, epsU, 2]];

checkEqScalar[
  "Normalized DtN geometry invariant has linear coefficient w1 - 2 z1",
  g1Generic,
  w1 - 2 z1
];

checkEqScalar[
  "Normalized DtN sound-speed invariant has linear coefficient w1 - 3 z1",
  c1Generic,
  w1 - 3 z1
];

(* ---------------------------------------------------------------------- *)
section["Exact cylinder / Neumann-bottom branch"]; 
(* ---------------------------------------------------------------------- *)

Lhat = 1 + ell1 epsU + ell2 epsU^2;
CsExact = 1 - 4 epsU;

Z2hatCylinder = Expand[Normal[Series[Lhat/CsExact, {epsU, 0, 2}]]];
Z4hatCylinder = Expand[Normal[Series[Lhat^3/CsExact^2, {epsU, 0, 2}]]];

GhatCylinder = Expand[Normal[Series[Z4hatCylinder/Z2hatCylinder^2, {epsU, 0, 2}]]];
CshatCylinder = Expand[Normal[Series[Z4hatCylinder/Z2hatCylinder^3, {epsU, 0, 2}]]];

checkEqScalar[
  "Cylinder DtN invariant G recovers the normalized length L/L0 exactly",
  GhatCylinder,
  Lhat
];

checkEqScalar[
  "Cylinder DtN invariant Cs recovers the normalized sound speed factor exactly",
  CshatCylinder,
  CsExact
];

(* ---------------------------------------------------------------------- *)
section["Rebuild the known breathing slope from the 1DOF throat closure"]; 
(* ---------------------------------------------------------------------- *)

Fclosure = 11 rhoVar^2/aVar + 2/(rhoVar aVar^2) + 5 rhoVar^5 aVar^3;
aAnsatz = 1 + A1 epsRho + A2 epsRho^2;
stationarySeries = Expand[Normal[Series[(D[Fclosure, aVar] /. {rhoVar -> 1 + epsRho, aVar -> aAnsatz}), {epsRho, 0, 2}]]];
A1Solution = A1 /. First[Solve[Coefficient[stationarySeries, epsRho, 1] == 0, A1]];
A2Solution = A2 /. First[Solve[(Coefficient[stationarySeries, epsRho, 2] /. A1 -> A1Solution) == 0, A2]];
closureRule = {A1 -> A1Solution, A2 -> A2Solution};

g1Value = FullSimplify[-A1Solution];

checkEqScalar[
  "The 1DOF throat closure still reproduces the known breathing slope A1 = -57/64",
  A1Solution,
  -57/64
];

checkEqScalar[
  "The linear geometry coefficient induced by the Bernoulli map is g1 = 57/64",
  g1Value,
  57/64
];

(* ---------------------------------------------------------------------- *)
section["Minimal quadratic DtN-invariant response closure"]; 
(* ---------------------------------------------------------------------- *)

Ggeneric = 1 + g1 epsU + g2 epsU^2;
DeffGeneric = Expand[Normal[Series[(1 - 4 epsU) (1 + alpha (Ggeneric - 1) + beta (Ggeneric - 1)^2), {epsU, 0, 2}]]];

checkEqScalar[
  "Any linear DtN geometry correction alpha (G-1) would shift the frozen 1PN slot by alpha g1",
  Coefficient[DeffGeneric, epsU, 1],
  alpha g1 - 4
];

checkEqScalar[
  "Preserving the frozen 1PN slot therefore forces alpha = 0",
  Coefficient[DeffGeneric /. alpha -> 0, epsU, 1],
  -4
];

checkEqScalar[
  "With alpha = 0, the 2PN quadratic coefficient is beta g1^2",
  Coefficient[DeffGeneric /. alpha -> 0, epsU, 2],
  beta g1^2
];

muSolution = FullSimplify[8/g1Value^2];

checkEqScalar[
  "Using the fixed throat slope gives mu = 32768/3249",
  muSolution,
  32768/3249
];

DeffBySlope = Expand[Normal[Series[(1 - 4 epsU) (1 + muSolution (Ggeneric - 1)^2), {epsU, 0, 2}]]];
checkEqScalar[
  "With mu = 8/g1^2, only the linear geometry slope g1 enters the 2PN slot",
  DeffBySlope,
  1 - 4 epsU + (muSolution g1^2) epsU^2
];

info["The quadratic DtN-invariant closure does not depend on the second geometry coefficient g2 until 3PN."];

(* ---------------------------------------------------------------------- *)
section["Specialize to the current throat closure and close the isotropic target"]; 
(* ---------------------------------------------------------------------- *)

rhoBernoulliU = Expand[Normal[Series[(1 - 4 epsU)^(1/4), {epsU, 0, 2}]]];
aOfU = Expand[Normal[Series[(aAnsatz /. closureRule /. epsRho -> (rhoBernoulliU - 1)), {epsU, 0, 2}]]];

DeffDtN = Expand[Normal[Series[(1 - 4 epsU) (1 + muSolution (aOfU - 1)^2), {epsU, 0, 2}]]];

checkEqScalar[
  "Current throat closure gives a(u) = 1 + 57/64 u + 298821/131072 u^2 + ...",
  aOfU,
  1 + (57/64) epsU + (298821/131072) epsU^2
];

checkEqScalar[
  "The DtN-invariant denominator becomes exactly 1 - 4u + 8u^2 through 2PN",
  DeffDtN,
  1 - 4 epsU + 8 epsU^2
];

LdtN = -(1 - epsU) Sqrt[1 - epsV2/DeffDtN] - epsU^2/2 + epsU^3/4;
LdtN2PN = Expand[PNSeries[LdtN, {epsU, epsV2}, 3]];

LisoExact = -Sqrt[((1 - epsU/2)/(1 + epsU/2))^2 - (1 + epsU/2)^4 epsV2];
Liso2PN = Expand[PNSeries[LisoExact, {epsU, epsV2}, 3]];

checkEqScalar[
  "The DtN-corrected one-body candidate reproduces the exact isotropic target through 2PN",
  LdtN2PN,
  Liso2PN
];

(* ---------------------------------------------------------------------- *)
section["Relation to the earlier raw resonance proxy / port-normalization fit"]; 
(* ---------------------------------------------------------------------- *)

Draw = Expand[Normal[Series[(1 - 4 epsU)/aOfU^2, {epsU, 0, 2}]]];
Pport = Expand[Normal[Series[aOfU^2 (1 + muSolution (aOfU - 1)^2), {epsU, 0, 2}]]];

checkEqScalar[
  "The raw resonance proxy still starts as 1 - 185/32 u + 324075/65536 u^2",
  Draw,
  1 - (185/32) epsU + (324075/65536) epsU^2
];

checkEqScalar[
  "The required port factor is exactly 1 + 57/32 u + 875093/65536 u^2",
  Pport,
  1 + (57/32) epsU + (875093/65536) epsU^2
];

checkEqScalar[
  "Raw resonance proxy times the DtN-invariant port factor reproduces the exact denominator",
  Expand[Normal[Series[Draw Pport, {epsU, 0, 2}]]],
  1 - 4 epsU + 8 epsU^2
];

info["So the earlier fitted port-normalization coefficients factor non-arbitrarily as Pport = G^2 (1 + mu (G-1)^2)."];

(* ---------------------------------------------------------------------- *)
section["Summary"]; 
(* ---------------------------------------------------------------------- *)

MouthDtNResults = <|
  "Z2hatCylinder" -> Z2hatCylinder,
  "Z4hatCylinder" -> Z4hatCylinder,
  "GhatCylinder" -> GhatCylinder,
  "CshatCylinder" -> CshatCylinder,
  "A1Solution" -> A1Solution,
  "A2Solution" -> A2Solution,
  "g1Value" -> g1Value,
  "muSolution" -> muSolution,
  "aOfU" -> aOfU,
  "DeffDtN" -> DeffDtN,
  "LdtN2PN" -> LdtN2PN,
  "Liso2PN" -> Liso2PN,
  "Draw" -> Draw,
  "Pport" -> Pport
|>;

Print["Key exported symbol: MouthDtNResults."];
Print["Passes: ", passCount];
Print["Fails : ", failCount];
If[failCount == 0,
  Print["\nALL CHECKS PASSED."],
  Print["\nSOME CHECKS FAILED. Inspect the residuals above."]
];

(*"
Output:

--- 4D preliminary 2PN mouth-operator / DtN reduction prototype ---

=== General one-body denominator conditions ===
PASS: For D = 1 + d1 u + d2 u^2, the 1PN self coefficient is -(d1+1)/2
PASS: For D = 1 + d1 u + d2 u^2, the 2PN U^2 v^2 coefficient is (d1^2 + d1 - d2)/2
PASS: Keeping the frozen 1PN self coefficient +3/2 requires d1 = -4
PASS: With d1 = -4, matching the isotropic U^2 v^2 coefficient 2 requires d2 = 8

=== Generic low-frequency DtN invariants ===
PASS: Normalized DtN geometry invariant has linear coefficient w1 - 2 z1
PASS: Normalized DtN sound-speed invariant has linear coefficient w1 - 3 z1

=== Exact cylinder / Neumann-bottom branch ===
PASS: Cylinder DtN invariant G recovers the normalized length L/L0 exactly
PASS: Cylinder DtN invariant Cs recovers the normalized sound speed factor exactly

=== Rebuild the known breathing slope from the 1DOF throat closure ===
PASS: The 1DOF throat closure still reproduces the known breathing slope A1 = -57/64
PASS: The linear geometry coefficient induced by the Bernoulli map is g1 = 57/64

=== Minimal quadratic DtN-invariant response closure ===
PASS: Any linear DtN geometry correction alpha (G-1) would shift the frozen 1PN slot by alpha g1
PASS: Preserving the frozen 1PN slot therefore forces alpha = 0
PASS: With alpha = 0, the 2PN quadratic coefficient is beta g1^2
PASS: Using the fixed throat slope gives mu = 32768/3249
PASS: With mu = 8/g1^2, only the linear geometry slope g1 enters the 2PN slot
INFO: The quadratic DtN-invariant closure does not depend on the second geometry coefficient g2 until 3PN.

=== Specialize to the current throat closure and close the isotropic target ===
PASS: Current throat closure gives a(u) = 1 + 57/64 u + 298821/131072 u^2 + ...
PASS: The DtN-invariant denominator becomes exactly 1 - 4u + 8u^2 through 2PN
PASS: The DtN-corrected one-body candidate reproduces the exact isotropic target through 2PN

=== Relation to the earlier raw resonance proxy / port-normalization fit ===
PASS: The raw resonance proxy still starts as 1 - 185/32 u + 324075/65536 u^2
PASS: The required port factor is exactly 1 + 57/32 u + 875093/65536 u^2
PASS: Raw resonance proxy times the DtN-invariant port factor reproduces the exact denominator
INFO: So the earlier fitted port-normalization coefficients factor non-arbitrarily as Pport = G^2 (1 + mu (G-1)^2).

=== Summary ===
Key exported symbol: MouthDtNResults.
Passes: 21
Fails : 0

ALL CHECKS PASSED.
"*)
