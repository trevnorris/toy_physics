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

banner["STAGE 105 — EXACT FIXING OF chi_Q"];

Clear[aThroat, cSound, omega, chiQ];
$Assumptions =
  Element[{aThroat, cSound, omega, chiQ}, Reals] &&
  aThroat > 0 && cSound > 0 && chiQ > 0;

omegaPole = 3*cSound/(2*aThroat);
sigmaQ = (9/8)/omegaPole^5;
Print["pole scale = ", fmt[omegaPole]];
Print["sigma_Q = ", fmt[sigmaQ]];
expectZero["sigma_Q - 4 a^5/(27 c_s^5)", sigmaQ - 4*aThroat^5/(27*cSound^5)];

retFactored =
  (omegaPole^2 + 3*(omegaPole^2 - omega^2 - I*chiQ*sigmaQ*omegaPole^2*omega^5))/
  (4*(omegaPole^2 - omega^2 - I*chiQ*sigmaQ*omegaPole^2*omega^5));
retWindow = Expand[ComplexExpand[Normal[Series[retFactored /. omega -> omega, {omega, 0, 5}]]]];
imagProjection = Expand[ComplexExpand[Im[retWindow]]];
imagPart5 = Coefficient[imagProjection, omega, 5];
realPart = Expand[ComplexExpand[Re[retWindow]]];
Print["retarded real projection = ", fmt[realPart]];
Print["retarded imaginary projection = ", fmt[imagProjection]];
expectZero["omega^2 coefficient", Coefficient[realPart, omega, 2] - aThroat^2/(9*cSound^2)];
expectZero["omega^4 coefficient", Coefficient[realPart, omega, 4] - 4*aThroat^4/(81*cSound^4)];
expectZero["imag omega^5 coefficient", imagPart5 - chiQ*aThroat^5/(27*cSound^5)];

chiBranch = FullSimplify[
  Reduce[imagPart5 - aThroat^5/(27*cSound^5) == 0, chiQ, Reals],
  Assumptions -> $Assumptions
];
chiWitness = chiQ /. ToRules[chiBranch];
Print["chi_Q condition from outgoing match (Reduce) = ", fmt[chiBranch]];
Print["chi_Q from exact outgoing match = ", fmt[chiWitness]];
expectZero["chi_Q - 1", chiWitness - 1];

Clear[z, xiQ];
$Assumptions = Element[{z, xiQ}, Reals];

(* Deformed-branch coefficients via polynomial inversion of the operator
   identity Lambda * Y = -3. *)
outgoingOp = -3 + z^2/3 + z^4/9 + I*xiQ*z^5/9;
trialBranch = b0 + b1*z + b2*z^2 + b3*z^3 + b4*z^4 + b5*z^5;
prodTrunc = Normal[Series[outgoingOp*trialBranch, {z, 0, 5}]] // Expand;
prodCoeffs = CoefficientList[prodTrunc, z];
While[Length[prodCoeffs] < 6, AppendTo[prodCoeffs, 0]];
coeffSys = Thread[(prodCoeffs + Join[{3}, ConstantArray[0, 5]]) == 0];
branchRules = First[Solve[coeffSys, {b0, b1, b2, b3, b4, b5}]];
normSeries = Expand[trialBranch /. branchRules];
Print["Deformed DtN normalized branch (polynomial-inverse) = ", fmt[normSeries]];
expectZero["deformed constant coefficient", (b0 /. branchRules) - 1];
expectZero["deformed z^2 coefficient", (b2 /. branchRules) - 1/9];
expectZero["deformed z^4 coefficient", (b4 /. branchRules) - 4/81];
expectZero["deformed imag z^5 coefficient", (b5 /. branchRules)/I - xiQ/27];

Print[""];
Print["Stage 105 Mathematica audit passed."];

Exit[0];
