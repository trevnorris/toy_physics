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

(* Independent derivation path: this engine derives chi_Q = 1 via Apart-based
   partial-fraction verification, SeriesCoefficient operator-form extraction,
   and Reduce over the reals, rather than the .py's manual 3/4 + (1/4)/(...)
   decomposition followed by Series/Normal/Coefficient/Solve. *)

polescl = 3*cSound/(2*aThroat);          (* pole scale (formerly named omegaQ) *)
sigmaQcan = (9/8)/polescl^5;
Print["pole scale = ", fmt[polescl]];
Print["sigma_Q^can = ", fmt[sigmaQcan]];
expectZero["sigma_Q^can - 4 a^5/(27 c_s^5)", sigmaQcan - 4*aThroat^5/(27*cSound^5)];

(* Retarded module written as a single ratio rather than as a partial-
   fraction decomposition. Apart should recover the canonical
   3/4 + 1/(4 denomRet) form; we assert that. *)
denomRet = 1 - omega^2/polescl^2 - I*chiQ*sigmaQcan*omega^5;
yQretRatio = (4 - 3*omega^2/polescl^2 - 3*I*chiQ*sigmaQcan*omega^5)/(4*denomRet);
expectZero[
  "Y_Q^ret partial-fraction form check (Apart of unfactored ratio)",
  Together[yQretRatio - (3/4 + 1/(4*denomRet))]
];

(* Coefficient extraction via SeriesCoefficient (operator form), distinct from
   Series/Normal/Coefficient. *)
cRet2 = SeriesCoefficient[yQretRatio, {omega, 0, 2}];
cRet4 = SeriesCoefficient[yQretRatio, {omega, 0, 4}];
cRet5 = SeriesCoefficient[yQretRatio, {omega, 0, 5}];
Print["Y_Q^ret omega^2 coeff = ", fmt[cRet2]];
Print["Y_Q^ret omega^4 coeff = ", fmt[cRet4]];
Print["Y_Q^ret omega^5 coeff = ", fmt[cRet5]];
expectZero["omega^2 coefficient", cRet2 - aThroat^2/(9*cSound^2)];
expectZero["omega^4 coefficient", cRet4 - 4*aThroat^4/(81*cSound^4)];
expectZero["imag omega^5 coefficient", cRet5/I - chiQ*aThroat^5/(27*cSound^5)];

(* chi_Q identification via Reduce over the reals (not Solve over chi_Q
   alone). The (cRet5/I) is a linear function of chiQ, so Reduce returns
   a single substitution rule. *)
chiReduce = Reduce[cRet5/I == aThroat^5/(27*cSound^5), chiQ, Reals];
chiVal = chiQ /. ToRules[chiReduce];
Print["chi_Q from exact outgoing match (Reduce) = ", fmt[chiVal]];
expectZero["chi_Q - 1", chiVal - 1];

Clear[z, xiQ];
$Assumptions = Element[{z, xiQ}, Reals];

(* Deformed-branch coefficients via polynomial inversion of Lambda * Y = -3
   (the operator identity), rather than Y := -3/Lambda. This avoids
   reproducing the .py's division-based normalization. *)
lamDeformed = -3 + z^2/3 + z^4/9 + I*xiQ*z^5/9;
yAnsatz = b0 + b1*z + b2*z^2 + b3*z^3 + b4*z^4 + b5*z^5;
prodTrunc = Normal[Series[lamDeformed*yAnsatz, {z, 0, 5}]] // Expand;
prodCoeffs = CoefficientList[prodTrunc, z];
While[Length[prodCoeffs] < 6, AppendTo[prodCoeffs, 0]];
coeffSys = Thread[(prodCoeffs + Join[{3}, ConstantArray[0, 5]]) == 0];
ySolved = First[Solve[coeffSys, {b0, b1, b2, b3, b4, b5}]];
yDeformedSeries = Expand[yAnsatz /. ySolved];
Print["Deformed DtN normalized branch (polynomial-inverse) = ", fmt[yDeformedSeries]];
expectZero["deformed constant coefficient", (b0 /. ySolved) - 1];
expectZero["deformed z^2 coefficient", (b2 /. ySolved) - 1/9];
expectZero["deformed z^4 coefficient", (b4 /. ySolved) - 4/81];
expectZero["deformed imag z^5 coefficient", (b5 /. ySolved)/I - xiQ/27];

Print[""];
Print["Stage 105 Mathematica audit passed."];

Exit[0];
