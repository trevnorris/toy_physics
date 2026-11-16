(* Independent amplitude and phase shooting calculation, in the same
   dimensionless units as the Cartesian boundary-value calculation.
   Both the initial amplitude slope and the conserved current are unknown. *)
ClearAll["Global`*"];
length = 1/2;
precision = 36;
endState = ParametricNDSolveValue[{
    amp''[w] == current^2/amp[w]^3 + 2 (amp[w]^8 - 1) amp[w],
    phase'[w] == current/amp[w]^2,
    amp[0] == 1, amp'[0] == slope, phase[0] == 0
  }, {amp[length], phase[length]}, {w, 0, length}, {slope, current},
  WorkingPrecision -> precision, AccuracyGoal -> 24, PrecisionGoal -> 24];

ampExit[s_?NumericQ, j_?NumericQ] := First[endState[s, j]];
phaseExit[s_?NumericQ, j_?NumericQ] := Last[endState[s, j]];
Do[
  roots = Check[FindRoot[{
      ampExit[s, j] == 9/10, phaseExit[s, j] == target
    }, {{s, N[-1/5, precision]}, {j, N[target/length, precision]}},
    WorkingPrecision -> precision, AccuracyGoal -> 18, PrecisionGoal -> 18],
    $Failed];
  If[roots === $Failed, Exit[1]];
  flux = j /. roots;
  mismatch = Max[Abs[(endState[s, j] /. roots) - {9/10, target}]];
  Print["phase=", NumberForm[N[target], {2, 1}],
    " current=", NumberForm[flux, {12, 10}, NumberPadding -> {"", "0"}],
    " endpoint_error_below_1e-16=", mismatch < 10^-16];
  If[! TrueQ[mismatch < 10^-16], Exit[1]],
  {target, {0, 1/10, 1/5}}];
