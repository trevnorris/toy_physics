(* Independent normal-coordinate solution of the field-inertia benchmark.
   z is time in relative-mode units, eta is the wave fraction of inertia. *)
ClearAll["Global`*"];
$HistoryLength = 0;
verify[condition_] := If[!TrueQ[FullSimplify[condition, 0 < eta < 1]], Exit[1]];
path[drive_] := Module[{gap, mean},
  gap = DSolveValue[{b''[z] + b[z] == drive, b[0] == 0, b'[0] == 0}, b[z], z];
  mean = DSolveValue[{c''[z] == (1 - eta) drive, c[0] == 0, c'[0] == 0}, c[z], z];
  {mean + eta gap, mean - (1 - eta) gap}
];
full = path[1];
Print["unit-drive {rim,wave}: ", InputForm[full]];
Do[
  {rim, wave} = path[drive];
  verify[D[rim, {z, 2}] == drive - eta (rim-wave)];
  verify[D[wave, {z, 2}] == (1 - eta) (rim-wave)];
  halfCheck = FullSimplify[{rim, wave} - drive full];
  verify[halfCheck == {0, 0}];
  Print["drive, gap at pi and 2pi, trajectory prediction residuals: ",
        InputForm[{drive, Simplify[(rim-wave) /. z -> Pi],
          Simplify[(rim-wave) /. z -> 2 Pi], halfCheck}]],
  {drive, {0, 1, -1, 1/2, -1/2}}];
(* Infer a restoring coefficient from displacement alone by integrating the
   relative displacement twice with zero initial data. *)
integratedGap = Integrate[(z-u) ((full[[1]]-full[[2]]) /. z->u), {u,0,z}];
verify[full[[2]] == (1-eta) integratedGap];
det = Factor[Det[{full /. z->Pi, full /. z->2 Pi}]];
Print["full-drive determinant: ", InputForm[det]];
verify[det == 4 Pi^2 (1-eta)];
Exit[0];
