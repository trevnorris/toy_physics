(* Direct shooting of the unit-source field and its running shell integral.
   Regular interior means zero derivative at a. The decaying exterior
   imposes h[b]+b h'[b]=0. No analytic escape formula enters the solve. *)
ClearAll["Global`*"];
oneMouth[a_] := Module[{b = a + 1/20, kappa = 7/8, km = 1, eta, shoot,
    base, unit, atEdge, mixed, z, self},
  eta = 3/(4 Pi (b^3-a^3));
  shoot[initial_] := Through[NDSolveValue[{
      h'[r] == v[r], kappa (v'[r] + 2 v[r]/r) == eta (km h[r] - 1),
      integral'[r] == 4 Pi r^2 eta h[r], h[a] == initial,
      v[a] == 0, integral[a] == 0}, {h, v, integral}, {r, a, b},
      WorkingPrecision -> 40, PrecisionGoal -> 25, AccuracyGoal -> 27][b]];
  base = shoot[0]; unit = shoot[1];
  atEdge[row_] := row[[1]] + b row[[2]];
  mixed = base - atEdge[base] (unit - base)/atEdge[unit - base];
  z = -4 Pi kappa b^2 mixed[[2]];
  self = mixed[[3]];
  If[! TrueQ[Abs[atEdge[mixed]] < 10^-20 && Abs[z + km self - 1] < 10^-19 &&
       z > 0 && self > 0], Exit[1]];
  Print["radius ", a, ", escape and self response = ", InputForm[N[{z, self}, 14]]];
  {z, self}
];
rows = oneMouth /@ {1, 2};
pairs[v_] := {v[[1]]^2, v[[1]] v[[2]], v[[2]]^2};
z = rows[[All, 1]]; self = rows[[All, 2]];
Print["held unit mean bend = ", InputForm[N[pairs[z/self], 14]]];
Print["unit integrated conormal, replacement = ", InputForm[N[pairs[z], 14]]];
Print["only the unit committed source = ", InputForm[N[-pairs[z], 14]]];
Print["required bend and conormal ratios = ",
  InputForm[N[{(self[[2]]/z[[2]])/(self[[1]]/z[[1]]), z[[1]]/z[[2]]}, 14]]];
Exit[0];
