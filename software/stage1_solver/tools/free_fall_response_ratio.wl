(* Independent test: integrate exterior carriage and pressure, then drive a
   two-element variational field with its far endpoint free to move. *)
ClearAll["Global`*"];
$HistoryLength = 0;
check[label_, value_, expected_] := Module[{result = Factor[value]},
  Print[label, " = ", ToString[result, InputForm]];
  If[!TrueQ[FullSimplify[result == expected]], Exit[1]]
];

(* n points from the target into the exterior. The inward load is minus
   the outward momentum flux. u is the uniform incident speed along x. *)
vx = u - w cosTheta;
vn = u cosTheta - w;
pressure = pInf - rho (u^2 + w^2 - 2 u w cosTheta)/2;
load = -2 Pi radius^2 Integrate[rho vx vn + pressure cosTheta, {cosTheta, -1, 1}];
intake = 4 Pi radius^2 rho w;
check["exterior pull per incident speed", D[load, u], intake];

(* Interpolate transverse displacements on two equal-depth elements.
   mu is the unscaled line inertia; the same lam multiplies kinetic and
   strain energies. The common other fields have slow response mb s^2. *)
cell = len/2;
element[left_, right_] := Module[{shape = left (1 - t) + right t},
  Integrate[lam mu cell (s^2 shape^2 + (c D[shape,t]/cell)^2)/2, {t, 0, 1}]
];
work = element[x, y] + element[y, z];
interior = First[Solve[{D[work, y] == 0, D[work, z] == 0}, {y, z}]];
reaction = Factor[(D[work, x] /. interior)/x];
momentumRate = lam mu s^2 cell ((x+y)/2 + (y+z)/2);
check["force minus integrated momentum rate per displacement",
      reaction - (momentumRate /. interior)/x, 0];
check["static unit-displacement field", Limit[({y,z} /. interior)/x, s -> 0], {1,1}];
check["held reaction per incident speed", -D[load,u], -intake];
slowGain = Factor[Limit[s^2 D[load,u]/(mb s^2 + reaction), s -> 0]];
check["scale 1 gain per exterior pull coefficient", (slowGain /. lam -> 1)/intake,
      1/(mb + mu len)];
check["scale 2 gain per exterior pull coefficient", (slowGain /. lam -> 2)/intake,
      1/(mb + 2 mu len)];
Exit[0];
