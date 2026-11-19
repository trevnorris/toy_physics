(* First variation of the mouth work at fixed physical data. eps includes
   1/(4 Pi kappa d). Retain arbitrary one-mouth responses at both radii. *)
ClearAll["Global`*"];
check[label_, actual_, expected_] := Module[{residual = Factor[actual - expected]},
  Print[label, " = ", InputForm[Factor[actual]]];
  If[residual =!= 0, Print["residual = ", InputForm[residual]]; Exit[1]]
];

isolated = DiagonalMatrix[{selfA, selfB}];
propagation = {{0, za zb}, {za zb, 0}};
held = {sa pa, sb pb};
y0 = LinearSolve[isolated, held];
yPrime = LinearSolve[isolated, -propagation.y0];
odd[expr_] := Expand[expr] /. {sa^2 -> 0, sb^2 -> 0, sa -> 1, sb -> 1};
check["held bend", odd[-held.yPrime/2], za zb pa pb/(selfA selfB)];

(* With fixed conormal data only the local field changes. Each source
   reservoir does work on its own local field, including its perturbation. *)
force[ya_, yb_, wa_, wb_] := Module[{y = {sa ya, sb yb}, w = {sa wa, sb wb}},
  odd[(y/2 - w).(propagation.y)]
];
load = force[qa + g, qb + g, g, g];
source = force[ja + g, jb + g, ja + g, jb + g];
mixed = force[qa + g, qb + g, g + lam qa, g + lam qb];
check["fixed integrated load", load, za zb (qa qb - g^2)];
check["fixed energy source", source, -za zb (ja + g) (jb + g)];
check["fractional source work", mixed,
  za zb ((1 - 2 lam) qa qb - lam g (qa + qb) - g^2)];
check["two-radius determinant",
  (mixed /. {qb -> qa, zb -> za}) (mixed /. {qa -> qb, za -> zb}) - mixed^2,
  -(1 - lam)^2 (za zb g (qa - qb))^2];

(* Direct solution outside two disjoint spheres with prescribed signed
   values. Start each image chain at its sphere center. Each reflection
   cancels the preceding source's potential on the other sphere. The sum
   of image strengths inside each sphere gives its integrated conormal
   in units of 4 Pi kappa. No response-matrix coefficients are used here. *)
dirichletEnergy[ra_, rb_, d_, va_, vb_, reflections_] := Module[
  {radii = {ra, rb}, centers = {0, d}, values = {va, vb}, totals = {0, 0},
   index, next, position, strength, delta},
  Do[
    index = seed; position = centers[[seed]]; strength = radii[[seed]] values[[seed]];
    totals[[index]] += strength;
    Do[
      next = 3 - index;
      delta = position - centers[[next]];
      strength *= -radii[[next]]/Abs[delta];
      position = centers[[next]] + radii[[next]]^2/delta;
      totals[[next]] += strength;
      index = next,
      {reflections}],
    {seed, 2}];
  -values.totals/2
];
oddEnergy[ra_, rb_, d_, depth_] :=
  (dirichletEnergy[ra, rb, d, 1, 1, depth] -
   dirichletEnergy[ra, rb, d, 1, -1, depth])/2;
scaledForce[ra_, rb_, dist_, depth_, step_] := Module[{d = N[dist, 65]},
  -d^2 (oddEnergy[ra, rb, d + step, depth] -
        oddEnergy[ra, rb, d - step, depth])/(2 step)
];
pairs = {{1, 1}, {1, 2}, {2, 2}};
limit = Times @@@ pairs;
previous = Infinity;
Do[
  values = (scaledForce[#[[1]], #[[2]], dist, 32, 10^-7] &) /@ pairs;
  control = (scaledForce[#[[1]], #[[2]], dist, 24, 10^-6] &) /@ pairs;
  error = Max[Abs[values/limit - 1]];
  Print["held unit bend, d/a = ", dist, ", d^2 Fodd/(4 Pi kappa) = ",
    InputForm[N[values, 14]]];
  If[! TrueQ[Max[Abs[values - control]] < 10^-10 && error < previous &&
      Min[values] > 0], Exit[1]];
  previous = error,
  {dist, {20, 40, 80}}];
If[! TrueQ[previous < 3/1000], Exit[1]];
Print["largest relative error from {1,2,4} at d/a=80 = ", InputForm[N[previous, 10]]];
Exit[0];
