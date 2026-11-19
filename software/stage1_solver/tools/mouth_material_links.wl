(* Round fixed boundary, geodesic triangle edges, ordinary exterior Maxwell
   cost gamma Integral[B.B]. A directed link is measured in the local site
   frames; both preferred and physical links transform by the same endpoint
   phases. The attachment energy is J Sum[1-Cos[relative edge angle]]. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label, " = ", ToString[value, InputForm]];
checks = {};
verify[condition_] := AppendTo[checks, TrueQ[condition]];

pauli = Table[PauliMatrix[k], {k, 3}];
project[v_] := (IdentityMatrix[2] + Sum[v[[k]] pauli[[k]], {k, 3}])/2;
{a, b, c} = Table[Array[var, 3], {var, {aa, bb, cc}}];
trace = Expand[Tr[project[a].project[b].project[c]]];
invariants = (1 + a.b + b.c + c.a + I Det[{a, b, c}])/4;
residual = Expand[trace - invariants];
emit["triangle projector trace, invariant expression, residual", {trace, invariants, residual}];
verify[residual === 0];
(* Tr[P1 P2 P3] equals <1|2><2|3><3|1>; normalization does not
   change its argument when the three overlaps are nonzero. *)

radial = {Sin[t] Cos[p], Sin[t] Sin[p], Cos[t]};
states = {{Cos[t/2], Exp[I p] Sin[t/2]}, {-Exp[-I p] Sin[t/2], Cos[t/2]}};
realSimplify[expr_] := FullSimplify[expr, Assumptions -> Element[{t, p}, Reals]];
norms = realSimplify[Conjugate[#].# & /@ states];
projectorDifferences = MapThread[
  realSimplify[Outer[Times, #1, Conjugate[#1]] - project[#2 radial]] &,
  {states, {1, -1}}];
connections = Table[realSimplify[-I Conjugate[z].D[z, u]], {z, states}, {u, {t, p}}];
curvatures = realSimplify[D[#[[2]], t] - D[#[[1]], p]] & /@ connections;
degrees = Integrate[#, {t, 0, Pi}, {p, 0, 2 Pi}]/(2 Pi) & /@ curvatures;
radialFields = realSimplify[curvatures/(r^2 Sin[t])];
exterior = FullSimplify[
  Integrate[4 Pi gamma r^2 #^2, {r, radius, Infinity}, Assumptions -> radius > 0],
  Assumptions -> radius > 0] & /@ radialFields;
emit["spinor norms, projector differences", {norms, projectorDifferences}];
emit["site connection {theta,phi}, curvature, integral/(2Pi), radial B, exterior energy",
  {connections, curvatures, degrees, radialFields, exterior}];
verify[norms == {1, 1} && projectorDifferences == ConstantArray[0, {2, 2, 2}]];
verify[degrees == {1, -1} && exterior == {Pi gamma/radius, Pi gamma/radius}];

(* No complex-linear unitary conjugation can reverse the sign of the trace
   of the Pauli triple. An antiunitary on the same two states can implement
   reversal and conjugates all overlaps. *)
anti = I pauli[[2]];
triple = {Tr[Dot @@ pauli], Tr[Dot @@ (-pauli)]};
antiResidual = (anti.Conjugate[#].ConjugateTranspose[anti] + #) & /@ pauli;
emit["trace before and after reversal, antiunitary reversal residual", {triple, antiResidual}];
verify[triple == {2 I, -2 I} && antiResidual == ConstantArray[0, {3, 2, 2}]];

(* Refine all oriented octants on the sphere. The flux comparator uses
   spherical excess, computed only from the three tangent-plane angles. *)
octants = (If[Det[#] > 0, #, Reverse[#]] &) /@
  (DiagonalMatrix /@ Tuples[{-1, 1}, 3]);
subdivide[{v_, w_, u_}] := Module[
  {vw = Normalize[v + w], wu = Normalize[w + u], uv = Normalize[u + v]},
  {{v, vw, uv}, {vw, w, wu}, {uv, wu, u}, {vw, wu, uv}}
];
meshes = NestList[Flatten[subdivide /@ #, 1] &, N[octants], 2];
interior[v_, w_, u_] := VectorAngle[w - (v.w) v, u - (v.u) v];
area[vertices_] := Total[Apply[interior, RotateLeft[vertices, #]] & /@ Range[0, 2]] - Pi;
(* Choose one of two local charts to keep denominators away from zero. *)
section[{x_, y_, z_}] := If[z >= 0,
  {1 + z, x + I y}/Sqrt[2 (1 + z)], {x - I y, 1 - z}/Sqrt[2 (1 - z)]];
link[z_, w_] := With[{overlap = Conjugate[z].w}, overlap/Abs[overlap]];
links[z_] := MapThread[link, {z, RotateLeft[z]}];
localFrame[v_] := Exp[I (4 v[[1]] - 9 v[[2]]^2 + 5 v[[1]] v[[3]])];
phases[polygon_, side_, rad_] := Module[{xyz, normals, z, u, ug, uc, un, flat},
  xyz = rad polygon;
  normals = side (Normalize /@ xyz);
  z = section /@ normals;
  u = links[z];
  ug = links[MapThread[Times, {localFrame /@ xyz, z}]];
  uc = links[Conjugate[z]];
  un = links[section /@ (-normals)];
  flat = links[ConstantArray[section[{0, 0, 1}], 3]];
  Arg[Times @@ #] & /@ {u, ug, uc, un, u^2, u uc, flat, ConstantArray[1, 3]}
];
tol = 10^-9;
meshRows = Flatten[Table[Module[
    {faces = meshes[[level + 1]], balls, fluxes, numbers, worstArea, worstGauge, cutoff},
    balls = area /@ faces;
    fluxes = phases[#, side, rad] & /@ faces;
    numbers = Total[fluxes]/(2 Pi);
    worstArea = Max[Abs[fluxes[[All, 1]] - side balls/2]];
    worstGauge = Max[Abs[fluxes[[All, 2]] - fluxes[[All, 1]]]];
    cutoff = Max[Abs[fluxes[[All, 1]]]];
    verify[Max[Abs[numbers - {side, side, -side, -side, 2 side, 0, 0, 0}]] < tol];
    verify[worstArea < tol && worstGauge < tol && cutoff < Pi];
    {Length[faces], rad, side, N[Total[balls]], Chop[numbers],
      {worstArea, worstGauge}, cutoff}],
    {level, 0, 2}, {rad, {1, 3}}, {side, {1, -1}}], 2];
emit["faces, radius, side, area, counts {native,gauged,conjugate,other line,square,opposite tensor channels,uniform,unlocked flat}, max residuals {half area,gauge}, max face angle",
  meshRows];

(* For every link U let delta=Arg[U/U*] at the candidate with a different
   integer. Each triangle angle changes by its oriented sum of deltas,
   modulo 2 Pi. All unreduced triangle sums of delta cancel on the closed
   triangulation. Some face must therefore reach/pass the principal cut;
   its three endpoint deviations have absolute sum at least Pi-max|phi*|.
   Positivity on every edge gives a global lower bound from that one edge.
   At the preferred links the exact radial exterior realizes every face's
   spherical excess flux. The flat trial A=0 costs at most 2 J per edge;
   any nonzero integer costs at least Pi gamma/R by spherical Cauchy. *)
certificate[rows_] := Module[{nface, rad, maxAngle, gapFactor, nedge, unitCost},
  {nface, rad, maxAngle} = rows[[{1, 2, 7}]];
  nedge = 3 nface/2;
  gapFactor = 1 - Cos[(Pi - maxAngle)/3];
  unitCost = Pi/rad; (* gamma=1 in common length and energy units *)
  {nface, rad, gapFactor, 16 gapFactor, N[unitCost], N[2 nedge/512],
    16 gapFactor > unitCost, 2 nedge/512 < unitCost}
];
bounds = certificate /@ Select[meshRows, #[[3]] == 1 &];
emit["faces, radius, wrong-sector attachment/J lower bound, wrong-sector lower (J=16), locked trial (gamma=1), flat upper (J=1/512), strong and weak comparisons", bounds];
verify[And @@ Flatten[bounds[[All, {7, 8}]]]];
emit["checks", checks];
If[!TrueQ[And @@ checks], Exit[1]];
