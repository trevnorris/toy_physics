(* Two complex two-state coordinates in a common internal SU(2) frame.
   Exterior assumptions: compact, independent face connections, identical
   positive kinetic and curl coefficients, and London stiffness only for the
   composites specified in each row. D = d - i q.A. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label <> " = " <> ToString[value, InputForm]];
checks = {};
verify[condition_] := AppendTo[checks, TrueQ[condition]];
eps = {{0, 1}, {-1, 0}};
singlet[x_, y_] := x.eps.y;
overlap[x_, y_] := Conjugate[x].y;
pair[x_, y_] := {singlet[x, y], overlap[x, y]};
top = {zt1, zt2};
bottom = {zb1, zb2};
base = pair[top, bottom];
swap = pair[bottom, top];
swapTarget = {-base[[1]], Conjugate[base[[2]]]};
emit["pair, exchange, exchange residual", {base, swap, Expand[swap - swapTarget]}];
verify[Expand[swap - swapTarget] == {0, 0}];

(* A generic common proper rotation; normalization det(u)=1 gives SU(2).
   These identities hold before imposing normalization. *)
u = {{aa, bb}, {-Conjugate[bb], Conjugate[aa]}};
suResidual = Expand[{Transpose[u].eps.u - Det[u] eps,
                     ConjugateTranspose[u].u - Det[u] IdentityMatrix[2]}];
emit["rotation determinant, rotated invariant tensors, residual",
  {Det[u], {Transpose[u].eps.u, ConjugateTranspose[u].u}, suResidual}];
verify[suResidual == ConstantArray[0, {2, 2, 2}]];

(* The SU(2) action is transitive on normalized rays; point the first
   north, keep the full relative polar angle and the second azimuth. *)
relative = {Cos[angle/2], Exp[I azimuth] Sin[angle/2]};
values = pair[{1, 0}, relative];
realAngles = Element[{angle, azimuth, xt, xb}, Reals];
strengths = FullSimplify[Conjugate[values] values, Assumptions -> realAngles];
angleTargets = {(1 - Cos[angle])/2, (1 + Cos[angle])/2};
angleResidual = FullSimplify[strengths - angleTargets, Assumptions -> realAngles];
endpoint = {strengths /. angle -> 0, strengths /. angle -> Pi};
emit["modulus squared, angle residual, parallel and antiparallel preparations",
  {strengths, angleResidual, endpoint}];
verify[angleResidual == {0, 0} && endpoint == {{0, 1}, {1, 0}}];

transformed = pair[Exp[I xt] top, Exp[I xb] bottom];
frameFactors = FullSimplify[transformed/base, Assumptions -> realAngles];
weights = Table[FullSimplify[D[factor, x]/(I factor), Assumptions -> realAngles],
  {factor, frameFactors}, {x, {xt, xb}}];
emit["endpoint frame factors and charge rows {top,bottom}", {frameFactors, weights}];
verify[weights == {{1, 1}, {-1, 1}}];
(* The directed phase link of the composite has the endpoint frame change
   q.(alpha_j-alpha_i), identical to q.(A_t,A_b) for directed links. *)
linkAngles = {at, ab};
deltas = {dxt, dxb};
covariant = Expand[weights.(linkAngles + deltas) - weights.deltas];
emit["transformed link with frame change removed, preferred combinations, residual",
  {covariant, weights.linkAngles, Expand[covariant - weights.linkAngles]}];
verify[Expand[covariant - weights.linkAngles] == {0, 0}];

(* Derive the two transverse inverse-propagator eigenvalues from the
   quadratic form. Phase gradients are longitudinal on a uniform patch. *)
potential = {at, ab};
elastic = Sum[{rv, rw}[[j]] (weights[[j]].potential)^2/2, {j, 2}];
hessian = Table[D[elastic, x, y], {x, potential}, {y, potential}];
basis = {{1, 1}, {1, -1}}/Sqrt[2];
normal = Simplify[basis.hessian.Transpose[basis]];
lambda = Diagonal[normal];
dispersion = (kb k^2 + lambda)/ke;
emit["London Hessian, even/odd Hessian, omega^2, inverse length squared",
  {hessian, normal, dispersion, lambda/kb}];
verify[normal == DiagonalMatrix[{2 rv, 2 rw}]];

(* A nowhere-zero section trivializes its line on S^2. A transition with
   integer winding (nt,nb) induces the following composite windings. *)
transition = frameFactors /. {xt -> nt p, xb -> nb p};
windings = FullSimplify[D[#, p]/(I #)] & /@ transition;
emit["composite transitions, windings, opposite-face preparation",
  {transition, windings, windings /. {nt -> side, nb -> -side}}];
verify[Expand[windings - weights.{nt, nb}] == {0, 0}];
verifiedMirror = Simplify[(windings /. {nt -> side, nb -> -side}) - {0, -2 side}];
verify[verifiedMirror == {0, 0}];

(* Nullspace is over the Lie algebra. Use primitive integer generators,
   plus finite stabilizer information from the compact angle constraints.
   A top-positive generator sets the common convention for probe charges. *)
primitive[v_] := Sign[First[Select[v, # != 0 &]]] v/(GCD @@ Abs[v]);
probes = {{1, 0}, {0, 1}, {1, 1}};
rows = Table[Module[{active, qs, generators, suscept, gaps, probeCharges},
    active = Flatten[Position[mask, 1]];
    qs = If[active == {}, {{0, 0}}, weights[[active]]];
    generators = primitive /@ NullSpace[qs];
    gaps = lambda /. {rv -> mask[[1]] rv, rw -> mask[[2]] rw};
    suscept = If[generators === {}, ConstantArray[0, {2, 2}],
      Transpose[generators].Inverse[generators.Transpose[generators]].generators];
    probeCharges = (probes.# &) /@ generators;
    {mask, qs.{nt, nb}, generators, gaps/ke, suscept,
      probeCharges, (ke #.# &) /@ generators}],
  {mask, {{0, 0}, {1, 0}, {0, 1}, {1, 1}}}];
emit["stiffness switches {v,w}, required zero windings, primitive generators, gap^2 even/odd, massless projector, signed probe charges {top,bottom,equal} per generator, unnormalized kinetic coefficients",
  rows];
verify[rows[[2, 3]] == {{1, -1}} && rows[[3, 3]] == {{1, 1}} && rows[[4, 3]] == {}];
verify[rows[[2, 6]] == {{1, -1, 0}} && rows[[3, 6]] == {{1, 1, 2}}];
verify[rows[[2, 5]] == {{1, -1}, {-1, 1}}/2 && rows[[2, 7]] == {2 ke}];
(* xt,xb in turns in [0,1). Both condensates require integral xt+xb and
   xb-xt; their integer range is exhausted by the following finite set. *)
compact = Flatten[Table[
  Solve[{xt + xb == j, xb - xt == h, 0 <= xt < 1, 0 <= xb < 1},
        {xt, xb}, Reals], {j, 0, 1}, {h, -1, 1}], 2];
compactPoints = {xt, xb} /. compact;
emit["determinant and surviving transformations in turns with both stiff",
  {Det[weights], compactPoints}];
verify[Sort[compactPoints] == {{0, 0}, {1/2, 1/2}}];

(* If an ordinary static source separately obeys a screened Poisson law,
   its radial gradient flux need not be constant. It differs from the
   first Chern flux of a connection on concentric spheres with no
   intervening magnetic current. The London transverse calculation alone
   does not specify a scalar-potential response. *)
green = Exp[-mu r]/(4 Pi r);
greenResidual = Simplify[-D[r^2 D[green, r], r]/r^2 + mu^2 green];
gradientFlux = Simplify[-4 Pi r^2 D[green, r]];
fluxLimits = {Limit[gradientFlux, r -> 0, Direction -> "FromAbove"],
  Limit[gradientFlux, r -> Infinity, Assumptions -> mu > 0]};
emit["static Helmholtz example {residual, gradient flux, limits}, closed-sector Chern flux",
  {{greenResidual, gradientFlux, fluxLimits}, 2 Pi windings}];
verify[greenResidual == 0 && fluxLimits == {1, 0}];
emit["checks", checks];
If[!TrueQ[And @@ checks], Exit[1]];
