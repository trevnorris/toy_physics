(* Differentiate along a branch with both integrated loads free to respond.
   Reconstruct the change of each local stored energy from its port work. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, value_] := Print[label, " = ", ToString[value, InputForm]];
check[label_, value_, expected_] := Module[{result = Factor[value]},
  emit[label, result]; If[Factor[result - expected] =!= 0, Exit[1]]];
g = {{self1, x z1 z2}, {x z1 z2, self2}};
y = {p1[x], p2[x]};
h = g.y;
stored = D[y.h/2, x];
local = -y.D[h, x];
check["stationary local work plus field storage",
  (stored + local) /. x -> 0, -z1 z2 p1[0] p2[0]];
heldLoads = {p1'[0] -> 0, p2'[0] -> 0};
emit["field and attachment slopes with both loads fixed",
  Factor /@ ({stored, local} /. x -> 0 /. heldLoads)];

(* A fixed-reference spring with force k (ref - h) stores y^2/(2 k).
   Solve for y at finite separation, then differentiate the full energy. *)
refs = {s1 a1, s2 a2};
compliance = DiagonalMatrix[{1/k1, 1/k2}];
springLoads = LinearSolve[g + compliance, refs];
total = springLoads.(g + compliance).springLoads/2;
springCoefficient = Factor[(D[total, x] /. x -> 0)/(s1 s2)];
check["two springs after elimination", springCoefficient,
  -z1 z2 k1 a1 k2 a2/((1 + k1 self1) (1 + k2 self2))];
check["stiff spring at both mouths", Limit[Limit[springCoefficient, k1 -> Infinity],
  k2 -> Infinity], -z1 z2 a1 a2/(self1 self2)];

(* Differentiate the straight segment energy with respect to its ends.
   The cylinder has two curved directions and one straight direction. *)
wallEnergy = 4 Pi r^2 length (sigma + bend (2/(3 r))^2);
axial = D[wallEnergy, length];
work = axial (rimGain - farGain) ell;
check["wall work vanishes for equal gains", work /. farGain -> rimGain, 0];
check["change of pull upon doubling radius", (axial /. r -> 2 r) - axial,
  12 Pi r^2 sigma];

archive = Import["docs/_measurements/unequal_mouth_electric.json", "RawJSON"];
lines = archive["runs"]["annular_mouth_radius_wolfram"]["stdout_claim_lines"];
(* Use the previous independent radial solve. Its first two lines store
   the escape and self response at radii one and two. *)
escapes = Table[Module[{line, match},
  line = SelectFirst[lines, StringStartsQ[#, "\"radius \"" <> ToString[i] <> ","] ||
    StringStartsQ[#, "\"radius \"" <> ToString[i] <> "\""] &];
  If[! StringQ[line], Exit[1]];
  match = StringCases[line, RegularExpression["InputForm\\[\\{([0-9]+\\.[0-9]+)"] -> "$1"];
  If[Length[match] != 1, Exit[1]]; ToExpression[First[match]]], {i, 2}];
If[! TrueQ[0 < First[escapes] < Last[escapes] < 1], Exit[1]];
required = Factor[(axial/(axial /. r -> 2 r)) First[escapes]/Last[escapes]];
emit["relative transfer, curvature and tension limits",
  N[{required /. sigma -> 0, required /. bend -> 0}, 12]];
emit["large/small strength for equal transfer in those limits",
  N[{1/required /. sigma -> 0, 1/required /. bend -> 0}, 12]];
Exit[0];
