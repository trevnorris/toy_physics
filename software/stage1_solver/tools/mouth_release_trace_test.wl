(* Work on the two sides of a phase-continuous carrier boundary, using the
   archived scalar density traces as common inputs with the virtual-work test. *)
ClearAll["Global`*"];
energy = k w[y] (n - n0[y])^2/2;
(* Variation of the two bulk gradient energies and of the real surface
   energy with respect to the common amplitude z = Sqrt[n]. *)
surfaceDerivative = D[energy /. n -> z^2, z];
amplitudeJump = jump /. First[Solve[-2 c jump + surfaceDerivative == 0, jump]];
logSlopeJump = Simplify[(2 amplitudeJump/z) /. z -> Sqrt[n], n > 0];
(* With the common normal into the drain, mixed momentum flux includes
   -(c/2) n d_x d_N log n. Subtract the outgoing flux from the incoming. *)
trace = n0[y] + u n1[y];
stressDifference = (c/2) trace D[logSlopeJump /. n -> trace, y];
driveSlope = Simplify[D[stressDifference, u] /. u -> 0];
boundaryDerivative = D[k n0[y] w[y] n1[y], y];
interiorWork = Simplify[driveSlope - boundaryDerivative];
Print["resting energy, extra normal log slope, extra mixed flux: ",
  InputForm[{energy /. n -> n0[y], logSlopeJump /. n -> n0[y], stressDifference /. u -> 0}]];
Print["linear flux before and after removing the compact boundary divergence: ",
  InputForm[{driveSlope, interiorWork}]];
If[!TrueQ[Simplify[interiorWork + k w[y] n1[y] n0'[y]] == 0], Exit[1]];

folder = Environment["MOUTH_RELEASE_FACE_DIR"];
count = 32; dx = 1/8;
readFace[name_] := Module[{path, values},
  path = FileNameJoin[{folder, name <> ".f64"}];
  If[!FileExistsQ[path] || FileByteCount[path] != 8 count^3, Exit[2]];
  values = BinaryReadList[path, "Real64", ByteOrdering -> -1];
  ArrayReshape[values, {count, count, count}]
];
rest = readFace["rest"]; tangent = readFace["tangent"];
coords = N[(Range[count] - (count+1)/2) dx];
weight = Table[Max[0., 1 - 4 (x^2+y^2+z^2)]^3, {x, coords}, {y, coords}, {z, coords}];
grad[a_] := (RotateRight[a, 2] - 8 RotateRight[a, 1] +
              8 RotateLeft[a, 1] - RotateLeft[a, 2])/(12 dx);
numFlux[n_, en_] := dx^3 Total[Flatten[n grad[en]]];
(* c and the slope-jump factor cancel in (c/2) n grad[2 e_n/c].
   Take k=1 to report the coefficient of this otherwise unselected term. *)
restFlux = numFlux[rest, weight (rest - rest)];
slopeFlux = numFlux[rest, weight tangent];
Print["rest and first-order net face difference per strength: ",
  InputForm[{restFlux, slopeFlux}]];
printPair[name_, speed_] := Module[{pos, neg, currents},
  pos = readFace[name]; neg = Reverse[pos];
  currents = {numFlux[pos, weight (pos-rest)], numFlux[neg, weight (neg-rest)]};
  Print["signed face differences and odd slope at |U|=", speed, ": ",
    InputForm[Append[currents, (First[currents]-Last[currents])/(2 speed)]]];
];
printPair["flow0.02", 0.02];
printPair["flow0.01", 0.01];
Exit[0];
