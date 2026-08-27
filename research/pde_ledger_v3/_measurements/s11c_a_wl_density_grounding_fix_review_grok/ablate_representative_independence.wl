(* Independence: grounded FACE_SHIFT density background is not a hand-written 0.
   Changing which §2b representative is used must change the emitted background
   face value. Uses the same post-fix traceSource structure. *)

$HistoryLength = 0;
ClearAll["Global`*"];

spatialCoordinates = {x1, x2, x3};
backgroundValue[expression_] := expression /. waveOrder -> 0;
shapeDerivative[expression_] := Together[Expand[
  D[expression, waveOrder] /. waveOrder -> 0]];

widthProfile[a_, b_, c_] := W0 (1 + etaBg w1Jet[0, 0, 0]);

heightLabMinus[] :=
  -widthProfile[x1, x2, x3]/2 +
    waveOrder (dofCenter zetaCenter[x1, x2, x3, time] -
      dofWidth deltaWidth[x1, x2, x3, time]/2);

traceSource[fieldZero_, fieldWave_, height_] := Module[
  {backgroundHeight, backgroundField, backgroundFace,
    backgroundNormalDerivative, backgroundAnsatz},
  backgroundHeight = backgroundValue[height];
  backgroundField = fieldZero[spatialCoordinates, normalCoordinate];
  backgroundFace = backgroundField /.
    normalCoordinate -> backgroundHeight;
  backgroundNormalDerivative = D[backgroundField, normalCoordinate] /.
    normalCoordinate -> backgroundHeight;
  backgroundAnsatz = backgroundFace +
    (normalCoordinate - backgroundHeight) backgroundNormalDerivative;
  (backgroundAnsatz +
      waveOrder fieldWave[spatialCoordinates, normalCoordinate]) /.
    normalCoordinate -> height
];

rho4Zero[coordinates_List, normal_] := rhoBr/W0;
rhoBrZero[coordinates_List, normal_] := rhoBr/widthProfile @@ coordinates;
rhoWave[coordinates_List, normal_] :=
  rhoBulkPerturbation @@ Append[coordinates, {normal, time}];

height = heightLabMinus[];
exactR4 = traceSource[rho4Zero, rhoWave, height];
exactRB = traceSource[rhoBrZero, rhoWave, height];

bgFace[exact_] := backgroundValue[exact];

Print["=== EXACT_TRACE background face (waveOrder -> 0) ==="];
Print["RHO4_CONSTANT face = ", bgFace[exactR4]];
Print["RHOBR_CONSTANT face = ", bgFace[exactRB]];
Print["faces_differ? ", Simplify[bgFace[exactR4] - bgFace[exactRB]] =!= 0];
Print["RHO4 - RHOBR = ", Simplify[bgFace[exactR4] - bgFace[exactRB]]];
