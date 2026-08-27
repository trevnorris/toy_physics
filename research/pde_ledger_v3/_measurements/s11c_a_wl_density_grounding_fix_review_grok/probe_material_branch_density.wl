(* Probe: does FACE_SHIFT density background use the §2c branched anchor?
   Engine FACE_SHIFT path: rho4Profile[density, spatialCoordinates] always.
   Engine elsewhere: branchedRho4Source uses inverseCoordinates for MATERIAL_ADVECTED.
   Compare SHAPE of density background alone for MATERIAL + RHOBR. *)

$HistoryLength = 0;
ClearAll["Global`*"];

spatialCoordinates = {x1, x2, x3};
backgroundValue[expression_] := expression /. waveOrder -> 0;
shapeDerivative[expression_] := Together[Expand[
  D[expression, waveOrder] /. waveOrder -> 0]];

widthProfile[a_, b_, c_] := W0 (1 + etaBg w1[a, b, c]);

inverseCoordinates[] := spatialCoordinates - waveOrder Through[
  {uOne, uTwo, uThree}[Sequence @@ Append[spatialCoordinates, time]]];

(* As FACE_SHIFT currently wires it: always spatial. *)
faceShiftBgMATERIALRHOBR = rhoBr/widthProfile @@ spatialCoordinates;

(* As branchedRho4Source wires it for MATERIAL_ADVECTED. *)
branchedBgMATERIALRHOBR = rhoBr/widthProfile @@ inverseCoordinates[];

Print["=== MATERIAL_ADVECTED RHOBR density background ==="];
Print["FACE_SHIFT path (spatial): ", faceShiftBgMATERIALRHOBR];
Print["branchedRho4 path (chi):   ", branchedBgMATERIALRHOBR];
Print["SHAPE FACE_SHIFT path: ", shapeDerivative[faceShiftBgMATERIALRHOBR]];
Print["SHAPE branched path:   ", shapeDerivative[branchedBgMATERIALRHOBR]];
Print["SHAPE residual (branched - FACE_SHIFT path) = ",
  Simplify[shapeDerivative[branchedBgMATERIALRHOBR] -
    shapeDerivative[faceShiftBgMATERIALRHOBR]]];
Print["material_transport_missing_from_FACE_SHIFT_path? ",
  Simplify[shapeDerivative[branchedBgMATERIALRHOBR] -
      shapeDerivative[faceShiftBgMATERIALRHOBR]] =!= 0];
