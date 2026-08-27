(* Focused FORM ablation of S11c-a traceSource density background.
   Copy of the post-fix trace-construction logic only — not a full engine run.
   BEFORE: §2b RHO4_CONSTANT representative (w-independent).
   AFTER:  replace background with an explicitly w-dependent test function.
   Report shapeDerivative of the density trace (LAB_HELD, FACE_MINUS, DELTA_W). *)

$HistoryLength = 0;
ClearAll["Global`*"];

spatialCoordinates = {x1, x2, x3};
waveOrder = waveOrder;
normalCoordinate = normalCoordinate;
time = time;
etaBg = etaBg;
W0 = W0;
rhoBr = rhoBr;
dofWidth = dofWidth;
dofCenter = dofCenter;

backgroundValue[expression_] := expression /. waveOrder -> 0;
shapeDerivative[expression_] := Together[Expand[
  D[expression, waveOrder] /. waveOrder -> 0]];

widthProfile[x1, x2, x3] := W0 (1 + etaBg w1Jet[0, 0, 0]);

(* Minimal LAB_HELD height (sign = -1). *)
heightLabMinus[] := Module[{inverse},
  inverse = spatialCoordinates - waveOrder Through[
    {uOne, uTwo, uThree}[Sequence @@ Append[spatialCoordinates, time]]];
  -widthProfile[x1, x2, x3]/2 +
    waveOrder (dofCenter zetaCenter @@ Append[inverse, time] +
      (-1) dofWidth (deltaWidth @@ Append[inverse, time])/2)
];

(* Post-fix traceSource (same structure as the engine). *)
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

rho4ConstantZero[coordinates_List, normal_] := rhoBr/W0;
rhoBulkWave[coordinates_List, normal_] :=
  rhoBulkPerturbation @@ Append[coordinates, {normal, time}];

(* FORM change: explicitly w-dependent background (not a coefficient rescale). *)
wDependentZero[coordinates_List, normal_] :=
  rhoBr/W0 + alphaTest * normal^2;

applyDeltaW[expression_] := expression /. {dofWidth -> 1, dofCenter -> 0};

height = heightLabMinus[];

beforeExact = traceSource[rho4ConstantZero, rhoBulkWave, height];
beforeShape = applyDeltaW[shapeDerivative[beforeExact]];

afterExact = traceSource[wDependentZero, rhoBulkWave, height];
afterShape = applyDeltaW[shapeDerivative[afterExact]];

Print["=== BEFORE (grounded RHO4_CONSTANT; w-independent) ==="];
Print["backgroundField = ", rho4ConstantZero[spatialCoordinates, normalCoordinate]];
Print["D[bg,w] = ", D[rho4ConstantZero[spatialCoordinates, normalCoordinate], normalCoordinate]];
Print["SHAPE_DERIVATIVE = ", beforeShape];

Print["=== AFTER (FORM ablation: bg -> rhoBr/W0 + alphaTest*w^2) ==="];
Print["backgroundField = ", wDependentZero[spatialCoordinates, normalCoordinate]];
Print["D[bg,w] = ", D[wDependentZero[spatialCoordinates, normalCoordinate], normalCoordinate]];
Print["SHAPE_DERIVATIVE = ", afterShape];

Print["=== residual (AFTER - BEFORE) — shift term must reappear ==="];
Print["AFTER - BEFORE = ", Together[Expand[afterShape - beforeShape]]];
Print["shift_reappeared? ", Simplify[afterShape - beforeShape] =!= 0];
