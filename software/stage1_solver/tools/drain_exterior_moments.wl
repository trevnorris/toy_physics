(* Independent 3D potential calculation. Represent the solution by a layer of
   sinks on a smaller concentric sphere. Fit its weights to the outer sphere's
   specified normal speed, then integrate the field at observation points.
   Unit length a=1 and reference integrated intake Qref=4 Pi. *)
ClearAll["Global`*"];
layerRadius = 1/2;
(* In layer coordinates with polar axis along the observing ray, z is the
   cosine of the angle to that ray. v is the cosine to the fixed intake axis.
   Averaging the angular weight first makes the remaining integral smooth. *)
v = u z + Sqrt[1-u^2] Sqrt[1-z^2] Cos[phi];
weightsOnLayer = {1, v, (3 v^2 - 1)/2};
angleWeights = FullSimplify[Integrate[weightsOnLayer, {phi, 0, 2 Pi}],
                            -1 <= u <= 1 && -1 <= z <= 1];
radialKernel = (layerRadius z - rad)/(4 Pi (rad^2 + layerRadius^2 - 2 rad layerRadius z)^(3/2));
responses[r0_, u0_] := responses[r0, u0] = Table[
  NIntegrate[Evaluate[(radialKernel angleWeights[[mode]]) /. {rad -> r0, u -> u0}],
             {z, -1, 1}, WorkingPrecision -> 35, AccuracyGoal -> 23,
             PrecisionGoal -> 20], {mode, 3}];
sampleAngles = {-1, 0, 1};
boundaryMatrix = responses[1, #] & /@ sampleAngles;
targetRate[angle_, p_, s_, rate_] := -rate (1 + p angle + s (3 angle^2 - 1)/2);
Do[
  {name, p, s, rate} = entry;
  weights = LinearSolve[boundaryMatrix, targetRate[#, p, s, rate] & /@ sampleAngles];
  testAngle = 1/3;
  error = responses[1, testAngle].weights - targetRate[testAngle, p, s, rate];
  Print[name, " [layer monopole, off-grid BC residual] = ",
        ToString[NumberForm[{First[weights], error}, {18,12}], OutputForm]];
  Do[
    vals = -distance^2 (responses[distance, #].weights & /@ {1, 0});
    Print[name, " d/a=", distance, " [scaled pole, equator] = ",
          ToString[NumberForm[vals, {18,12}], OutputForm]],
    {distance, {2,5,10}}];
  If[Abs[error] > 10^-14, Exit[1]],
  {entry, {{"round", 0, 0, 1}, {"dipole", 1/2, 0, 1},
           {"quadrupole", 0, 1/2, 1}, {"double", 0, 0, 2}}}];
Exit[0];
