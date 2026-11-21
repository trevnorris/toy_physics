(* Collective momentum of one thin inertial order front per sidewall area.
   q is the logistic order fraction. Its normal derivative is q(1-q)/delta.
   x=cos(theta) projects a local normal displacement onto the bend axis. *)
ClearAll["Global`*"];
$HistoryLength = 0;
$Assumptions = inertia > 0 && kappa > 0 && delta > 0 && half > 0 && radius > 0;
emit[label_, actual_, expected_] := Module[{value = FullSimplify[actual]},
  Print[label <> " = " <> ToString[value, InputForm]];
  If[FullSimplify[value == expected] =!= True, Exit[1]]];

slope = q (1-q)/delta;
surfaceNorm = Integrate[slope, {q, 0, 1}];
emit["squared normal gradient integrated through one front", surfaceNorm, 1/(6 delta)];
tension = kappa surfaceNorm;
massArea = inertia surfaceNorm;
emit["front mass per area divided by tension", massArea/tension, inertia/kappa];
bend = Sin[Pi u]^2;
normalVelocity = peakSpeed bend x;
kinetic = 2 half radius^2 (2 Pi) Integrate[massArea normalVelocity^2/2, {x,-1,1}, {u,0,1}];
labMomentum = 2 half radius^2 (2 Pi) Integrate[massArea normalVelocity x, {x,-1,1}, {u,0,1}];
radial = 8 Pi radius^2 half massArea;
emit["total sideways momentum per peak speed and radial inertia", labMomentum/(peakSpeed radial), 1/6];
emit["bend-conjugate momentum per peak speed and radial inertia",
     D[kinetic, peakSpeed]/(peakSpeed radial), 1/8];
(* Verify the positive-excess integral controlling the optimal contrast
   using its antiderivative evaluated on the central interval. *)
primitive = 3 u/8 - Sin[2 Pi u]/(4 Pi) + Sin[4 Pi u]/(32 Pi);
emit["primitive differentiated minus squared bend shape", D[primitive,u]-bend^2, 0];
Exit[0];
