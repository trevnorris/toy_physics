(* Independent volume and surface checks on a finite shell. Set the radial
   current from Gauss's law on the whole ball (reference) or around the
   prescribed additional point intake. Sample local 4D continuity away
   from the origin and integrate both faces and both spherical cuts. *)
ClearAll["Global`*"];
halfDepth = 1;
numberDensity = 3 + (3/100) time;
verticalCurrent[w_] := -(3/25) w;
coordinates = {x, y, z};
radialNorm = Sqrt[coordinates.coordinates];
ballVolume[rad_] := 4 Pi rad^3/3;
sphereArea[rad_] := 4 Pi rad^2;
inflowPerVolume = verticalCurrent[-halfDepth] - verticalCurrent[halfDepth];
storedPerVolume = 2 halfDepth D[numberDensity, time];

vectorCurrent[extra_] := With[
  {flux = ballVolume[radialNorm] (inflowPerVolume - storedPerVolume) - extra},
  coordinates/radialNorm * flux/(2 halfDepth sphereArea[radialNorm])
];

radialFlux[extra_, radius_] := Module[{current, normal, rules},
  current = vectorCurrent[extra];
  normal = {Sin[theta] Cos[phi], Sin[theta] Sin[phi], Cos[theta]};
  rules = Thread[coordinates -> radius normal];
  Integrate[FullSimplify[(current.normal /. rules) radius^2 Sin[theta],
              0 < theta < Pi && 0 < phi < 2 Pi],
            {theta, 0, Pi}, {phi, 0, 2 Pi}, {w, -halfDepth, halfDepth}]
];
shellBudget[extra_, radius_] := Module[{core, faces, storage},
  core = -radialFlux[extra, 1];
  faces = Integrate[4 Pi s^2 (verticalCurrent[halfDepth] - verticalCurrent[-halfDepth]),
                    {s, 1, radius}];
  storage = Integrate[4 Pi s^2 storedPerVolume, {s, 1, radius}];
  {core, faces, storage}
];

Do[
  residual = FullSimplify[D[numberDensity, time] + Div[vectorCurrent[extra], coordinates]
                           + D[verticalCurrent[w], w], radialNorm > 0];
  Print["extra=", InputForm[extra], " local continuity residual=", InputForm[residual]];
  If[residual =!= 0, Exit[1]],
  {extra, {0, 24 Pi}}
];
Do[
  fluxes = radialFlux[#, radius] & /@ {0, 24 Pi};
  budgets = shellBudget[#, radius] & /@ {0, 24 Pi};
  speeds = fluxes/(sphereArea[radius] 2 halfDepth (numberDensity /. time -> 0));
  Print["R=", radius, " [background speed, throat speed, difference, extra inward flux/pi]=",
        InputForm[{First[speeds], Last[speeds], Last[speeds] - First[speeds],
                   (First[fluxes] - Last[fluxes])/Pi}]];
  Print["R=", radius, " budgets/pi [core, faces, storage]=", InputForm[budgets/Pi],
        " outward fluxes/pi=", InputForm[fluxes/Pi]];
  If[Simplify[fluxes + (Total /@ budgets)] =!= {0, 0}, Exit[1]],
  {radius, {1, 5, 10}}
];
Exit[0];
