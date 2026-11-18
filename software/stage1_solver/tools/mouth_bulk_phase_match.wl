(* Fourier and radial checks of the ideal half-space receiver, m=hbar=nBulk=1
   in the numerical example. +w points from the feeding slab into the bulk. *)
ClearAll["Global`*"];
$HistoryLength = 0;
mode = coeff Exp[-k z];
amplitude = coeff /. First[Solve[(D[mode, z] /. z -> 0) == j/nBulk, coeff]];
(* Inverse 3D transform of the unit Neumann datum at positive z.
   After angular integration the radial k^2 measure cancels both 1/k factors. *)
kernel = -Integrate[Exp[-k z] Sin[k r], {k, 0, Infinity},
                    Assumptions -> r > 0 && z > 0]/(2 Pi^2 r);
fluxNorm = Integrate[4 Pi r^2 D[kernel,z], {r,0,Infinity}, Assumptions -> z>0];
contrast[expr_] := (expr /. y -> 0) - ((expr /. y -> 1/2)+(expr /. y -> -1/2))/2;
affine = Simplify[contrast[c0 + wx x + wy y + wz z]];
Print["decaying mode amplitude, inverse kernel, unit flux, affine contrast: ",
      InputForm[{amplitude,kernel,fluxNorm,affine}]];

(* Evaluate the unforced, spherically symmetric intake without Cartesian
   gridding. The radial integral is split at its integrable logarithm. *)
areaWeight = 4 Pi Integrate[a^2 (1-a^2)^3,{a,0,1}];
intake[a_] := (1/4) (1-a^2)^3/areaWeight;
phiAt[r0_] := -NIntegrate[a intake[a] Log[(r0+a)/Abs[r0-a]],
                          {a,0,r0,1},WorkingPrecision->40, AccuracyGoal->26]/(Pi r0);
vals = phiAt /@ {1/2, 1/Sqrt[2]};
Print["continuum radial receiver potentials, contrast: ", InputForm[N[Append[vals,Subtract@@vals],16]]];
If[!TrueQ[Simplify[amplitude == -j/(k nBulk)] && fluxNorm == 1 && affine == 0], Exit[1]];
Exit[0];
