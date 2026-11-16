(* Fixed radius and intake label the state. Compare conjugate momentum for
   a potential-only boundary with a possible residual kinetic term, after
   accounting for the brane. The coefficient is a function of radius and
   intake; there is no assigned particle mass or fit to the Newtonian law. *)
stateCoefficient = b[a, j];
candidate = stateCoefficient v[t]^2/2;
momentum = D[candidate, v[t]];
drive = D[momentum, t];
(* This term is independent of all carrier field data for the imposed
   motion, with global a and j held fixed. A scalar static surface energy
   alone has no velocity derivative at fixed surface geometry. *)
coords = {density, phase, dn, dphase};
fieldVariations = D[candidate, #] & /@ coords;
potentialMomentum = D[energy[a,j,density,dn],v[t]];

(* One possible matched surface source is the kinetic energy of normal
   displacement. For a translated sphere the normal speed is v cos(theta).
   Integrate the dipole to get the state coefficient for a uniform areal
   inertia density. This is a conditional example, not an input value. *)
rimKinetic = Integrate[s[a,j] (v[t] mu)^2/2 * 2 Pi a^2,{mu,-1,1}];
rimCoefficient = D[rimKinetic, {v[t],2}];
(* Matched-speed crossings remove forces depending on speed alone. *)
opposed[force_] := ((force /. v'[t] -> acc) - (force /. v'[t] -> -acc))/(2 acc);
mouthLoad = steady[v[t]] - carrier accVar;
totalDrive = drive - (mouthLoad /. accVar -> v'[t]);
inferred = Simplify[opposed[totalDrive] + opposed[mouthLoad /. accVar -> v'[t]]];
Print["static surface momentum = ", InputForm[potentialMomentum]];
Print["candidate boundary momentum, drive = ", InputForm[{momentum,drive}]];
Print["prescribed-state carrier variations = ", InputForm[fieldVariations]];
Print["normal-only isotropic rim coefficient = ", InputForm[rimCoefficient]];
Print["total drive contrast per acceleration contrast = ",InputForm[Simplify[opposed[totalDrive]]]];
Print["total drive plus inward mouth load contrasts per acceleration contrast = ", InputForm[inferred]];
If[fieldVariations =!= {0,0,0,0} || potentialMomentum =!= 0 ||
   Simplify[inferred-stateCoefficient] =!= 0 || Simplify[drive-stateCoefficient v'[t]] =!= 0, Exit[1]];
