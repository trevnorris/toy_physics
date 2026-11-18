(* Build the carrier momentum tensor in four Cartesian components. Integrate
   on the in-brane sphere before choosing any profile along the bulk axis. *)
ClearAll["Global`*"];
$HistoryLength = 0;
rhat = {Cos[alpha], Sin[alpha] Cos[phi], Sin[alpha] Sin[phi], 0};
ex = {1, 0, 0, 0};
ew = {0, 0, 0, 1};
mu = ex.rhat;
normal = sn rhat + cs ew;
base = vr rhat + vw ew;
perturbation = ar mu rhat + at (ex - mu rhat) + aw mu ew;
vel = base + u perturbation;
rho = rho0 + u rho1 mu;
pressure = p0 + u p1 mu;
tensor = rho Outer[Times, vel, vel] + pressure IdentityMatrix[4];
angularSlope[field_] := Module[{slope},
  slope = TrigReduce[Expand[D[field, u] /. u -> 0]];
  Simplify[Integrate[slope Sin[alpha], {alpha, 0, Pi}, {phi, 0, 2 Pi}]/(4 Pi)]
];
mass = angularSlope[rho vel.normal];
flux = angularSlope[ex.tensor.normal];
vn = sn vr + cs vw;
jBase = rho0 vn;
jSlope = rho1 vn + rho0 (sn ar + cs aw);
Print["mass angular mean = ", InputForm[mass]];
Print["three times momentum angular mean = ", InputForm[Expand[3 flux]]];
If[mass =!= 0 || Expand[3 flux - vr jSlope - jBase (ar + 2 at) - sn p1] =!= 0, Exit[1]];

(* Independent ideal limit. On a uniform radial discharge, differentiate the
   kinetic tensor and Bernoulli pressure directly, averaging on S^3 about x. *)
n = {Cos[gamma], Sin[gamma], 0, 0};
potential = drift rad + dip/rad^3;
speed = q/(2 Pi^2 rad^3 density);
v = speed n + u (D[potential, rad] n n[[1]] + potential/rad (ex - n n[[1]]));
idealTensor = density (Outer[Times, v, v] - IdentityMatrix[4] (v.v)/2);
idealFlux = Simplify[4 Pi rad^3 Integrate[
    TrigReduce[D[ex.idealTensor.n, u] /. u -> 0] Sin[gamma]^2, {gamma, 0, Pi}]];
Print["uniform plus decaying dipole = ", InputForm[idealFlux]];
If[Expand[idealFlux - q drift] =!= 0, Exit[1]];
Exit[0];
