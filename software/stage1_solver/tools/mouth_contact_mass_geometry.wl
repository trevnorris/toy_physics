(* Linearize about a marked, nonrotating sphere.  Two small angular
   displacements on each patch are identified by a measured registration
   with the physical tangent rotation vector lengthMap s.  Face f rotates
   by f times that vector.  All s values refer to the same marked face.
   A bearing translates with the observed sphere and carries two signed
   radial arms.  mu and carrier are actual masses per reference area.
   Positive mu, arm, lengthMap, radius, nonnegative carrier, |eta| <= 1. *)
ClearAll["Global`*"];
emit[label_, x_] := WriteString["stdout", label, " = ", ToString[x, InputForm], "\n"];
checks = {};
compare[label_, actual_, expected_] := Module[{residual},
  residual = Factor[actual - expected];
  emit[label, {actual, expected, residual}];
  AppendTo[checks, AllTrue[Flatten[{residual}], # === 0 &]]];

signs = {-1, 1};
weight[f_, e_] := mu (1 + f e eta)/4;
normal = {0, 0, 1};
velocity = {vx, vy, vz};
spin = {sx, sy, 0};
pointVelocity[f_, e_, v_, w_, n_] := v + Cross[f lengthMap w, e arm n];
energy[v_, w_, n_] := Expand[carrier v.v/2 +
  Sum[With[{uvel = pointVelocity[f, e, v, w, n]},
    weight[f, e] uvel.uvel/2], {f, signs}, {e, signs}]];
localEnergy = energy[velocity, spin, normal];
localMass = D[localEnergy, {vx, 2}];
localInertia = D[localEnergy, {sx, 2}];
crossMass = -D[localEnergy, sx, vy];
survey = {
  carrier + Sum[weight[f, e], {f, signs}, {e, signs}],
  Sum[weight[f, e] (lengthMap e arm)^2, {f, signs}, {e, signs}],
  Sum[weight[f, e] f lengthMap e arm, {f, signs}, {e, signs}]};
compare["velocity Hessian and face-resolved mass survey {mass,inertia,cross}",
  {localMass, localInertia, crossMass}, Factor /@ survey];
emit["{cross/inertia, cross/inertia under face exchange, balanced cross}",
  Factor /@ {crossMass/localInertia, (crossMass/localInertia /. eta -> -eta),
             crossMass /. eta -> 0}];
armMasses = Table[Factor[Sum[weight[f, e], {f, signs}]], {e, signs}];
faceMasses = Table[Factor[Sum[weight[f, e], {e, signs}]], {f, signs}];
emit["face totals, arm totals, angular inertia at eta={-1,0,1}",
  {faceMasses, armMasses, Table[localInertia /. eta -> h, {h, -1, 1}]}];
relativeSpin = spin - (crossMass/localInertia) Cross[velocity, normal];
tangent = velocity - (velocity.normal) normal;
transverseRemainder = Factor[localMass - crossMass^2/localInertia];
emit["local mass remaining transverse to normal after square completion", transverseRemainder];
completeEnergy = localInertia relativeSpin.relativeSpin/2 +
   transverseRemainder tangent.tangent/2 + localMass (velocity.normal)^2/2;
compare["point energy and angular square with its translating mass terms",
   localEnergy, Expand[completeEnergy]];
boost = {bx, by, bz};
pointMomentum = Expand[carrier velocity + Sum[weight[f, e]
   pointVelocity[f, e, velocity, spin, normal], {f, signs}, {e, signs}]];
compare["boost increment and measured total point momentum term",
   Expand[energy[velocity + boost, spin, normal] - localEnergy],
   Expand[boost.pointMomentum + localMass boost.boost/2]];

(* Orient the translation axis along z.  The two tangential harmonics of
   degree l are grad_S P_l and its rotation by n.  Only the latter has an
   azimuthal component for m=0.  Other m integrate to zero.  Integration
   by parts reduces the uniform azimuthal overlap for every l to
   2 Integrate[u P_l(u)], the inner product with P_1. *)
overlap[ell_, angularProfile_] := 2 Pi Integrate[
    angularProfile (1 - u^2) D[LegendreP[ell, u], u], {u, -1, 1}];
compare["uniform toroidal overlaps, degrees 1 through 6, and Legendre inner products",
   Table[overlap[l, 1], {l, 6}],
   Table[4 Pi Integrate[u LegendreP[l, u], {u, -1, 1}], {l, 6}]];
emit["toroidal overlaps for cross profile 1+gradient u, degrees 1 through 3",
   Table[overlap[l, 1 + gradient u], {l, 3}]];
ez = {0, 0, 1};
sphereNormal = {Sqrt[1 - u^2], 0, u};
meridian = {u, 0, -Sqrt[1 - u^2]};
compare["translation drive projected onto meridian",
   meridian.Cross[ez, sphereNormal], 0];
integrateSurface[x_] := Factor[2 Pi radius^2 Integrate[Expand[x], {u, -1, 1}]];

(* s = H (zhat cross n), geometric z velocity = v.  The integrated
   energy uses velocities of the individual points at each latitude. *)
dipoleEnergy = integrateSurface[energy[v ez, hvel Cross[ez, sphereNormal], sphereNormal]];
shellMass = Factor[D[dipoleEnergy, {v, 2}]];
shellInertia = Factor[D[dipoleEnergy, {hvel, 2}]];
shellCross = Factor[-D[dipoleEnergy, v, hvel]];
area = 4 Pi radius^2;
compare["dipole Hessian {M,I,K} and geometric moments",
   {shellMass, shellInertia, shellCross},
   {area localMass, 2 area localInertia/3, 2 area crossMass/3}];
linearMomentum = Factor[D[dipoleEnergy, v]];
angularMomentum = Factor[D[dipoleEnergy, hvel]];
(* Finite rotation about the local azimuthal axis gives a second route to
   true material center at linear order. zoffset is measured from the
   translating geometric center. *)
zoffset[f_, e_] := With[{phi = f lengthMap hpos Sqrt[1 - u^2]},
   radius u + e arm (u Cos[phi] - Sqrt[1 - u^2] Sin[phi])];
zMoment = carrier radius u + Sum[weight[f, e] zoffset[f, e], {f, signs}, {e, signs}];
center0 = integrateSurface[zMoment /. hpos -> 0]/shellMass;
centerSlope = Factor[integrateSurface[D[zMoment, hpos] /. hpos -> 0]/shellMass];
compare["geometric-to-material center {offset,slope} and kinetic prediction",
   {center0, centerSlope}, {0, Factor[-shellCross/shellMass]}];
compare["constituent momentum and mass times material center velocity",
   linearMomentum, Factor[shellMass (v + centerSlope hvel)]];
centerEnergy = Expand[dipoleEnergy /. v -> cvel - centerSlope hvel];
compare["material-center energy and separated translation and relative rotation",
   centerEnergy,
   Expand[shellMass cvel^2/2 + (shellInertia - shellCross^2/shellMass) hvel^2/2]];
clampForce = Factor[D[linearMomentum, hvel] hacc];
holdingTorque = Factor[D[angularMomentum, v] vacc];
emit["force on held geometric center during angular drive, reciprocal torque holding H fixed",
   {clampForce, holdingTorque}];
compare["clamp-inferred, reciprocal and surveyed cross moment",
   {-clampForce/hacc, -holdingTorque/vacc}, {2 area crossMass/3, 2 area crossMass/3}];
emit["free marked-center velocity at zero momentum, relative dipole inertia",
   {Factor[v /. First[Solve[linearMomentum == 0, v]]],
    Factor[D[centerEnergy, {hvel, 2}]]}];
emit["all identity checks", And @@ checks];
Exit[If[And @@ checks, 0, 1]];
