(* Cartesian potential gradients and finite-amplitude Bernoulli tractions on
   two external spheres. The receiver is an idealized four-dimensional channel. *)
ClearAll["Global`*"];
$HistoryLength = 0;
oddFlux[d_, rho_, q_, radius_, speed_, dipole_] := Module[
  {coords, r, omega, normal, potential, vel, tensor, axial, halfDifference, integrate},
  coords = Table[Unique["coord"], d];
  r = Sqrt[coords.coords];
  omega = 2 Pi^(d/2)/Gamma[d/2];
  normal = PadRight[{Cos[theta], Sin[theta]}, d];
  potential = q r^(2-d)/((2-d) rho omega) +
              eps First[coords] (speed + dipole/r^d);
  vel = FullSimplify[Grad[potential, coords] /. Thread[coords -> radius normal],
                     0 < theta < Pi];
  tensor = rho (Outer[Times, vel, vel] - IdentityMatrix[d] (vel.vel)/2);
  axial = First[tensor.normal];
  halfDifference = Expand[((axial /. eps -> 1) - (axial /. eps -> -1))/2];
  integrate = 2 Pi^((d-1)/2)/Gamma[(d-1)/2] radius^(d-1);
  FullSimplify[integrate Integrate[halfDifference Sin[theta]^(d-2), {theta, 0, Pi}]]
];
(* Distinct densities, radii and dipoles test cancellation independently of
   a shared volume current or vanishing local receiver velocity. *)
Do[
  feed = oddFlux[3, 2, -1/4, radius, u, 3 u];
  bulk = oddFlux[4, 3, 1/4, radius+1, v, -5 u];
  forces = Table[Expand[feed + (bulk /. v -> fraction u)], {fraction, {0, 1/2, 1}}];
  Print["radius pair ", {radius, radius+1}, " brane and bulk ", InputForm[{feed, bulk}],
        " holder for W/U=0,1/2,1 ", InputForm[forces]];
  If[!TrueQ[FullSimplify[feed == -u/4 && bulk == v/4]], Exit[1]],
  {radius, {2, 5}}];
Exit[0];
