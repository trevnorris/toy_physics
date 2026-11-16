(* Shell-wise solver source and Cartesian stationary carrier on a flat axial segment. *)
ClearAll["Global`*"];
verify[label_, value_, expected_] := Module[{residual = Factor[value - expected]},
  Print[label, " = ", InputForm[Factor[value]], "   residual = ", InputForm[residual]];
  If[residual =!= 0, Exit[1]]
];

vWall = k (r/radius[w])^4;
k1 = -D[vWall, radius[w]];
source = -k1 den;
verify["Path-A local density kernel", source, -4 k r^4 den/radius[w]^5];
verify["local axial flare force", -den D[vWall, w], -radius'[w] source];
(* Shell weights are fixed in the solver. Integrating r^4 den gives moment B. *)
shellSource = (source /. den -> 1) /. r^4 -> (b + eps db);
shellSource = shellSource /. radius[w] -> rad + eps eta;
verify["increment of integrated Path-A RHS", D[shellSource, eps] /. eps -> 0,
  -4 k (rad db - 5 b eta)/rad^6];

(* Cartesian field x+i y, no gauge. The real and imaginary fields obey the same GNLS. *)
dens = x[w]^2 + y[w]^2;
grad2 = x'[w]^2 + y'[w]^2;
numberFlux = hb/m (x[w] y'[w] - y[w] x'[w]);
pressure = ke dens^5;
enthalpy = 5 ke dens^4/4;
curvatureFactor = 2 m (pot[w] + enthalpy - mu)/hb^2;
ode = {x''[w] -> curvatureFactor x[w], y''[w] -> curvatureFactor y[w]};
noether = hb^2 grad2/(2 m) + (mu - pot[w]) dens - ke dens^5/4;
stress = m numberFlux^2/dens + pressure + hb^2/(4 m) (D[dens, w]^2/dens - D[dens, {w, 2}]);
verify["Cartesian axial flux minus hydrodynamic stress", (stress /. ode) - noether, 0];
verify["Cartesian flux gradient plus density times potential gradient",
  (D[noether, w] /. ode) + dens pot'[w], 0];

(* Fix local phase to zero. At unit amplitude in these units x'=a' and y'=j. *)
local = {hb -> 1, m -> 1, mu -> 1, ke -> 4/5, pot[w] -> 0, x[w] -> 1, y[w] -> 0};
excess = (stress /. ode /. local) - 4/5;
Do[
  verify["local excess at {j, a'} = " <> ToString[point[[1 ;; 2]], InputForm],
    excess /. {x'[w] -> point[[2]], y'[w] -> point[[1]]}, point[[3]]],
  {point, {{0, 0, 0}, {1/5, 0, 1/50}, {0, 1/5, 1/50}}}];
(* The advective part is j^2/n. With a resting base its first derivative vanishes. *)
advective = m (eps dj)^2/(nb + eps dn);
verify["linear advective response at rest", D[advective, eps] /. eps -> 0, 0];
verify["quadratic advective response at rest", SeriesCoefficient[advective, {eps, 0, 2}], m dj^2/nb];
