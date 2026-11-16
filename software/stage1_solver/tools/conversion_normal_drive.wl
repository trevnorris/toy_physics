(* Reversible normal load from a parcel stretch, followed by a separate
   travelling-front balance. All quantities are per unit transverse area. *)
ClearAll["Global`*"];
well[q_] := a q^2 (1-q)^2;
(* Unit reference length dilates by b: n -> n/b, grad chi -> grad chi/b.
   Order is carried for this virtual test, and external core terms are absent. *)
energyPerReference = b (K0 (n/b)^5/4 + well[q] + k (s/b)^2/2);
momentumFlux = Expand[-D[energyPerReference, b] /. b -> 1];
orderFlux = Simplify[momentumFlux - K0 n^5];
normalForce = -D[momentumFlux /. {n -> dens[w], q -> chi[w], s -> chi'[w]}, w];
orderDerivative = D[well[q], q] /. q -> chi[w];
forceResidual = Simplify[normalForce + dens[w] D[5 K0 dens[w]^4/4, w]
    - (orderDerivative - k chi''[w]) chi'[w]];
endJump = Simplify[(orderFlux /. {q -> 0, s -> 0}) - (orderFlux /. {q -> 1, s -> 0})];
(* On the decreasing stage-006 kink, s = -q(1-q)/ell. *)
kinkFlux = Simplify[orderFlux /. {k -> 2 a ell^2, s -> -q (1-q)/ell}];
Print["material force identity residual = ", forceResidual];
Print["bare order [flat-end flux jump, kink flux] = ", ToString[{endJump, kinkFlux}, InputForm]];

(* Integrate n Gamma across a front in the front frame. The integration
   variable is its order fraction; no lab current or front speed is input. *)
orderRate = -lam q (1-q);
slope = -q (1-q)/ell;
integratedSource = Integrate[Cancel[n0 orderRate/slope], {q, 1, 0}];
relativeVelocity = Simplify[(u - v) /. First[Solve[-n0 (u-v) == integratedSource, u]]];
localResidual = Simplify[relativeVelocity slope - orderRate];
stationaryFrontCurrent = n0 relativeVelocity;
restFrontSpeed = v /. First[Solve[-v == relativeVelocity, v]];
restMaterialCurrent = Simplify[n0 (restFrontSpeed + relativeVelocity)];
Print["u - V = ", ToString[relativeVelocity, InputForm]];
Print["fixed-cut number currents [stationary front, stationary material] = ",
    ToString[{stationaryFrontCurrent, restMaterialCurrent}, InputForm]];
If[{forceResidual, endJump, kinkFlux, localResidual} =!= {0, 0, 0, 0}, Exit[1]];
