(* Project the companion potential before differentiating. Integrate the
   pole at distance rsrc > r by changing the angular variable to distance. *)
ClearAll["Global`*"];
mu[l_] := (rsrc^2 + r^2 - l^2)/(2 r rsrc);
projectedPotential = Factor[(3/2) Integrate[
    Expand[mu[l] (strength/(4 Pi l)) l/(r rsrc)], {l, rsrc-r, rsrc+r}]];
normalDipole = Factor[D[projectedPotential, r]];
centerVelocity = D[strength/(4 Pi (rsrc-x)), x] /. x -> 0;
Print["dipole of potential, radial derivative, center velocity = ",
  ToString[{projectedPotential, normalDipole, centerVelocity}, InputForm]];
If[Together[normalDipole-centerVelocity] =!= 0, Exit[1]];

(* Both aperture cuts are oriented into the drain and carry the same
   positive inward number rate. The transverse density stress is invariant. *)
heldFeed = m j vfeed + stressFeed;
heldExit = m j vexit + stressExit;
movingFeed = m j (vfeed-u) + stressFeed;
movingExit = m j (vexit-u) + stressExit;
Print["changes of feed, receiving and feed-minus-receiving flux = ",
  ToString[Simplify[{movingFeed-heldFeed, movingExit-heldExit,
    (movingFeed-movingExit)-(heldFeed-heldExit)}], InputForm]];
If[Expand[(movingFeed-movingExit)-(heldFeed-heldExit)] =!= 0, Exit[1]];

(* Independently solve the projected relaxation law for a step of each
   signed slip from a common, unstrained state. *)
ang = Integrate[Cos[th]^2 Sin[th], {th,0,Pi}, {az,0,2 Pi}];
load[slip_] := k a^2 ang DSolveValue[
  {e'[t] + e[t]/tau == slip, e[0] == 0}, e[t], t];
heldLoad = FullSimplify[load[u], Assumptions -> tau > 0];
movingLoad = FullSimplify[load[0-(-u)], Assumptions -> tau > 0];
power = FullSimplify[v Limit[load[-v], t -> Infinity, Assumptions -> tau > 0],
  Assumptions -> tau > 0];
Print["held step, translating step at V=-U, work rate for drift at V = ",
  ToString[{heldLoad, movingLoad, power}, InputForm]];
If[Simplify[heldLoad-movingLoad] =!= 0 ||
   Simplify[power+4 Pi a^2 k tau v^2/3] =!= 0, Exit[1]];
