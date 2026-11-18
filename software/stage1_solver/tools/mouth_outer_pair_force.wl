(* Paired-minus-isolated exterior momentum on a sphere enclosing one target.
   The source is at +x. Mass flux is outward-positive in each ideal exterior.
   Each angular integral keeps the full, spatially varying companion velocity. *)
ClearAll["Global`*"];
$HistoryLength = 0;
outwardDifference[d_, rho_, qt_, qs_, b_, sep_] := Module[
  {omega, measure, distance, targetR, compX, compR, comp2, mu, adv, pressure},
  omega = 2 Pi^(d/2)/Gamma[d/2];
  measure = 2 Pi^((d-1)/2)/Gamma[(d-1)/2] b^(d-1) Sin[theta]^(d-2);
  mu = Cos[theta];
  targetR = qt/(rho omega b^(d-1));
  distance = Sqrt[b^2 + sep^2 - 2 b sep mu];
  compX = qs (b mu - sep)/(rho omega distance^d);
  compR = qs (b - sep mu)/(rho omega distance^d);
  comp2 = (qs/(rho omega distance^(d-1)))^2;
  adv = rho ((targetR mu + compX) (targetR + compR) - targetR^2 mu);
  pressure = -rho (2 targetR compR + comp2) mu/2;
  NIntegrate[Evaluate[measure (adv + pressure)], {theta, 0, Pi},
    WorkingPrecision -> 42, AccuracyGoal -> 30, PrecisionGoal -> 26]
];
Do[
  feed = outwardDifference[3, 2, -1/4, -3/8, cut, distance];
  bulk = outwardDifference[4, 3, 1/4, 3/8, cut, distance];
  Print["distance, cut radius, outward brane and bulk flux, total holder: ",
    InputForm[N[{distance, cut, feed, bulk, feed+bulk}, 15]]];
  (* Scale each independently integrated reading using its own measure. *)
  scaled = {feed (4 Pi) 2 distance^2, bulk (2 Pi^2) 3 distance^3};
  Print["scaled fluxes: ", InputForm[N[scaled, 15]]];
  If[Max[Abs[scaled + 3/32]] > 10^-23, Exit[1]],
  {distance, {8, 16}}, {cut, {1, 3}}];
Exit[0];
