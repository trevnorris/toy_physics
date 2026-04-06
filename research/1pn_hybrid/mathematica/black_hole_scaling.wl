(* ---------------------------------------------------------------------- *)
(* Script: black_hole_scaling.wl                                           *)
(* Purpose: Turn the representative impedance fit into an approximate      *)
(* horizon-radius and compactness model.                                   *)
(* ---------------------------------------------------------------------- *)

ClearAll["Global`*"]

lambdaStar = N[(Sqrt[2] * Pi)/BesselJZero[0, 1], 16];
sigmaFS = N[1/(1 + lambdaStar^2), 16];
kappa = 2.0;

Mmodel[a_] := a * (1 + sigmaFS/(a^2 + 1));
rH[a_] := kappa * a;
compactness[a_] := Mmodel[a]/rH[a];

sampleRadii = {1.0, 2.0, 4.0, 8.0};
sampleMasses = N[Mmodel /@ sampleRadii, 16];
sampleHorizons = N[rH /@ sampleRadii, 16];
fitExpr = Normal[Fit[Transpose[{Log[sampleMasses], Log[sampleHorizons]}], {1, x}, x]];
logFitExponent = N[Coefficient[fitExpr, x], 16];

Print["--- Representative Horizon Scaling ---"];
Print["lambda_* = ", lambdaStar];
Print["sigma_fs = ", sigmaFS];
Print["log-fit exponent for r_H(M) = ", logFitExponent];
Print["a=1  -> {M, r_H, C} = ", N[{Mmodel[1.0], rH[1.0], compactness[1.0]}, 16]];
Print["a=2  -> {M, r_H, C} = ", N[{Mmodel[2.0], rH[2.0], compactness[2.0]}, 16]];
Print["a=4  -> {M, r_H, C} = ", N[{Mmodel[4.0], rH[4.0], compactness[4.0]}, 16]];
Print["a=8  -> {M, r_H, C} = ", N[{Mmodel[8.0], rH[8.0], compactness[8.0]}, 16]];
Print["Compactness limit = ", Limit[compactness[a], a -> Infinity]];

(*"
Output:

--- Representative Horizon Scaling ---
lambda_* = 1.84748657712012805107642087596939904276`16.
sigma_fs = 0.2265926068524372168300616723354750574`15.81056168507603
log-fit exponent for r_H(M) = 1.0516533277732314
a=1  -> {M, r_H, C} = {1.1132963034262187, 2., 0.5566481517131093}
a=2  -> {M, r_H, C} = {2.090637042740975, 4., 0.5226592606852437}
a=4  -> {M, r_H, C} = {4.053315907494691, 8., 0.5066644884368364}
a=8  -> {M, r_H, C} = {8.027888320843378, 16., 0.5017430200527111}
Compactness limit = 0.5
"*)
