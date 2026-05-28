ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 166 — EXACT BUNDLE INVERSION OF THE LAST FOUR DRIFTS"];

Clear[dTheta, dKs, dKq, dP, drho, da, dcs, dZ, dN0, dD0];
$Assumptions = Element[{dTheta, dKs, dKq, dP, drho, da, dcs, dZ, dN0, dD0}, Reals];

eq1 = dTheta == 2*drho;
eq2 = dKs == 2*da + drho;
eq3 = dKq == dZ + 2*dcs - 2*da;
eq4 = dP == 5*(dcs - da);
sol = Solve[{eq1, eq2, eq3, eq4}, {drho, da, dcs, dZ}, Reals][[1]];

drhoSol = FullSimplify[drho /. sol, Assumptions -> $Assumptions];
daSol = FullSimplify[da /. sol, Assumptions -> $Assumptions];
dcsSol = FullSimplify[dcs /. sol, Assumptions -> $Assumptions];
dZSol = FullSimplify[dZ /. sol, Assumptions -> $Assumptions];

Print["drho = ", fmt[drhoSol]];
Print["da   = ", fmt[daSol]];
Print["dcs  = ", fmt[dcsSol]];
Print["dZ   = ", fmt[dZSol]];

banner["General inversion forms (paper Sec. 2)"];
expectZero["drho general", drhoSol - dTheta/2];
expectZero["da general", daSol - (dKs/2 - dTheta/4)];
expectZero["dcs general", dcsSol - (dKs/2 - dTheta/4 + dP/5)];
expectZero["dZ general", dZSol - (dKq - 2*dP/5)];

banner["Forward verification"];
expectZero["Theta law", (dTheta - 2*drho) /. sol];
expectZero["Ks law", (dKs - 2*da - drho) /. sol];
expectZero["Kq law", (dKq - dZ - 2*dcs + 2*da) /. sol];
expectZero["P0 law", (dP - 5*(dcs - da)) /. sol];

banner["Independent matrix-inverse cross-check"];
(* Forward map M: (drho, da, dcs, dZ) -> (dTheta, dKs, dKq, dP) from eq1..eq4. *)
Mmat = {
  {2, 0, 0, 0},   (* dTheta = 2 drho *)
  {1, 2, 0, 0},   (* dKs    = drho + 2 da *)
  {0, -2, 2, 1},  (* dKq    = -2 da + 2 dcs + dZ *)
  {0, -5, 5, 0}   (* dP     = -5 da + 5 dcs *)
};
inv = Inverse[Mmat];
solVec = inv . {dTheta, dKs, dKq, dP};
expectZero["matrix drho", solVec[[1]] - dTheta/2];
expectZero["matrix da", solVec[[2]] - (dKs/2 - dTheta/4)];
expectZero["matrix dcs", solVec[[3]] - (dKs/2 - dTheta/4 + dP/5)];
expectZero["matrix dZ", solVec[[4]] - (dKq - 2*dP/5)];
(* Round-trip: forward map of the matrix solution recovers the observables.   *)
(* Sum-of-squares scalarization (zero iff every component residual is zero;    *)
(* expectZero tests res === 0, which is False for a length-4 list).            *)
expectZero["matrix round-trip", Total[(Mmat . solVec - {dTheta, dKs, dKq, dP})^2]];

banner["Equivalent full-bundle form with P_0 = N_0 / D_0"];
dcsBundle = FullSimplify[dcsSol /. dP -> dN0 - dD0, Assumptions -> $Assumptions];
dZBundle = FullSimplify[dZSol /. dP -> dN0 - dD0, Assumptions -> $Assumptions];
Print["dcs(bundle) = ", fmt[dcsBundle]];
Print["dZ(bundle)  = ", fmt[dZBundle]];
expectZero[
  "bundle identity for dcs",
  dcsBundle - (dKs/2 - dTheta/4 + (dN0 - dD0)/5)
];
expectZero[
  "bundle identity for dZ",
  dZBundle - (dKq - 2*(dN0 - dD0)/5)
];

banner["Frozen-wall corollary"];
drhoFrozen = FullSimplify[drhoSol /. dTheta -> 0, Assumptions -> $Assumptions];
daFrozen = FullSimplify[daSol /. dTheta -> 0, Assumptions -> $Assumptions];
dcsFrozen = FullSimplify[dcsSol /. dTheta -> 0, Assumptions -> $Assumptions];
dZFrozen = FullSimplify[dZSol /. dTheta -> 0, Assumptions -> $Assumptions];
Print["drho|frozen = ", fmt[drhoFrozen]];
Print["da|frozen   = ", fmt[daFrozen]];
Print["dcs|frozen  = ", fmt[dcsFrozen]];
Print["dZ|frozen   = ", fmt[dZFrozen]];
expectZero["frozen drho", drhoFrozen];
expectZero["frozen da", daFrozen - dKs/2];
expectZero["frozen dcs", dcsFrozen - (dKs/2 + dP/5)];
expectZero["frozen dZ", dZFrozen - (dKq - 2*dP/5)];

banner["Explicit Family-1 wall density"];
thetaChi = SetPrecision[4.06863235008162, 30];
rhoChi = Sqrt[thetaChi/25];
Print["rho_w^(chi) = ", fmt[N[rhoChi, 18]], " * lambda_mu^(-1)"];

Print[""];
Print["Carry-forward formulas:"];
Print["  delta ln rho_w = 1/2 delta ln Theta_w"];
Print["  delta ln a     = 1/2 delta ln K_s - 1/4 delta ln Theta_w"];
Print["  delta ln c_s   = 1/2 delta ln K_s - 1/4 delta ln Theta_w + 1/5 delta ln P_0"];
Print["  delta ln Z_q   = delta ln K_q - 2/5 delta ln P_0"];
Print["  with delta ln P_0 = delta ln N_0 - delta ln D_0 on the isotropic bundle."];

Print[""];
Print["Stage 166 Mathematica audit passed."];

Exit[0];
