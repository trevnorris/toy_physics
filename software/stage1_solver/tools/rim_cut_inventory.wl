(* Independent angular momentum projection and frame check for the signed receiver. *)
ClearAll["Global`*"];
angularMean[f_] := Integrate[f, {u, -1, 1}]/2;
(* u=x/r. Cartesian x derivative of f(r,u) at fixed other coordinates. *)
dx[f_] := u D[f, r] + (1-u^2)/r D[f, u];
field = f[r] + I h[r] + e u (p[r] + I q[r]);
conjField = f[r] - I h[r] + e u (p[r] - I q[r]);
px = hb (conjField dx[field] - field dx[conjField])/(2 I);
projected = Simplify[angularMean[Coefficient[Expand[px], e, 1]]];
rho = f[r]^2 + h[r]^2;
drho = 2(f[r] p[r] + h[r] q[r]);
g = hb/m (f[r] q[r] - h[r] p[r])/rho;
v0 = hb/m (f[r] h'[r] - h[r] f'[r])/rho;
hydro = m (rho (D[g,r] + 2 g/r) + drho v0)/3;
angularResidual = Simplify[Together[projected - hydro]];
Print["angular projection = ", InputForm[projected], " hydro residual = ", angularResidual];

dir = DirectoryName[$InputFileName];
log = FileNameJoin[{dir, "..", "..", "..", "_scratch", "s11c", "rim_cut_inventory01", "stdout"}];
rows = ImportString[#, "RawJSON"]& /@ Import[log, "Lines"];
check[row_] := Module[{boost, mass, dp, relative, lab, imbalance, fraction},
  boost = row["boost_flux"];
  mass = row["mass"];
  relative = row["cut_flow_comoving"] + row["cut_stress"];
  lab = relative + boost;
  dp = row["s"] row["p_relative"] + mass row["A"];
  imbalance = dp + relative + row["side_to_support"] - (row["feed_lab"] - boost);
  fraction = Abs[imbalance/dp];
  If[row["A"] > 0, Print["level=", row["level"], " s=", row["s"],
    " relative_cut=", InputForm[relative], " storage_rate=", InputForm[dp],
    " closure_error_over_storage_rate=", InputForm[fraction]]];
  Max[Abs[{lab - row["cut_lab"], dp - row["storage_rate_lab"],
    imbalance - row["balance_residual"], boost - row["rate"] row["A"]/row["s"]}]]
];
residuals = check /@ rows;
cases = Sort[({#["level"], #["s"], #["A"]}& /@ rows)];
expected = Sort[Tuples[{{64,128},{0.01,0.02},{-0.0001,0.0001}}]];
Print["samples=", Length[rows], " frame check max residual=", InputForm[Max[residuals]]];
If[!TrueQ[cases == expected && angularResidual == 0 && Max[residuals] < 10^-14], Exit[1]];
