(* Local dimensionless stationary GNLS in carrier form. The background is
   Exp[I w/hbar], with m=1, n0=1, V=1 and chemical potential set by the EOS.
   Test the isolated and driven families directly in the complex equation. *)
ClearAll["Global`*"];
hbar = 1/Sqrt[6];
enthalpy[n_] := 35 n^4/96;
chemicalPotential = enthalpy[1] + 1/2;
laplacian[field_] := D[field, {x, 2}] + D[field, {w, 2}];
base = Exp[I w/hbar];
carrierResidual[densityChange_, phasePotential_] := Module[{variation, residual},
  variation = base (densityChange/2 + I phasePotential/hbar);
  residual = -hbar^2 laplacian[variation]/2 +
              (enthalpy[1] - chemicalPotential) variation +
              enthalpy'[1] densityChange base;
  FullSimplify[residual/base]
];
localReadings[densityChange_, phasePotential_] := Module[{envelope, current, density},
  envelope = densityChange/2 + I phasePotential/hbar;
  density = ComplexExpand[2 Re[envelope]];
  (* From j = hbar Im[conj(psi) grad(psi)] and v=j/n at positive density. *)
  current = ComplexExpand[hbar Im[{D[envelope, x],
             D[envelope, w] + I (envelope + Conjugate[envelope])/hbar}]];
  FullSimplify[{density, current[[1]], current[[2]]-density, current[[2]]}]
];
branches = {
  {"isolated", a Exp[-Sqrt[11] w], a Exp[-Sqrt[11] w]/Sqrt[11]},
  {"paired", b Sin[x] (8 Exp[-3 w] - 3 Exp[-2 w]),
             u x + b Sin[x] (3 Exp[-3 w] - 2 Exp[-2 w])}
};
Do[
  {name, densityChange, phasePotential} = branch;
  residual = carrierResidual[densityChange, phasePotential];
  readings = localReadings[densityChange, phasePotential];
  Print[name, " linearized complex GNLS residual = ", InputForm[residual]];
  Print[name, " outlet changes [n,v_x,v_w,j_w] = ",
        InputForm[FullSimplify[readings /. w -> 0]]];
  Print[name, " downstream changes = ", InputForm[Limit[readings, w -> Infinity]]];
  If[residual =!= 0, Exit[1]],
  {branch, branches}];
