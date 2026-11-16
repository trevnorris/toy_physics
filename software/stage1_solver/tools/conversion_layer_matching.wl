(* Steady normal tube balance and moving-sheet momentum, one throat, zero sideways flow. *)
ClearAll["Global`*"];
verify[label_, expr_] := Module[{v = Factor[expr]},
  Print[label, " = ", ToString[v, InputForm]];
  If[v =!= 0, Exit[1]]];

(* The moving Heaviside contributes [flux]-W[density] to each conservation law. *)
numberSheet = np up - nm um - vel (np - nm);
momentumSheet = tp + m np up^2 - tm - m nm um^2 - vel m (np up - nm um);
downstream = First[Solve[numberSheet == 0, up]];
reaction = Factor[-momentumSheet /. downstream /. um -> vel + j/nm];
Print["momentum delivered to layer per area = ", ToString[Expand[reaction], InputForm]];
verify["reaction minus normal stress and recoil", reaction - tm + tp - m j^2 (1/nm - 1/np)];
verify["derivative of reaction with respect to surface speed at fixed j", D[reaction, vel]];

(* Start with d/dz of planar normal order stress. Use the stationary kinetic law
   to eliminate the second derivative, with curved Laplacian (A chi')'/A. *)
slope = chi'[z];
mu = f'[chi[z]] - kap D[area[z] slope, z]/area[z] + muExtra[z];
kinetics = flow slope/(area[z] den[z]) == -mob[z] mu;
second = First[Solve[kinetics, chi''[z]]];
stress = kap slope^2/2 - f[chi[z]];
stressGradient = D[stress, z] /. second;
verify["normal order stress gradient identity",
  stressGradient - flow slope^2/(area[z] den[z] mob[z]) +
  kap area'[z]/area[z] slope^2 - muExtra[z] slope];

(* Explicit planar profile in the independent coordinate x=exp(z/ell). *)
profile = 1/(1+x);
qz = x D[profile, x]/ell;
profileNorm = Integrate[qz^2 ell/x, {x,0,Infinity}, Assumptions -> ell > 0];
Print["planar kink integral of slope squared = ", ToString[profileNorm, InputForm]];
verify["positive planar weight normalization", ell profileNorm - 1/6];
pureWell[q_] := a q^2 (1-q)^2;
verify["planar bulk order stress jump", pureWell[1] - pureWell[0]];
