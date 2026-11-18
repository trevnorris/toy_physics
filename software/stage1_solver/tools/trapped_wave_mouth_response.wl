(* Complex forced-mode response, with axial force from virtual translation
   work on rings of a sphere. Frequencies refer to the fast shear oscillation. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, expr_, check_] := Module[{value = Factor[expr]},
  Print[label, " = ", InputForm[value]];
  If[!TrueQ[FullSimplify[value == check]], Exit[1]]
];
phasor = h amp/(det + I loss);
bc = ComplexExpand[Re[phasor]];
bs = -ComplexExpand[Im[phasor]];
emit["real displacement susceptibility times drive", bc, h amp det/(det^2 + loss^2)];
emit["sine displacement susceptibility times drive", bs, h amp loss/(det^2 + loss^2)];
(* A ring at polar angle th moves normally by dx Cos[th] under translation.
   sIn and sOut are the in-phase and quadrature strain changes per unit speed. *)
p[th_] := kap ComplexExpand[Abs[amp + u Cos[th] (sIn - I sOut)]^2]/2;
work = Integrate[2 Pi rad^2 Sin[th] dx Cos[th] p[th], {th, 0, Pi}];
load = D[work, dx];
emit["pressure without drive", p[th] /. u -> 0, kap amp^2/2];
emit["force without drive", load /. u -> 0, 0];
emit["force with purely quadrature response", load /. sIn -> 0, 0];
emit["force slope with both phases", D[load, u], 4 Pi rad^2 kap amp sIn/3];
emit["slope of the forced mode", D[load, u] /. {sIn -> bc, sOut -> bs},
     4 Pi rad^2 kap h amp^2 det/(3 (det^2 + loss^2))];
(* Lab velocity amplitude of f(x-X) exp(i omega t), differentiated for
   arbitrary in-plane translation. The Hessian is a tensor locally. *)
vel = I om f - {vx, vy, vz}.{fx, fy, fz};
energy = rho ComplexExpand[Abs[vel]^2]/4;
emit["translational xx kinetic Hessian per volume", D[energy, {vx, 2}], rho fx^2/2];
emit["translational xy kinetic Hessian per volume", D[energy, vx, vy], rho fx fy/2];
Exit[0];
