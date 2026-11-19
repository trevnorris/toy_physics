(* Independent section-degree and weak transport calculation in units R=1. *)
ClearAll["Global`*"];
check[label_, value_, expected_] := Module[{v = Factor[value]},
  Print[label, " = ", InputForm[v]];
  If[!TrueQ[FullSimplify[v == expected]], Exit[1]]
];

(* Integrating the angular area around a section gives half the change
   in the axial component of its unit normal between the two poles. *)
axial[c_] := (r c - b)/Sqrt[r^2 + b^2 - 2 r b c];
degree = (axial[1] - axial[-1])/2;
check["section outside the zero", FullSimplify[degree, r > b > 0], 1];
check["section inside the zero", FullSimplify[degree, b > r > 0], 0];

(* A radial section has fixed chi throughout its time history. The radial
   flux of -dchi ^ omega is its pullback to (t,theta,phi). *)
chiFace = chi[r0];
check["three pulled-back components of dchi on a cut",
      D[chiFace, t]^2 + D[chiFace, theta]^2 + D[chiFace, phi]^2, 0];
weight[r_] := 30 (r-2)^2 (3-r)^2;
left[b_] = Integrate[weight[r], {r,b,3}];
check["remaining layer count at entry, halfway, exit",
      {left[2], left[5/2], left[3]}, {1,1/2,0}];

(* Test local continuity with an arbitrary time-independent radial test
   function f. The derivative of the moving lower integration limit gives
   the loss term. The delta carried at the zero weighs A[b] f[b]; its flux
   against f' follows from its actual coordinate speed. *)
lowerTerm = D[Integrate[weight[r] f[r], {r,b,3}, GenerateConditions->False], b];
check["layer weak source per unit tilt", lowerTerm, -weight[b] f[b]];
corrected = lowerTerm + D[A[b] f[b], b] - A[b] f'[b];
check["weak balance coefficient after adding the moving zero",
      Cancel[corrected/f[b]], A'[b] - weight[b]];
carried[b_] = Integrate[weight[r], {r,2,b}];
check["local residual with initially unweighted zero",
      corrected /. Derivative[1][A][b] -> carried'[b], 0];
check["layer plus assigned zero weight", left[b]+carried[b], 1];
check["uncompensated weak source at halfway, unit speed and f=1",
      lowerTerm /. {b->5/2, f[_]->1}, -15/8];
check["signed change across upper and reflected layer",
      Table[s (left[3]-left[2]), {s,{1,-1}}], {-1,1}];
