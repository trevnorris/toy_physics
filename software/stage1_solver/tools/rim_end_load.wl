(* Independent calculation via endpoint virtual work of the two fields.
   q'' = (s/c)^2 q. For positive real s, the quadratic boundary form is
   half the integral of tension q'^2 + density s^2 q^2 on the solution. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, expr_, target_] := Module[{value=Factor[expr]},
  Print[label, " = ", InputForm[value]];
  If[!TrueQ[FullSimplify[value == target, s>0 && tau>0 && z1>0 && z2>0 && 0<h<1]],Exit[1]]
];
profile = (left Sinh[s (tau-u)] + right Sinh[s u])/Sinh[s tau];
(* u is travel time from the left end. Force is +/- Z dq/du. *)
boundaryWork = z (right (D[profile,u] /. u->tau) - left (D[profile,u] /. u->0))/2;
quad = FullSimplify[boundaryWork, s>0 && tau>0];
layer[a_,b_,imp_] := quad /. {left->a, right->b, z->imp};
response[za_,zb_] := Module[{action,stationary,force},
  action = layer[x,y,za] + layer[y,v,zb];
  stationary = First[Solve[{D[action,y]==0,D[action,v]==0},{y,v}]];
  force = FullSimplify[(D[action,x] /. stationary)/x,
       s>0 && tau>0 && za>0 && zb>0];
  force
];
k12=response[z1,z2]; k21=response[z2,z1];
hyp=Tanh[s tau];
emit["two free-end field layers",k12, s hyp (z1+z2)/(1+hyp^2 z2/z1)];
emit["reversed free-end layers",k21, s hyp (z1+z2)/(1+hyp^2 z1/z2)];
emit["slow inertia in each order", Limit[{k12,k21}/s^2,s->0],{1,1} tau (z1+z2)];
emit["leading dynamic difference", SeriesCoefficient[k12-k21,{s,0,4}],
     tau^3 (z1+z2)(z1/z2-z2/z1)];
near = (a x^2+2 b x y+d y^2)/2;
rim[end_] := D[near,x]/x /. First[Solve[D[near+end y^2/2,y]==0,y]];
emit["load sensitivity of a common reciprocal near region",rim[g2]-rim[g1],
     b^2 (g2-g1)/((d+g1)(d+g2))];
Exit[0];
