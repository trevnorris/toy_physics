(* Derive the endpoint response from virtual work and condensed dynamic
   stiffness. Independently use the relative equation for the release. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, v_, target_] := Module[{value=Factor[v]},
  Print[label," = ", InputForm[value]];
  If[!TrueQ[FullSimplify[value == target,
      mr>0 && mw>0 && k>0 && s>0 && nu>0 && eps>0 && t>0 && g>0 && 1/4<r<1]], Exit[1]]
];
form = (k (x-y)^2 + s^2 (mr x^2 + mw y^2))/2;
heldWave = y /. First[Solve[D[form,y]==0,y]];
response = Factor[(D[form,x] /. y->heldWave)/(s^2 x)];
emit["rim response at positive Laplace rate", response, mr+mw k/(k+mw s^2)];
z = Factor[(heldWave/x /. s->I nu) /. k->mw g];
ineq = Reduce[g>nu^2 && z-1<=eps && nu>0 && eps>0, g, Reals];
Print["subresonant tracking region: ",InputForm[ineq]];
If[!TrueQ[Resolve[ForAll[{nu,eps,g}, Implies[nu>0 && eps>0,
     Equivalent[g>nu^2 && z-1<=eps,g>=nu^2 (1+1/eps)]]],Reals]],Exit[1]];
(* Solve the three independent rim readings for g without fixing either mass. *)
data={m0==mr+mw, m1==(response /. k->mw g),
      m2==((response /. k->mw g) /. s->2 s)};
fromData=Factor[g /. First[Solve[data,{mr,mw,g}]]];
emit["three-readout inversion", fromData /. m1->m0-r*(m0-m2),4*s^2*(1-r)/(4*r-1)];
omega2=k*(1/mr+1/mw);
slip=DSolveValue[{q''[t]+omega2 q[t]==f/mr,q[0]==0,q'[0]==0},q[t],t];
cm=f*t^2/(2 (mr+mw));
emit["normalized step gap",slip/cm,2*(1-Cos[Sqrt[omega2]t])/(t^2*k/mw)];
emit["one percent rate factor",Sqrt[1+1/(1/100)],Sqrt[101]];
Exit[0];
