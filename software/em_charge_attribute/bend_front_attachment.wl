(* Integrate the order balance in a window moving with the rim. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, val_] := Print[label, " = ", ToString[val, InputForm]];
check[label_, a_, b_] := Module[{res = FullSimplify[a-b, Assumptions -> ell>0 && c>0 && rho>0]},
  emit[label, {a,b,res}]; If[res =!= 0, Exit[1]]];
frontRel = len + (follow-1) z[t];
smooth = 1-10 u^3+15 u^4-6 u^5;
(* Both endpoints lie outside the layer. Width is measured from the center
   to an edge, and u ranges from 0 to 1 across the transition. *)
inventory = n (frontRel-width-left + 2 width Integrate[smooth,{u,0,1}]);
inFlux = n v + jIn; outFlux = jOut;
netProduction = D[inventory,t]+outFlux-inFlux;
effectiveLoss = Simplify[-netProduction+jOut-jIn];
emit["order inventory in window", inventory];
emit["signed production", netProduction];
emit["corrected conversion for bulk-fixed and comoving fronts",
     {effectiveLoss/.follow->0,effectiveLoss/.follow->1}];
check["conversion excess per density", Expand[effectiveLoss/n-v], -D[frontRel,t]];

(* Obtain endpoint loads from the axial translation current of a uniform
   transverse wave, independently of differentiating a mode frequency. *)
momentum = -rho D[field[x,t],t] D[field[x,t],x];
flux = rho (D[field[x,t],t]^2+c^2 D[field[x,t],x]^2)/2;
check["translation current divergence", Expand[D[momentum,t]+D[flux,x]],
  rho D[field[x,t],x](c^2 D[field[x,t],{x,2}]-D[field[x,t],{t,2}])];
wave[x_,t_] := amp Sin[Pi x/(2 ell)] Cos[c Pi t/(2 ell)];
period = 4 ell/c;
mean[g_] := Integrate[TrigReduce[g],{t,0,period},
                      Assumptions->ell>0 && c>0 && Element[x,Reals]]/period;
energyDensity = rho (D[wave[x,t],t]^2 + c^2 D[wave[x,t],x]^2)/2;
localFlux = FullSimplify[mean[flux /. field -> wave],Assumptions->ell>0 && c>0];
stored = Simplify[Integrate[mean[energyDensity],{x,0,ell}, Assumptions->ell>0 && c>0],
                  Assumptions->ell>0 && c>0];
forces = Factor[{-1,1} localFlux/stored];
emit["endpoint loads per wave energy", forces];
shear = mean[rho c^2 D[wave[x,t],x]^2/2];
moment = Integrate[x shear,{x,0,ell},Assumptions->ell>0 && c>0];
norm = Integrate[shear,{x,0,ell},Assumptions->ell>0 && c>0];
centroid = Factor[moment/norm/ell];
emit["strain energy centroid fraction from rim", centroid];
emit["strain centroid response", Expand[D[z[t]+centroid frontRel,t]/z'[t]]];
emit["work per energy and rim shift", Simplify[forces . {1,follow}]];
freq = (2 Pi/period) /. ell->frontRel;
check["frequency drift and independently integrated conversion",
  Factor[D[freq,t]/freq], Factor[(effectiveLoss/n-v)/frontRel]];
Exit[0];
