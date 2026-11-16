(* Check the acceleration-odd extraction for a periodic linear response,
   and the relation between the whole-plane mixed flux and the boundary
   form used in the independent carrier momentum inventory. *)
vel = v0 + eps Cos[w t];
fieldObservable = f0 + eps (fc Cos[w t] + fs Sin[w t]);
storedMomentum = p0 + eps (pc Cos[w t] + ps Sin[w t]);
atCrossings[h_] := FullSimplify[{h /. t -> 3 Pi/(2 w), h /. t -> Pi/(2 w)},
                                  Assumptions -> w > 0];
oddPerA[h_] := Simplify[(Subtract @@ atCrossings[h])/(Subtract @@ atCrossings[D[vel,t]])];
Print["crossing speeds = ", InputForm[atCrossings[vel]]];
Print["crossing accelerations = ", InputForm[atCrossings[D[vel,t]]]];
Print["flux contrast / acceleration contrast = ", InputForm[oddPerA[fieldObservable]]];
Print["stored momentum rate contrast / acceleration contrast = ",
      InputForm[oddPerA[D[storedMomentum,t]]]];

(* m=hbar=1, density exp[L], phase T. The derivative remainder integrates
   to zero at a periodic boundary or around a rim with zero normal slope. *)
density = Exp[logn[x,z]];
adv = density D[phase[x,z],x] D[phase[x,z],z];
stress = -density D[logn[x,z],x,z]/4;
gradientProduct = adv + D[density,x] D[logn[x,z],z]/4;
surfaceIdentity = Simplify[adv+stress-gradientProduct+D[density D[logn[x,z],z],x]/4];
Print["mixed momentum flux identity residual = ", InputForm[surfaceIdentity]];
If[Simplify[oddPerA[fieldObservable]+fs/w] =!= 0 ||
   Simplify[oddPerA[D[storedMomentum,t]]-pc] =!= 0 || surfaceIdentity =!= 0, Exit[1]];
