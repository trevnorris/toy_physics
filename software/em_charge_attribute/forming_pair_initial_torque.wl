(* Initial angular impulse rate about the energy centroid from the force on
   the confining potential's wave field and the stress through the midplane. *)
ClearAll["Global`*"];
shape = x + tilt (w^2-1/2);
pot = x^2 + w^4/4;
dx = D[shape,x]-x shape; dw = D[shape,w]-w shape;
en = (dx^2+dw^2+pot shape^2)/2;
(* Strip the common Gaussian and first integrate the full x direction. *)
alongX[p_] := With[{coeffs=CoefficientList[Expand[p],x]},
  Sum[coeffs[[j+1]] Gamma[(j+1)/2],{j,0,Length[coeffs]-1,2}]];
overHalf[p_,side_] := With[{coeffs=CoefficientList[Expand[alongX[p]/.w->side w],w]},
  Sum[coeffs[[j+1]] Gamma[(j+1)/2]/2,{j,0,Length[coeffs]-1}]];
Do[
  energy=overHalf[en,side]; xc=overHalf[x en,side]/energy; wc=overHalf[w en,side]/energy;
  trap=-overHalf[shape^2 ((x-xc)D[pot,w]-(w-wc)D[pot,x]),side]/2;
  cut=side alongX[Expand[((x-xc)(dw^2-en)+wc dw dx)/.w->0]];
  spinTorque=Factor[trap+cut];
  Print["side ",side," linear coefficient = ", ToString[Factor[D[spinTorque,tilt]/.tilt->0],InputForm]];
  Do[Print["side ",side," tilt ",val," initial intrinsic torque = ",
           ToString[N[spinTorque/.tilt->val,16],InputForm]],{val,{1/100,-1/100,0}}],
  {side,{1,-1}}];
