(* Integrate the static midplane Green function over each spherical load.
   Substitution from polar cosine to chord length gives the radial bounds. *)
ClearAll["Global`*"];
check[label_, expr_, target_, assumptions_] := Module[{value=FullSimplify[expr,assumptions]},
  Print[label," = ",ToString[value,InputForm]];
  If[!TrueQ[FullSimplify[value==target,assumptions]],Exit[1]]
];
shell[r_,radius_] := Integrate[1/(2 r radius),{chord,Abs[r-radius],r+radius},
                               Assumptions->r>0 && radius>0];
h=(push shell[r,a]-back shell[r,b])/(4 Pi kap);
inner=FullSimplify[h,0<a<r<b && kap>0];
outer=FullSimplify[h,0<a<b<r && kap>0];
flux[q_]:=-4 Pi kap r^2 D[q,r];
check["inside return surface",inner,(push/r-back/b)/(4 Pi kap),0<a<r<b && kap>0];
check["outside return surface",outer,(push-back)/(4 Pi kap r),0<a<b<r && kap>0];
check["oriented slope fluxes",{flux[inner],flux[outer]},{push,push-back},r>0 && kap>0];
check["complete return far bend",outer /. back->push,0,0<a<b<r && kap>0];
(* Cartesian carrier trace and its polar density-stress form. *)
n=xx[w]^2+yy[w]^2;
j=xx[w] yy'[w]-yy[w] xx'[w];
cart=(xx'[w]^2+yy'[w]^2-xx[w] xx''[w]-yy[w] yy''[w])/2;
check["normal matter-stress residual",cart-j^2/n+n D[Log[n],{w,2}]/4,0,n>0];
Exit[0];
