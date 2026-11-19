(* Construct the rim compliance in face coordinates from the exterior energy
   for specified endpoint heights, then solve for the heights under load. *)
ClearAll["Global`*"];
assume = ri > 0 && ell > 0 && cm > 0 && ct > 0 && r >= ri;
check[name_, expr_, target_] := Module[{v=FullSimplify[expr,assume]},
  Print[name," = ",ToString[v,InputForm]];
  If[!TrueQ[FullSimplify[v == target,assume]],Exit[1]]
];
(* Spherical Laplace/Helmholtz fields with arbitrary rim heights x,y. *)
mid = (x+y) ri/(2 r);
half = (x-y) ri Exp[-(r-ri)/ell]/(2 r);
(* On shell, integration by parts of the positive quadratic energy leaves
   the inner boundary term. Both exterior fields decay at infinity. *)
lap[q_] := D[r^2 D[q,r],r]/r^2;
check["bulk equations",{lap[mid],lap[half]-half/ell^2},{0,0}];
energy = -2 Pi ri^2 ((cm mid D[mid,r] + ct half D[half,r]) /. r->ri);
rim = First[Solve[{D[energy,x]==fu,D[energy,y]==fl},{x,y}]];
top = {fu->f,fl->0};
mirror = {fu->0,fl->-f};
faces = Factor /@ ({mid+half,mid-half} /. rim);
blend = (cm/ct) Exp[-(r-ri)/ell]/(1+ri/ell);
check["one-face displacements per force", (faces /. top)/f,
      {1+blend,1-blend}/(4 Pi cm r)];
check["reflection residual",(faces /. mirror) + Reverse[faces /. top],{0,0}];
check["weak half-load residual",(faces /. {fu->-f/2,fl->0})+(faces /. top)/2,{0,0}];
check["asymptotic slope times r squared",
      r^2 D[Mean[faces /. top],r],-f/(4 Pi cm)];
check["independent faces at zero width curvature, equal tensions",
      Limit[faces /. top /. ct->cm,ell->Infinity],{f/(2 Pi cm r),0}];
Print["illustrative opposite/loaded face ratios = ",
      ToString[N[Table[(Last[faces]/First[faces]) /. top /. ct->cm /. ell->ri /. r->n ri,
                        {n,1,3}],12],InputForm]];
Exit[0];
