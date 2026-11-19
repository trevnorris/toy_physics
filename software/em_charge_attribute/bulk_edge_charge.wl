(* Real-space displacement and a compact packet in the Gaussian Maxwell theory. *)
ClearAll["Global`*"];
assume = a > 0 && sig > 0 && lam > 0 && eps > 0 && b > a && Element[qq,Reals];
moment[p_] := Integrate[(1-z^2)^p,{z,-1,1}];
h[w_,b_] := (1-((w-b)/a)^2)^2/(a moment[2]); (* inside its support *)
zn = Integrate[Exp[-(w/lam)^2],{w,-Infinity,Infinity},Assumptions->lam>0];
coupled = Integrate[qq h[w,b]/Sqrt[zn],{w,b-a,b+a},Assumptions->assume];
pulse = Integrate[qq h[cut,u],{u,cut-a,cut+a},Assumptions->a>0];
h2 = Integrate[h[w,b]^2,{w,b-a,b+a},Assumptions->assume];
(* Radial displacement follows from enclosed 3D number. Integrate its energy
   by parts to avoid the squared error function. Surface terms vanish at both ends. *)
f[r_] := Exp[-r^2/(2 sig^2)]/(2 Pi sig^2)^(3/2);
enc = Integrate[4 Pi r^2 f[r],{r,0,s},Assumptions->s>0 && sig>0];
shell = FullSimplify[2 enc D[enc,s]/s, s>0 && sig>0];
fieldIntegral = Integrate[shell,{s,0,Infinity},Assumptions->sig>0]/(4 Pi);
scale = qq^2 fieldIntegral/(2 eps);
energies = FullSimplify[{scale/zn, scale h2 Exp[((b+a)/lam)^2]},assume];
ends = FullSimplify[{Limit[enc^2/s,s->0,Direction->"FromAbove",Assumptions->sig>0],
                     Limit[enc^2/s,s->Infinity,Assumptions->sig>0]},assume];
Print["normalized mode source and integrated passage = ",InputForm[{coupled,pulse}]];
Print["squared depth norm and radial surface terms = ",InputForm[{h2,ends}]];
Print["real-space field integral and energy bounds at positive depth = ",
      InputForm[{fieldIntegral,energies}]];
If[!TrueQ[FullSimplify[coupled Sqrt[zn]==pulse && pulse==qq &&
      ends=={0,0} && fieldIntegral>0,assume]],Exit[1]];
