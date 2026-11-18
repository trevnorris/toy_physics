(* Independent end-driven wall and Gaussian transverse packet calculations. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, expr_, ref_, assumptions_:True] := Module[{d},
  d=FullSimplify[expr-ref,assumptions];
  Print[label," computed=",InputForm[expr]," comparison=",InputForm[ref]," difference=",InputForm[d]];
  If[!TrueQ[d==0],Exit[1]]
];
(* Start at the cut with its force and displacement and integrate back to the rim.
   Here k=s/c, with mu=T/c^2. Normalize by the acceleration s^2 X. *)
cond= c>0 && b>0 && s>0 && ten>0 && x0>0;
k=s/c;
profile=end Cosh[k (b-w)] + flux Sinh[k (b-w)]/(ten k);
Do[
  load=If[j==0,0,ten k Tanh[k b]];
  rules=First[Solve[{(profile/.w->0)==x0,flux==load end},{end,flux}]];
  stored=Integrate[(ten/c^2) s profile,{w,0,b},Assumptions->cond] /. rules;
  cut=(flux/.rules);
  incoming=((-ten D[profile,w] /. w->0) /. rules);
  emit["wall balance", incoming-cut, s stored,cond];
  emit["wall slow storage",Limit[s stored/(s^2 x0),s->0],ten b/c^2,cond];
  emit["wall slow export",Limit[cut/(s^2 x0),s->0],j ten b/c^2,cond],
  {j,0,1}];

(* Project the field's translational current directly onto a Gaussian with
   center Y. Integrals are over one complete transverse section, per unit
   spectator area. At the sample the origin is at the packet center. *)
packet=f[z] Exp[-(x-y[z,t])^2/(2 a^2)] Cos[omega t];
ut=D[packet,t]; ux=D[packet,x]; uz=D[packet,z];
point={Derivative[0,1][y][z,t]->v,Derivative[1,0][y][z,t]->slope,y[z,t]->0};
mean[expr_]:=Expand[expr] /. {Cos[omega t]^2->1/2,Sin[omega t] Cos[omega t]->0};
integral[expr_]:=Integrate[expr,{x,-Infinity,Infinity},Assumptions->a>0];
norm=integral[D[Exp[-x^2/(2 a^2)],x]^2];
mom=integral[mean[(-rho ut ux)/.point]];
out=integral[mean[(mod uz ux)/.point]];
emit["wave momentum integrated across section",mom,rho norm f[z]^2 v/2,a>0];
emit["wave flux integrated across section",out,-mod norm f[z]^2 slope/2,a>0];
Do[
  prim[u_]:=Integrate[1-Cos[(2 n+1) Pi v0],{v0,0,u}];
  weight=(prim[1/4]-prim[0])/(prim[1]-prim[0]);
  Print["mode=",n," first-quarter ratio=",InputForm[weight]," numerical=",InputForm[N[weight,14]]],
  {n,0,1}];
Exit[0];
