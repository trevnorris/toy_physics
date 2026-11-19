(* Independently integrate the volume swept by the continuing unit normal. *)
ClearAll["Global`*"];
check[name_,actual_,expected_] := Module[{v=FullSimplify[actual]},
  Print[name," = ",ToString[v,InputForm]];
  If[!TrueQ[FullSimplify[v==expected]],Exit[1]]
];
er={Sin[h]Cos[f],Sin[h]Sin[f],Cos[h]};
normal[s_]:=Append[s Sin[p] er,Cos[p]];
volume[s_]:=Module[{n=normal[s]},
  FullSimplify[Det[{D[n,p],D[n,h],D[n,f],n}]]];
Do[
  rate[s]=FullSimplify[Integrate[volume[s],{h,0,Pi},{f,0,2Pi}]/Pi^2];
  check["signed normal-volume rate "<>ToString[s],rate[s], 4s Sin[p]^2/Pi];
  polar[s,alpha_]=Integrate[rate[s],{p,0,alpha}];
  check["continued cap, tube, flared end "<>ToString[s],
    polar[s,#]& /@ {0,Pi/2,Pi},{0,s,2s}],
  {s,{1,-1}}];
depth[b_,u_]:=3(1-u^2)(1+b u)/4;
weighted[b_,angle_] := Integrate[depth[b,u] polar[1,angle],{u,-1,1}];
check["vertical layer across different loss profiles",weighted[b,Pi/2],1];
check["symmetrically placed conversion",weighted[0,Pi/2+Pi u/6],1];
check["sensitivity at vertical",rate[1]/.p->Pi/2,4/Pi];
Do[
  reading=N[weighted[k,Pi/2+Pi u/6],24];
  Print["ordered turn at skew ",ToString[k,InputForm]," = ",
    ToString[NumberForm[reading,{15,12},NumberPadding->{"","0"}]]],
  {k,{-1/2,0,1/2}}];

(* A conserved 1+1 reduced current from a moving signed primitive. The
   pointwise conservation check tests the exchange caused by weighting. *)
potential[xx_,tt_]:=Exp[-(xx-2tt)^2];
ord=1/(1+Exp[x-t]);
cur={D[potential[x,t],x],-D[potential[x,t],t]};
check["unweighted moving current divergence",D[cur[[1]],t]+D[cur[[2]],x],0];
check["weighted divergence evaluated at x=t=1",
  (D[ord cur[[1]],t]+D[ord cur[[2]],x])/.{x->1,t->1},-1/(2E)];

(* Matter time-derivative term with a real phase. Probe modes enter only
   through the gauge potential. Compare a flat and a depth-weighted probe. *)
timeTerm=-n (hbar phaseRate+q0 test[x] eps);
stored=orderPot[chi]+stiff chi strain^2/2;
check["matter, order and shear variation", -D[timeTerm-stored,eps],n q0 test[x]];
(* Use smooth, equally normalized, displaced charge profiles rather than
   point samples. The constant Maxwell zero mode probes the unweighted
   integral of J^0. Translation in depth leaves it invariant. *)
g[w_,center_]:=Exp[-(w-center)^2]/Sqrt[Pi];
zero=Integrate[g[w,0]-g[w,l],{w,-Infinity,Infinity},Assumptions->Element[l,Reals]];
check["constant-mode projection of balanced depths",zero,0];
local=Integrate[Exp[-w^2](g[w,0]-g[w,l]),{w,-Infinity,Infinity},
  Assumptions->Element[l,Reals]];
check["depth-weighted diagnostic for the same profiles",local,(1-Exp[-l^2/2])/Sqrt[2]];
