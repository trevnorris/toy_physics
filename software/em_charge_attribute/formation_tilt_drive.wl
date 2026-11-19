(* Mixed stress flux over a round cut, then an independent Hermite-basis
   calculation of the confinement forces and instantaneous wave acceleration. *)
ClearAll["Global`*"];
emit[name_, v_] := Print[name, " = ", ToString[v, InputForm]];
angular[f_] := Integrate[Integrate[f Sin[th],{ph,0,2 Pi}],{th,0,Pi}]/(4 Pi);
ex=Sin[th] Cos[ph];
(* First assemble the traction from the stress tensor in a radial-depth frame. *)
traction={{rr,rw},{rw,ww}}.{nr,nw}-vc {pr,pw};
local=Det[{{r ex,wrel},{ex traction[[1]],traction[[2]]}}];
emit["full angular mean of moving-cut flux",angular[local]];
emit["angular first moment coefficient",Factor[angular[ex local]/(local/ex)]];
emit["carrier and parent-wave lateral force factors",angular /@ {ex,ex^3}];
If[angular[local]=!=0,Exit[1]];

(* The highest application here is w^2 on |2>, so the oscillator basis
   covers every intermediate state exactly. Derivative norms give kinetics. *)
dim=7;
lower=SparseArray[Table[{j,j+1}->Sqrt[j],{j,1,dim-1}],{dim,dim}];
ident=IdentityMatrix[dim,SparseArray];
coord=(lower+Transpose[lower])/Sqrt[2];
deriv=(lower-Transpose[lower])/Sqrt[2];
xx=KroneckerProduct[coord,ident]; wwop=KroneckerProduct[ident,coord];
dx=KroneckerProduct[deriv,ident]; dw=KroneckerProduct[ident,deriv];
state[n_,m_]:=UnitVector[dim^2,n dim+m+1];
parent=state[1,0]; target=state[0,2];
components[v_] := Module[{wv=wwop.v,w2=wwop.(wwop.v)},
  {dx.v,dw.v,xx.v-b0 v-b1 wv-b2(w2-v/2),(w2-depth^2 v)/2}];
overlap[a_,b_] := Expand[Total[MapThread[Dot,{a,b}]]];
cp=components[parent]; ct=components[target];
mix=overlap[ct,cp]; emit["target acceleration for general centerline",-mix];
trial=Cos[beta] cp+Sin[beta] ct;
energy=overlap[trial,trial]/2+stiff b2^2/2 /. {b0->0,b1->0};
orig={b2->0,beta->0};
emit["energy slopes at the symmetric configuration",(D[energy,#]& /@ {b2,beta}) /. orig];
hess=Simplify[Table[D[energy,i,j],{i,{b2,beta}},{j,{b2,beta}}]/.orig];
emit["projected energy Hessian",hess];
tiltRelaxed=Factor[Det[hess]/hess[[1,1]]];
emit["relaxed curvature",tiltRelaxed];
emit["held curvatures at first and last geometry",Table[hess[[2,2]]/.depth->value,{value,{0,2}}]];
contrast=Simplify[energy-(energy/.{b2->-b2,beta->-beta})];
emit["energy contrast after reflection",contrast];
If[!TrueQ[Simplify[-mix==b2 && contrast==0]],Exit[1]];
