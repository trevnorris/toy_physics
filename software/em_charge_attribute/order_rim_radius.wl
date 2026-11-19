(* Independent parametric area variation. Both coordinates vary before
   specializing to a segment tangent to the w axis. H is mean curvature. *)
ClearAll["Global`*"];
$HistoryLength=0;
emit[label_,v_]:=Print[label," = ",ToString[v,InputForm]];
check[label_,err_]:=Module[{res=FullSimplify[err,Assumptions->a>0 && len>0]},
  emit[label,res]; If[res=!=0,Exit[1]]];
speed=Sqrt[r1^2+z1^2];
meridian=(r1 z2-z1 r2)/speed^3;
angular=z1/(rad speed);
integrand=4 Pi rad^2 speed (sig+kb ((meridian+2 angular)/3)^2);
dt[f_]:=D[f,rad] r1+D[f,r1] r2+D[f,z1] z2+D[f,r2] r3+D[f,z2] z3;
pz=D[integrand,z1]-dt[D[integrand,z2]];
atLip[orientation_]:=Factor[pz /. {rad->a,r1->0,z1->orientation,
                                    z2->0,z3->0}];
plus=atLip[1]; minus=atLip[-1];
emit["pull on near end toward +w at axial tangent",plus];
check["mirror sum",plus+minus];
(* Translate the entire parametrized surface and both its ends together. *)
shifted=integrand/.{z1->D[height[t]+offset,t],z2->D[height[t]+offset,{t,2}]};
check["whole-sheet rigid translation work",D[shifted,offset]];
straight=Simplify[plus /. r2->0];
emit["straight segment desired separation force",straight];
check["cylinder end work vs direct extension",
      straight-D[4 Pi a^2 len(sig+kb (2/(3 a))^2),len]];
segment=4 Pi a^2 (far-near)(sig+kb (2/(3 a))^2);
check["cylinder near and far reactions including signs",
      (D[segment,near]+straight)^2+(D[segment,far]-straight)^2];
(* Historical reference radius varies as 1+smoothstep/2 in a length
   37/20. Compare only the bending contribution at the mouth. *)
rtest=a(1+(3 (z/len)^2-2 (z/len)^3)/2);
qb=D[rtest,{z,2}]/.z->0;
flare=Factor[(plus/.sig->0)/.r2->qb];
emit["flare fixed-length end bending",flare];
emit["similar flare / cylinder bending",
     Factor[(flare/(straight/.sig->0))/.len->37 a/20]];
(* Direct spatial quadrature of a finite symmetric portion of the actual
   curved layer. kappa=ell=1, planar tension=1/6. Subtract the leading
   rim area term. The omitted kink tails shrink exponentially with a. *)
weight[d_]:=1/(16 Cosh[d/2]^4);
second=6 NIntegrate[d^2 weight[d],{d,-Infinity,Infinity},
                    WorkingPrecision->40,AccuracyGoal->25];
emit["second planar moment divided by tension and width squared",N[second,18]];
check["planar norm",Integrate[weight[d],{d,-Infinity,Infinity}]-1/6];
errors=Table[
  excess=6 NIntegrate[(aa+d)^2 weight[d],{d,-aa/2,0,aa/2},
                      WorkingPrecision->40,AccuracyGoal->24]-aa^2;
  emit["curved trial after subtracting area, radius "<>ToString[aa],N[excess,16]];
  Abs[excess-second], {aa,{20,40,80}}];
If[!TrueQ[errors[[1]]>errors[[2]]>errors[[3]] && errors[[3]]<10^-20 &&
           Abs[second-(Pi^2/3-2)]<10^-20],Exit[1]];
Exit[0];
