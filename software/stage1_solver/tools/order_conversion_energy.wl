(* Time-path and virtual-variation checks of the same supplied historical
   order terms. No throat profile or new physical kinetics is introduced. *)
ClearAll["Global`*"];
$HistoryLength=0;
Print["Wolfram version = ",$Version];
emit[label_,value_] := Print[label," = ",ToString[Factor[value],InputForm]];
check[label_,lhs_,rhs_] := Module[{residual=FullSimplify[lhs-rhs]},
 Print[label," | ",ToString[lhs,InputForm]," | ",ToString[rhs,InputForm],
       " | residual ",ToString[residual,InputForm]];
 If[residual=!=0,Exit[1]]];

well[c_]:=aa c^2 (1-c)^2;
(* Homogeneous prescribed path, one full unit of order lost in time 1/gam.
   This is a kinematic witness, not a solution of the optional relaxation law. *)
path=1-gam t;
timePower=Factor[D[well[path],t]];
emit["same negative balance rate, power at chi=3/4",timePower/.t->1/(4 gam)];
emit["same negative balance rate, power at chi=1/4",timePower/.t->3/(4 gam)];
emit["bare well barrier above either endpoint",well[path]/.t->1/(2 gam)];
pathChange=Integrate[timePower,{t,0,1/gam},Assumptions->gam>0];
check["integrated bare order change",pathChange,well[0]-well[1]];

(* Use independently time-dependent order and shear variables before relating
   a path derivative to the partial order variation. *)
energy=well[c[t]]+c[t] ss[t]+gg[c[t]];
totalRate=D[energy,t];
orderRate=(D[well[z],z]/.z->c[t])+ss[t]+gg'[c[t]];
check["additional live shear-path work",totalRate-orderRate c'[t],c[t] ss'[t]];
emit["complete homogeneous endpoint drop with live shear and extra energy",
     (well[1]+sIn+gIn)-(well[0]+gOut)];

(* Virtual compression of a homogeneous parcel of fixed constituent count
   derives pressure; retain the common C_ref on both sides. *)
parcelEnergy=volume (kk (count/volume)^5/4+cref count/volume);
parcelPressure=-D[parcelEnergy,volume];
head=FullSimplify[(parcelEnergy+parcelPressure volume)/count];
headDensity=FullSimplify[head/.count->nn volume];
headDifference=FullSimplify[(headDensity/.nn->nf)-(headDensity/.nn->nb)];
emit["n=5 constituent enthalpy",headDensity];
emit["n=5 constituent enthalpy difference",headDifference];
check["common conserved-number reference cancellation",D[headDifference,cref],0];
check["equal-density bare constituent head difference",headDifference/.nb->nf,0];

(* Required net per-mass change, with the exterior cross energy as supplied
   benchmark input. No new value of the conversion energy is chosen. *)
kCross=da db/(4 Pi rho3 distance)+da db/(4 Pi^2 rho4 distance^2);
power=D[kCross,distance](vb-va);
yieldRules=First@Solve[{da dl1==power,db dl2==power},{dl1,dl2}];
emit["required pair-minus-alone net energy per mass, mouth 1",dl1/.yieldRules];
emit["required pair-minus-alone net energy per mass, mouth 2",dl2/.yieldRules];
variation=D[rate[t] yield[t],t];
emit["variation of conversion power",variation];
check["fixed rate and unchanged net energy yield",variation/.{rate'[t]->0,yield'[t]->0},0];
check["fixed rate, live net energy yield",variation/.rate'[t]->0,rate[t] yield'[t]];
