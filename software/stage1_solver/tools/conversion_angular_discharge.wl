(* Independently check the S12 trace freedom, then measure held finite-cut
   flow directly for two choices of the decaying exterior response. *)
ClearAll["Global`*"];
$HistoryLength=0;
Print["Wolfram version = ",$Version];
exact[label_,actual_,expected_] := Module[{a=FullSimplify[actual],r},
  r=FullSimplify[a-expected];
  Print[label," calculated=",ToString[a,InputForm]," expected=",ToString[expected,InputForm],
    " residual=",ToString[r,InputForm]];
  If[r=!=0,Exit[1]]];
(* Upper graph face: currents per projected area, oriented from slab to bulk.
   C is a coefficient of the order source, while pIn/pOut are constituent rates. *)
upperMaterial = pOut-pIn;
upperOrder = -f pIn+bOut-bIn;
change={pIn->pIn+eta,pOut->pOut+eta,bIn->bIn-f eta};
exact["upper material jump change",(upperMaterial/.change)-upperMaterial,0];
exact["upper order source change",(upperOrder/.change)-upperOrder,0];
exact["outward material change",(pOut/.change)-pOut,eta];
exact["order source change if bIn fixed",
    (upperOrder/.pIn->pIn+eta)-upperOrder,-f eta];

(* Both unit mouths emit the same mass, rho=1 and q_i=2 Pi^2.
   d is a preparation label. All fields are time independent in each solve. *)
d0=4;
amp=k (d-d0) (d0/d)^2;
sq=x^2+y^2;
phi=-1/(2 sq)-1/(2 ((x-d)^2+y^2))+amp x/sq^2;
vn=(D[phi,x] Cos[t]+D[phi,y] Sin[t])/.{x->r Cos[t],y->r Sin[t]};
surface[expr_,lo_,hi_] := NIntegrate[Evaluate[4 Pi Sin[t]^2 expr],
    {t,lo,hi},WorkingPrecision->35,AccuracyGoal->22,PrecisionGoal->22];
moment[ds_,rad_,slope_] := 4/(2 Pi^2) surface[
   Cos[t] (vn/.{r->rad,d->ds,k->slope}),0,Pi];
extract[ds_,slope_] := With[{ri=6/5,ro=8/5},
   (moment[ds,ro,slope]-moment[ds,ri,slope])/(3(ri^-4-ro^-4))];
fmt[v_] := ToString[N[v,12],InputForm];
Do[
  value=extract[ds,kk];
  expected=amp/.{d->ds,k->kk};
  dv=(vn/.{r->6/5,d->ds,k->kk})-(vn/.{r->6/5,d->ds,k->0});
  halves=(6/5)^3 {surface[dv,0,Pi/2],surface[dv,Pi/2,Pi]};
  Print["d=",ds," k=",kk," two shell B=",fmt[value],
    " changes in facing and away mass rates=",fmt[halves]];
  If[Max[Abs[{value-expected,Total[halves],First[halves]+4 Pi expected/(6/5)}]]>10^-19,
    Exit[1]],
 {kk,{0,1/20}},{ds,{39/10,4,41/10}}];
exact["midpoint derivative of admissible family",D[amp,d]/.d->d0,k];
