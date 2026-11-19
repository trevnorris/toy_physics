(* Calculate the closing-plane pressure independently in cylindrical coordinates. *)
ClearAll["Global`*"];
$HistoryLength=0;
emit[label_,v_]:=Print[label," = ",ToString[v,InputForm]];
verify[label_,x_]:=Module[{res=FullSimplify[x, Assumptions->
    rho>0 && rad>0 && a>0 && xrad>0 && n0>0 && m>0 && k>0 && xi>0]},
  emit[label,res]; If[res=!=0,Exit[1]]];
rsq=xrad^2+w^2;
pot=-rate/(2 Pi^2 rho rsq)-eps rate a^2 (4 w^2/rsq-1)/(4 Pi^2 rho rsq^2);
vx=D[pot,xrad]; vw=D[pot,w];
verify["4D source-free equation",D[vx,xrad]+2 vx/xrad+D[vw,w]];
verify["flat-plane no-through condition",vw/.w->0];
point={xrad->a Sin[t],w->a Cos[t]};
normal=FullSimplify[(vx Sin[t]+vw Cos[t])/.point,Assumptions->0<t<Pi/2&&a>0];
emit["source normal velocity",normal];
verify["mass rate on emitting cut",Integrate[4 Pi a^3 rho Sin[t]^2 normal,{t,0,Pi/2}]-rate];
(* No momentum arrives from infinity for these decaying harmonic modes.
   Momentum across the curved cut therefore equals the pressure force on
   the exterior part of the intact plane, with the uniform bath removed. *)
planeP=FullSimplify[-rho (vx^2+vw^2)/2/.w->0];
plane=Integrate[-4 Pi xrad^2 planeP,{xrad,rad,Infinity},Assumptions->rad>0];
emit["export inferred from exterior plane",Collect[Expand[plane],eps]];
verify["far cut",Limit[plane,rad->Infinity]];
emit["shaped to uniform inner-cut ratio",
  Factor[(plane/.{rad->a,eps->1/4})/(plane/.{rad->a,eps->0})]];

(* Now keep the neutral carrier amplitude and its gradient energy. Order lambda
   counts the weak speed. The density and pressure first change at order two. *)
amp=Sqrt[n0](1+lambda^2 dn[rad]/(2 n0)); den=amp^2;
vel=lambda rate/(m Pi^2 rad^3 den);
coef[expr_]:=FullSimplify[SeriesCoefficient[expr,{lambda,0,2}],Assumptions->n0>0];
head=coef[hp (den-n0)+m vel^2/2-2 k (D[amp,{rad,2}]+3 D[amp,rad]/rad)/amp];
rr=coef[m den vel^2+n0 hp (den-n0)+2 k(D[amp,rad]^2-amp D[amp,{rad,2}])];
curv=First[Solve[head==0,dn''[rad]]];
projected=Integrate[4 Pi rad^3 Sin[t]^2 Cos[t],{t,0,Pi/2}];
flux=Expand[projected (rr/.curv)];
emit["amplitude equation hemisphere flux",flux];
mode=A BesselK[1,rad/xi]/rad;
verify["equal-rate density mode",
  FunctionExpand[D[mode,{rad,2}]+3 D[mode,rad]/rad-mode/xi^2]];
change=Coefficient[flux,dn'[rad]]D[mode,rad];
verify["equal-rate flux shift",FunctionExpand[change+4 Pi k A rad BesselK[2,rad/xi]/xi]];
(* Finite-cut witness in units k=A=xi=1. A scales with lambda^2 in the density. *)
emit["normalized shift at R/xi=1",N[-4 Pi BesselK[2,1],14]];
Exit[0];
