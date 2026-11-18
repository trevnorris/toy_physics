(* Measure carriage and Bernoulli pressure separately on finite target spheres.
   The target and source dipoles both use the laboratory +x axis. *)
ClearAll["Global`*"];
$HistoryLength=0;
Print["Wolfram version = ",$Version];
area[n_] := 2 Pi^(n/2)/Gamma[n/2];
fmt[v_] := ToString[N[v,12],InputForm];
vel[n_,dens_,q_,b_,center_] := Module[{ss=(x-center)^2+y^2,pot},
  pot = -q ss^((2-n)/2)/(dens area[n] (n-2)) + b (x-center) ss^(-n/2);
  {D[pot,x],D[pot,y]}];
measure[n_,dens_,q1_,q2_,targ_,src_,sep_,rad_] := Module[
  {v,v0,sphere,terms,normal={Cos[theta],Sin[theta]}},
  v = vel[n,dens,q1,targ,0]+vel[n,dens,q2,src,sep];
  v0 = vel[n,dens,q1,0,0];
  terms = dens {v[[1]] (v.normal)-v0[[1]] (v0.normal),
                -(v.v-v0.v0) Cos[theta]/2};
  sphere[expr_] := NIntegrate[Evaluate[
    area[n-1] rad^(n-1) Sin[theta]^(n-2) (expr /. {x->rad Cos[theta],y->rad Sin[theta]})],
    {theta,0,Pi},WorkingPrecision->43,PrecisionGoal->27,AccuracyGoal->30];
  sphere /@ terms];
Do[
  dens=n-1;
  sign=If[n==3,-1,1]; q1=sign 3/4; q2=sign 5/4;
  kt=1/20; ks=-1/30;
  Do[
    bt=kt/sep^2; bs=ks/sep^2;
    base=measure[n,dens,q1,q2,0,0,sep,rad];
    pair=measure[n,dens,q1,q2,bt,bs,sep,rad];
    difference=pair-base;
    monopole=-q1 q2/(dens area[n] sep^(n-1));
    finite=(n-1)(q2 bt-q1 bs)/sep^n+n(n-1) dens area[n] bt bs/sep^(n+1);
    err=Max[Abs[Total[base]-monopole],Abs[Total[difference]-finite]];
    Print["n=",n," d=",sep," R=",rad," full holder=",fmt[Total[pair]],
      " redistribution {carriage,pressure,total}=",fmt[Append[difference,Total[difference]]],
      " d^2 * carriage change=",fmt[sep^2 First[difference]],
      " d^(n+2) * total change=",fmt[sep^(n+2) Total[difference]],
      " max formula error=",fmt[err]];
    If[err>10^-25 || Abs[First[difference]-finite]<10^-11,Exit[1]],
    {sep,{8,16,32}},{rad,{1,3/2}}],
  {n,{3,4}}];
