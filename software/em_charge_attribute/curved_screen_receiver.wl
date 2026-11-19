(* The one-receiver calculation integrated from the source side. The
   finite probe is a uniform ball with unit total height-source strength. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, val_] := Print[label, " = ", ToString[val, InputForm]];
check[label_, actual_, expected_] := Module[{v=Factor[actual]},
  emit[label,v]; If[! TrueQ[FullSimplify[v == expected,Assumptions->d>b>0 && a>0]],Exit[1]]];

(* Upper-side oriented area element of the wall, outward from the fluid.
   Keep both other slopes and the unknown pressure variation in the test. *)
normal = {eps fx, eps fy, eps fz, -1};
tangent = {1, 0, 0, eps fx0};
check["wall pressure flux, derivative at zero",D[(p0+eps p1)normal.tangent,eps]/.eps->0,
  p0 (fx-fx0)];
(* The exact radial pressure tail also follows by the primitive's
   derivative and decay, independent of the Python integral. *)
weight[r_] := 2 md^2/(Pi^3 rho) (1/(3 r^3)-2 e a^2/(5 r^5)+e^2 a^4/(7 r^7));
pPlane[r_] := -md^2 (1-e a^2/r^2)^2/(2 Pi^4 rho r^6);
check["pressure tail derivative", D[weight[r],r],4 Pi r^2 pPlane[r]];
check["no boundary constant at infinity",Limit[weight[r],r->Infinity],0];
check["inverse-square coefficient",Limit[d^2 amp weight[d]/d^2,d->Infinity],0];
check["leading finite-cut coefficient",Limit[d^5 amp weight[d]/d^2,d->Infinity],
  2 amp md^2/(3 Pi^3 rho)];

(* At target-centric radius r, only a cap of the source ball is present.
   Weight each point charge by its center-gradient component times the
   intact-screen pressure tail outside its own distance from the target. *)
capCos[r_,d_,b_] := (r^2+d^2-b^2)/(2 r d);
capMeasure[r_,d_,b_] := 3 (1-capCos[r,d,b]^2)/(4 b^3);
integral = Integrate[Expand[capMeasure[r,d,b]],{r,d-b,d+b},Assumptions->d>b>0];
check["source cap weights sum to exact center gradient",integral,1/d^2];

Do[
  prev = 0;
  Do[
    measured = NIntegrate[Evaluate[capMeasure[r,dd,1] weight[r] /.
       {md->1,rho->1,a->1,e->shape}],{r,dd-1,dd+1},WorkingPrecision->35,
       AccuracyGoal->27,PrecisionGoal->22];
    scaled = measured dd^5 3 Pi^3/2;
    bounds = {weight[dd+1],weight[dd-1]} /. {md->1,rho->1,a->1,e->shape};
    If[!TrueQ[First[bounds] <= measured dd^2 <= Last[bounds]],Exit[1]];
    emit["shape " <> ToString[shape] <> ", d " <> ToString[dd] <> ", normalized d^5 J",N[scaled,14]];
    If[prev>0,emit["doubling exponent",N[Log[2,prev/measured],12]]]; prev=measured,
    {dd,{8,16,32,64}}];
  If[!TrueQ[Abs[scaled-1]<1/1000],Exit[1]],
  {shape,{0,1/4}}];
Exit[0];
