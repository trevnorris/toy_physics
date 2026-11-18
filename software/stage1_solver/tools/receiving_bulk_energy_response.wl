(* Independent finite-radius quadrature from primitive Cartesian fields. *)
ClearAll["Global`*"];
$HistoryLength = 0;
Print["Wolfram version = ", $Version];
rho = 1; q1 = 2; q2 = 3; d0 = 4; omega = 2 Pi^2;
k1 = q1/(rho omega);
sq1 = x^2+y^2; sq2 = (x-sep)^2+y^2;
(* The pair-only amplitude goes to zero in the isolated limit. *)
amplitude = kk (sep-d0) (d0/sep)^2;
potential = -i1 q1/(2 rho omega sq1)-i2 q2/(2 rho omega sq2) +
  i1 i2 amplitude x/sq1^2;
{vx, vy} = D[potential, #]& /@ {x,y};
dt = -speed D[potential, sep];
energy = rho (vx^2+vy^2)/2;
pressure = 2/5-energy-rho dt;
vn = vx Cos[th]+vy Sin[th];
(* Normal points out of the region enclosed by the cut. *)
flux = (energy+pressure+1/7) vn;
restForce = -(rho vx vn+pressure Cos[th]);
fourStates = {{1,1,1},{-1,1,0},{-1,0,1},{1,0,0}};
interaction[expr_] := Total[Map[
  Function[row, row[[1]] (expr /. {i1->row[[2]],i2->row[[3]]})], fourStates]];
intPower = Factor[interaction[flux] /. sep->d0];
intForce = Factor[interaction[restForce] /. {sep->d0,speed->0}];
surface[expr_, center_, radius_] := NIntegrate[
  Evaluate[4 Pi radius^3 Sin[th]^2 (expr /.
    {x->center+radius Cos[th],y->radius Sin[th]})],
  {th,0,Pi}, WorkingPrecision->36,AccuracyGoal->23,PrecisionGoal->23];
fmt[v_] := ToString[N[v,12],InputForm];
require[ok_] := If[!TrueQ[ok], Exit[1]];
expected = q1 q2/(rho omega d0^3);
Do[
  f = surface[intForce /. kk->k, 0, rad];
  powers = Table[surface[intPower /. {kk->k,speed->v}, 0, rad],
                 {v,{0,1/10,-1/10,1/20,-1/20}}];
  full = (powers[[2]]-powers[[3]])/(2/10);
  half = (powers[[4]]-powers[[5]])/(2/20);
  Print["R=",rad," k/k1=",k/k1," static force=",fmt[f],
    " powers(0,+,-,+half,-half)=",fmt[powers],
    " odd slopes=",fmt[{full,half}]," slope / static force=",fmt[full/f]];
  target = expected-k q2/(4 d0^3);
  require[Max[Abs[{f-expected,full-target,half-target,powers[[1]],
    powers[[2]]+powers[[3]],powers[[4]]+powers[[5]]}]] < 10^-21],
  {rad,{1/2,1}},{k,{0,k1}}];
(* Verify conserved rates away from the comparison event, with active dipole. *)
Do[
  massFlux = rho vn /. {i1->1,i2->1,sep->ds,kk->k1};
  rates = {surface[massFlux,0,1],surface[massFlux,ds,3/5]};
  Print["d=",ds," two outward mass rates=",fmt[rates]];
  require[Max[Abs[rates-{q1,q2}]] < 10^-21],{ds,{39/10,41/10}}];
(* Fixed offset in pressure or internal energy cancels in the original four states.
   Removing the time-dependent pressure on the fixed cut kills the odd slope. *)
steadyFlux = Factor[interaction[(energy+pressure+rho dt+1/7) vn] /. sep->d0];
removed = surface[steadyFlux /. {kk->k1,speed->1/10},0,1];
Print["interaction power after removing unsteady pressure=",fmt[removed]];
require[Abs[removed] < 10^-21];
