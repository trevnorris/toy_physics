(* Quasistatic translation work, bulk stress and two physical mouth sizes. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[label_, value_] := Print[label, " = ", ToString[value, InputForm]];
check[label_, actual_, expected_] := Module[{v = FullSimplify[actual]},
  emit[label, v]; If[! TrueQ[FullSimplify[v == expected]], Exit[1]]];

(* Exterior traction for an applied field with a general direction. *)
normal = {Sin[t] Cos[p],Sin[t] Sin[p],Cos[t]};
incident = {gx,gy,gz};
own = -mass normal/(4 Pi k r^2);
cross[n_] := k (own (incident.n) + incident (own.n) - (own.incident)n);
inward = -r^2 Integrate[TrigReduce[Expand[cross[normal] Sin[t]]],
   {t,0,Pi},{p,0,2 Pi}];
check["three-component scalar input",inward,mass incident];

(* In the material tangent frame the signed angle between the flow and
   face normal is (1-b) times the local tilt, to first order. The two
   cuts have the same directed thrust upon subtracting ambient pressure. *)
turn = sign thrust (Sin[(1-bi) angle]-Sin[(1-bo) angle]);
check["intrinsic tangent-force slope per angle",D[turn,angle]/.angle->0,
      sign thrust (bo-bi)];
vGuide = {ux,uw}; vi = {vix,viw}; vo = {vox,vow};
arrives = mdot ((vGuide+vi).(vGuide+vi) - (vGuide+vo).(vGuide+vo))/2;
check["moving control, heat-free residual allowing unequal relative speeds",
  arrives-mdot (vi-vo).vGuide-mdot (vi.vi-vo.vo)/2,0];

(* Direct average of the incident gradient over a three-dimensional
   sampling shell. Change variables to distance t from the source. *)
avg = Integrate[-(d^2-b^2+t^2)/(4 d^2 b t^2),{t,d-b,d+b},
    Assumptions-> d>b>0];
check["sampled gradient on a concentric shell",avg,-1/d^2];
check["asymmetric sampling leading differential",Limit[d^3(-1/(d+c)^2+1/d^2),d->Infinity],2 c];

(* Hemispherical bulk cut. w is the polar axis, x is tangent. Integrate
   the full 4D tensor angular weights instead of supplying their moments. *)
mu = Sin[theta] Cos[phi];
hem[expr_] := Integrate[2 Pi Sin[theta]^2 Sin[phi] expr,
  {theta,0,Pi/2},{phi,0,Pi}];
phiPert = dip x/(x^2+y^2+z^2+w^2)^2;
onCut = {x->r mu,y->r Sin[theta] Sin[phi],z->0,w->r Cos[theta]};
vx = Simplify[D[phiPert,x] /. onCut,Assumptions->r>0];
vr = D[dip mu/r^3,r]; baseV = q/r^3;
ad = hem[r^3 rho baseV (vx+mu vr)];
pr = hem[-r^3 rho baseV mu vr];
check["integrated carried horizontal momentum",ad,-3 Pi^2 rho q dip/(4 r^4)];
check["integrated pressure momentum",pr,3 Pi^2 rho q dip/(4 r^4)];
check["total decaying-dipole flux",ad+pr,0];
check["uniform bulk flow per unit incident speed",hem[r^3 rho baseV],Pi^2 rho q];
check["hypothetical independent pressure coefficient",hem[mu^2],Pi^2/4];

(* Three pairings of physical radii, including both measurement orderings
   of the mixed pair. The target responds independently of its source. *)
ordered = -Outer[Times,{b1,b2},{l1,l2}];
check["equal diagonal forces, mixed-force ratio",
  (ordered[[1,2]]/ordered[[2,1]]) /. b2 -> b1 l1/l2,l2^2/l1^2];
check["reciprocity residual with common ratio", (ordered[[1,2]]-ordered[[2,1]]) /.
  {b1->c l1,b2->c l2},0];
Exit[0];
