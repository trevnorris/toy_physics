(* Held controls keep beta fixed. t labels any monotonic unwinding history. *)
ClearAll["Global`*"];
check[label_, expr_, expected_] := Module[{v = Factor[expr]},
  Print[label, " = ", InputForm[v]];
  If[!TrueQ[FullSimplify[v == expected]], Exit[1]]
];

(* Physical density and one spatial current component, varied before the
   probe is removed. Lother contains the displayed probe-independent terms. *)
matter = -n (hbar phaseT + q a0) - n (hbar phaseN-q an)^2/(2 m) + lOther;
check["electric density from temporal probe", -D[matter,a0], q n];
flow = D[matter,an] /. an->0;
check["electric current at zero spatial probe", flow, q n hbar phaseN/m];
check["current relative to cut at speed vc",
      (flow + vc D[matter,a0]) /. phaseN->m vn/hbar, q n (vn-vc)];

(* Count constituents directly from the material chart and their crossing
   time. A fixed-r sphere moves in w as beta varies. At the instant of
   comparison the held and moving embeddings coincide. *)
height = eps (r^2/(2 radius)-b[t] r Cos[theta]);
tangents = {{1,0,0,D[height,r]},{0,r,0,D[height,theta]},{0,0,r Sin[theta],0}};
metricDet = Factor[Det[tangents.Transpose[tangents]]];
areaFactor2 = 1 + D[height,r]^2 + (D[height,theta]/r)^2;
check["volume Jacobian squared per r^4 sin^2(theta)",
      metricDet/(r^4 Sin[theta]^2), areaFactor2];
jac = r^2 Sin[theta] Sqrt[areaFactor2];
physicalDensity = n0/Sqrt[areaFactor2];
chargePerR = Integrate[q physicalDensity jac, {theta,0,Pi}]*2 Pi;
chartSpeed = intake/(4 Pi n0 r^2);
check["moving charge per unit radius", chargePerR, 4 Pi q n0 r^2];
check["collar inventory difference at the same tilt",
      Integrate[chargePerR-(chargePerR /. b[t]->b0),{r,2 radius,3 radius}],0];
check["charge crossing the bulk-side edge in the moving chart",
      (chargePerR chartSpeed /. r->3 radius),q intake];

(* Measure the independently completed geometric current by inventory
   upstream of an interior cut c before and after the zero passes it.
   A held shape has constant upstream inventory and no crossing dose. *)
h[y_] := 30 (y-2)^2 (3-y)^2;
layerBefore = Integrate[h[y], {y,b,c}];
pointBefore = Integrate[h[y], {y,2,b}];
(* After b>c there are no regular degrees or point in the upstream slab. *)
upBefore = Factor[layerBefore+pointBefore];
check["completed upstream inventory independent of point position",
      D[upBefore,b],0];
check["oriented geometric dose at two interior cuts and in edge limit",
      s {upBefore /. c->5/2,upBefore /. c->8/3,Limit[upBefore,c->3]},
      s {1/2,64/81,1}];

(* Sliding held snapshots have the same instantaneous geometric inventory,
   yet the rate of each individually held member is zero. *)
storedLayer = s Integrate[h[y],{y,b[t],3}];
check["matched spatial count difference",
      storedLayer - ((storedLayer /. b[t]->b0) /. b0->b[t]),0];
check["rate difference at mid-layer",
      (D[storedLayer,t]-D[storedLayer /. b[t]->b0,t]) /. b[t]->5/2,
      -15 s b'[t]/8];
