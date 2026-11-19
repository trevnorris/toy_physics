(* Independent virtual-work calculation using face locations and a 4D cap. *)
ClearAll["Global`*"];
$HistoryLength=0;
emit[label_,v_]:=Print[label," = ",ToString[v,InputForm]];
check[label_,x_]:=Module[{v=FullSimplify[x,Assumptions->a>0 && ell>0 && k>0]},
  emit[label,v]; If[v=!=0,Exit[1]]];
facePlus=center+width/2+shift;
faceMinus=center-width/2+shift;
check["common displacement: thickness variation", D[facePlus-faceMinus,shift]];
check["common displacement: midpoint variation minus one",D[(facePlus+faceMinus)/2,shift]-1];
wallPull=4 Pi sig a^2+16 Pi bend/9;
near=facePlus;
far=far0+beta shift;
(* The near and far faces of the straight segment exert opposite tractions. *)
work=FullSimplify[wallPull D[near,shift]-wallPull D[far,shift]];
emit["positive throat wall force on the chosen motion",Factor[work]];
check["translation of a fixed shape",work /. beta->1];
check["fixed receiver point", (work /. beta->0)-wallPull];
(* At fixed radius a uniform bulk pressure jump pushes the hemispherical
   cap. Integrate its projected area in polar angle in four-space. *)
polarArea=4 Pi a^3 Sin[theta]^2;
cap[p_] := Integrate[polarArea Cos[theta] p,{theta,0,Pi/2}];
area=cap[1];
emit["cap axial area",area];
emit["slab-normal pressure needed to support the straight end",Expand[wallPull/area]];
check["cap and opposite flat cut of uniform pressure",cap[dp]-4 Pi a^3 dp/3];
(* A nonuniform pressure samples its projected mean, not its area mean. *)
emit["cap axial force, pressure p0+p1 cos(theta)",Expand[cap[p0+p1 Cos[theta]]]];
(* Both supports act on the cap. Balance their forces first; use the
   opposite force on the slab for its attachment, plus the near wall pull. *)
tSupport=k (y-z-rest);
ysol=y/.First[Solve[bulk-s P-tSupport==0,y]];
emit["displacement response of the far attachment",D[ysol,z]];
emit["complete slab reaction",Simplify[(s P+tSupport)/.y->ysol]];
check["equal positions with compensating reference",
      (ysol /. rest->(s P-bulk)/k)-z];
(* The scalar is sampled at the moving mouth plane. A jointly carried
   scalar profile and evaluation point give a constant, even with a slope. *)
sample=prof[w-shift];
check["jointly carried parent scalar value variation",D[sample/.w->shift,shift]];
e0=kM h0^2/2-s g h0;
e1=kM (h0+c delta/ell)^2/2-s g (h0+c delta/ell);
emit["mouth load under specified overlap",Factor[-Coefficient[Expand[e1-e0],delta,1]]];
Exit[0];
