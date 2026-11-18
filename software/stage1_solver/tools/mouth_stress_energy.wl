(* Noether energy and current from the real/imaginary-field action using
   independent derivative jets, followed by the Euler-Lagrange equations.
   The calculation is the known neutral carrier sector only. *)
ClearAll["Global`*"];
$HistoryLength=0;
Print["Wolfram version = ",$Version];
emit[label_,value_] := Print[label," = ",ToString[Factor[value],InputForm]];
check[label_,lhs_,rhs_] := Module[{residual=FullSimplify[lhs-rhs]},
 Print[label," | ",ToString[lhs,InputForm]," | ",ToString[rhs,InputForm],
       " | residual ",ToString[residual,InputForm]];
 If[residual=!=0,Exit[1]]];

c=hb^2/(2 mass);
ax=Table[Unique["ax"],{4}]; bx=Table[Unique["bx"],{4}];
atx=Table[Unique["atx"],{4}]; btx=Table[Unique["btx"],{4}];
number=a^2+b^2;
lagrangian=hb (b at-a bt)-c (ax.ax+bx.bx)-vv number-uu[number];
energy=Expand[at D[lagrangian,at]+bt D[lagrangian,bt]-lagrangian];
current=Table[D[lagrangian,ax[[i]]] at+D[lagrangian,bx[[i]]] bt,{i,4}];
emit["neutral carrier energy density",energy];
emit["neutral carrier energy current, component 0",First[current]];
(* EL = partial L / partial field - dt(partial L / partial field_t)
   - sum_i di(partial L / partial field_i), evaluated on derivative jets. *)
eqa=D[lagrangian,a]-hb bt+2 c lapA;
eqb=D[lagrangian,b]+hb at+2 c lapB;
rates=First@Solve[{eqa==0,eqb==0},{at,bt}];
energyRate=D[energy,a] at+D[energy,b] bt+D[energy,vv] vt+
 Total[Table[D[energy,ax[[i]]] atx[[i]]+D[energy,bx[[i]]] btx[[i]],{i,4}]];
currentDivergence=-2 c (ax.atx+bx.btx+at lapA+bt lapB);
source=FullSimplify[(energyRate+currentDivergence)/.rates];
check["neutral carrier local energy source",source,number vt];

(* Transform the real/imaginary jets by the amplitude-angle Jacobian. *)
coordinates={Sqrt[nn] Cos[theta],Sqrt[nn] Sin[theta]};
timeJets=D[coordinates,nn] nt+D[coordinates,theta] thetat;
spaceJets=D[coordinates,nn] nx+D[coordinates,theta] thetax;
polarCurrent=FullSimplify[-2 c timeJets.spaceJets,Assumptions->nn>0];
emit["polar energy current, one component",polarCurrent];
check["stationary energy per number current",
 polarCurrent/.{nt->0,thetat->-mu/hb,thetax->mass vel/hb},mu nn vel];
emit["moving normal-cut energy current",polarCurrent-en wn];
emit["implemented held-feed power increment",jn (mu0+incident^2/2)-jn mu0];

(* Two material patches with different application velocities. *)
force=tone+ttwo;
power=tone uone+ttwo utwo;
check["stress power split",power,vv0 force+tone (uone-vv0)+ttwo (utwo-vv0)];
emit["zero-resultant witness, force",force/.{tone->tr,ttwo->-tr}];
emit["zero-resultant witness, power",power/.{tone->tr,ttwo->-tr,uone->us+sl,utwo->us-sl}];

(* Virtual translation differentiates energy with respect to the profile's
   coordinate; explicit time translation must give the opposite field power. *)
coupling=nn potential[x-xc];
reaction=-D[coupling,xc];
sourceMoving=D[coupling,xc] vcenter;
check["translated-profile energy transfer",sourceMoving,-reaction vcenter];
