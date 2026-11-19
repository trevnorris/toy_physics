(* Independent pullback computation from the unit normal in four dimensions. *)
ClearAll["Global`*"];
check[name_,actual_,expected_] := Module[{result=FullSimplify[actual]},
  Print[name, " = ", ToString[result, InputForm]];
  If[!TrueQ[FullSimplify[result==expected]], Exit[1]]
];
v={v1,v2,v3};
dir=Append[v,1]/Sqrt[1+v.v];
(* Use slope coordinates with positive density for an increasing gradient.
   Rotation in the three brane directions allows evaluation along one axis. *)
frame = Append[(D[dir,#]& /@ v),dir] /. {v1->0,v2->0,v3->u};
check["oriented normal volume Jacobian at axial slope",
  FullSimplify[Det[frame],Element[u,Reals]],1/(1+u^2)^2];
(* The unit normal at the lip traces the latitude t=arctan|v| in S^3.
   Angular derivatives of the latitude drop out of the boundary two-form. *)
n={Sin[t] Sin[b] Cos[f],Sin[t] Sin[b] Sin[f],Sin[t] Cos[b],Cos[t]};
tangents=(D[n,#]& /@ {t,b,f}) /. f->0;
area=FullSimplify[Sqrt[Det[tangents.Transpose[tangents]]],0<b<Pi && 0<t<Pi/2];
cap[angle_]=FullSimplify[Integrate[area,{b,0,Pi},{t,0,angle}] 2/Pi,
  0<=angle<=Pi/2];
check["full vertical normal cap",cap[Pi/2],1];
check["finite 45 and 60 degree caps",{cap[Pi/4],cap[Pi/3]},
  {1/2-1/Pi,2/3-Sqrt[3]/(2Pi)}];
er={Sin[b] Cos[f],Sin[b] Sin[f],Cos[b]};
lip[s_] := Module[{radial=-s er,solid},
  solid=TrigExpand[radial.Cross[D[radial,b],D[radial,f]]]//FullSimplify;
  -cap[Pi/2] Integrate[solid,{b,0,Pi},{f,0,2Pi}]/(4Pi)];
Do[check["signed exterior on lips "<>ToString[signs,InputForm],
   Total[lip /@ signs],Total[signs]],
   {signs,{{1},{-1},{1,1},{1,-1},{-1,1},{-1,-1}}}];
check["upper exterior sphere tail times r^6",
  Limit[-r^6 cap[ArcTan[c/r^2]],r->Infinity,Assumptions->c>0],-4c^3/(3Pi)];
(* A regularized closed dimple has opposite changes of inclination on its
   inner and outer portions. A finite test angle keeps both terms nonzero. *)
check["smooth cap and outer dimple integrals at 60 degrees",
  {-cap[Pi/3],cap[Pi/3]}, {-2/3+Sqrt[3]/(2Pi),2/3-Sqrt[3]/(2Pi)}];

(* A changing, nonradial graph with nonzero mixed derivatives. The four
   signed minors give the full spacetime current, including shape changes. *)
coordinates={time,x,y,z};
height=(1+time)x^2/2 + y^2 + z^2/2 + x y z + time y z;
slopes=D[height,#]& /@ {x,y,z};
jac=Table[D[slopes[[i]],coordinates[[j]]],{i,3},{j,4}];
minors=Table[(-1)^(j-1)Det[Transpose[Delete[Transpose[jac],j]]],{j,4}];
check["current on a changing nonspherical test graph at (1,1,2,3)",
  minors /. {time->1,x->1,y->2,z->3},{3,-5,1,2}];
check["spacetime divergence of unweighted minors",
  Expand[Sum[D[minors[[i]],coordinates[[i]]],{i,4}]],0];
check["current contracted with derivative of slope norm",
  Expand[minors.(D[slopes.slopes,#]& /@ coordinates)],0];

(* Mechanical control with bulk-fixed rear reference and the same unit
   tension, wave-speed factor, rear stiffness and 4 Pi kappa. Vary action
   at assigned radii. Solve for the unique positive reflecting length. *)
Do[
  {a,j}=prep;
  positive=Select[len /. Solve[len+(1+1/a)(1-j/len^2)==3,len,Reals],#>0&];
  If[Length[positive]!=1,Exit[1]];
  ell=First[positive];
  Print["radius, action, wave energy, scalar flux = ",
    ToString[{a,j,j/ell,1-j/ell^2},InputForm]];
  If[!TrueQ[FullSimplify[ell+(1+1/a)(1-j/ell^2)==3]],Exit[1]],
  {prep,{{1,2},{2,25/6}}}];
