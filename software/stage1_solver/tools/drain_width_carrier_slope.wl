(* Use y=n^2. The fixed-current, local Bernoulli row is the cubic
   y^3-y+2 I^2/W^2=0. Select its subsonic root tending to one at I=0. *)
ClearAll["Global`*"];
row = y^3-y+2 flow^2/width^2;
dpdy = D[(y^(5/2)-1)/5,y];
slope = Simplify[-dpdy D[row,width]/D[row,y]];
Print["pressure slope in y=n^2 = ",InputForm[slope]];
weak = Limit[(slope /. y->1)/flow^2,flow->0];
Print["leading coefficient divided by squared areal current = ",InputForm[weak]];
If[Simplify[weak-1/width^3] =!= 0,Exit[1]];
Do[
  local = row /. {flow->1/(16 Pi radius^2),width->1};
  roots = y /. NSolve[local==0,y,Reals,WorkingPrecision->45];
  base = First[Select[roots,#>1/Sqrt[3]&]];
  value = N[slope /. {y->base,flow->1/(16 Pi radius^2),width->1},16];
  Print["radius = ",radius," slope = ",InputForm[value]];
  If[!TrueQ[value>0 && Abs[local /. y->base]<10^-35],Exit[1]],
  {radius,{8,16}}];
Exit[0];
