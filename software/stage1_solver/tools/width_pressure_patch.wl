(* Solve the radial traction experiment by matching the interior and exterior
   displacements and gradient tractions. Compute the zero-stiffness case
   separately, retaining the boundary condition at infinity. *)
ClearAll["Global`*"];
$HistoryLength = 0;
assume = rad > 0 && length > 0 && grad > 0 && mass > 0;
lap[f_] := D[rad^2 D[f,rad],rad]/rad^2;
match[inner_,outer_] := First[Solve[
    ({inner == outer, D[inner,rad] == D[outer,rad]} /. rad -> length),{a,b}]];
inner = press/(grad mass^2) + a Sinh[mass rad]/rad;
outer = b Exp[-mass rad]/rad;
eq = FullSimplify[{grad (mass^2 inner-lap[inner])-press, grad (mass^2 outer-lap[outer])},assume];
coeff = match[inner,outer];
center = FullSimplify[Limit[inner /. coeff,rad->0]/press,assume && press!=0];
Print["radial residuals = ",InputForm[eq]," ; center per unit pressure = ",InputForm[center]];
If[eq =!= {0,0} || !TrueQ[FullSimplify[
   center == (1-(1+mass length) Exp[-mass length])/(grad mass^2),assume]],Exit[1]];
inSoft = a-press rad^2/(6 grad);
outSoft = b/rad;
soft = FullSimplify[Limit[inSoft /. match[inSoft,outSoft],rad->0],assume];
Print["zero-stiffness center width = ",InputForm[soft]];
Print["zero-stiffness exterior width = ",InputForm[FullSimplify[outSoft /. match[inSoft,outSoft],assume]]];
If[!TrueQ[FullSimplify[soft == press length^2/(2 grad),assume]],Exit[1]];
Exit[0];
