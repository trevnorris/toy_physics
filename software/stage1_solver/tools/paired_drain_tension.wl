(* Derive the stress in Cartesian coordinates, then obtain hemisphere
   tractions by a closing disk through the target centre. Check each
   hemisphere directly on the curved surface as well. *)
ClearAll["Global`*"];
xyz = {x,y,z};
ex = {1,0,0};
potential[center_,rate_] := rate/(4 Pi Sqrt[(xyz-center).(xyz-center)]);
velocity[center_,rate_] := Grad[potential[center,rate],xyz];
strain[v_] := Module[{j = D[v,{xyz}]},(j + Transpose[j])/2];

atMiddle = FullSimplify[strain[velocity[{0,0,0},qs]+velocity[d ex,qs]] /.
                         {x->d/2,y->0,z->0},d>0];
Print["midpoint strain eigenvalues in the axial and transverse directions = ",
      ToString[Diagonal[atMiddle],InputForm]];
If[!TrueQ[FullSimplify[atMiddle == qs DiagonalMatrix[{8,-4,-4}]/(Pi d^3),d>0]],Exit[1]];
incident = 2 eta strain[velocity[d ex,qs]];
divergence = FullSimplify[Div[incident,xyz],(x-d)^2+y^2+z^2>0];
diskIntegrand = FullSimplify[incident[[1,1]] /. {x->0,y->s,z->0},d>0 && s>=0];
halfFlux = FullSimplify[2 Pi Integrate[diskIntegrand*s,{s,0,b},
                        Assumptions->d>b>0],d>b>0];
Print["source stress divergence = ",ToString[divergence,InputForm],
      " ; source-facing closed-by-disk traction = ",ToString[halfFlux,InputForm]];
If[divergence =!= {0,0,0} || !TrueQ[FullSimplify[
       halfFlux == eta qs b^2/(d^2+b^2)^(3/2),d>b>0]],Exit[1]];

n = {m,Sqrt[1-m^2],0};
curved = FullSimplify[(ex.incident.n) /. Thread[xyz -> b n],-1<m<1 && d>b>0];
Do[
  rules = {d->sep,b->radius,eta->1,qs->1};
  front = NIntegrate[Evaluate[2 Pi b^2 curved /. rules],{m,0,1},WorkingPrecision->40,AccuracyGoal->30];
  back = NIntegrate[Evaluate[2 Pi b^2 curved /. rules],{m,-1,0},WorkingPrecision->40,AccuracyGoal->30];
  wanted = halfFlux /. rules;
  Print["d,b = ",{sep,radius}," front/back = ",ToString[N[{front,back},12],InputForm]];
  If[!TrueQ[Max[Abs[{front-wanted,back+wanted}]] < 10^-26],Exit[1]],
  {sep,{8,16}},{radius,{1,3}}];

targetStrain = strain[velocity[{0,0,0},qt]];
dipole = FullSimplify[(ex.(2 beta m targetStrain).n) /. Thread[xyz -> b n],b>0 && -1<m<1];
viscSlope = FullSimplify[2 Pi b^2 Integrate[dipole,{m,-1,1}],b>0];
Print["target radial strain times viscosity dipole, axial traction per U = ",ToString[viscSlope,InputForm]];
If[!TrueQ[FullSimplify[viscSlope == 4 qt beta/(3 b),b>0]],Exit[1]];
Exit[0];
