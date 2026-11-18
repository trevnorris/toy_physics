(* Use the Stokes stream function for the two velocity components. Derive
   the change in stress between spheres from the convective volume force. *)
ClearAll["Global`*"];
surface[expr_] := Factor[2 Pi r^2 Integrate[Expand[expr Sin[th]], {th, 0, Pi}]];
check[label_, lhs_, rhs_] := Module[{diff = Factor[lhs-rhs]},
  Print[label, " | ", ToString[lhs, InputForm], " | ", ToString[rhs, InputForm],
        " | residual ", ToString[diff, InputForm]];
  If[diff =!= 0, Exit[1]]];

stream = (w r^2/2+c r+h/r) Sin[th]^2;
radial = Simplify[D[stream,th]/(r^2 Sin[th])];
polar = Simplify[-D[stream,r]/(r Sin[th])];
traction = eta (r D[polar/r,r]+D[radial,th]/r);
bc = First@Solve[{(radial/.{r->a,th->0})==b,
                  (traction/.{r->a,th->Pi/2})==0},{c,h}];
radial = Factor[radial/.bc];
polar = Factor[polar/.bc];
base = -q/r^2;
rate = 4 Pi dens q;
axial = radial Cos[th]-polar Sin[th];
carried = surface[-dens (base (v+axial)+radial base Cos[th])];
check["intake perturbation",surface[-dens radial],0];
check["held rim carriage per mass",Factor[(carried/.{r->a,v->0})/rate],w/3+b];
check["moving rim carriage per mass",Factor[(carried/.r->a)/rate],v+w/3+b];
check["held normal locking limit",(carried/.{v->0,b->0,r->a}),rate w/3];

(* In the mouth frame, div(sigma) = dens times the convective acceleration.
   Its volume integral gives the signed change in stress between surfaces.
   This uses velocity alone, independently of any supplied pressure profile. *)
convRadial = base D[radial,r]+radial D[base,r];
convPolar = base (D[polar,r]+polar/r);
shellSource = surface[dens (convRadial Cos[th]-convPolar Sin[th])];
stressChange = Factor[Integrate[shellSource,{r,a,radius},Assumptions->radius>a>0]];
carryChange = Factor[(carried/.r->radius)-(carried/.r->a)];
check["stress change vs inward carriage lost",stressChange,-carryChange];

(* Obtain the rim stress from the zeroth-order stress and the reciprocal
   volume integral using a unit radial dipole in a resting far medium. *)
normalStress = 2 eta D[radial,r]-(2 eta c Cos[th]/r^2/.bc);
leading = surface[normalStress Cos[th]-(traction/.bc) Sin[th]]/.r->a;
testFlow = {radial,polar}/.{w->0,b->1};
workDensity = surface[dens testFlow . {convRadial,convPolar}];
correction = Factor[-Integrate[workDensity,{r,a,Infinity},Assumptions->a>0]];
relativeTotal = Factor[leading+correction+(carried/.{r->a,v->0})];
check["relative total at rim",relativeTotal,4 Pi eta a (w-b)+rate (w+b)/2];
