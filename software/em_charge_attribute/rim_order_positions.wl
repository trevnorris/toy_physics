(* Independent variation of the order functional at a regular initial level,
   using the first moment of a ball and its perturbed volume. *)
ClearAll["Global`*"];
emit[key_, v_] := Print[key, " = ", ToString[v,InputForm]];
ass = r>0 && k>0 && p>0 && Element[{amp,phase,mob,mu},Reals];
strainPotential[phi_] := mu amp^2 k^2 Sin[k r n+phi]^2/4;
(* r(n)=R+h*dR(n), and d chi / dh = -mob*dF/dchi. *)
dr[phi_] := mob strainPotential[phi]/p;
(* First moment variation divided by base volume. First-moment numerator
   vanishes at rest, so the differentiated volume contributes zero. *)
cmRate=FullSimplify[2 Pi r^3 Integrate[n dr[phase],{n,-1,1},
   Assumptions->ass]/(4 Pi r^3/3),Assumptions->ass];
emit["front centroid rate for general phase",cmRate];
reflected=FullSimplify[((cmRate/.phase->Pi/4)-(cmRate/.phase->-Pi/4))/2,
   Assumptions->ass];
emit["centered crest, reflection residual and half-amplitude factor",
   FullSimplify[{cmRate/.phase->0,cmRate+(cmRate/.phase->-phase),
   (reflected/.amp->amp/2)/reflected},Assumptions->ass && amp>0 && mu>0 && mob>0]];
speed=reflected/.{mob->1,mu->1,amp->1/10,k->1/2,p->1};
shape[s_]:=1+(3s^2-2s^3)/2;
emit["midplane and end front speeds and ratio",
   N[{speed/.r->1,speed/.r->3/2,(speed/.r->3/2)/(speed/.r->1)},16]];
bow=2 NIntegrate[(speed/.r->shape[s]) Cos[Pi s],{s,0,1},WorkingPrecision->32];
emit["projected front bow velocity",bow];

(* Unit step in each normalized mode. Integrate the first-order system with
   a matrix exponential, retaining the drive as a constant fourth state. *)
matrix[c_] := {{0,1,0,0},{-1-c,0,c,0},{c,0,-1-c,1},{0,0,0,0}};
initial={0,0,0,1};
emit["onset wall velocity, acceleration, jerk and front velocity",
   { (matrix[c].initial)[[1]], (MatrixPower[matrix[c],2].initial)[[1]],
     (MatrixPower[matrix[c],3].initial)[[1]], (matrix[c].initial)[[3]] }];
Do[emit["contact "<>ToString[c]<>" normalized wall/order at t=1",
   (MatrixExp[N[matrix[c],28]].initial)[[{1,3}]]],{c,{0,1}}];
(* Passive exchange: derivative of the stored energy plus drag loss
   minus input power, evaluated from the equations. *)
energy=(v^2+x^2+y^2+c (x-y)^2)/2;
yd=f-y-c (y-x);
balance=Expand[D[energy,x]v+D[energy,v](-x+c (y-x))+D[energy,y]yd+yd^2-f yd];
emit["power balance residual",balance];
If[!TrueQ[balance==0] || !NumberQ[bow],Exit[1]];
