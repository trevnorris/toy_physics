(* Compute the actual graph normal, then intersect its translated gradient
   sphere with a regular ray. Laboratory directions are identical in all trials. *)
ClearAll["Global`*"];
emit[label_, expr_, assumptions_: True] := Module[{v=FullSimplify[expr,assumptions]},
  Print[label," = ",ToString[v,InputForm]]; v];
verify[a_,b_] := If[!TrueQ[FullSimplify[a==b]],Exit[1]];
vars={x,y,z};
height=e ((x^2+y^2+z^2)/(2 rad)-b z);
p=Grad[height,vars];
emit["projected normal numerator",p];
point=First[Solve[Thread[p=={0,0,0}],vars]];
emit["critical point",vars/.point];
metric=IdentityMatrix[3]+Outer[Times,p,p];
verify[emit["metric determinant at critical point", Det[metric/.point]],1];
second=Table[D[height,vars[[i]],vars[[j]]],{i,3},{j,3}];
verify[emit["three principal curvature magnitudes for e,rad>0",
  Eigenvalues[(Inverse[metric].second/Sqrt[1+p.p])/.point],e>0 && rad>0],
  {e/rad,e/rad,e/rad}];

(* Positive x ray in gradient space, in units with e=1. Each positive
   intersection contributes sign(dot(outward normal, ray direction)). *)
degree[u_,beta_] := Total[Function[solution,Sign[(ray/.solution)/u]] /@
   Solve[ray^2+beta^2==u^2 && ray>0,ray,Reals]];
verify[emit["regular-ray counts at radii 2 and 3, tilt 5/2",
            {degree[2,5/2],degree[3,5/2]}],{0,1}];
order[v_] := 1-10(v-2)^3+15(v-2)^4-6(v-2)^5;
count[beta_] := Module[{cuts=Sort[DeleteDuplicates[Join[{2,3},Select[{beta},2<#<3&]]]]},
  Total[Function[pair,degree[Mean[pair],beta]*(order[First[pair]]-order[Last[pair]])] /@
    Partition[cuts,2,1]]];
betas={0,9/4,5/2,11/4,4};
verify[emit["upper layer counts",count /@ betas],{1,459/512,1/2,53/512,0}];
verify[emit["mirror counts",(-count[#]& /@ betas)],{-1,-459/512,-1/2,-53/512,0}];

(* Check continuity on the time-dependent induced metric in Cartesian
   graph coordinates. The laboratory w velocity contains the wall speed. *)
h=height /. b->bfun[t];
grad=Grad[h,vars];
jac=Sqrt[1+grad.grad];
vxyz=flux vars/(4 Pi n0 (vars.vars)^(3/2));
vw=D[h,t]+grad.vxyz;
dens=n0/jac;
verify[emit["moving-sheet continuity",
  D[jac dens,t]+Div[jac dens vxyz,vars],vars.vars>0],0];
verify[emit["flow across sphere",4 Pi (vars.vars) jac dens (vxyz.vars)/Sqrt[vars.vars],
  vars.vars>0],flux];
pole={x->0,y->0,z->rad bfun[t]};
emit["finite moving w velocity at the crossing",vw/.pole,rad>0 && bfun[t]>0];
(* Average the radial kinetic term over the sphere; evaluate the remaining
   elementary integrals by moments in units of radius. *)
ke = Integrate[(1+e^2(u-b c)^2)/(2 u^2),{c,-1,1},{u,2,3}];
verify[emit["dimensionless static kinetic energy",ke],1/6+e^2(1+b^2/18)];
