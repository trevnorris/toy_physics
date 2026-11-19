(* Obtain the wall energy for arbitrary monotone profiles from the square
   completion, then compare a relaxed column with an unmodified wall trial. *)
ClearAll["Global`*"];
$Assumptions = And @@ Thread[{a,kap,aa,nu0,q,len,c,r,tau,ell}>0];
v[u_] := u^2 (1-u)^2;
wall[coeff_] := Integrate[Sqrt[2 kap coeff] u (1-u),{u,0,1}];
(* The profile that saturates the bound has chi' = sqrt(2 coeff/kap) chi(1-chi). *)
neutral = wall[a];
width[coeff_] := Sqrt[kap/(2 coeff)];
coeff = a/(1+q/(aa nu0))^2;
relaxed = FullSimplify[aa (wall[coeff]-neutral)/neutral];
(* In the neutral minimizer the gradient and potential contribute equally. *)
frozen = FullSimplify[aa ((1+coeff/a)/2-1)];
Print["width ratio and relaxed excess / neutral area energy = ",
      InputForm[FullSimplify[{width[coeff]/width[a],relaxed}]]];
Print["area derivative, relaxed and frozen limits = ",
      InputForm[FullSimplify[{D[relaxed,aa],Limit[relaxed,aa->Infinity],
                             Limit[frozen,aa->Infinity]}]]];
Print["relaxed excess above its dilute limit = ",
      InputForm[FullSimplify[relaxed-Limit[relaxed,aa->Infinity]]]];
(* At a prescribed width, integrate the square left over the neutral bound. *)
ds = Sqrt[kap/(2 a)];
gap = Integrate[kap/2 (u (1-u) (1/ell-1/ds))^2 ell/(u (1-u)),{u,0,1}];
Print["fixed-width excess / [a (ell-delta)^2/(6 ell)] = ",
      InputForm[FullSimplify[gap/(a (ell-ds)^2/(6 ell)),ell!=ds && $Assumptions]]];
(* An areal excess is a separate constrained branch. Work in the area variable. *)
disk = c q^2 Sqrt[Pi/aa]+tau aa;
pressureBalance = FullSimplify[2 t (D[disk,aa] /. aa->t^(2/3)),t>0 && $Assumptions];
areaPower = t /. First[Solve[pressureBalance==0,t]];
areaStar = FullSimplify[areaPower^(2/3)];
rStar = FullSimplify[Sqrt[areaStar/Pi]];
curvature = FullSimplify[D[disk /. aa->Pi r^2,{r,2}] /. r->rStar];
Print["constrained radius and radial curvature = ",InputForm[{rStar,curvature}]];
If[!TrueQ[FullSimplify[D[relaxed,aa]<0 && curvature>0 &&
        Limit[relaxed-frozen,aa->Infinity]==0]],Exit[1]];
