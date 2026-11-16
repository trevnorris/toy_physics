(* Build the mixed current directly from the complex carrier. The
   stress representative is hb^2 Re[(d_i psi)^* d_j psi] minus
   hb^2/4 d_i d_j |psi|^2, plus EOS pressure on the diagonal (m=1). *)
ClearAll["Global`*"];
hb = 1/Sqrt[6];
z2 = Sin[x] Exp[-2 w];
z3 = Sin[x] Exp[-3 w];
rho1 = b z2 + d (8 z3 - 3 z2);
pot1 = u x + 2 b z2/3 + d (3 z3 - 2 z2);
(* Work with the envelope after removing the background phase Exp[I w/hb]. *)
envelope = 1 + ep (rho1/2 + I pot1/hb);
dx = D[envelope,x];
dw = D[envelope,w] + I envelope/hb;
(* Expand in real coordinates before taking the real drive derivative.
   Differentiating Re or Conjugate as a complex function loses terms. *)
first[quantity_] := Simplify[D[ComplexExpand[quantity],ep] /. ep -> 0];
jx1 = first[hb Im[Conjugate[envelope] dx]];
jw1 = first[hb Im[Conjugate[envelope] dw]];
dens = ComplexExpand[Conjugate[envelope] envelope];
mixed = first[hb^2 Re[Conjugate[dx] dw] - hb^2 D[dens,x,w]/4];
diagonal = first[hb^2 Re[Conjugate[dx] dx] - hb^2 D[dens,{x,2}]/4 + 7 dens^5/24];
(* At the base state, j_w/n=1 and j_x=0. *)
carried = jx1;
capillary = Simplify[mixed-carried];
balance = Simplify[D[diagonal,x]+D[mixed,w]];
Print["bulk derivatives and their sum = ",InputForm[Simplify[
    {D[diagonal,x],D[mixed,w],balance}]]];
cut = Simplify[{jw1,carried,capillary,mixed} /. w -> 0];
Print["entrance [number, carried, density gradient, total] = ",InputForm[cut]];
Print["deep [number, total] = ",InputForm[Limit[{jw1,mixed},w->Infinity]]];
strip = Integrate[cut,{x,-a,a},Assumptions->a>0];
Print["full strip = ",InputForm[Simplify[strip /. a->Pi]]];
(* Integrate the ball on radial shells; x/r is the angular cosine. *)
ball = Integrate[2 Pi r^2 (cut /. x -> r mu), {r,0,a},{mu,-1,1},Assumptions->a>0];
Print["ball [number, carried, density gradient, total] = ",InputForm[Simplify[ball]]];
If[balance =!= 0 || !FreeQ[cut,_Complex], Exit[1]];
