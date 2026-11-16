(* Neutral GNLS, m = hbar = 1. Check the moving-aperture equation and
   lab-frame momentum in one brane direction. The other coordinates are
   spectators to the translation. Normal flux laws move with the opening. *)
residual[p_] := I D[p,t] + (D[p,{x,2}] + D[p,{z,2}])/2 - h p;
moving = Exp[-I mu t] f[x-v t,z];
held = Exp[I u x - I (mu+u^2/2) t] f[x,z];
movingEq = Simplify[(residual[moving]/Exp[-I mu t]) /. x -> y+v t];
heldEq = Simplify[(residual[held]/Exp[I u x-I (mu+u^2/2) t]) /. x -> y];
frameDiff = Simplify[movingEq - (heldEq /. u -> -v)];
Print["moving residual = ", InputForm[movingEq]];
Print["held residual = ", InputForm[heldEq]];
Print["difference at u=-v = ", InputForm[frameDiff]];

(* Use two real envelope components to avoid assumptions about conjugates. *)
momentum[p_,q_] := p D[q,x] - q D[p,x];
rot = {{Cos[u x],-Sin[u x]}, {Sin[u x],Cos[u x]}};
f0 = {a[x],b[x]};
heldMomentum = Simplify[Apply[momentum,rot.f0]];
boostCheck = Simplify[heldMomentum - Apply[momentum,f0] - u f0.f0];
Print["carried momentum from lateral phase factor = ",
      InputForm[Simplify[heldMomentum-Apply[momentum,f0]]]];
Print["boost identity residual = ", InputForm[boostCheck]];

(* First derivative of a periodic carrier's integrated momentum. *)
f1 = {c[x],d[x]};
derivative = D[Apply[momentum,f0+v f1],v] /. v -> 0;
boundaryTerm = a[x] d[x]-b[x] c[x];
partsCheck = Simplify[derivative - 2*(c[x] b'[x]-d[x] a'[x]) - D[boundaryTerm,x]];
Print["first momentum variation minus integrated form is d_x of ", InputForm[boundaryTerm]];
Print["local integration by parts residual = ", InputForm[partsCheck]];
If[{frameDiff,boostCheck,partsCheck} =!= {0,0,0}, Exit[1]];
