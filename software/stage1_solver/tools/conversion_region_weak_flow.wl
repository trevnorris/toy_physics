(* Direct complex carrier, and the original fraction balance. Use the
   signed axial coordinate w (mouth zero, downstream positive). *)
ClearAll["Global`*"];
hb = 1/Sqrt[6];
h[n_] := 35 n^4/96;
mu = h[1] + 1/2;
psi0 = Exp[I w/hb];
rho1 = amp Sin[x] Exp[-2 w];
pot1 = u x + p Sin[x] Exp[-2 w];
psi1 = psi0 (rho1/2 + I pot1/hb);
lap[f_] := D[f, {w, 2}] + D[f, {x, 2}];
residual = FullSimplify[(-hb^2 lap[psi1]/2 + (h[1]-mu) psi1 + h'[1] rho1 psi0)/psi0];
pRule = First[Solve[(residual /. {x -> Pi/2, w -> 0}) == 0, p]];
Print["carrier residual before mode matching = ", InputForm[residual]];
Print["mode and carrier residual = ", InputForm[{p, residual} /. pRule // Simplify]];
pot = pot1 /. pRule;
assume = Element[{x,w,u,amp,b}, Reals];
normalFlux1 = FullSimplify[hb Im[Conjugate[psi1] D[psi0,w] +
                        Conjugate[psi0] D[psi1,w]] /. pRule, assume];
Print["current derivative at mouth = ", InputForm[normalFlux1 /. w -> 0]];
Print["paired opposite-point current sum = ",
      InputForm[Simplify[normalFlux1 + (normalFlux1 /. x -> -x)]]];

(* In the stationary response of this dynamical local test, linearize
   div(n chi u) + n chi (1-chi)=0 in conservative form. Try a translated
   order profile with the same decaying exponential as the imposed mode. *)
chi0 = 1/(1 + Exp[w-b]);
chi1 = s Sin[x] Exp[-2 w] D[chi0,w];
dens = 1 + e rho1;
frac = chi0 + e chi1;
ux = e D[pot,x];
uw = 1 + e D[pot,w];
ord = D[dens frac ux,x] + D[dens frac uw,w] + dens frac (1-frac);
r0 = Simplify[ord /. e -> 0];
r1 = Simplify[D[ord,e] /. e -> 0];
sRule = First[Solve[(r1 /. {x -> Pi/2, w -> b}) == 0, s]];
Print["order base, linear coefficient and fitted mode = ", InputForm[{r0,r1,s /. sRule}]];
Print["order residual after matching = ", InputForm[Simplify[r1 /. sRule]]];
loss = dens frac (1-frac);
(* Read level and peak by evaluating both conditions on w=b+e d. *)
positionShift[condition_] := Module[{condition1},
  condition1 = Simplify[(D[condition /. w -> b+e d,e] /. e -> 0) /. sRule];
  d /. First[Solve[condition1 == 0,d]]
];
Print["[isolated level, loss slope, loss curvature] at b = ",
      InputForm[Simplify[{frac,D[loss,w],D[loss,{w,2}]} /. {e->0,w->b}]]];
Print["[half-order shift, net-loss peak shift] = ",
      InputForm[Simplify[{positionShift[frac-1/2],positionShift[D[loss,w]]}]]];
