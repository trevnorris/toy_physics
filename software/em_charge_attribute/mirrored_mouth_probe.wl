(* Independent weak-probe calculation from real density and phase. *)
ClearAll["Global`*"];
emit[label_, value_, expected_] := Module[{v = FullSimplify[value]},
  Print[label, " = ", ToString[v, InputForm]];
  If[!TrueQ[FullSimplify[v == expected]], Exit[1]]];

psi = Sqrt[density] Exp[I th];
dt = I (thdot + charge eps probe/hb) psi;
dtbar = -I (thdot + charge eps probe/hb) Sqrt[density] Exp[-I th];
timeTerm = I hb (Sqrt[density] Exp[-I th] dt - psi dtbar)/2;
emit["linear action reading", -D[timeTerm, eps], charge density probe];

(* Transport on two paths of unequal residence at one common inward rate. *)
flow = {frac rate, (1-frac) rate};
occupancy = flow {tau1, tau2};
depth = {w1, w2};
mode[k_, x_] := HermiteH[k, x/lam];
reading[k_, sign_] := charge occupancy . (mode[k, sign #]& /@ depth);
Do[emit["half mirror difference mode " <> ToString[k],
   (reading[k, 1] - reading[k, -1])/2,
   If[k == 1, 2 charge (occupancy.depth)/lam, 0]], {k, 0, 2}];
op[f_] := -Exp[(w/lam)^2] D[Exp[-(w/lam)^2] D[f,w],w];
emit["flat transverse operator", op[mode[0,w]], 0];
emit["odd mode mass squared", op[mode[1,w]]/mode[1,w], 2/lam^2];
emit["centered receiver odd overlap", mode[1,0], 0];

(* Round layer integrated on a physical sphere. Matter dwell times may
   differ at the two radii. A separate bound source can be recovered. *)
order = (1+Cos[Pi u])/2;
Do[
 area = rr^2 Sin[theta];
 along = -D[order,u]/(4 Pi rr^2);
 count = Integrate[area along,{theta,0,Pi},{phi,0,2 Pi},{u,0,1}];
 Do[
   total = reading[0,s] + bound s count;
   reconstructed = charge Total[flow {tau1,tau2}];
   emit["residual {R, orientation} " <> ToString[{rr,s}], total-reconstructed,s bound],
   {s,{-1,1}}],
 {rr,{1,2}}];
