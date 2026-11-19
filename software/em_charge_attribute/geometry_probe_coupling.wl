(* Independent weak-probe action on the material coordinates of the collar. *)
ClearAll["Global`*"];
check[label_, expr_, expected_] := Module[{v=FullSimplify[expr]},
  Print[label, " = ", ToString[v, InputForm]];
  If[!TrueQ[FullSimplify[v == expected]], Exit[1]]
];
nav = {Sin[theta] Cos[az],Sin[theta] Sin[az],Cos[theta]};
degree[vec_] := Integrate[vec.Cross[D[vec,theta],D[vec,az]]/(4 Pi),
                          {theta,0,Pi},{az,0,2 Pi}];
check["degree of oriented sections", {degree[nav],degree[-nav]}, {1,-1}];
(* Shape and width are time dependent. A nonpolynomial order profile is
   used independently of the polynomial profile in the other script. *)
chi[u_] := (1+Cos[Pi u])/2;
path = 2 radius[t]+width[t] u;
shellCharge = s Integrate[-chi'[u],{u,0,1}];
check["canonical charge from trial potential energy",
      D[phi g0 shellCharge/Sqrt[zint],phi], s g0/Sqrt[zint]];
(* Gauge transformation of -g(phi-v A_r) is g times the total time
   derivative of the gauge function along each chart label. *)
localShift = g[radius[t],act[t]] (Derivative[0,1][lam][path,t]+
             D[path,t] Derivative[1,0][lam][path,t]);
endpoint = g[radius[t],act[t]] lam[path,t];
check["gauge variation after endpoint subtraction, per gauge function",
      Cancel[(localShift-D[endpoint,t])/lam[path,t]],
      -D[g[radius[t],act[t]],t]];

(* On a mirror pair, every neutral even term has the same value. Allow
   a fully shape-dependent carrier inventory and asymmetric well depths.
   Compare arbitrary reflected profiles, each partitioned into two samples. *)
cell[height_,ph_] := brane[height]+act0 freq[height]+
                      q ph (pop[height]+bath[height]);
upper = cell[h1,phi]+cell[h2,phi];
mirror = cell[h2,phi]+cell[h1,phi];
check["common static mode source per mouth", D[upper,phi],
      q (pop[h1]+pop[h2]+bath[h1]+bath[h2])];
check["native orientation-odd response at fixed carrier coupling",
      D[(upper-mirror)/2,phi],0];
check["orientation-odd response if a signed potential coupling is inserted",
      D[((upper+g0 phi)-(mirror-g0 phi))/2,phi],g0];
