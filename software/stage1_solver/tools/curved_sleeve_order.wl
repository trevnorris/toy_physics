(* Shape variation of an S^2 sleeve in four spatial dimensions. Distinguish
   the thin-wall load from a finite-thickness order change sampled at a point. *)
ClearAll["Global`*"];
well[q_] := a q^2 (1-q)^2;
muN = D[K0 n^5/4 + well[c] + k grad^2/2, n];
Print["constituent chemical potential, direct order derivative = ", ToString[{muN,D[muN,c]},InputForm]];

(* Use the divergence of the compressive order stress at z=0. For a profile
   of signed distance, div Pi along N = PiNN' + H (PiNN-PiTT). *)
profile = (1 - Tanh[z/(2 ell)])/2 - eps Exp[-z^2/b^2];
normalStress = -well[profile] + k D[profile,z]^2/2;
tangentStress = -well[profile] - k D[profile,z]^2/2;
normalForce = (-D[normalStress,z] - hcurv (normalStress-tangentStress)) /. z -> 0;
initialSlope = D[profile,z] /. {eps->0,z->0};
orderAtSurface = Simplify[(normalForce /. eps->0)/initialSlope];
response = Simplify[(D[normalForce,eps] /. eps->0)/initialSlope];
controls = Simplify[Table[response/a /. {b->size,k->2 a ell^2},{size,{ell,3 ell}}]];
Print["order derivative at surface, first order force change per profile slope = ",
      ToString[{orderAtSurface,response},InputForm]];
Print["normalized change at b=ell and 3ell = ", ToString[controls,InputForm]];

(* Radial EL derivative of the actual surface area. A radial variation has
   normal component delta R / sqrt(1+R'^2). *)
area = 4 Pi rad^2 Sqrt[1+p^2];
rules = {rad -> radius[w], p -> radius'[w]};
scale = Sqrt[1+radius'[w]^2];
el = (D[area,rad] /. rules) - D[D[area,p] /. rules,w];
signedCurvature = Simplify[el/(4 Pi radius[w]^2)];
axial = signedCurvature radius'[w]/scale;
Print["curvature sum and axial traction per tension = ",
      ToString[{signedCurvature,Simplify[axial]},InputForm]];
boundaryAxial = (area - p D[area,p]) /. rules;
residual = Simplify[axial (area /. rules) - D[boundaryAxial,w]];
Print["axial force minus rim derivative per tension = ",residual];
file = FileNameJoin[{"software","stage1_solver","frozen","m1c",
    "834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8","wp1_background_10x8.json"}];
cfg = Import[file,"RawJSON"]["config"]["branch"];
If[cfg["geometry_profile"] =!= "cubic_smoothstep", Exit[1]];
{lo,hi,r0,r1} = Rationalize[Lookup[cfg,{"w_min","w_max","r_mouth","r_exit"}],10^-12];
profileR[x_] := With[{s=(x-lo)/(hi-lo)}, r0+(r1-r0)(3 s^2-2 s^3)];
axialMid = Simplify[axial /. radius->profileR /. w->(lo+hi)/2];
rim = boundaryAxial /. radius -> profileR;
total = Simplify[(rim /. w->hi) - (rim /. w->lo)];
Print["archived [midpoint axial traction, integrated axial side force] per tension = ",
      ToString[{axialMid,total},InputForm]];
If[residual =!= 0 || !VectorQ[N[{axialMid,total}],NumberQ],Exit[1]];
