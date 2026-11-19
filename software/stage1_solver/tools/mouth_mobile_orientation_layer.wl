(* A retained surface of orientable sites on one translating spherical mouth.
   Coordinates X, u and s are the marked center, tangential dipole shift along
   the projected translation axis, and excess geometric-frame angular link. *)
ClearAll["Global`*"];
emit[label_, v_] := WriteString["stdout", label, " = ", ToString[v, InputForm], "\n"];
thetaNew = theta - u Sin[theta]/r;
jac = Sin[thetaNew] D[thetaNew,theta]/Sin[theta];
densitySlope = Simplify[D[1/jac,u] /. u->0];
emit["fractional density change per tangential displacement",densitySlope];
rate = rateDipole Cos[theta]; (* net rate of site addition per physical area *)
areaMean[f_] := Integrate[f Sin[theta],{theta,0,Pi}]/2;
source = rate /. rateDipole-> -2 sigma ud/r;
emit["source holding density fixed and its surface average",{source,areaMean[source]}];
(* Normalize the energy by two thirds of the initial area. *)
surfaceAverage[f_] := 3 areaMean[f]/2;
tag = {r Sin[thetaNew], x+r Cos[thetaNew]};
vel = (xd D[tag,x]+ud D[tag,u]) /. u->0;
kin = surfaceAverage[sigma vel.vel/2]+ m (sd-b ud)^2/2;
pot = surfaceAverage[ka (u densitySlope)^2/2]+h s^2/2+bond(y-s)^2/2;
mass = Table[D[kin,i,j],{i,{xd,ud,sd}},{j,{xd,ud,sd}}];
spring = Table[D[pot,i,j],{i,{x,u,s,y}},{j,{x,u,s,y}}];
emit["translation-tangential-angular mass matrix",mass];
emit["tangential restoring coefficient",spring[[2,2]]];
dyn = spring - w^2 PadRight[mass,{4,4}];
ports={1,4}; internal={2,3};
(* Keep both reciprocal reactions; virtual X holds the exterior slip y fixed. *)
portMatrix = dyn[[ports,ports]]-dyn[[ports,internal]] .
    Inverse[dyn[[internal,internal]]] . dyn[[internal,ports]];
cross = Factor[portMatrix[[2,1]]];
emit["cross port reciprocity residual",Factor[cross-portMatrix[[1,2]]]];
slow = Assuming[ka>0 && bond>0 && h>0, FullSimplify[Limit[cross/w^4,w->0]]];
emit["slow twist from translation per omega^4",slow];
massSlow = Factor[-Limit[(portMatrix[[1,1]]-Limit[portMatrix[[1,1]],w->0])/w^2,w->0]];
emit["slow translating layer mass per two-thirds area",massSlow];
checks = {densitySlope===2 Cos[theta]/r, areaMean[source]===0,
    Simplify[spring[[2,2]]-2 ka/r^2]===0,
    Factor[cross-portMatrix[[1,2]]]===0,
    Factor[slow-b bond m sigma r^2/(2 ka (bond+h))]===0,
    massSlow===3 sigma/2};
emit["identity checks",checks];
If[!And@@checks,Exit[1]];
