(* Construct the net-force Green field from equilibrium and surface traction.
   Compute bend work from virtual changes of material lengths and supports. *)
ClearAll["Global`*"];
say[label_, value_] := Print[label, " = ", ToString[value, InputForm]];
xyz = {x, y, z};
rad = Sqrt[xyz.xyz];
ax = {0, 0, 1};
u = aa ax/rad + dd z xyz/rad^3;
stress = Table[
   mu (D[u[[i]], xyz[[j]]] - D[u[[j]], xyz[[i]]]) +
     KroneckerDelta[i,j] bulk Div[u, xyz], {i,3}, {j,3}];
eq = Simplify[Div[stress, xyz], rad > 0];
(* Use u at a unit sphere, boundary stress times its outward normal.
   Axial symmetry reduces the independent angular variable to cos(theta). *)
tractionZ = Simplify[(stress.xyz)[[3]] /. {x->Sqrt[1-t^2], y->0, z->t},
                     -1<t<1];
net = -2 Pi Integrate[tractionZ, {t,-1,1}];
coeff = First[Solve[{(eq[[1]] /. {x->1,y->1,z->2})==0, net==1}, {aa,dd}]];
bulkResidual = Simplify[eq /. coeff, rad>0 && bulk>0 && mu>0];
say["coefficients of isotropic and radial tensors", {aa,dd}/.coeff];
say["bulk equation and net force", {bulkResidual, Simplify[net/.coeff]}];
(* Signed sources at d +/- b/2 and +/- b/2, all unit magnitudes. *)
locations = {-b/2,b/2};
weights = {-1,1};
g[r_] = Simplify[(aa+dd)/r /. coeff];
pair = Factor[Sum[weights[[i]] weights[[j]] g[d+locations[[j]]-locations[[i]]],
                   {i,2},{j,2}]];
say["work for two finite force pairs", pair];
power[expr_] := Limit[d D[expr,d]/expr, d->Infinity];
powers = power /@ {g[d], D[g[d],d], pair, D[pair,d]};
say["four far field powers", powers];

(* Virtual work at fixed physical imposed offsets. The near trace zN relaxes.
   p is the force entering the node and zN the near displacement. *)
sol = First[Solve[{p cap[h]+zN==drive, p==(ext[h]+contact[h])zN}, {p,zN}]];
store = p^2 cap[h]/2 + (ext[h]+contact[h]) zN^2/2;
onShell = Together[store/.sol];
qBend = (p^2 cap'[h] - zN^2 (ext'[h]+contact'[h]))/2;
workResiduals = Together /@ {
   D[onShell,h]+(qBend/.sol), D[onShell,drive]-(p/.sol),
   (p-contact[h] zN-ext[h] zN)/.sol};
ambiguity = Simplify[{-D[-sg charge h,drive], -D[-sg charge h,h]}];
say["bend work, drive work and force residuals", workResiduals];
say["effect of independent signed baseline on lateral and bend readings", ambiguity];

(* Compute the displacement along the distributed layer from its strain.
   Integrate its local material virtual work at fixed force for this variation. *)
len = len0-alpha sg ell h;
tension = t0 (1+v)^2;
slope = p/tension;
drop = Integrate[len slope, {v,0,1}];
qLocal = Integrate[tension slope^2 D[len,h]/2,{v,0,1}];
trialSol = First[Solve[{drop+zN==drive, p==4 zN}, {p,zN}]];
values = {p,2 zN,-2 zN,qLocal} /. trialSol /. {h->0,t0->2,len0->3,ell->1};
cases = {{1,1,1/10},{1,1,-1/10},{-1,1,1/10},{1,1,1/20},
         {1,1,-1/20},{1,1,0},{1,0,1/10}};
rows = (values /. Thread[{sg,alpha,drive}->#])& /@ cases;
say["control rows P S H Q", rows];
Print["audit ", StringRiffle[ToString[#,InputForm]& /@ Flatten[{powers,rows}], " "]];
If[bulkResidual =!= {0,0,0} || Simplify[(net/.coeff)-1] =!= 0 ||
    powers =!= {-1,-2,-3,-4} || workResiduals =!= {0,0,0} ||
    ambiguity =!= {0,sg charge} || rows[[1,4]] == 0 ||
    rows[[1,4]] =!= -rows[[3,4]] || rows[[4,4]] =!= rows[[1,4]]/4,
   Exit[1]];
