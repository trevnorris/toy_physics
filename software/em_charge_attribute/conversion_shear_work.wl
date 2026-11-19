(* Work in the translated shear profile's coordinates. u=tanh[(q-X)/(2 delta)].
   Translating the front relative to that profile changes the order coefficient. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[name_, value_] := Print[name, " = ", ToString[value, InputForm]];
check[name_, left_, right_] := Module[{error = Factor[left - right]},
  emit[name, {left, right, error}];
  If[error =!= 0, Exit[1]]];
order = (1 - u)/((1 - u) + (1 + u) Exp[(x - y)/width]);
jac = 2 width/(1 - u^2);
origin = {x -> 0, y -> 0};
work[integrand_] := Factor[Integrate[Cancel[integrand /. origin], {u, -1, 1}]];
check["front weight per unit stored energy",
  Factor[(D[order, y]/order) /. origin], (1 - (order /. origin))/width];

(* Independent check from a displacement variation of u_t(q) at fixed domain:
   the stress conjugate to grad u_t has only the tangential component on this cut.
   A normal displacement u_n(q) adds no curl in this planar geometry. *)
grad = Array[Unique["grad"] &, {4, 4}];
curlSquare = Sum[(grad[[i,j]] - grad[[j,i]])^2, {i,1,3}, {j,i+1,4}];
energyDensity = chi mu curlSquare/2;
pureShear = Map[# -> If[# === grad[[4,1]], gamma, 0] &, Flatten[grad]];
cutLoad = Table[D[energyDensity, grad[[4,j]]] /. pureShear, {j,1,4}];
emit["displacement traction on normal cut for tangent 1 shear", cutLoad];
check["direct normal displacement traction", Last[cutLoad], 0];

Do[
  integrand = jac order peak (1 - u^2)^p;
  stored = work[integrand];
  front = work[-D[integrand, y]];
  profile = work[-D[integrand, x]];
  ratio[p] = Factor[-front/stored];
  emit["{p, E/A, F_front/A, F_profile/A, -F_front/E}",
    {p, stored, front, profile, ratio[p]}];
  check["common translation", front + profile, 0];
  check["zero excitation", front /. peak -> 0, 0],
  {p, {1, 2}}];
emit["equal energy ratio", Factor[ratio[2]/ratio[1]]];
Exit[0];
