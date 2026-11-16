(* Independent lab-frame calculation. Integrate momentum at an exterior
   sphere and subtract momentum stored between that sphere and the intake. *)
ClearAll["Global`*"];
$HistoryLength = 0;
xyz = {x, y, z};
positive = a > 0 && r > a && rho > 0 && mass > 0 && intake > 0;

(* Release from rest with constant trial acceleration. Slip and prescribed
   uniform suction give the relative normal fluid speed at the moving mouth. *)
dx = x - acc*t^2/2;
radius = Sqrt[dx^2 + y^2 + z^2];
potential = flow*x + c/radius + (d0 + d1*t)*dx/radius^3;
velocity = D[potential, #] & /@ xyz;
onSphere[expr_, rr_] := FullSimplify[
  expr /. {x -> rr*u, y -> rr*Sqrt[1-u^2], z -> 0},
  Assumptions -> rr > 0 && -1 < u < 1];
normal = {u, Sqrt[1-u^2], 0};
(* At rest the first time derivative of the surface position and normal
   vanishes. The velocity boundary condition itself changes with acc. *)
bcValue = onSphere[(velocity /. t -> 0).normal, a] + intake/(4*Pi*a^2);
bcRate = onSphere[(D[velocity, t] /. t -> 0).normal, a] - acc*u;
conditions = Join[CoefficientList[bcValue, u], CoefficientList[bcRate, u]];
coeffs = First[Solve[Thread[conditions == 0], {c, d0, d1}]];
flowField = velocity /. coeffs;
pressure = -rho*(D[potential, t] + (velocity.velocity - flow^2)/2) /. coeffs;
radialCurrent = ((rho*First[flowField]*flowField + {pressure, 0, 0}).normal) /. t -> 0;
inward = Factor[-2*Pi*r^2*Integrate[onSphere[radialCurrent, r], {u, -1, 1}]];
shellDensityRate = onSphere[rho*(D[First[flowField], t] /. t -> 0), r];
storedRate = Integrate[2*Pi*r^2*Integrate[shellDensityRate, {u, -1, 1}],
                       {r, a, outer}, Assumptions -> outer > a > 0];

(* Balance the whole instrument. The outlet passes mass at the input rate
   with mean lateral speed fraction*flow at the observation, without an
   additional lateral traction. No accumulation other than mass*acc. *)
residual = inward - storedRate - rho*intake*fraction*flow + support + probe - mass*acc;
holding = Factor[support /. First[Solve[(residual /. {probe -> 0, acc -> 0}) == 0, support]]];
response = Factor[acc /. First[Solve[(residual /. support -> 0) == 0, acc]]];
Print["normal_bc_residuals = ", InputForm[Simplify[conditions /. coeffs]]];
Print["inward_exterior = ", InputForm[inward]];
Print["shell_storage_rate = ", InputForm[storedRate]];
Print["external_hold = ", InputForm[holding]];
Print["released_acceleration = ", InputForm[response]];

units = {a -> 1, rho -> 3/(4*Pi), intake -> 4*Pi, mass -> 1};
values = Association[];
Do[
  state = Join[units, {fraction -> kValue, flow -> uValue}];
  name = StringJoin["k=", ToString[kValue, InputForm], ",U=", ToString[uValue, InputForm]];
  row = Prepend[Table[response /. state /. probe -> applied, {applied, {-3/200, 0, 3/200}}],
                holding /. state];
  AssociateTo[values, name -> (ToString[#, InputForm] & /@ row)],
  {kValue, {0, 1}}, {uValue, {1/100, 0, -1/100}}];
Print["measurements = ", ExportString[values, "RawJSON", "Compact" -> True]];
Exit[0];
