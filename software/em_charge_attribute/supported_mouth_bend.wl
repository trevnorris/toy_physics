(* Integrate the driver reaction using the transmitted tension as a parameter. *)
ClearAll["Global`*"];
check[name_, value_, expected_] := Module[{v = Factor[value]},
  Print[name, " = ", ToString[v, InputForm]];
  If[!TrueQ[FullSimplify[v == expected]], Exit[1]]
];

(* The harmonic midplane carries constant radial load in the shell. The
   regular blank has zero radial load inside the same driven sphere. *)
$Assumptions = outer > inner > 0 && mu > 0 && ell > 0 && J > 0 && nu > 0 &&
  gap >= 0 && k > 0 && kb > 0 && kc > 0;
shellDrop = Integrate[q/(4 Pi mu radius^2), {radius, inner, outer}];
shellLoad = q /. First[Solve[shellDrop == z-b, q]];
check["radial rim stiffness", shellLoad/(z-b), 4 Pi mu/(1/inner - 1/outer)];
check["regular blank enclosed load", Limit[4 Pi mu radius^2 finiteSlope, radius -> 0], 0];

(* At equilibrium z=b+p/k, Y=rest+eta b-p/kb and p=V'(ell). *)
rim = b + p/k;
reflector = rest + eta b - p/kb;
span = reflector - rim;
dpdb = d /. First[Solve[d == kc (D[span, b] + D[span, p] d), d]];
ks = 1/(1/k + 1/kb + 1/kc);
check["transmitted tension derivative", dpdb, (eta-1) ks];
along[expr_] := D[expr,b] + D[expr,p] dpdb;
check["rim and reflecting end response", {along[rim], along[reflector]},
  {1-(1-eta)ks/k, eta+(1-eta)ks/kb}];
(* On the stationary branch the two spring energies are p^2/(2 k_i).
   Work is integrated along this branch; the core gains p d(ell). *)
workSlope = along[p^2/(2 k) + p^2/(2 kb)] + p along[span];
check["total excess work slope", workSlope, (eta-1)p];
check["excess work curvature", along[Expand[workSlope] // Simplify], (1-eta)^2 ks];
check["comoving reference work", workSlope /. eta -> 1, 0];
freq = Sqrt[gap^2 + nu^2/ell^2];
wavePressure = -J D[freq,ell];
check["wave incremental stiffness", FullSimplify[-D[wavePressure,ell]],
  J nu^2 (3 gap^2 ell^2 + 2 nu^2)/(ell^3 (gap^2 ell^2+nu^2)^(3/2))];
check["wave storage slope for fixed rear reference",
  -wavePressure (along[span] /. eta->0), wavePressure ks/kc];

(* Independent finite work: change variables from b to the positive
   reflecting length and integrate the driver reaction, keeping the rear
   reference fixed. This avoids subtraction of separately stored energies. *)
load[u_] := 2 - 1/u^2;
drive[u_] := 3 - u - 2 load[u];
workPrimitive = Integrate[-load[u] D[drive[u],u],u];
Print["b       rim_shift      reflector_shift  wave_change    excess_work"];
Do[
  roots = u /. NSolve[Numerator[Together[drive[u]-driveValue]] == 0, u, Reals,
    WorkingPrecision -> 45];
  positive = Select[roots, # > 0 &];
  If[Length[positive] != 1, Exit[1]];
  len = First[positive];
  tension = load[len];
  zNow = driveValue + tension;
  endNow = 3 - tension;
  work = (workPrimitive /. u -> len) - (workPrimitive /. u -> 1);
  row = {driveValue, zNow-1, endNow-2, 1/len-1, work};
  Print[StringRiffle[(ToString[NumberForm[N[#,16],{14,10}, NumberPadding->{"","0"},
    ExponentFunction->(Null&)]]&) /@ row," "]];
  If[Abs[endNow-zNow-len] > 10^-35, Exit[1]],
  {driveValue,{1/50,-1/50,1/100,-1/100}}];
