(* Solve the radial field jump, then integrate the force on the second shell. *)
ClearAll["Global`*"];
check[name_, actual_, expected_] := Module[{value = FullSimplify[actual]},
  Print[name, " = ", ToString[value, InputForm]];
  If[!TrueQ[FullSimplify[TrigToExp[value-expected] == 0]], Exit[1]]
];
$Assumptions = a > 0 && mu > 0 && r > 0 && dist > 2 a;
inside = ii Sinh[mu r]/r;
outside = oo Exp[-mu r]/r;
matching = {(inside-outside /. r->a)==0,
  (r^2 (D[inside,r]-D[outside,r]) /. r->a)==1};
coeff = First[Solve[matching,{ii,oo}]];
amp = oo /. coeff;
check["transmission of one thickness shell", amp, Sinh[mu a]/(mu a)];
check["self response from the jump solution", outside /. coeff /. r->a,
  (1-Exp[-2 mu a])/(2 mu a^2)];
(* Reciprocity of the source and receiver shell averages gives two factors. *)
cross = amp^2 Exp[-mu dist]/dist;
check["opposite-face response with the thickness restoration removed",
  Limit[cross-1/dist,mu->0],0];
check["far difference divided by twice squared isolated load",
  Limit[dist (-((1/dist+cross)-(-1/dist+cross)))/2,dist->Infinity],-1];
check["squared perpendicular gap to quadratic order",
  SeriesCoefficient[((h+2 ep t)/Sqrt[1+ep^2 v^2]-h)^2,{ep,0,2}],4t^2];

(* Numerical response and force on each extended receiver, integrated from
   the external solution, without using the closed double-shell kernel.
   Here 4 Pi k_m=4 Pi k_t=a=lambda=k_b=A=1. c is cos(theta). *)
form[m_] := If[m==0,1,Sinh[m]/m];
radius[c_,dd_] := Sqrt[dd^2+1+2 dd c];
mean[m_,dd_] := form[m]/2 NIntegrate[
  Exp[-m radius[c,dd]]/radius[c,dd],{c,-1,1},WorkingPrecision->55];
grad[m_,dd_] := -form[m]/2 NIntegrate[
  Exp[-m radius[c,dd]] (1+m radius[c,dd]) (dd+c)/radius[c,dd]^3,
  {c,-1,1},WorkingPrecision->55];
selfValue = 1+(1-Exp[-2])/2;
p0 = 1/20;
rear = 1+(1+selfValue)p0;
load[u_] := 1+p0-1/u^2;
mutualAtLength[u_] := (rear-u)/load[u]-1-selfValue;
(* Vary the shared compliance in the relaxed slab. At fixed rim depths the
   field contributes -(load1 load2) dg. The relaxed local terms are stationary. *)
workDensity = Cancel[-load[u]^2 D[mutualAtLength[u],u]];
primitive = Integrate[workDensity,u];
Print["s1 s2 distance excess_energy/P0^2 outward_force/P0^2 load/P0"];
Do[
  {s1,s2} = signs;
  mutual = s1 s2 mean[0,dd] + mean[1,dd];
  ell = u /. FindRoot[u+(1+selfValue+mutual) load[u]==rear,{u,1},
    WorkingPrecision->45];
  pp = load[ell];
  energyChange = ((primitive /. u->ell)-(primitive /. u->1))/p0^2;
  force = (pp/p0)^2 (s1 s2 grad[0,dd]+grad[1,dd]);
  values = {energyChange, force, pp/p0};
  Print[StringRiffle[Join[ToString /@ {s1,s2,dd},
    (ToString[NumberForm[N[#,16],{15,10},NumberPadding->{"","0"},
      ExponentFunction->(Null&)]]& /@ values)]," "]];
  If[ell<=0 || Abs[mutualAtLength[ell]-mutual]>10^-34 ||
    Abs[force - (pp/p0)^2 (-s1 s2/dd^2 + (D[cross,dist] /.
      {a->1,mu->1,dist->dd}))]>10^-34,Exit[1]],
  {signs,{{1,1},{1,-1},{-1,1},{-1,-1}}},{dd,{8,4}}];
