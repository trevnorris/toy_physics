(* Independent pressure-inclusive control-volume and finite-annulus checks. *)
ClearAll["Global`*"];
$HistoryLength=0;
emit[label_, value_]:=Print[label," = ",ToString[value,InputForm]];
check[label_,actual_,expected_]:=Module[{v=Factor[actual]},emit[label,v];
  If[Factor[v-expected]=!=0,Exit[1]]];
(* Use exit speed v. A_i=A, A_o=m A and uniform static head. *)
md=rho v m area;
pin = rho (v^2-(m v)^2)/2;
ti = md m v + area pin;
to = md v;
check["positive difference of directed thrusts",ti-to,rho area v^2(m-1)^2/2];
tangent={1,sl};
inward={-sl,1}/Sqrt[1+sl^2]; outlet={0,1};
result=ti inward-to outlet;
check["tangent derivative, outlet follows independent reference",
  D[tangent.result,sl]/.sl->0,-to];
check["bent scalar strength in units of ell",result[[2]]/.sl->0,ti-to];
check["co-rotating outlet and inlet",tangent.((ti-to)inward),0];
check["add upstream elbow carried by same mouth",
  D[tangent.(result+ti (outlet-inward)),sl]/.sl->0,ti-to];
check["pressure included in transported head difference",md ((m v)^2/2+pin/rho-v^2/2),0];

(* Solve r^2 f''+2r f'-(2+V r^2) f=0 as a Taylor series in (r-a).
   Both edges are within its convergence radius. Check order refinement
   and the differential residual at the midpoint and the far edge. *)
seriesGain[aa_,order_]:=Module[{bb=aa+1/20,span,pot,c,n,term,edge,norm,res},
  span=bb^3-aa^3; pot=N[3/(4 Pi span (7/8)),55];
  c[-2]=0; c[-1]=0; c[0]=aa; c[1]=1;
  Do[c[n+2]=(-2 aa (n+1)^2 c[n+1]+(2+pot aa^2-n(n+1)) c[n]+
    2 pot aa c[n-1]+pot c[n-2])/(aa^2 (n+1)(n+2)),{n,0,order-2}];
  term=Sum[c[n] t^n,{n,0,order}]; edge=(term/.t->bb-aa);
  res=(aa+t)^2 D[term,{t,2}]+2 (aa+t)D[term,t]-(2+pot (aa+t)^2)term;
  If[!TrueQ[Max[Abs[(res/.t->#)]& /@ {1/40,1/20}]<10^-10],Exit[1]];
  norm=3/((D[term,t]/.t->bb-aa)+2 edge/bb);
  norm (bb^2 edge-aa^3)/span];
gain[aa_]:=Module[{u=seriesGain[aa,32],coarse=seriesGain[aa,18]},
  emit["radius "<>ToString[aa]<>", annulus slope gain",N[u,14]];
  If[!TrueQ[Abs[u-coarse]<10^-15],Exit[1]];u];
(* Previous guarded isolated-load solves, with the same annulus and
   global normalization. Read only the two SymPy radius output lines. *)
archive=Import["docs/_measurements/unequal_mouth_electric.json","RawJSON"];
lines=archive["runs"]["annular_mouth_radius_sympy"]["stdout_claim_lines"];
z=Table[Module[{line,match},line=SelectFirst[lines,StringStartsQ[#,"radius "<>ToString[i]<>", escape = "]&];
  match=StringCases[line,RegularExpression["escape = ([0-9.]+),"]->"$1"];
  If[Length[match]!=1,Exit[1]];ToExpression[First[match]]],{i,2}];
u=gain/@{1,2};
ratios=1-u/z ((ti/(ti-to))/.m->2//Simplify);
emit["area ratio 2, signed response/emission",N[ratios,14]];
emit["required mass-rate ratio to equalize emission",Sqrt[8 First[z]/Last[z]]];
emit["equal emission, mixed-force ratio for two targets",First[ratios]/Last[ratios]];
If[!TrueQ[0<Last[ratios]/First[ratios]<1 && Max[ratios]<0],Exit[1]];
Exit[0];
