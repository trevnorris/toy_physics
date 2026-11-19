(* Transverse plane waves use Exp[I (k z - w t)] and q = k^2.
   Each calculation uses one transverse polarization. All inertias and
   bare elastic/Maxwell coefficients are positive. The couplings are real. *)
ClearAll["Global`*"];
emit[label_, value_] := Print[label, " = ", ToString[value, InputForm]];
checks = {};
compare[label_, lhs_, rhs_] := Module[{residual},
  residual = Factor[Together[lhs - rhs]];
  emit["check " <> label <> " {left,right,difference}", {lhs, rhs, residual}];
  AppendTo[checks, residual === 0]
];
linearRoot[polynomial_, variable_] :=
  Factor[-Coefficient[polynomial, variable, 0]/Coefficient[polynomial, variable, 1]];

(* Order at zero and its coefficient for a rational expression. *)
leading[expr_, variable_] := Module[{r = Cancel[expr], num, den, n, d},
  If[r === 0, Return[{Infinity, 0}]];
  num = Expand[Numerator[r]];
  den = Expand[Denominator[r]];
  n = Exponent[num, variable, Min];
  d = Exponent[den, variable, Min];
  {n - d, Factor[Coefficient[num, variable, n]/Coefficient[den, variable, d]]}
];
at[poly_, variable_, order_] := Factor[Coefficient[Expand[poly], variable, order]];

(* Hamilton's equations with the velocities and canonical momenta eliminated.
   Source terms are retained in the Hamiltonian and in the returned RHS. *)
phaseSystem[coords_, momenta_, velocities_, ham_] := Module[
  {velocityLaw, pv, dotp, eqs, mat, rhs, n = Length[coords]},
  velocityLaw = D[ham, #] & /@ momenta;
  pv = momenta /. First[Solve[Thread[velocityLaw == velocities], momenta]];
  dotp = Table[Sum[D[pv[[i]], coords[[j]]] (-I w coords[[j]]) +
    D[pv[[i]], velocities[[j]]] (-w^2 coords[[j]]), {j, n}], {i, n}];
  eqs = Expand[dotp + ((D[ham, #] & /@ coords) /. Thread[momenta -> pv] /.
      Thread[velocities -> -I w coords])];
  mat = Table[Coefficient[eqs[[i]], coords[[j]]], {i,n}, {j,n}];
  rhs = -(eqs /. Thread[coords -> ConstantArray[0,n]]);
  {pv, mat, rhs}
];

Print["Definitions: q=k^2>0; y=w^2; E_T=I w A; sources {j,F} couple to {A,u}."];
Print["The spectral weight Z of a diagonal coordinate response is the coefficient of 1/(yPole-y)."];
Print["Electric-current response E/j carries the additional factor I w compared with A/j."];

(* At fixed longitudinal electric displacement disp the Gauss constraint
   gives phi / freeCharge = (E_L / disp) / q. *)
longH = (disp - coupling u)^2/(2 eps) + (s + bulk q) u^2/2;
uStatic = linearRoot[D[longH,u], u];
electricStatic = Factor[D[longH,disp] /. u -> uStatic];
greenStatic = Factor[electricStatic/(disp q)];
emit["longitudinal {fixed-D energy,u,E,phi/freeCharge}", {longH,uStatic,electricStatic,greenStatic}];
emit["coupling(q)",g0+f q];
emit["general phi/freeCharge",greenStatic /. coupling->g0+f q];
staticCases = {
   {"unpinned local",{s->0,coupling->g0}},
   {"unpinned gradient",{s->0,coupling->f q}},
   {"pinned local",{coupling->g0}}
};
Do[Module[{kernel, apart, ir},
  kernel = Factor[greenStatic /. item[[2]]];
  apart = Apart[kernel,q];
  ir = Factor[Cancel[q kernel] /. q->0];
  emit[item[[1]] <> " {static kernel,partial fractions,Coulomb coefficient}",{kernel,apart,ir}];
  compare[item[[1]] <> " static Gauss and polarization",
    Together[(eps electricStatic + coupling uStatic) /. item[[2]]], disp]
], {item,staticCases}];
zeroStiffH = longH /. {s->0,bulk->0};
zeroStiffU = linearRoot[D[zeroStiffH,u],u];
emit["strict s=bulk=0, coupling(q)!=0 at q>0: {energy,minimizer,E,phi/freeCharge}",
  {zeroStiffH,zeroStiffU,Factor[D[zeroStiffH,disp] /. u->zeroStiffU],
    Factor[(D[zeroStiffH,disp] /. u->zeroStiffU)/(disp q)]}];
zeroStiffDecoupled = zeroStiffH /. coupling->0;
emit["strict s=bulk=coupling(q)=0 at q>0: {energy,u stationarity,E,phi/freeCharge}",
  {zeroStiffDecoupled,D[zeroStiffDecoupled,u],D[zeroStiffDecoupled,disp],
    Factor[D[zeroStiffDecoupled,disp]/(disp q)]}];

ham = (pa + coupling u)^2/(2 eps) + pu^2/(2 rho) +
    nu q a^2/2 + (s + mu q) u^2/2 - j a - force u;
{mom, m, source} = phaseSystem[{a,u},{pa,pu},{va,vu},ham];
det = Expand[Det[m]] /. {w^4->y^2,w^2->y};
adj = Map[Factor, Det[m] Inverse[m], {2}];
emit["transverse {Hamiltonian,momenta,Fourier matrix,sources}", {ham,mom,m,source}];
emit["transverse determinant in y",det];
emit["response numerators, rows {E,u}, columns {j,F}, shared denominator det",
  Map[Factor,DiagonalMatrix[{I w,1}].adj,{2}]];

(* Mass-normalized squared-frequency data. a0 and b0 are the two uncoupled
   squared frequencies at this q; eta is the squared normalized coupling. *)
poly = Expand[(det/(eps rho)) /. {nu -> alpha eps,mu -> beta rho,
  s -> tau rho,coupling^2 -> eta eps rho}];
pm = Expand[poly /. {alpha->a0/q,tau->b0-beta q}];
ne = Factor[adj[[1,1]]/rho /. w^2->y /. {mu->beta rho,s->rho (b0-beta q)}];
nm = Factor[adj[[2,2]]/eps /. w^2->y /. nu->a0 eps/q];
emit["{alpha,beta,tau,eta} definitions",{nu/eps,mu/rho,s/rho,coupling^2/(eps rho)}];
emit["{monic dispersion, eps A/j numerator, rho u/F numerator}",{pm,ne,nm}];
roots = y /. Solve[pm==0,y];
emit["squared-frequency roots",roots];
weights[yr_, rules_, variable_] := Map[
  leading[(# /. rules /. y->yr),variable]&,
  {-ne/D[pm,y],-nm/D[pm,y]}];
emit["{electric,mechanical} normalized weights at yPole",
  {-ne/D[pm,y],-nm/D[pm,y]} /. y->yPole];
emit["u/E on a source-free eigenmode",Factor[-m[[2,1]]/(I w m[[2,2]]) /. w^2->y]];
emit["P(0),P(a0),P(b0)",Factor[pm /. y->#]& /@ {0,a0,b0}];
compare["canonical matrix inverse 11",Factor[(m.adj)[[1,1]]],Factor[Det[m]]];
compare["canonical matrix inverse 12",Factor[(m.adj)[[1,2]]],0];
compare["monic dispersion physical normalization",Expand[pm /. {a0->nu q/eps,
    b0->(s+mu q)/rho,eta->coupling^2/(eps rho)}],det/(eps rho)];

Print["Asymptotic pairs {power,coefficient} below are for q=k^2 unless k is explicitly given."];
Print["Weights are ordered {eps times electric coordinate weight, rho times mechanical coordinate weight}."];
(* Nonzero gap at q=0: continue the simple zero and obtain the other
   branch from the trace, so tiny oscillator strengths remain available. *)
Do[Module[{name=item[[1]], r, p, t1, t2, low, high, yl, yh},
  r = {a0->alpha q,b0->item[[2]]+beta q,eta->gam};
  p = pm /. r;
  t1 = linearRoot[at[p /. y->v1 q,q,1],v1];
  t2 = linearRoot[at[p /. y->t1 q+v2 q^2,q,2],v2];
  low = t1 q+t2 q^2;
  high = Expand[-Coefficient[p,y,1]-low];
  yl = leading[low,q]; yh = leading[high,q];
  emit[name <> " {low y leading,high y leading,high q correction}",{yl,yh,at[high,q,1]}];
  emit[name <> " {low weights,high weights}",{weights[low,r,q],weights[high,r,q]}];
  emit[name <> " low u/E leading",leading[g0/(rho (item[[2]]+beta q-low)),q]];
  emit[name <> " high u/E leading",leading[g0/(rho (item[[2]]+beta q-high)),q]];
  compare[name <> " low characteristic to q^2",Normal[Series[p /. y->low,{q,0,2}]],0];
  compare[name <> " high characteristic to q^2",Normal[Series[p /. y->high,{q,0,2}]],0]
], {item,{{"unpinned local",0},{"pinned local",tau}}}];
emit["gam in local cases",g0^2/(eps rho)];

(* With pure gradient coupling eta=d q^2. Strata of unequal and equal
   bare speeds must be expanded independently. *)
gradRules = {a0->alpha q,b0->beta q,eta->d q^2};
gradP = pm /. gradRules;
gradAtSpeed[speed_] := Module[{coef,corr,branch},
  coef = at[gradP /. y->speed q+v2 q^2,q,3];
  corr = linearRoot[coef,v2];
  branch = speed q+corr q^2;
  emit["gradient unequal {bare speed squared,y expansion,weights,u/E leading}",
    {speed,branch,weights[branch,gradRules,q], leading[f q/(rho (beta q-branch)),q]}];
  compare["gradient unequal dispersion through q^3 at " <> ToString[speed],
    at[gradP /. y->branch,q,3],0]
];
emit["d in gradient case",f^2/(eps rho)];
Scan[gradAtSpeed,{alpha,beta}];
equalRules = {a0->alpha k^2,b0->alpha k^2,eta->d k^4};
equalP = pm /. equalRules;
equalCoef = at[equalP /. y->alpha k^2+split k^3,k,6];
equalSquare = linearRoot[Expand[equalCoef] /. split^2->split2,split2];
emit["gradient equal {y/k^2 limit, (k^3 splitting coefficient)^2}",{alpha,equalSquare}];
Do[Module[{branch = alpha k^2 + sign Sqrt[equalSquare] k^3},
  emit["gradient equal {sign,k-based weights,k-based u/E leading}",
    {sign,weights[branch,equalRules,k],leading[f k^2/(rho (alpha k^2-branch)),k]}];
  compare["gradient equal lowest secular coefficient sign " <> ToString[sign],
    at[equalP /. y->branch,k,6],0]
], {sign,{-1,1}}];

(* A physical reference field: common displacement uc is neutral, relative
   displacement pr carries the polarization. Both have canonical momenta. *)
ham2 = (pa+g0 pr)^2/(2 eps) + pc^2/(2 mc) + pp^2/(2 mp) +
  nu q a^2/2 + muc q uc^2/2 + h q uc pr + (s+mup q) pr^2/2 -
  j a - fc uc - fp pr;
{mom2,m2,source2} = phaseSystem[{a,uc,pr},{pa,pc,pp},{va,vc,vp},ham2];
det2 = Collect[Expand[Det[m2]],w,Factor] /. {w^6->y^3,w^4->y^2,w^2->y};
num2 = Expand[Det[m2[[2;;3,2;;3]]]] /. {w^4->y^2,w^2->y};
emit["reference {Hamiltonian,momenta,Fourier matrix,sources}",{ham2,mom2,m2,source2}];
emit["reference {A/j numerator,denominator}",{num2,det2}];
emit["reference q=0 denominator",Factor[det2 /. q->0]];
gap2 = linearRoot[Cancel[(det2 /. q->0)/y^2],y];
speeds2 = at[det2 /. y->v q,q,2];
emit["reference {gap squared,long-wave secular coefficient,squared speeds}",
  {gap2,speeds2,Factor[v /. Solve[speeds2==0,v]]}];
emit["reference h=0 {factored numerator,factored denominator,reduced electric-current coordinate response}",
  {Factor[num2 /. h->0],Factor[det2 /. h->0],Factor[(num2/det2) /. h->0]}];
emit["reference resultant in y of current numerator and denominator",Factor[Resultant[num2,det2,y]]];
emit["reference electric numerator at y=0",Factor[num2 /. y->0]];

(* Resolve how strongly the common acoustic mode registers in A/j.
   Degenerate speeds require a separate limit because the two low poles
   can share the leading electric weight at their crossing. *)
cCommon = muc/mc;
cElectric = nu s/(eps s+g0^2);
zweight2[branch_, variable_, rules_:{}] :=
  leading[(-num2/D[det2,y] /. rules /. y->branch),variable];
Do[Module[{coef,corr,branch},
  coef = at[det2 /. y->speed q + vv q^2,q,3];
  corr = linearRoot[coef,vv];
  branch = speed q+corr q^2;
  emit["reference unequal {bare acoustic squared speed, electric A/j weight leading}",
    {speed,zweight2[branch,q]}];
  compare["reference unequal acoustic characteristic at q^3",
    at[det2 /. y->branch,q,3],0]
],{speed,{cCommon,cElectric}}];
emit["reference gapped A/j weight at q=0",Factor[(-num2/D[det2,y]) /. q->0 /. y->gap2]];
(* Use v as the common squared speed at exact tuning, imposing it on
   both bare parameters before expanding. *)
tuned = {muc->mc v,nu->v (eps s+g0^2)/s,q->k^2};
tunedSec = at[det2 /. tuned /. y->v k^2+split k^3,k,6];
tunedSquare = linearRoot[Expand[tunedSec] /. split^2->split2,split2];
emit["reference tuned h!=0 {squared speed,(k^3 splitting coefficient)^2}",
  {v,tunedSquare}];
Do[Module[{branch = v k^2 + sign Sqrt[tunedSquare] k^3},
  emit["reference tuned {sign,A/j weight leading in k}",{sign,zweight2[branch,k,tuned]}];
  compare["reference tuned lowest secular coefficient sign " <> ToString[sign],
    at[det2 /. tuned /. y->branch,k,6],0]
], {sign,{-1,1}}];
compare["reference electric minor via decoupled common cancellation",
  Factor[num2/det2 /. h->0],(s+mup q-mp y)/((nu q-eps y)(s+mup q-mp y)-g0^2 y)];
compare["reference response satisfies the common-mode dark test",
  Factor[num2 /. y->cCommon q],-h^2 q^2];
emit["all computed consistency checks",checks];
Exit[If[And@@checks,0,1]];
