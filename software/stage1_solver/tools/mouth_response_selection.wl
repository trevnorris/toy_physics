(* Work directly in a basis of the two distinct normal decay rates. *)
ClearAll["Global`*"];
lap[f_] := D[f, {x, 2}] + D[f, {w, 2}];
dn = Sin[x] (a2 Exp[-2 w] + a3 Exp[-3 w]);
df = u x + Sin[x] (p2 Exp[-2 w] + p3 Exp[-3 w]);
continuity = D[dn, w] + lap[df];
(* Extract coefficients with y=exp(-w). *)
coefficients[f_] := CoefficientList[Expand[(f/Sin[x] // Simplify) /. w -> -Log[y]]
                                        // PowerExpand // Expand, y];
phaseRules = First[Solve[Thread[coefficients[continuity] == 0], {p2, p3}]];
bern = D[df, w] + (35 dn - lap[dn])/24;
bulk = Simplify[{continuity, bern} /. phaseRules];
Print["potential amplitudes = ", InputForm[phaseRules]];
Print["bulk continuity and phase = ", InputForm[bulk]];

(* Solve steady material advection of the log-odds order coordinate with
   the undisturbed far profile. Match its two decaying coefficients. *)
dy = Sin[x] (s2 Exp[-2 w] + s3 Exp[-3 w]);
advective = D[df + dy, w] /. phaseRules;
orderRules = First[Solve[Thread[coefficients[advective] == 0], {s2, s3}]];
orderResidual = Simplify[advective /. orderRules];
Print["order coefficients and residual = ", InputForm[{orderRules, orderResidual}]];

normal = (dn + D[df, w]) /. phaseRules;
boundary = Simplify[{normal, dn} /. {w -> 0, x -> Pi/2}];
amplitudes = First[Solve[Thread[boundary == {jread, hread}], {a2, a3}]];
Print["[current,density] sine coefficients = ", InputForm[boundary]];
Print["selected density modes given both readings = ", InputForm[amplitudes]];
Print["unchanged-current modes = ", InputForm[Reduce[First[boundary] == 0, a2, Reals]]];
(* Full carrier Cauchy traction on the flat cut, positive normal into the
   sleeve. The mixed pressure component vanishes and the density Hessian
   contributes +(1/24) n d_x d_w log(n) in test units. *)
density = 1 + eps dn;
tx = D[density D[Log[density], x, w]/24, eps] /. eps -> 0;
ux = D[df, x] /. phaseRules;
side = FullSimplify[{ux, tx} /. w -> 0 /. amplitudes];
Print["mouth [lateral velocity,carrier Cauchy traction] = ", InputForm[side]];
Print["at zero local current change = ", InputForm[side /. jread -> 0]];
chi[y_] := 1/(1 + Exp[y]);
field = chi[w - b + eps dy];
orderFlux = D[kappa D[field, x] D[field, w], eps] /. eps -> 0;
normalizedOrderFlux = FullSimplify[(orderFlux /. w -> 0) / chi'[-b]^2
                                    /. orderRules /. amplitudes];
Print["conditional order flux divided by isolated order slope squared = ",
      InputForm[normalizedOrderFlux]];
residual = Simplify[(boundary /. amplitudes) - {jread, hread}];
Print["matched current and density residuals = ", InputForm[residual]];
If[bulk =!= {0, 0} || orderResidual =!= 0 || residual =!= {0, 0}, Exit[1]];
