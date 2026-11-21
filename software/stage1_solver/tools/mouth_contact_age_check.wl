(* Scalar elastic contacts of stiffness mu per unit ordered fraction, in a
   material parcel. A birth records the current shear as its rest shear.
   Removal samples every live contact with the same instantaneous hazard h.
   q,m1,m2 are the zeroth, first and second moments of the rest shear,
   measured per unit original parcel population. A tag marks all entry contacts. *)
ClearAll["Global`*"];
$HistoryLength = 0;
emit[key_, value_] := Print[key, " = ", ToString[value, InputForm]];

state = {q, m1, m2, s, mark};
flow = {birth - h q, birth s - h m1, birth s^2 - h m2, shearRate, -h mark};
parcelRate[expr_] := Expand[(D[expr, #] & /@ state) . flow];
stress = mu (q s - m1);
storage = mu (q s^2 - 2 s m1 + m2)/2;
normal = stress/(mu q);
oldShare = mark/q;
emit["cohort traction, elastic energy", {stress, storage}];
emit["stress rate, energy rate from rest-shear moments", parcelRate /@ {stress, storage}];
emit["rates inferred from old count and old share",
     Factor /@ {-parcelRate[mark]/mark, -parcelRate[oldShare]/oldShare}];
emit["energy above mean-strain spring", Factor[storage - stress^2/(2 mu q)]];
identities = Factor /@ {parcelRate[stress] - (mu q shearRate - h stress),
    parcelRate[storage] - (stress shearRate - h storage),
    parcelRate[normal] - (shearRate - birth normal/q),
    parcelRate[oldShare] + birth oldShare/q};
emit["moment identities: traction, energy, mean strain, tag dilution", identities];

(* Direct crossing test in arbitrary time units. chi and the imposed two-tone
   shear are identical in both cases. Scale the small shear by its amplitude,
   traction by mu times that amplitude, energies by mu times amplitude^2.
   Entry contacts all have rest shear zero for this energy experiment.
   B=a chi and D=(1+a)chi. The nonconstant chi reaches 1/4 at flight end. *)
chi[t_] := 3 Exp[-t]/4;
flight = Log[3];
shear[t_] := Cos[2 t] + Sin[5 t]/3;
survival[t_, u_, a_] := Exp[-(1 + a) (t - u)];
newWeight[t_, u_, a_] := a chi[u] survival[t, u, a];
entryWeight[t_, a_] := chi[0] survival[t, 0, a];
quadrature[expr_, domain__] := NIntegrate[expr, domain, WorkingPrecision -> 35,
    AccuracyGoal -> 23, PrecisionGoal -> 23, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}];
moment[t_?NumericQ, a_, j_] := entryWeight[t, a] shear[t]^j +
    quadrature[newWeight[t, u, a] (shear[t] - shear[u])^j, {u, 0, t}];
(* For exp(i omega t) input, solve g' + a g = s', g(0)=s(0).
   The exponential term retains the transient from the entry condition. *)
tone[t_, a_, w_] := (I w Exp[I w t] + a Exp[-a t])/(a + I w);
predicted[t_, a_] := chi[t] (Re[tone[t, a, 2]] + Im[tone[t, a, 5]]/3);
toneResidual = FullSimplify[D[tone[t, a, w], t] + a tone[t, a, w] - I w Exp[I w t]];
emit["single-tone ODE and entry residuals", {toneResidual, Factor[tone[0, a, w] - 1]}];

check[a_] := Module[{points, traction, directError, finalMass, finalStorage,
    work, removal, residual, old, renewed, transfer},
  points = flight {1/4, 1/2, 3/4, 1};
  traction = moment[#, a, 1] & /@ points;
  directError = Max[Abs[traction - (predicted[#, a] & /@ points)]];
  finalMass = entryWeight[flight, a] + quadrature[newWeight[flight, u, a], {u, 0, flight}];
  finalStorage = moment[flight, a, 2]/2;
  work = quadrature[predicted[t, a] shear'[t], {t, 0, flight}];
  (* Sum energy removed over the entry cohort and birth cohorts. Parameterize
     birth time by u=v t so both numerical integration limits are fixed. *)
  removal = (1 + a)/2 (quadrature[entryWeight[t, a] shear[t]^2, {t, 0, flight}] +
      quadrature[t newWeight[t, v t, a] (shear[t] - shear[v t])^2,
        {t, 0, flight}, {v, 0, 1}]);
  residual = work - (finalStorage - chi[0] shear[0]^2/2 + removal);
  old = Simplify[entryWeight[flight, a]];
  renewed = 1 - old/chi[flight];
  transfer = Table[tone[flight, a, w] Exp[-I w flight], {w, {2, 5}}];
  emit["renewal parameter a", a];
  emit["exit order, old amount, old share, new share", {chi[flight], old, old/chi[flight], renewed}];
  emit["exit total mass from birth cohorts", N[finalMass, 14]];
  emit["two-tone cohort tractions at 1/4, 1/2, 3/4, 1 flight", N[traction, 14]];
  emit["exit mean-strain / complex unit tone at each frequency including entry transient", N[transfer, 14]];
  emit["maximum absolute cohort / tone superposition discrepancy", N[directError, 4]];
  emit["entry storage, exit storage, imposed work, energy removed, budget residual",
      N[{chi[0] shear[0]^2/2, finalStorage, work, removal, residual}, 14]];
  And[TrueQ[directError < 10^-18], TrueQ[Abs[finalMass - chi[flight]] < 10^-18],
      TrueQ[Abs[residual] < 10^-18], TrueQ[removal >= 0]]];
cases = check /@ {0, 1};
If[!TrueQ[identities === {0, 0, 0, 0} && toneResidual === 0 && And @@ cases], Exit[1]];
