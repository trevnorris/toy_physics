(* Use the local complex carrier and its normal derivatives directly.
   At an intake of density one remove the common tangential/time phase
   factor. The normal derivative is pure imaginary when density slope is zero.
   Determine the next normal derivative from the stationary GNLS. *)
ClearAll["Global`*"];
mouthField = 1;
normalDerivative = I s;
secondDerivative = c + I d;
(* hbar=m=1, h(1)=1, V=0 and mu=3/2. *)
residual = -(secondDerivative - k^2 mouthField)/2 + (1 - 3/2) mouthField;
normalCurvature = First[Solve[Thread[ComplexExpand[{Re[residual], Im[residual]}] == {0, 0}],
                              {c, d}]];
current = ComplexExpand[Im[Conjugate[mouthField] normalDerivative]];
densityCurvature = ComplexExpand[
  2 Re[Conjugate[mouthField] secondDerivative + Conjugate[normalDerivative] normalDerivative]
] /. normalCurvature;
currentDerivative = ComplexExpand[Im[Conjugate[normalDerivative] normalDerivative +
                    Conjugate[mouthField] secondDerivative]] /. normalCurvature;
Print["normal carrier curvature = ", InputForm[secondDerivative /. normalCurvature]];
Print["normal current and its derivative = ", InputForm[{current, currentDerivative}]];
Do[
  {label, tangentialK, slope} = entry;
  Print[label, " [tangential k, inward current, density curvature] = ",
        InputForm[{k, current, densityCurvature} /. {k -> tangentialK, s -> slope}]],
  {entry, {{"isolated", 0, 1/2}, {"incident, continuation A", 1/4, 1/2},
           {"incident, continuation B", 1/4, 9/16}}}];
