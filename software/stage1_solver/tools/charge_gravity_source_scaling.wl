(* Compare the source ratios by eliminating the two mouth amplitudes.
   All Q and electric inputs are positive magnitudes in a common normalization. *)
ClearAll["Global`*"];
$HistoryLength = 0;
assumptions = a > 0 && b > 0 && k > 0 && unit > 0 && m1 > 0 && m2 > 0;
drainMatch = First[Solve[{a x1 == k m1, a x2 == k m2}, {x1, x2}]];
electricRatio = FullSimplify[(b x2)/(b x1) /. drainMatch, assumptions];
matchingMasses = FullSimplify[Reduce[electricRatio == 1 && m2 > 0, m2, Reals],
  assumptions];
gainRatio = FullSimplify[(bSecond/b /. First[Solve[
  (bSecond x2 == b x1) /. drainMatch, bSecond]]), assumptions];
(* With a fixed first electric reading, predict the second without a new fit. *)
firstMatch = First[Solve[(b x1 == unit) /. drainMatch, b]];
secondReading = FullSimplify[b x2 /. drainMatch /. firstMatch, assumptions];
Print["electric magnitude ratio using one common amplitude = " <> ToString[electricRatio, InputForm]];
Print["equal units require = " <> ToString[matchingMasses, InputForm]];
Print["electric gain ratio for equal units = " <> ToString[gainRatio, InputForm]];
Print["second magnitude after calibrating the first = " <> ToString[secondReading, InputForm]];
If[!TrueQ[FullSimplify[electricRatio == m2/m1 && Equivalent[matchingMasses,m2 == m1] &&
    gainRatio == m1/m2 && secondReading == unit m2/m1, assumptions]], Exit[1]];
Exit[0];
