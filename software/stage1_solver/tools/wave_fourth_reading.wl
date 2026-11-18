(* Independent fourth reading from the archived benchmark's time-domain
   displacement. Integrate directly at the withheld rate, using a unit push. *)
ClearAll["Global`*"];
$HistoryLength = 0;
archive = Import["docs/_measurements/trapped_wave_push.json", "RawJSON"];
line = First[archive["runs"]["wolfram"]["stdout_claim_lines"]];
archived = "(z^2 - eta*z^2)/2 + eta*(1 - Cos[z])";
If[!StringContainsQ[line, "InputForm[{" <> archived <> ","], Exit[1]];
rim[z_] := (z^2-eta*z^2)/2 + eta*(1-Cos[z]);
(* For force 1 and a zero-displacement, zero-velocity start, m(p) is
   the reciprocal of p^3 times the Laplace transform of the displacement. *)
integral = Integrate[Exp[-3 u z] rim[z], {z,0,Infinity},
                     Assumptions -> u>0 && 0<eta<1];
value = Factor[1/((3 u)^3 integral)];
Print["independent fourth from the Wolfram time history: ", InputForm[value]];
If[!TrueQ[FullSimplify[value == (1+9 u^2)/(1-eta+9 u^2),
   u>0 && 0<eta<1]],Exit[1]];
Exit[0];
