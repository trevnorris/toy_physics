(* Define parameters *)
(* R is the instantaneous distance, v is velocity, cs is sound speed *)
(* nDotV is the projection of velocity onto the line of sight *)

(* 1. The Retarded Potential (Liénard-Wiechert) *)
(* This represents the TRUE total solution to the wave equation for a moving point source *)
PhiRetarded = -mu / (R * (1 - nDotV / cs));

(* 2. The Poisson Potential (Instantaneous) *)
(* This represents the 0PN Newtonian constraint *)
PhiPoisson = -mu / R;

(* 3. Series Expansion of the Retarded Potential *)
(* We assume v is small compared to cs *)
PhiRetardedExpanded = Series[PhiRetarded, {cs, Infinity, 2}];

Print["--- Expansion of the True Retarded Potential (PhiTotal) ---"];
Print[Normal[PhiRetardedExpanded]];

(* 4. The 'Double Counting' Problem in the current text *)
(* The text defines PhiTotal = PhiPoisson + PhiL *)
(* But later derives PhiL as the full Retarded Potential *)
PhiTotalCurrentText = PhiPoisson + Normal[PhiRetardedExpanded];

Print["\n--- PhiTotal in Current Text (Static Limit Check) ---"];
(* Check the term that does not depend on cs (the static limit) *)
StaticLimitCurrent = Select[PhiTotalCurrentText, FreeQ[#, cs] &];
Print["Static Limit (should be -mu/R): ", StaticLimitCurrent];

(* 5. The Fix: Defining PhiL as the Difference *)
(* PhiL should be the difference between the Full Retarded and the Poisson potentials *)
PhiLagCorrected = Normal[PhiRetardedExpanded] - PhiPoisson;

Print["\n--- Corrected PhiL (Lag Potential) ---"];
Print[PhiLagCorrected];

(* 6. Physical Interpretation of the Corrected Lag *)
(* Notice the leading order term is now Order(1/cs), not Order(1) *)
Print["\n--- Leading order of Corrected PhiL ---"];
Print["Does it vanish in static limit (cs->infinity)? ", Limit[PhiLagCorrected, cs -> Infinity] == 0];
