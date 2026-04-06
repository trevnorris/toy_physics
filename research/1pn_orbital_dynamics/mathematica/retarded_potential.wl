(* Define parameters *)
(* R is the instantaneous distance, v is velocity, cs is sound speed *)
(* nDotV is the projection of velocity onto the line of sight *)

ClearAll["Global`*"];

(* 1. The Retarded Potential (Liénard-Wiechert) *)
(* This is the full scalar retarded solution for a moving point source. *)
PhiRetarded = -mu / (R * (1 - nDotV / cs));

(* 2. The Poisson Potential (Instantaneous) *)
(* This represents the 0PN Newtonian constraint *)
PhiPoisson = -mu / R;

(* In the central-field test-mass limit, the source is static:
   nDotV -> 0, so PhiRetarded == PhiPoisson exactly. *)
PhiStatic = Simplify[PhiRetarded /. nDotV -> 0];
Print["Static-source check: PhiRetarded(static) = ", PhiStatic];
Print["(Should match PhiPoisson = ", PhiPoisson, ")"];

(* 3. Series Expansion of the Retarded Potential *)
(* We assume v is small compared to cs *)
PhiRetardedExpanded = Normal[Series[PhiRetarded, {cs, Infinity, 2}]];

Print["--- Expansion of the full retarded scalar potential PhiRetarded ---"];
Print[PhiRetardedExpanded];

(* 4. Paper-consistent lag definition *)
(* The paper defines PhiL as the residual between the retarded and *)
(* Poisson solutions, so that Phi = PhiPoisson + PhiL = PhiRetarded. *)
PhiLag = Simplify[PhiRetardedExpanded - PhiPoisson];

Print["\n--- Paper definition: PhiL = PhiRetarded - PhiPoisson ---"];
Print[PhiLag];

PhiTotalPaper = Simplify[PhiPoisson + PhiLag];

Print["\n--- Reconstructed total potential Phi = PhiPoisson + PhiL ---"];
Print[PhiTotalPaper];

StaticLimitPaper = Simplify[PhiTotalPaper /. nDotV -> 0];
Print["Static-source limit from reconstructed Phi: ", StaticLimitPaper];

Print["\n--- Static-limit checks ---"];
Print["PhiL vanishes for static source? ", Simplify[PhiLag /. nDotV -> 0] == 0];
Print["PhiL vanishes as cs->infinity? ", Limit[PhiLag, cs -> Infinity] == 0];

(*"
Output:
Static-source check: PhiRetarded(static) = -(mu/R)
(Should match PhiPoisson = -(mu/R))
--- Expansion of the full retarded scalar potential PhiRetarded ---
-(mu/R) - (mu*nDotV)/(cs*R) - (mu*nDotV^2)/(cs^2*R)

--- Paper definition: PhiL = PhiRetarded - PhiPoisson ---
-((mu*nDotV*(cs + nDotV))/(cs^2*R))

--- Reconstructed total potential Phi = PhiPoisson + PhiL ---
-((mu*(cs^2 + cs*nDotV + nDotV^2))/(cs^2*R))
Static-source limit from reconstructed Phi: -(mu/R)

--- Static-limit checks ---
PhiL vanishes for static source? True
PhiL vanishes as cs->infinity? True
"*)
