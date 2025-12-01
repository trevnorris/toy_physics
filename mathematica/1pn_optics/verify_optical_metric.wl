(* Clear existing definitions to avoid namespace conflicts *)
ClearAll["Global`*"]

(* ================================================================= *)
(* PART 1: DEFINITIONS AND SETUP *)
(* ================================================================= *)
Print["--- PART 1: DEFINITIONS ---"]

(* Define constants *)
(* G: Gravitational Constant, M: Mass, c: Speed of Light *)
(* We use epsilon as a bookkeeping parameter for 1PN order expansion *)
epsilon = (G * M) / (c^2 * r);

(* Define the Refractive Index N(r) for the n=5 stiff vacuum *)
(* From Section 4 of your paper *)
nIndex[r_] := 1 + 2 * (G * M) / (c^2 * r)

(* ================================================================= *)
(* PART 2: THE OPTICAL METRIC SIGN CORRECTION *)
(* ================================================================= *)
Print["\n--- PART 2: OPTICAL METRIC & LIGHT SPEED ---"]

(* The Reviewer suggests the spatial part must be N^2, not 1/N^2.
   Let's verify the coordinate speed of light for both cases.
   We assume a radial null geodesic: ds^2 = 0.
*)

(* Case A: Your original Draft (1/N^2) *)
metricSpatialOriginal = 1 / nIndex[r]^2;
(* ds^2 = -c^2 dt^2 + (1/N^2) dr^2 = 0 *)
vLightOriginal = Solve[-c^2 + metricSpatialOriginal * v^2 == 0, v][[2]];
(* We take the positive root *)

(* Case B: The Corrected Metric (N^2) *)
metricSpatialCorrected = nIndex[r]^2;
(* ds^2 = -c^2 dt^2 + (N^2) dr^2 = 0 *)
vLightCorrected = Solve[-c^2 + metricSpatialCorrected * v^2 == 0, v][[2]];

Print["Coordinate speed of light (dr/dt):"]
Print["Original (1/N^2): ", Simplify[Series[vLightOriginal[[1, 2]], {G, 0, 1}]]]
Print["Corrected (N^2):  ", Simplify[Series[vLightCorrected[[1, 2]], {G, 0, 1}]]]

(* RESULT INTERPRETATION:
   Original -> c * (1 + 2GM/rc^2) -> c * N -> Speed INCREASES near mass. (Incorrect)
   Corrected -> c * (1 - 2GM/rc^2) -> c / N -> Speed DECREASES near mass. (Correct)
*)


(* ================================================================= *)
(* PART 3: EXTRACTING PPN GAMMA FROM OPTICAL METRIC *)
(* ================================================================= *)
Print["\n--- PART 3: PPN GAMMA EXTRACTION ---"]

(* Standard PPN Isotropic Metric Expansion:
   g_rr = 1 + 2 * gamma * (GM / rc^2)
*)

(* Expand the Corrected Spatial Metric to 1st order *)
expandedSpatial = Normal[Series[metricSpatialCorrected, {G, 0, 1}]];

Print["Expanded Spatial Metric (g_rr): ", expandedSpatial]

(* Extract Gamma *)
(* We equate the coefficient of (GM/rc^2) to 2*gamma *)
term1PN = Coefficient[expandedSpatial, (G * M)/(c^2 * r)];
gammaOptical = term1PN / 2;

Print["Implied PPN Gamma for Optical Metric: ", gammaOptical]

(* RESULT EXPECTATION:
   N^2 ~ (1 + 2e)^2 ~ 1 + 4e.
   4e = 2 * gamma * e  =>  gamma = 2.
*)


(* ================================================================= *)
(* PART 4: THE '10 vs 6' CALCULATION *)
(* ================================================================= *)
Print["\n--- PART 4: PRECESSION COEFFICIENT (10 vs 6) ---"]

(* Standard Formula for Perihelion Precession in PPN:
   dPhi = (6 * Pi * GM / (a(1-e^2)c^2)) * ((2 - beta + 2*gamma) / 3)

   For the 'Naive Optical Test Particle' Model:
   1. We insert Newtonian Potential manually -> beta = 1
   2. We use the Optical Spatial Metric -> gamma = 2 (calculated above)
*)

betaNaive = 1;
gammaNaive = 2;

precessionFactor = (2 - betaNaive + 2 * gammaNaive) / 3;
totalPrecessionCoeff = 6 * precessionFactor; (* relative to Pi * epsilon *)

Print["Naive Optical Parameters: beta = ", betaNaive, ", gamma = ", gammaNaive]
Print["Precession Factor ((2-beta+2gamma)/3): ", precessionFactor]
Print["Total Precession Coefficient (Factor * 6): ", totalPrecessionCoeff]

(* RESULT EXPECTATION:
   (2 - 1 + 4) / 3 = 5/3
   5/3 * 6 = 10.
   The coefficient is 10.
*)


(* ================================================================= *)
(* PART 5: PSI RATIOS (SOLITON HYPOTHESIS) *)
(* ================================================================= *)
Print["\n--- PART 5: PSI RATIOS ---"]

(* Psi is defined via: g_rr = 1 + 2*Psi/c^2

   Psi_opt comes from Optical Metric (gamma = 2)
   Psi_orb comes from Orbital beta (beta = 3/2)
*)

(* Optical Psi *)
(* g_rr = 1 + 4 epsilon *)
(* 2 * Psi / c^2 = 4 * GM / rc^2 *)
(* Psi = 2 * GM / r *)
psiOpticalCoeff = 2;

(* Orbital Psi *)
(* From Paper I: Sigma = 3/2 * epsilon *)
(* 1 + Sigma = 1 + 2 * Psi / c^2 *)
(* Psi = 1/2 * Sigma * c^2 = 1/2 * (3/2) * GM/r = 3/4 * GM/r *)
psiOrbitalCoeff = 3/4;

ratio = psiOrbitalCoeff / psiOpticalCoeff;

Print["Psi_Optical Coefficient: ", psiOpticalCoeff]
Print["Psi_Orbital Coefficient: ", psiOrbitalCoeff]
Print["Ratio (Orbital / Optical): ", ratio]

(*"
Output:

--- PART 2: OPTICAL METRIC & LIGHT SPEED ---
Coordinate speed of light (dr/dt):
Original (1/N^2): SeriesData[G, 0, {c, (2*M)/(c*r)}, 0, 2, 1]
Corrected (N^2):  SeriesData[G, 0, {c, (-2*M)/(c*r)}, 0, 2, 1]

--- PART 3: PPN GAMMA EXTRACTION ---
Expanded Spatial Metric (g_rr): 1 + (4*G*M)/(c^2*r)
Implied PPN Gamma for Optical Metric: 2

--- PART 4: PRECESSION COEFFICIENT (10 vs 6) ---
Naive Optical Parameters: beta = 1, gamma = 2
Precession Factor ((2-beta+2gamma)/3): 5/3
Total Precession Coefficient (Factor * 6): 10

--- PART 5: PSI RATIOS ---
Psi_Optical Coefficient: 2
Psi_Orbital Coefficient: 3/4
Ratio (Orbital / Optical): 3/8
"*)

