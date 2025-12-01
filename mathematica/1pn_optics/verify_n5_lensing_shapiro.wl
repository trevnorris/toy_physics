(* Verification Script for Paper II: Lensing, Shapiro, and Redshift 
   Author: AI Assistant / Trevor Norris
   Purpose: Verify the n=5 superfluid vacuum produces 1PN-compliant optics.
*)

(* Clear previous definitions to prevent namespace conflicts *)
ClearAll["Global`*"];

(* --------------------------------------------------------- *)
(* 1. SETUP: FLUID DYNAMICS AND EQUATION OF STATE (EOS)      *)
(* --------------------------------------------------------- *)

Print["--- 1. Deriving Refractive Index from EOS ---"];

(* Define variables - CamelCase used to avoid underscores *)
(* P: Pressure, rho: Density, n: Polytropic index, K: Constant *)
(* c0: Vacuum sound speed, G: Gravitational constant, M: Mass *)
(* r: Radial coordinate *)

(* The General Polytropic EOS *)
eqnState = P == K * rho^n;

(* Squared sound speed definition: cs^2 = dP/drho *)
(* We use a symbol cs0sq for the background sound speed to keep algebra clean *)
(* cs0sq = K * n * rho0^(n-1) *)

(* The Flux Tube Condition: Pressure gradient balances Newtonian tension *)
(* dP/dr = rho * G * M / r^2 *)
(* Delta P = - G * M * rho0 / r *)
pressureDeficit = - (G * M * rho0) / r;

(* Relate Pressure perturbation to Density perturbation: dP = cs0^2 * drho *)
(* deltaRho = Delta P / cs0sq *)
(* We use the symbol cs0sq here to represent c_s0^2 *)
densityPerturbation = pressureDeficit / cs0sq;
fractionalDensityChange = Simplify[densityPerturbation / rho0];

Print["Fractional Density Perturbation (delta rho / rho):"];
Print[fractionalDensityChange];
(* Expected: -GM / (r cs0^2) *)


(* --------------------------------------------------------- *)
(* 2. REFRACTIVE INDEX DERIVATION FOR GENERAL n              *)
(* --------------------------------------------------------- *)

(* Sound speed scales as rho^((n-1)/2) *)
(* fracDeltaCs = (n-1)/2 * fracDeltaRho *)
fractionalSoundSpeedChange = ((n - 1) / 2) * fractionalDensityChange;

(* Refractive Index N(r) = c0 / cs(r) ~ 1 - deltaCs/c0 *)
(* N(r) = 1 - fractionalSoundSpeedChange *)
refractiveIndexPerturbation = - fractionalSoundSpeedChange;

(* We want to find alpha such that N(r) = 1 + alpha * GM / (r cs0^2) *)
(* Solve for alpha by dividing the perturbation by the target term *)
termOfInterest = (G * M) / (r * cs0sq);
alphaCoeff = Simplify[refractiveIndexPerturbation / termOfInterest];

Print["Derived Refractive Index N(r) perturbation term:"];
Print[refractiveIndexPerturbation];

Print["Coefficient alpha (in terms of n):"];
Print[alphaCoeff];
(* Expected: (n-1)/2 *)


(* --------------------------------------------------------- *)
(* 3. LENSING INTEGRAL (LIGHT BENDING)                       *)
(* --------------------------------------------------------- *)

Print["\n--- 2. Lensing Integral ---"];

(* Geometric optics geometry *)
(* Path z from -infinity to +infinity *)
(* Impact parameter b *)
(* r = Sqrt[b^2 + z^2] *)

(* The gradient of ln(N) transverse to the path *)
(* grad_perp(N) ~ dN/db *)
(* N(b, z) ~ 1 + alpha * GM / (c0^2 * Sqrt[b^2 + z^2]) *)
(* dN/db = - alpha * GM * b / (c0^2 * (b^2 + z^2)^(3/2)) *)

(* We substitute alphaCoeff into the integrand *)
integrandLensing = - alphaCoeff * (G * M * b) / (cs0sq * (b^2 + z^2)^(3/2));

(* Perform the integral over z from -infinity to infinity *)
deflectionAngle = Integrate[Abs[integrandLensing], {z, -Infinity, Infinity}, Assumptions -> {b > 0, G > 0, M > 0, cs0sq > 0}];

Print["Calculated Deflection Angle (Delta Theta) for general n:"];
Print[deflectionAngle];

(* Standard GR result: 4GM / (b c^2) *)
(* We equate our result to GR to find n *)
targetGR = (4 * G * M) / (b * cs0sq);

(* Updated Solve to avoid inverse function warnings by restricting n to physical values > 1 *)
solutionN = Solve[{deflectionAngle == targetGR, n > 1}, n];

Print["\n--- RESULT: Value of n required to match GR Lensing ---"];
Print[solutionN];

(* Verify n=5 explicitly *)
n5Check = deflectionAngle /. n -> 5;
Print["Deflection Angle for n=5:"];
Print[n5Check];
Print["Does it match GR? "];
Print[Simplify[n5Check == targetGR]];


(* --------------------------------------------------------- *)
(* 4. SHAPIRO DELAY (TIME OF FLIGHT)                         *)
(* --------------------------------------------------------- *)

Print["\n--- 3. Shapiro Delay Check ---"];

(* Shapiro delay is integral of (N - 1) dz *)
(* N - 1 = alpha * GM / (r c0^2) *)
(* For n=5, alpha should be 2 *)

alphaN5 = alphaCoeff /. n -> 5;
integrandShapiro = alphaN5 * (G * M) / (cs0sq * Sqrt[b^2 + z^2]);

(* Integrate from -L to L (finite distance to source/receiver) *)
(* We assume L >> b *)
timeDelay = Integrate[integrandShapiro, {z, -L, L}, Assumptions -> {L > 0, b > 0, L > b}];

Print["Shapiro Time Delay for n=5 (integrated from -L to L):"];
Print[timeDelay];

(* Check for the characteristic Log term *)
(* The dominant term for L >> b is Log[4 L^2 / b^2] = 2 Log[2L/b] *)
(* Or roughly 2 * Log[L/b] depending on exact limits *)
(* GR expression usually has 4 GM/c^3 * Log[...]. *)
(* Our coefficient is 2 * alpha * GM/c^3. *)

Print["Coefficient of Log term (GR expects 4):"];
Print[2 * alphaN5];


(* --------------------------------------------------------- *)
(* 5. REDSHIFT / TRILEMMA CHECK                              *)
(* --------------------------------------------------------- *)

Print["\n--- 4. Redshift & Mass Scaling ---"];

(* Redshift z = delta_omega / omega *)
(* If omega ~ m (internal clock mass) and m ~ rho (cavitation) *)
(* Then z ~ delta_rho / rho *)

redshiftFraction = fractionalDensityChange;
Print["Redshift (delta rho / rho):"];
Print[redshiftFraction];

(* GR Expectation: z = - Delta Phi / c^2 = - (-GM/r) / c^2 = GM / (r c^2) *)
(* Note: The sign convention in the derivation: *)
(* Delta P = - GM rho / r. *)
(* Pressure drops near mass. Density drops near mass. *)
(* Clock rate (mass) decreases near mass. *)
(* Lower mass => Lower frequency => Redshift (relative to infinity). *)

(* Check coefficient *)
(* We found delta_rho/rho = - GM / (r cs^2) earlier. *)
(* Magnitude matches GR weak field potential exactly. *)

Print["Does magnitude match GR potential GM/(r c^2)? "];
(* Added Assumptions to allow Simplify to process the Abs[] correctly *)
Print[Simplify[Abs[redshiftFraction] == (G * M) / (r * cs0sq), Assumptions -> {G > 0, M > 0, r > 0, cs0sq > 0}]];

Print["\n--- CONCLUSION ---"];
Print["1. n=5 is uniquely required to match the factor of 4 in Lensing."];
Print["2. This sets alpha=2 in the refractive index."];
Print["3. alpha=2 yields the correct factor of 4 in Shapiro delay."];
Print["4. The density drop matches the gravitational redshift potential."];

(*"
Output:

--- 1. Deriving Refractive Index from EOS ---
Fractional Density Perturbation (delta rho / rho):
-((G*M)/(cs0sq*r))
Derived Refractive Index N(r) perturbation term:
(G*M*(-1 + n))/(2*cs0sq*r)
Coefficient alpha (in terms of n):
(-1 + n)/2

--- 2. Lensing Integral ---
Calculated Deflection Angle (Delta Theta) for general n:
(G*M*Abs[1 - n])/(b*cs0sq)

--- RESULT: Value of n required to match GR Lensing ---
{{n -> 5}}
Deflection Angle for n=5:
(4*G*M)/(b*cs0sq)
Does it match GR?
True

--- 3. Shapiro Delay Check ---
Shapiro Time Delay for n=5 (integrated from -L to L):
(4*G*M*Log[(L + Sqrt[b^2 + L^2])/b])/cs0sq
Coefficient of Log term (GR expects 4):
4

--- 4. Redshift & Mass Scaling ---
Redshift (delta rho / rho):
-((G*M)/(cs0sq*r))
Does magnitude match GR potential GM/(r c^2)?
True
"*)
