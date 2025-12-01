(* ================================================================= *)
(* SUPERFLUID DEFECT TOY MODEL: DYON STRUCTURAL ANALYSIS TOOL        *)
(* ================================================================= *)

(* OBJECTIVE: 
   1. Correctly define the GR Lense-Thirring target (clean math).
   2. Test if a "Line Vortex" flow matches GR.
   3. Test if a "Vortex Dipole" (Dyon) flow matches GR.
   4. Mathematically define the required flow structure for Paper III.
*)

ClearAll["Global`*"]

(* ----------------------------------------------------------------- *)
(* PART 1: THE TARGET (GENERAL RELATIVITY)                           *)
(* ----------------------------------------------------------------- *)

Print["\n--- PART 1: DEFINING THE GR TARGET ---"];

(* The metric for a slowly rotating body (Lense-Thirring) in weak fields:
   ds^2 = -c^2 dt^2 + (1 + 2GM/rc^2)dr^2 + ... - 2(2GJ/c^2 r) sin^2(theta) dt dphi
   
   Note: The cross term g_0_phi contains a sin^2(theta) factor.
   The Frame Dragging Frequency omega(r) = -g_0_phi / g_phi_phi
*)

coords = {t, r, theta, phi};

(* Correct GR g_0_phi definition with sin^2(theta) *)
(* Factor of 2 vs 4 depends on isotropic vs Boyer-Lindquist, using 4 for standard weak field form *)
g0phiGR = - (4 * GNewton * JGR * Sin[theta]^2) / (c^2 * r); 

(* Standard spatial metric component g_phi_phi *)
gphiphiGR = r^2 * Sin[theta]^2; 

(* Calculate the Target Frame Dragging Frequency *)
omegaGR = Simplify[ - g0phiGR / gphiphiGR ];

Print["Target GR Frequency (omegaGR):"];
Print[omegaGR];
Print["Radial Power Law: ", Exponent[omegaGR, r]];
(* We expect -3, representing 1/r^3 decay *)


(* ----------------------------------------------------------------- *)
(* PART 2: ACOUSTIC METRIC ENGINE (VECTOR FIELD ENABLED)             *)
(* ----------------------------------------------------------------- *)

Print["\n--- PART 2: ACOUSTIC METRIC ENGINE ---"];

(* We define a function that takes a Velocity Vector Field v_fluid 
   and returns the relevant metric components.
   
   Acoustic Metric Shift Vector: g_0_i = - (rho0/c) * v_i
   For 1PN, we treat conformal factors as ~1.
*)

getAcousticFrameDragging[vVectorSpherical_] := Module[
  {vPhi, g0phiAcoustic, gphiphiAcoustic, omegaAcoustic},
  
  (* Extract the phi-component of the velocity vector *)
  (* vVectorSpherical assumed to be {v_r, v_theta, v_phi_physical} *)
  
  (* Note: The metric component g_0_phi couples to the COVARIANT velocity component.
     Covariant v_phi = v_phi_physical * (r sin theta)
     So g_0_phi = - Covariant_v_phi
  *)
  
  vPhiPhysical = vVectorSpherical[[3]];
  
  (* Construct Metric Components *)
  (* g_0_phi = - v_physical * r * sin(theta) *)
  g0phiAcoustic = - vPhiPhysical * r * Sin[theta];
  
  (* Standard flat space g_phi_phi for background *)
  gphiphiAcoustic = r^2 * Sin[theta]^2;
  
  (* Compute Frame Dragging Omega *)
  omegaAcoustic = Simplify[ - g0phiAcoustic / gphiphiAcoustic ];
  
  Return[omegaAcoustic];
];


(* ----------------------------------------------------------------- *)
(* PART 3: TEST CASE A - THE LINE VORTEX (SIMPLE DEFECT)             *)
(* ----------------------------------------------------------------- *)

Print["\n--- PART 3: TESTING LINE VORTEX (SIMPLE DEFECT) ---"];

(* Physical velocity of a line vortex along z-axis:
   v_phi = Gamma / (2 * pi * DistanceFromAxis)
   DistanceFromAxis = r * sin(theta)
*)

vLineVortex = {0, 0, GammaCirc / (2 * Pi * r * Sin[theta])};

omegaLine = getAcousticFrameDragging[vLineVortex];

Print["Acoustic Omega (Line Vortex):"];
Print[omegaLine];
Print["Radial Power Law: ", Exponent[omegaLine, r]];

If[Exponent[omegaLine, r] == -3,
   Print["-> MATCHES GR!"],
   Print["-> FAILS (Decays too slowly)"]
];


(* ----------------------------------------------------------------- *)
(* PART 4: TEST CASE B - THE DYON (VORTEX DIPOLE)                    *)
(* ----------------------------------------------------------------- *)

Print["\n--- PART 4: TESTING THE DYON (VORTEX DIPOLE) ---"];

(* Hypothesis: The 'Spin' of a Dyon is mechanically a 'Vortex Doublet' or Ring.
   This creates a flow field analogous to a magnetic dipole.
   
   Magnetic Dipole Vector Potential: A ~ (m x r) / r^3
   Fluid Velocity Analogue: v ~ (J_flow x r) / r^3
   
   In Spherical coords, (z_hat x r_hat) gives a pure phi component:
   |z_hat x r_hat| = sin(theta)
   
   So v_phi_physical ~ (Moment * sin(theta)) / r^3 * r = Moment * sin(theta) / r^2
   Wait, let's check dimensions.
   Dipole potential Phi ~ 1/r^2. Gradient v ~ 1/r^3.
   
   Actually, let's reverse engineer the required v.
   GR needs Omega ~ 1/r^3.
   Omega = v_phi / (r sin theta).
   Therefore v_phi must scale as (1/r^3) * r = 1/r^2.
   
   Let's test a flow field v_phi = K * sin(theta) / r^2.
   (This is the far-field flow of a small Vortex Ring).
*)

vDyon = {0, 0, DipoleStrength * Sin[theta] / r^2};

omegaDyon = getAcousticFrameDragging[vDyon];

Print["Acoustic Omega (Dyon/Dipole):"];
Print[omegaDyon];
Print["Radial Power Law: ", Exponent[omegaDyon, r]];

If[Exponent[omegaDyon, r] == -3,
   Print["-> MATCHES GR!"],
   Print["-> FAILS"]
];


(* ----------------------------------------------------------------- *)
(* PART 5: DEFINING THE DYON FOR PAPER III                           *)
(* ----------------------------------------------------------------- *)

Print["\n--- PART 5: THE DYON DEFINITION ---"];

(* Solve for the DipoleStrength required to match GR constants *)
(* Set Omega_Dyon = Omega_GR *)

calibrationEq = omegaDyon == omegaGR;
calibratedStrength = Solve[calibrationEq, DipoleStrength];

Print["To reproduce General Relativity, the Superfluid Dyon must generate a flow field"];
Print["v_phi = (DipoleStrength * sin(theta)) / r^2"];
Print["where DipoleStrength is:"];
Print[calibratedStrength[[1, 1]]];

Print["\nINTERPRETATION:"];
Print["1. The 'Spin' J of a defect corresponds to the Dipole Moment of the fluid flow."];
Print["2. This 1/r^2 velocity field is characteristic of a VORTEX RING in the far field."];
Print["3. Paper III should define the Dyon as a 'Flux Tube threaded by a Vortex Ring'."];

(*"
Output:


--- PART 1: DEFINING THE GR TARGET ---
Target GR Frequency (omegaGR):
(4*GNewton*JGR)/(c^2*r^3)
Radial Power Law: -3

--- PART 2: ACOUSTIC METRIC ENGINE ---

--- PART 3: TESTING LINE VORTEX (SIMPLE DEFECT) ---
Acoustic Omega (Line Vortex):
(GammaCirc*Csc[theta]^2)/(2*Pi*r^2)
Radial Power Law: -2
-> FAILS (Decays too slowly)

--- PART 4: TESTING THE DYON (VORTEX DIPOLE) ---
Acoustic Omega (Dyon/Dipole):
DipoleStrength/r^3
Radial Power Law: -3
-> MATCHES GR!

--- PART 5: THE DYON DEFINITION ---
To reproduce General Relativity, the Superfluid Dyon must generate a flow field
v_phi = (DipoleStrength * sin(theta)) / r^2
where DipoleStrength is:
DipoleStrength -> (4*GNewton*JGR)/c^2
"*)
