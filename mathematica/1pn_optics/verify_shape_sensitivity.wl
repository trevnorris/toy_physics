(* Verification Script for Paper II: Shape Sensitivity and WEP
   Author: AI Assistant / Trevor Norris
   Purpose:
     1. Quantify how sensitive the model is to the 'Sphere' assumption.
     2. Show that a cylindrical (flux-tube-like) hydrodynamic shape fails.
     3. Demonstrate that WEP is naturally preserved if the vacuum permeates
        matter, but generically violated if it treats a planet as a solid body
        unless densities are fine-tuned.
*)

ClearAll["Global`*"];

(* ================================================ *)
Print["--- 1. SHAPE SENSITIVITY (THE STIFFNESS CONSTRAINT) ---"];
(* ================================================ *)

(* We model the defect as an Ellipsoid of revolution with semi-axes (a, b, c).
   Motion is along the 'a' axis.
   If a = b = c, it is a Sphere.
   If a < b = c, it is Prolate (Cigar).
   If a > b = c, it is Oblate (Pancake).

   Hydrodynamic Added Mass Coefficient k for motion along axis 'a':
   k = alpha / (2 - alpha)
   where alpha is a geometric integral depending on eccentricity.
*)

(* Define the geometric integral alpha for an ellipsoid.
   ratio = b/a
   Note: Formulas differ for Prolate vs Oblate. We implement both. *)

kappaEllipsoid[ratioBA_] := Module[{e, alpha, k},
  If[ratioBA == 1,
    Return[0.5]; (* Sphere *)
  ];

  If[ratioBA < 1,
    (* Prolate (Cigar): b < a *)
    e = Sqrt[1 - ratioBA^2];
    alpha = (2 (1 - e^2)/e^3) (Log[(1 + e)/(1 - e)]/2 - e);,
    (* Oblate (Pancake): b > a *)
    e = Sqrt[1 - 1/ratioBA^2];
    alpha = (2/e^3) (e - Sqrt[1 - e^2] ArcSin[e]);
  ];

  k = alpha/(2 - alpha);
  k
];

(* Generate data for shapes ranging from Needle (0.1) to Sphere (1.0) to Pancake (10.0) *)
ratios = {0.1, 0.5, 0.9, 1.0, 1.1, 2.0, 10.0};
kappas = kappaEllipsoid /@ ratios;

Print[""];
Print["Table: Added Mass Coefficient (kappa) vs Shape (b/a ratio)"];
Print["Ratio < 1: Prolate (Tube-like)"];
Print["Ratio = 1: Sphere"];
Print["Ratio > 1: Oblate (Flat)"];
Print["-------------------------------------------------"];
Do[
  Print[
    "b/a = ",
    N[ratios[[i]], {3, 1}],
    "  ->  kappa = ",
    N[kappas[[i]], {3, 3}]
  ],
  {i, Length[ratios]}
];

(* Check Precession Deviation for 10% deformation *)
kappaSphere   = 0.5;
kappaDeformed = kappaEllipsoid[1.1]; (* 10% oblate *)
betaSphere    = 1 + kappaSphere;     (* beta = 1.5 *)
betaDeformed  = 1 + kappaDeformed;
precessionFactorSphere   = 3 + 2 betaSphere;   (* = 6.0 *)
precessionFactorDeformed = 3 + 2 betaDeformed;

Print[""];
Print["--- SENSITIVITY RESULT ---"];
Print["GR Target Factor: 6.000"];
Print["Spherical Defect Factor: ", N[precessionFactorSphere, {4, 3}]];
Print["10% Oblate Defect Factor: ", N[precessionFactorDeformed, {4, 3}]];
Print["Percent Error in Precession: ",
  N[100 (precessionFactorDeformed - 6)/6, {3, 2}],
  "%"
];


(* ================================================ *)
Print[""];
Print["--- 2. FLUX TUBE LIMIT (CYLINDER ANALYSIS) ---"];
(* ================================================ *)

(* Analytic limits for a Cylinder moving perpendicular to axis (Broadside)
   and parallel to axis (End-on). These are standard added-mass limits. *)

kappaCylBroad = 1.0;
kappaCylEnd   = 0.0;

betaBroad = 1 + kappaCylBroad;
betaEnd   = 1 + kappaCylEnd;

precessionBroad = 3 + 2 betaBroad; (* = 7 *)
precessionEnd   = 3 + 2 betaEnd;   (* = 5 *)

Print["If the defect behaves like a Cylinder segment in 3D:"];
Print["Broadside Motion (kappa = 1.0) -> beta = 2.0 -> Precession Factor = ",
  precessionBroad, " (Target 6)"];
Print["End-on Motion    (kappa = 0.0) -> beta = 1.0 -> Precession Factor = ",
  precessionEnd, " (Target 6)"];
Print["Conclusion: A cylindrical hydrodynamic profile cannot match the 1PN precession;"];
Print["the defect must be effectively spherical on the brane."];


(* ================================================ *)
Print[""];
Print["--- 3. WEAK EQUIVALENCE PRINCIPLE (WEP) TEST ---"];
(* ================================================ *)

(* Goal:
   Compare three cases:

   (1) Single defect in the superfluid.
   (2) Planet composed of N defects, but PERMEABLE: superfluid flows around
       each defect individually (sum of micro added masses).
   (3) Planet composed of N defects, but IMPERMEABLE: the superfluid treats
       the planet as one solid body of macroscopic density rhoMatter.

   We show:
     - Case (2) has exactly the same acceleration as a single defect.
     - Case (3) has a different acceleration unless rhoMatter is tuned
       to a special value, i.e., WEP is generically violated.
*)

ClearAll[Ndefects, mDefect, rhoVac, rhoMatter, vDefect,
         kappaSingle, kappaMacro, gField];

(* Parameters (kept symbolic for transparency) *)
Ndefects    = Symbol["Ndefects"];     (* number of defects in the planet *)
mDefect     = Symbol["mDefect"];      (* bare mass of one defect *)
rhoVac      = Symbol["rhoVac"];       (* superfluid vacuum density *)
rhoMatter   = Symbol["rhoMatter"];    (* bulk macroscopic density of the planet *)
vDefect     = Symbol["vDefect"];      (* effective void volume per defect *)
kappaSingle = Symbol["kappaSingle"];  (* added-mass coefficient for single defect *)
kappaMacro  = Symbol["kappaMacro"];   (* added-mass coefficient for solid planet *)
gField      = Symbol["gField"];       (* external gravitational field strength *)

assumptions = {
  Ndefects > 0, mDefect > 0, rhoVac > 0, rhoMatter > 0,
  vDefect > 0, kappaSingle > 0, kappaMacro > 0, gField > 0
};

(* (1) Single defect *)
accelSingle =
  (mDefect*gField) / (mDefect + kappaSingle*rhoVac*vDefect);

(* (2) Permeable planet: sum of micro added masses *)
accelPermeable =
  Simplify[
    (Ndefects*mDefect*gField) /
    (Ndefects*mDefect + Ndefects*kappaSingle*rhoVac*vDefect),
    Assumptions -> assumptions
  ];

(* (3) Solid planet: one macroscopic body of density rhoMatter
   Total bare mass M = Ndefects*mDefect
   Macroscopic volume Vplanet = M / rhoMatter = Ndefects*mDefect / rhoMatter
   Added mass = kappaMacro*rhoVac*Vplanet
*)
accelSolid =
  Simplify[
    (Ndefects*mDefect*gField) /
    (Ndefects*mDefect + kappaMacro*rhoVac*(Ndefects*mDefect/rhoMatter)),
    Assumptions -> assumptions
  ];

(* Normalize by gField *)
accelSingleNorm =
  Simplify[accelSingle/gField, Assumptions -> assumptions];
accelPermeableNorm =
  Simplify[accelPermeable/gField, Assumptions -> assumptions];
accelSolidNorm =
  Simplify[accelSolid/gField, Assumptions -> assumptions];

Print["Acceleration Ratios (normalized to gField):"];
Print["Single Defect:      a/g = ", accelSingleNorm];
Print["Permeable Planet:   a/g = ", accelPermeableNorm];
Print["Solid Planet:       a/g = ", accelSolidNorm];

(* Check equality conditions *)
wepPermeable =
  Simplify[accelSingleNorm == accelPermeableNorm, Assumptions -> assumptions];

rhoMatterSolution =
  Simplify[
    Solve[accelSingleNorm == accelSolidNorm, rhoMatter],
    Assumptions -> assumptions
  ];

Print[""];
Print["--- WEP ANALYSIS ---"];
Print["Do Single Defect and Permeable Planet match exactly? ", wepPermeable];

Print["Condition for Solid Planet to match Single Defect:"];
Print[rhoMatterSolution];

(* Express the tuning condition in terms of the effective defect density *)
rhoDefectEff =
  Simplify[mDefect/vDefect, Assumptions -> assumptions];

rhoMatterRequired =
  Simplify[
    rhoMatter /. First[rhoMatterSolution],
    Assumptions -> assumptions
  ];

ratioRequired =
  Simplify[rhoMatterRequired/rhoDefectEff, Assumptions -> assumptions];

Print["Effective defect density: rho_defect = mDefect / vDefect = ", rhoDefectEff];
Print["Required rhoMatter for exact equality: rhoMatter = ", rhoMatterRequired];
Print["So rhoMatter / rho_defect must equal: ", ratioRequired];

Print[""];
Print["Interpretation:"];
Print["  - If the vacuum permeates matter (case 2), a composite planet"];
Print["    falls with the same acceleration as a single defect (WEP OK)."];
Print["  - If the vacuum treats the planet as a solid body (case 3),"];
Print["    the acceleration depends on the macroscopic density rhoMatter."];
Print["    Exact agreement with the single-defect acceleration requires a"];
Print["    fine-tuned relation between rhoMatter and the defect density."];
Print["  - Thus, in this hydrodynamic model, WEP favors a permeating"];
Print["    (ghost-like) superfluid vacuum over a solid-obstacle picture."];

(*"
Output:

--- 1. SHAPE SENSITIVITY (THE STIFFNESS CONSTRAINT) ---

Table: Added Mass Coefficient (kappa) vs Shape (b/a ratio)
Ratio < 1: Prolate (Tube-like)
Ratio = 1: Sphere
Ratio > 1: Oblate (Flat)
-------------------------------------------------
b/a = 0.1  ->  kappa = 0.02070591807721217
b/a = 0.5  ->  kappa = 0.2100150489766415
b/a = 0.9  ->  kappa = 0.4402768490684319
b/a = 1.  ->  kappa = 0.5
b/a = 1.1  ->  kappa = 0.5602400845244107
b/a = 2.  ->  kappa = 1.115060485695698
b/a = 10.  ->  kappa = 6.184128758091638

--- SENSITIVITY RESULT ---
GR Target Factor: 6.000
Spherical Defect Factor: 6.
10% Oblate Defect Factor: 6.120480169048822
Percent Error in Precession: 2.0080028174803624%

--- 2. FLUX TUBE LIMIT (CYLINDER ANALYSIS) ---
If the defect behaves like a Cylinder segment in 3D:
Broadside Motion (kappa = 1.0) -> beta = 2.0 -> Precession Factor = 7. (Target 6)
End-on Motion    (kappa = 0.0) -> beta = 1.0 -> Precession Factor = 5. (Target 6)
Conclusion: A cylindrical hydrodynamic profile cannot match the 1PN precession;
the defect must be effectively spherical on the brane.

--- 3. WEAK EQUIVALENCE PRINCIPLE (WEP) TEST ---
Acceleration Ratios (normalized to gField):
Single Defect:      a/g = mDefect/(mDefect + kappaSingle*rhoVac*vDefect)
Permeable Planet:   a/g = mDefect/(mDefect + kappaSingle*rhoVac*vDefect)
Solid Planet:       a/g = rhoMatter/(rhoMatter + kappaMacro*rhoVac)

--- WEP ANALYSIS ---
Do Single Defect and Permeable Planet match exactly? True
Condition for Solid Planet to match Single Defect:
{{rhoMatter -> (kappaMacro*mDefect)/(kappaSingle*vDefect)}}
Effective defect density: rho_defect = mDefect / vDefect = mDefect/vDefect
Required rhoMatter for exact equality: rhoMatter = (kappaMacro*mDefect)/(kappaSingle*vDefect)
So rhoMatter / rho_defect must equal: kappaMacro/kappaSingle
"*)
