(* throat_spec_derivation.wl
   PURPOSE:
   1) Import the frozen Paper 7 operational definitions from definitions.wl (single source of truth).
   2) Derive the constrained stationary state equation (fixed N) in operator form.
   3) Derive the exact force-balance (Hellmann–Feynman) conditions for (a, L).
   4) Run pointwise + slice diagnostics for localization / “orthogonality” of the force kernels.
   5) Define the standardized leakage flux operator J_w and the measurement protocol text.

   NOTE:
   - This script enforces the 6-argument Paper 7 hard-mode signature:
       Vconf[x,y,z,w,a,L]
     and will Abort[] if that does not evaluate to an expression containing a and L.
*)

ClearAll["Global`*"];

Get[FileNameJoin[{DirectoryName[$InputFileName], "definitions.wl"}]];

debug = False; (* set True when diagnosing definition drift *)

(* ------------------------------------------------------------------------- *)
(* Helpers *)
(* ------------------------------------------------------------------------- *)

IsDefinedSym[s_Symbol] := (OwnValues[s] =!= {} || DownValues[s] =!= {} || UpValues[s] =!= {});

RequireDef[name_String] := Module[{sym = ToExpression[name]},
  If[Head[sym] =!= Symbol, Print["[ERROR] Bad symbol name: ", name]; Abort[]];
  If[IsDefinedSym[sym], Print["[OK] ", name], Print["[ERROR] Missing required definition: ", name]; Abort[]];
];

Print["=== Paper 7: Spec Derivation & Diagnostics (uses frozen definitions.wl) ==="];

(* Print freeze summary/signature if available *)
If[IsDefinedSym[Paper7FreezeSummary], Paper7FreezeSummary[], Print["[WARN] Paper7FreezeSummary[] not found."]];
If[debug,
  If[IsDefinedSym[Paper7FreezeSignature], Print["Freeze signature: ", Paper7FreezeSignature[]], Print["[NOTE] Paper7FreezeSignature[] not found (optional)."]];
];

(* Debug: confirm Vconf signature *)
If[debug && IsDefinedSym[Vconf],
  Print["Debug: Vconf DownValues (first): ", Short[DownValues[Vconf], 4]]
];

(* Required core objects *)
Print["-- 1. Verifying required frozen definitions --"];
RequireDef["Vconf"];
RequireDef["Egeom"];
Print["[DONE] Frozen definitions loaded."];

(* Quick Γ sanity: ensure frozen quadrature helpers are wired *)
If[IsDefinedSym[GammaIntegrate],
  Block[{rPort = 1, Anum, Aexact, relErr},
    Anum = GammaIntegrate[Function[{th, ph}, 1]];
    Aexact = 4 Pi rPort^2;
    relErr = N[Anum/Aexact - 1];
    Print["[Gamma] area check rel error: ", relErr];
  ],
  Print["[WARN] GammaIntegrate not found; Γ helpers may be missing."]
];

(* ========================================================================= *)
(* 2. CONSTRAINED STATIONARY STATE & FORCE BALANCE *)
(* ========================================================================= *)
Print["-- 2. Constrained Stationary State & Force Balance --"];

coords = {x, y, z, w};

(* Physical assumptions (safe even if some symbols are unused) *)
physAssumps = {
  a > 0, L > 0, dr > 0, dw > 0, p > 0, epsW > 0,
  m > 0, hbar > 0, K > 0,
  V0 > 0, OmOut > 0, OmIn > 0
};

(* Stationary fields: use independent conjugate symbol for symbolic safety *)
P0  = psi0[x, y, z, w];
P0b = psi0b[x, y, z, w];
rho0 = P0*P0b;

Lap4D = Laplacian[P0, coords];

(* Enforce hard-mode Vconf signature and evaluation *)
Vexplicit = Evaluate @ Vconf[x, y, z, w, a, L];

If[debug,
  Print["Debug: Vexplicit contains a? ", !FreeQ[Vexplicit, a]];
  Print["Debug: Vexplicit contains L? ", !FreeQ[Vexplicit, L]];
  Print["Debug: Vexplicit preview: ", Short[Vexplicit, 6]];
];

If[FreeQ[Vexplicit, a] || FreeQ[Vexplicit, L],
  Print["[ERROR] Vconf[x,y,z,w,a,L] did not evaluate to an expression containing a and L."];
  Print["        This would force kernA/kernL = 0 (invalid for Paper 7)."];
  Print["        Check definitions.wl and avoid NumericQ restrictions on a/L patterns."];
  Abort[];
];

(* GNLS operator (n=5 EOS baseline: (5K/4) rho^4 psi ) *)
Hpsi = -(hbar^2/(2*m))*Lap4D + Vexplicit*P0 + (5*K/4)*rho0^4*P0;

EqStationary = (Hpsi - mu*P0 == 0);

Print["Constrained stationary equation (fixed N), operator form:"];
Print["  -(hbar^2/(2 m)) ∇^2 ψ0 + Vconf(...) ψ0 + (5K/4) |ψ0|^8 ψ0 = μ ψ0"];

(* Force kernels: fluid force F_q = -∫ rho * ∂_q Vconf dX *)
kernA = -D[Vexplicit, a];
kernL = -D[Vexplicit, L];

(* Geometry energy is frozen as Egeom[a,L] in definitions.wl *)
EgeomExpr = Egeom[a, L];
forceGeomA = -D[EgeomExpr, a];
forceGeomL = -D[EgeomExpr, L];

Print["Force balance (radius a):"];
Print["  Integral[ |psi|^2 * (", Short[kernA], ") dX ] + (", forceGeomA, ") == 0"];

Print["Force balance (length L):"];
Print["  Integral[ |psi|^2 * (", Short[kernL], ") dX ] + (", forceGeomL, ") == 0"];

Print["Sanity checks (PossibleZeroQ after simplification):"];
Print["  kernA zero?      ", PossibleZeroQ[Simplify[kernA, Assumptions -> physAssumps]]];
Print["  kernL zero?      ", PossibleZeroQ[Simplify[kernL, Assumptions -> physAssumps]]];
Print["  forceGeomA zero? ", PossibleZeroQ[Simplify[forceGeomA, Assumptions -> physAssumps]]];
Print["  forceGeomL zero? ", PossibleZeroQ[Simplify[forceGeomL, Assumptions -> physAssumps]]];

(* ========================================================================= *)
(* 3. LOCALIZATION / ORTHOGONALITY DIAGNOSTICS *)
(* ========================================================================= *)
Print["-- 3. Localization Diagnostics (pointwise + slice) --"];

EvaluateAtPoint[expr_, rVal_, wVal_] := Module[{subbed},
  subbed = expr /. {x -> rVal, y -> 0, z -> 0, w -> wVal};
  Simplify[subbed /. {Tanh[0] -> 0, Sech[0] -> 1}, Assumptions -> physAssumps]
];

(* Pointwise sensitivity test points *)
rWall = a;   wWall = 0;
rCap  = 0;   wCap  = L/2;

valKaWall = EvaluateAtPoint[kernA, rWall, wWall];
valKLWall = EvaluateAtPoint[kernL, rWall, wWall];

valKaCap  = EvaluateAtPoint[kernA, rCap, wCap];
valKLCap  = EvaluateAtPoint[kernL, rCap, wCap];

Print["Pointwise sensitivity:"];
Print["  At lateral wall   (r=a, w=0)   -> dV/dL ~ ", Short[valKLWall]];
Print["  At endcap region  (r=0, w=L/2) -> dV/da ~ ", Short[valKaCap]];

(* Slice localization: compare center vs cap for L-kernel at r=0 *)
sliceLForce = kernL /. {x -> 0, y -> 0, z -> 0};

valCenterL = Simplify[sliceLForce /. w -> 0, Assumptions -> physAssumps];
valCapL    = Simplify[sliceLForce /. w -> (L/2), Assumptions -> physAssumps];

Print["Slice localization check (r=0):"];
Print["  L-kernel at center (w=0):   ", Short[valCenterL, 6]];
Print["  L-kernel at cap    (w=L/2): ", Short[valCapL, 6]];

ratioExpr = FullSimplify[valCenterL/valCapL, Assumptions -> physAssumps];
Print["  Ratio (center/cap): ", Short[ratioExpr, 6]];
Print["  Note: Expect ratio << 1 when (L/2 - epsW) >> dw (tight axial localization)."];

(* Optional numeric pass/fail if parameters are numeric in the current session *)
If[And @@ (NumericQ /@ {L, epsW, dw, V0, p, m, OmIn, OmOut, a}),
  ratioNum = Quiet @ N[ratioExpr];
  Print["  Ratio numeric: ", ratioNum];
  If[NumericQ[ratioNum] && Abs[ratioNum] < 10^-6,
    Print["[PASS] Tight axial localization (|center/cap| < 1e-6)."],
    Print["[NOTE] Axial localization not tight under current numeric parameters."]
  ];
];

(* ========================================================================= *)
(* 4. LEAKAGE FLUX OPERATOR (STANDARDIZED) *)
(* ========================================================================= *)
Print["-- 4. Leakage Flux Definition --"];

(* Generic dynamic fields (matches master_checks convention) *)
P  = psi[x, y, z, w];
Pb = psib[x, y, z, w];

JwOp = (hbar/(2*I*m)) * (Pb*D[P, w] - P*D[Pb, w]);

Print["Leakage operator JwOp (standard form):"];
Print["  ", JwOp];

Print["Measurement protocol (standardized):"];
Print["  Evaluate JwOp on the two slices w = ±Wcut, then integrate over the 3-ball R3 < Rmeasure."];
Print["  (Wcut and Rmeasure are part of the Paper 7 operational definitions / run config.)"];

Print["=== throat_spec_derivation.wl complete ==="];

(*"
Output:

=== Paper 7: Spec Derivation & Diagnostics (uses frozen definitions.wl) ===
=== Paper 7 Operational Definitions (Executable) ===
Version: definitions.wl v1.0
Primary Vconf: F1_modulated_brane_trap
Geometry DOFs: {a,L} only (baseline)
Axis conventions: Throat axis = w, brane plane = w=0
Constraint: FixedNormN
Projection weight: W(w)=|χ0(w)|^2 (Ω_out ground state)
Projection window: Numerics: |w| <= WprojCut = wProjNSig * ellW
Projection renorm: Use WwProj(w)=W(w)/WprojNorm on |w|<=WprojCut; tail mass = 1-WprojNorm
Gamma: R3=rPort, w=0
Gamma measure: dμ = r_port^2 sinθ dθ dφ
Gamma normal: outward radial on Γ
Gamma quadrature: θ: Gauss-Legendre (nTheta=32), φ: uniform (nPhi=64)
Measurement region: Γ: sphere R3=rPort at w=0
Flux primary: leakage, outward-oriented: Jw_out = sign(w) Jw at w=±Wcut over R3<Rmeasure
Ports: Y_lm complex basis, ordering = PortsList[lmax], lmax baseline = 2
Drive: Vconf -> Vconf + eps cos(ωt) exp(-w^2/wDrive^2) exp(-(R3-rPort)^2/rDrive^2) Y_lm(theta(x,y,z),phi(x,y))
Response operator: j_mouth(ω)=Zeff(ω) u(ω); leakage j_leak reported separately
Response extraction: Lock-in: discard 5 periods, measure 10 periods; nw=32; Zeff via EstimateZeffRobust (CondThresh=10000000000 RidgeLambda=0.0 SVTol=0.000000000001)
Zeff cond thresh: 10000000000
Zeff ridge lambda: 0.
Zeff sv tol: 1/1000000000000
Freeze Signature (SHA256): 39121031452442260063295664310682675073183234823340385570560837303905828468514
-- 1. Verifying required frozen definitions --
[OK] Vconf
[OK] Egeom
[DONE] Frozen definitions loaded.
[Gamma] area check rel error: -5.003900183098153*^-9
-- 2. Constrained Stationary State & Force Balance --
Constrained stationary equation (fixed N), operator form:
  -(hbar^2/(2 m)) ∇^2 ψ0 + Vconf(...) ψ0 + (5K/4) |ψ0|^8 ψ0 = μ ψ0
Force balance (radius a):
  Integral[ |psi|^2 * (Short[(m*(-OmIn^2 + OmOut^2)*w^2*Sech[(-a + Sqrt[x^2 + y^2 + z^2])/dr]^2*(1 + (-1 - Tanh[(-1/2*L + Sqrt[epsW^2 + w^2])/dw])/2))/(4*dr) + (p*V0*Sech[(-a + Sqrt[x^2 + y^2 + z^2])/dr]^2*(1 + Tanh[(-a + Sqrt[x^2 + y^2 + z^2])/dr])^(-1 + p))/(2^p*dr)]) dX ] + (-4*a^2*L*Pi*Pvac - (8*a^2*Pi + 8*a*L*Pi)*Sig) == 0
Force balance (length L):
  Integral[ |psi|^2 * (Short[(2^(-1 - p)*p*V0*Sech[(-1/2*L + Sqrt[epsW^2 + w^2])/dw]^2*(1 + Tanh[(-1/2*L + Sqrt[epsW^2 + w^2])/dw])^(-1 + p))/dw + (m*(-OmIn^2 + OmOut^2)*w^2*Sech[(-1/2*L + Sqrt[epsW^2 + w^2])/dw]^2*(1 + (-1 - Tanh[(-a + Sqrt[x^2 + y^2 + z^2])/dr])/2))/(8*dw)]) dX ] + ((-4*a^3*Pi*Pvac)/3 - 4*a^2*Pi*Sig) == 0
Sanity checks (PossibleZeroQ after simplification):
  kernA zero?      False
  kernL zero?      False
  forceGeomA zero? False
  forceGeomL zero? False
-- 3. Localization Diagnostics (pointwise + slice) --
Pointwise sensitivity:
  At lateral wall   (r=a, w=0)   -> dV/dL ~ Short[(2^(-1 - p)*p*V0*Sech[(-2*epsW + L)/(2*dw)]^2*(1 + Tanh[(2*epsW - L)/(2*dw)])^(-1 + p))/dw]
  At endcap region  (r=0, w=L/2) -> dV/da ~ Short[(Sech[a/dr]^2*(2^(4 - p)*p*V0*(1 - Tanh[a/dr])^(-1 + p) + (L^2*m*(OmIn^2 - OmOut^2)*(-1 + Tanh[(-L + Sqrt[4*epsW^2 + L^2])/(2*dw)]))/2))/(16*dr)]
Slice localization check (r=0):
  L-kernel at center (w=0):   Short[(2^(-1 - p)*p*V0*Sech[(-2*epsW + L)/(2*dw)]^2*(1 + Tanh[(2*epsW - L)/(2*dw)])^(-1 + p))/dw, 6]
  L-kernel at cap    (w=L/2): Short[(Sech[(L - Sqrt[4*epsW^2 + L^2])/(2*dw)]^2*(-1/2*(L^2*m*(OmIn^2 - OmOut^2)*(1 + Tanh[a/dr])) + 2^(4 - p)*p*V0*(1 + Tanh[(-L + Sqrt[4*epsW^2 + L^2])/(2*dw)])^(-1 + p)))/(32*dw), 6]
  Ratio (center/cap): Short[(2^(4 - p)*p*V0*Cosh[(L - Sqrt[4*epsW^2 + L^2])/(2*dw)]^2*Sech[(-2*epsW + L)/(2*dw)]^2*(1 + Tanh[(2*epsW - L)/(2*dw)])^(-1 + p))/(-1/2*(L^2*m*(OmIn - OmOut)*(OmIn + OmOut)*(1 + Tanh[a/dr])) + 2^(4 - p)*p*V0*(1 + Tanh[(-L + Sqrt[4*epsW^2 + L^2])/(2*dw)])^(-1 + p)), 6]
  Note: Expect ratio << 1 when (L/2 - epsW) >> dw (tight axial localization).
-- 4. Leakage Flux Definition --
Leakage operator JwOp (standard form):
  ((-1/2*I)*hbar*(psib[x, y, z, w]*Derivative[0, 0, 0, 1][psi][x, y, z, w] - psi[x, y, z, w]*Derivative[0, 0, 0, 1][psib][x, y, z, w]))/m
Measurement protocol (standardized):
  Evaluate JwOp on the two slices w = ±Wcut, then integrate over the 3-ball R3 < Rmeasure.
  (Wcut and Rmeasure are part of the Paper 7 operational definitions / run config.)
=== throat_spec_derivation.wl complete ===
"*)
