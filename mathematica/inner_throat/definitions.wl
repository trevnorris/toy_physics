(* ::Package:: *)

(* ============================================================================
   definitions.wl  (Paper 7 Hard-Mode 4D)
   ----------------------------------------------------------------------------
   PURPOSE
     Canonical, machine-checkable operational definitions for Paper 7.

   DESIGN GOAL
     - Executable mirror of inner_throat_freeze_sheet.md (Sections A–K).
     - Defines baseline confinement potential family (F1), wave confinement mu2A,
       geometry measures/energy, brane projection weight W(w), measurement region Γ,
       port basis, drive protocol, and response-extraction helpers.

   IMPORTANT PRACTICE
     Many project scripts begin with:
       ClearAll["Global`*"];
     That will erase any previously loaded definitions. Therefore:
       Add  Get["./definitions.wl"];
     immediately AFTER the ClearAll line in each script that uses Vconf/mu2A.

   HOW TO USE
     In any WL script (e.g., master_checks.wl), add near the top:
       Get["./definitions.wl"];

     Then you may call:
       Vconf[x,y,z,w,a,L]
       mu2A[x,y,z,w,a,L]
       Egeom[a,L]
       Ww[w]
       etc.

   PARAMETER CONVENTION
     This file NEVER overwrites numeric assignments. Symbols such as
       hbar, m, K, c, dr, dw, p, epsW, V0, OmOut, OmIn, muIn, muOut,
       Pvac, Sig, rPort, Rmeasure, Wcut, wDrive, rDrive, epsAng
     are referenced but NOT assigned here. If you want numbers, set them
     in your driver script after loading this file.

   ========================================================================== *)

ClearAll[
  SmoothStep, SmoothAbs,
  R3, ThetaXYZ, PhiXYZ,
  GateR, GateW, GateThroat,
  OmegaSq, Vbrane, VwallRadial, VcapAxial, Vconf,
  mu2A,
  V4Tube, A4Tube, Egeom, dEgeomda, dEgeomdL,
  ellW, Chi0, Ww,
  WprojCut, WprojNorm, WwProj, WprojTailMass,
  WQuadRule, BraneProjectW,
  RhoSym, J4DSym, JwSym,
  J4DFromPsi, JxyzFromPsi, JwFromPsi,
  JbraneXYZFromPsi,
  jMouthOnGammaFromPsi, jPortsFromPsiAtTime,
  GeometryDOFs, ThroatAxis, BranePlane,
  OutwardSignW, JwOutwardSym,
  GammaParam, GammaMeasure, GammaNormal,
  nThetaDefault, nPhiDefault, GammaQuadParams, GammaThetaQuadRule, GammaQuadGrid,
  GammaIntegrate, PortCoeffOnGamma,
  PortsList, PortKey, PortP, PortOnXYZ,
  DriveEnvelope, Vdrive, VconfDriven,
  EffortLinearized,
  ComplexAmplitudeLockIn, EstimateZeff,
  nWDefault, nDiscardPeriodsDefault, nMeasurePeriodsDefault,
  WindowIndicesByPeriods, ComplexAmplitudeLockInWindow,
  BuildPortAmplitudeMatrices, EstimateZeffRobust,
  Paper7Config, Paper7FreezeSignature, Paper7FreezeSummary
];

(* ---------------------------
   C) Smooth step conventions
   --------------------------- *)
SmoothStep[u_] := (1 + Tanh[u])/2;

(* SmoothAbs with explicit epsilon to avoid singular d/dw at w=0.
   Uses epsW by default (symbol or value supplied by user). *)
SmoothAbs[u_, eps_: epsW] := Sqrt[u^2 + eps^2];

(* ---------------------------
   D) Gates & coordinate helpers
   --------------------------- *)
R3[x_, y_, z_] := Sqrt[x^2 + y^2 + z^2];

GateR[R_, a_, dr_] := 1 - SmoothStep[(R - a)/dr];

GateW[w_, L_, dw_, eps_: epsW] :=
  1 - SmoothStep[(SmoothAbs[w, eps] - L/2)/dw];

GateThroat[x_, y_, z_, w_, a_, L_, dr_, dw_, eps_: epsW] :=
  GateR[R3[x,y,z], a, dr] * GateW[w, L, dw, eps];

(* ---------------------------
   D3) Family 1 explicit Vconf(x,y,z,w;a,L)
   (Modulated w-trap + radial wall + axial endcaps)
   --------------------------- *)
OmegaSq[x_, y_, z_, w_, a_, L_] := Module[{G},
  G = GateThroat[x,y,z,w,a,L,dr,dw,epsW];
  OmOut^2 - (OmOut^2 - OmIn^2) * G
];

Vbrane[x_, y_, z_, w_, a_, L_] :=
  (1/2) * m * OmegaSq[x,y,z,w,a,L] * w^2;

VwallRadial[x_, y_, z_, a_] := Module[{R = R3[x,y,z]},
  V0 * (SmoothStep[(R - a)/dr])^p
];

VcapAxial[w_, L_] :=
  V0 * (SmoothStep[(SmoothAbs[w, epsW] - L/2)/dw])^p;

Vconf[x_, y_, z_, w_, a_, L_] :=
  Vbrane[x,y,z,w,a,L] + VwallRadial[x,y,z,a] + VcapAxial[w,L];

(* ---------------------------
   E) Wave confinement mu2A (optional track)
   Frozen: mu2A = mu_in^2 + (mu_out^2 - mu_in^2)*(1 - G)
   --------------------------- *)
mu2A[x_, y_, z_, w_, a_, L_] := Module[{G},
  G = GateThroat[x,y,z,w,a,L,dr,dw,epsW];
  muIn^2 + (muOut^2 - muIn^2) * (1 - G)
];

(* ---------------------------
   G) Geometry measures & energy (baseline 4D tube convention)
   Tube = B^3(a) × [0, L]
   --------------------------- *)
V4Tube[a_, L_] := (4 Pi/3) * a^3 * L;

A4Tube[a_, L_] := (4 Pi a^2) * L + 2 * (4 Pi/3) * a^3;

Egeom[a_, L_] := Pvac * V4Tube[a, L] + Sig * A4Tube[a, L];

dEgeomda[a_, L_] := D[Egeom[a, L], a];
dEgeomdL[a_, L_] := D[Egeom[a, L], L];

(* ---------------------------
   H) Brane observable map W(w) = |χ0(w)|^2
   χ0 is far-field harmonic oscillator ground state for Ω_out
   --------------------------- *)
ellW[] := Sqrt[hbar/(m * OmOut)];

Chi0[w_] := (1/(Pi * ellW[]^2))^(1/4) * Exp[-w^2/(2 * ellW[]^2)];

(* Weight (real, normalized): ∫ W(w) dw = 1 *)
Ww[w_] := (1/(Pi * ellW[]^2))^(1/2) * Exp[-w^2/(ellW[]^2)];

(* --- H2) Projection cutoff convention (FROZEN for numerics) ---
   We compute brane projections by integrating over |w| <= WprojCut.
   Default: WprojCut = wProjNSig * ellW[], with wProjNSig set in driver scripts.
*)
WprojCut[] := wProjNSig * ellW[];

(* Mass captured by truncated window (analytic, since W is Gaussian) *)
WprojNorm[] := Erf[WprojCut[]/ellW[]];  (* = Erf[wProjNSig] *)

(* Renormalized truncated weight *)
WwProj[w_] := (Ww[w]/WprojNorm[]) * Boole[Abs[w] <= WprojCut[]];

(* Tail mass (diagnostic) *)
WprojTailMass[] := 1 - WprojNorm[];

(* ===========================
   Step 5: Brane projection helper (numerical)
   =========================== *)

(* Default quadrature nodes in w for brane projection *)
nWDefault = 32;

(* Quadrature on [-WprojCut, WprojCut] *)
WQuadRule[nw_Integer : nWDefault] := Module[
  {data, a, b, eps = 10^-12, c1Neg, c2Neg, sum1, sum2, target, pick},
  a = -WprojCut[]; b = WprojCut[];
  data = Quiet @ Check[GaussianQuadratureWeights[nw, {a, b}], $Failed];
  If[data === $Failed || !MatchQ[data, {{_, _} ..}], Return[$Failed]];

  c1Neg = Min[data[[All, 1]]] < -eps;
  c2Neg = Min[data[[All, 2]]] < -eps;

  Which[
    c1Neg && !c2Neg, data,
    c2Neg && !c1Neg, Transpose[{data[[All, 2]], data[[All, 1]]}],
    True,
    target = b - a;
    sum1 = Abs[Total[data[[All, 2]]] - target];
    sum2 = Abs[Total[data[[All, 1]]] - target];
    pick = If[sum2 < sum1, Transpose[{data[[All, 2]], data[[All, 1]]}], data];
    pick
  ]
];

(* Weighted projection: ∫_{-WprojCut}^{WprojCut} (Ww/WprojNorm) f(w) dw *)
BraneProjectW[f_, nw_Integer : nWDefault] := Module[{rule, nodes, weights, Wvals, fvals},
  rule = WQuadRule[nw];
  If[rule === $Failed, Return[$Failed]];
  nodes = rule[[All, 1]];
  weights = rule[[All, 2]];
  Wvals = (Ww /@ nodes)/WprojNorm[];
  fvals = f /@ nodes;
  Total @ MapThread[#1*#2*#3 &, {weights, Wvals, fvals}]
];

(* Effort variable (linearized enthalpy perturbation) *)
(* h(ρ) = (5K/4) ρ^4  ==>  δh ≈ 5K ρ0^3 δρ *)
EffortLinearized[rhoBrane_, rho0_] := 5*K*rho0^3*(rhoBrane - rho0);

(* ---------------------------
   B) Symbolic current operators for the standard symbols psi/psib
   (Used by master_checks and symbolic manipulations)
   --------------------------- *)
RhoSym[x_, y_, z_, w_, t_] := psib[x,y,z,w,t] * psi[x,y,z,w,t];

J4DSym[x_, y_, z_, w_, t_] := (hbar/(2 I m)) * {
  psib[x,y,z,w,t]*D[psi[x,y,z,w,t], x] - psi[x,y,z,w,t]*D[psib[x,y,z,w,t], x],
  psib[x,y,z,w,t]*D[psi[x,y,z,w,t], y] - psi[x,y,z,w,t]*D[psib[x,y,z,w,t], y],
  psib[x,y,z,w,t]*D[psi[x,y,z,w,t], z] - psi[x,y,z,w,t]*D[psib[x,y,z,w,t], z],
  psib[x,y,z,w,t]*D[psi[x,y,z,w,t], w] - psi[x,y,z,w,t]*D[psib[x,y,z,w,t], w]
};

JwSym[x_, y_, z_, w_, t_] := (hbar/(2 I m)) * (
  psib[x,y,z,w,t]*D[psi[x,y,z,w,t], w] - psi[x,y,z,w,t]*D[psib[x,y,z,w,t], w]
);

(* ===========================
   Step 5: Numeric current from psi(x,y,z,w,t)
   Assumption: psiFun supports Derivative[...] (e.g., InterpolatingFunction).
   =========================== *)

J4DFromPsi[psiFun_][x_?NumericQ, y_?NumericQ, z_?NumericQ, w_?NumericQ, t_?NumericQ] := Module[
  {psi, dpx, dpy, dpz, dpw, psib, dpsibx, dpsiby, dpsibz, dpsibw},
  psi = psiFun[x,y,z,w,t];
  dpx = (Derivative[1,0,0,0,0][psiFun])[x,y,z,w,t];
  dpy = (Derivative[0,1,0,0,0][psiFun])[x,y,z,w,t];
  dpz = (Derivative[0,0,1,0,0][psiFun])[x,y,z,w,t];
  dpw = (Derivative[0,0,0,1,0][psiFun])[x,y,z,w,t];

  psib = Conjugate[psi];
  dpsibx = Conjugate[dpx];
  dpsiby = Conjugate[dpy];
  dpsibz = Conjugate[dpz];
  dpsibw = Conjugate[dpw];

  (hbar/(2 I m)) * {
    psib*dpx - psi*dpsibx,
    psib*dpy - psi*dpsiby,
    psib*dpz - psi*dpsibz,
    psib*dpw - psi*dpsibw
  }
];

JxyzFromPsi[psiFun_][x_?NumericQ, y_?NumericQ, z_?NumericQ, w_?NumericQ, t_?NumericQ] :=
  Take[J4DFromPsi[psiFun][x,y,z,w,t], 3];

JwFromPsi[psiFun_][x_?NumericQ, y_?NumericQ, z_?NumericQ, w_?NumericQ, t_?NumericQ] :=
  J4DFromPsi[psiFun][x,y,z,w,t][[4]];

(* ---------------------------
   C1) Geometry DOFs & axis conventions (FROZEN)
   --------------------------- *)
GeometryDOFs[] := {a, L};
ThroatAxis[] := "w";
BranePlane[] := "w==0";

(* ---------------------------
   I3) Leakage sign conventions (FROZEN)
   Outward from the slab |w|<Wcut:
     + cap has outward normal +ŵ
     - cap has outward normal -ŵ
   So outward-oriented current is Jw_out = sign(w) * Jw.
   --------------------------- *)
OutwardSignW[w_] := Piecewise[{{1, w > 0}, {-1, w < 0}}, 0];

JwOutwardSym[x_, y_, z_, w_, t_] := OutwardSignW[w] * JwSym[x,y,z,w,t];

(* ---------------------------
   I) Measurement region Γ: 2-sphere on brane (R3=rPort, w=0)
   --------------------------- *)
GammaParam[th_, ph_] := {
  rPort * Sin[th] * Cos[ph],
  rPort * Sin[th] * Sin[ph],
  rPort * Cos[th],
  0
};

GammaMeasure[th_, ph_] := rPort^2 * Sin[th];  (* dμ = r^2 sinθ dθ dφ *)

GammaNormal[th_, ph_] := {  (* outward radial in brane 3-space at w=0 *)
  Sin[th] * Cos[ph],
  Sin[th] * Sin[ph],
  Cos[th],
  0
};

(* ---------------------------
   I1) Γ quadrature conventions (FROZEN)
   --------------------------- *)
nThetaDefault = 32;
nPhiDefault = 64;

GammaQuadParams[] := <|"nTheta" -> nThetaDefault, "nPhi" -> nPhiDefault|>;

GammaThetaQuadRule[n_Integer] := Module[
  {data, nodes, weights, roots, a = 0, b = Pi, x},
  data = Quiet @ Check[GaussianQuadratureWeights[n, {a, b}], $Failed];
  If[data =!= $Failed && MatchQ[data, {{_, _} ..}],
    If[OrderedQ[data[[All, 2]]],
      weights = data[[All, 1]];
      nodes = data[[All, 2]],
      weights = data[[All, 2]];
      nodes = data[[All, 1]]
    ];
    Return[Transpose[{nodes, weights}]];
  ];

  roots = x /. NSolve[LegendreP[n, x] == 0, x, Reals];
  roots = SortBy[roots, N];
  weights = 2/((1 - roots^2) * (D[LegendreP[n, x], x] /. x -> roots)^2);
  nodes = (b + a)/2 + (b - a)/2 * roots;
  weights = weights * (b - a)/2;
  Transpose[{nodes, weights}]
];

GammaQuadGrid[nTheta_Integer, nPhi_Integer] := GammaQuadGrid[nTheta, nPhi] = Module[
  {thetaRule, dphi, phis},
  thetaRule = GammaThetaQuadRule[nTheta];
  dphi = (2 Pi)/nPhi;
  phis = Table[(k - 1) * dphi, {k, 1, nPhi}];
  Flatten[
    Map[
      Function[{pair},
        With[{th = pair[[1]], wth = pair[[2]]},
          Table[{th, ph, wth * dphi * rPort^2 * Sin[th]}, {ph, phis}]
        ]
      ],
      thetaRule
    ],
    1
  ]
];

Options[GammaIntegrate] = {"nTheta" -> Automatic, "nPhi" -> Automatic};

GammaIntegrate[f_, opts : OptionsPattern[]] := Module[
  {params, nTheta, nPhi, grid},
  params = GammaQuadParams[];
  nTheta = OptionValue["nTheta"]; nPhi = OptionValue["nPhi"];
  If[nTheta === Automatic, nTheta = params["nTheta"]];
  If[nPhi === Automatic, nPhi = params["nPhi"]];
  grid = GammaQuadGrid[nTheta, nPhi];
  Total[Map[#[[3]] * f[#[[1]], #[[2]]] &, grid]]
];

Options[PortCoeffOnGamma] = Options[GammaIntegrate];

PortCoeffOnGamma[fOnGamma_, lmax_Integer, opts : OptionsPattern[]] := Module[
  {ports = PortsList[lmax]},
  Table[
    GammaIntegrate[Function[{th, ph}, Conjugate[PortP[key][th, ph]] * fOnGamma[th, ph]], opts],
    {key, ports}
  ]
];

(* ===========================
   Step 5: Mouth flux observable on Γ (for Zeff)
   j(th,ph,t) = J_brane(x(th,ph),y,z,t) · n_hat(th,ph)
   =========================== *)

JbraneXYZFromPsi[psiFun_][x_?NumericQ, y_?NumericQ, z_?NumericQ, t_?NumericQ, nw_Integer : nWDefault] :=
  BraneProjectW[Function[{ww}, JxyzFromPsi[psiFun][x,y,z,ww,t]], nw];

jMouthOnGammaFromPsi[psiFun_][th_?NumericQ, ph_?NumericQ, t_?NumericQ, nw_Integer : nWDefault] := Module[
  {X, x, y, z, nHat, Jb},
  X = GammaParam[th, ph];
  x = X[[1]]; y = X[[2]]; z = X[[3]];
  nHat = Take[GammaNormal[th, ph], 3];
  Jb = JbraneXYZFromPsi[psiFun][x,y,z,t,nw];
  Jb . nHat
];

Options[jPortsFromPsiAtTime] = Join[Options[GammaIntegrate], {"lmax" -> 2, "nW" -> nWDefault}];

jPortsFromPsiAtTime[psiFun_, t_?NumericQ, opts : OptionsPattern[]] := Module[
  {lmax, nw, fOnGamma},
  lmax = OptionValue["lmax"];
  nw = OptionValue["nW"];
  fOnGamma = Function[{th, ph}, jMouthOnGammaFromPsi[psiFun][th, ph, t, nw]];
  PortCoeffOnGamma[fOnGamma, lmax, FilterRules[{opts}, Options[GammaIntegrate]]]
];

(* Ports: spherical harmonics up to lmax *)
PortKey[l_Integer, mm_Integer] := {l, mm};

PortsList[lmax_Integer] := Flatten[
  Table[PortKey[l, mm], {l, 0, lmax}, {mm, -l, l}],
  1
];

PortP[{l_Integer, mm_Integer}][th_, ph_] := SphericalHarmonicY[l, mm, th, ph];

(* XYZ -> angles on sphere (with epsAng regularization) *)
ThetaXYZ[x_, y_, z_] := Module[{R = Sqrt[x^2 + y^2 + z^2 + epsAng^2]},
  ArcCos[z/R]
];

PhiXYZ[x_, y_] := ArcTan[x, y]; (* atan2(y,x) *)

PortOnXYZ[key:{l_Integer, mm_Integer}][x_, y_, z_] :=
  PortP[key][ThetaXYZ[x,y,z], PhiXYZ[x,y]];

(* ---------------------------
   J) Drive protocol: potential modulation localized near Γ
      Vconf -> Vconf + eps cos(ω t) f(X) P_k(s(X))
   --------------------------- *)
DriveEnvelope[x_, y_, z_, w_] := Module[{R = R3[x,y,z]},
  Exp[-(w^2)/(wDrive^2)] * Exp[-(R - rPort)^2/(rDrive^2)]
];

Vdrive[x_, y_, z_, w_, t_, key:{l_Integer, mm_Integer}, eps_, omega_] :=
  eps * Cos[omega * t] * DriveEnvelope[x,y,z,w] * PortOnXYZ[key][x,y,z];

VconfDriven[x_, y_, z_, w_, a_, L_, t_, key:{l_Integer, mm_Integer}, eps_, omega_] :=
  Vconf[x,y,z,w,a,L] + Vdrive[x,y,z,w,t,key,eps,omega];

(* ---------------------------
   K) Effective operator definition Z^eff(ω)
   Minimal helpers for extraction from time series (lock-in style).
   --------------------------- *)

(* Complex amplitude estimate at frequency ω using demodulation.
   For uniform sampling, this approximates ∫ s(t) e^{-iωt} dt / T. *)
ComplexAmplitudeLockIn[tList_List, sList_List, omega_] := Module[
  {n = Length[tList], dt, T, phasor},
  If[n != Length[sList] || n < 4, Return[$Failed]];
  dt = (tList[[-1]] - tList[[1]])/(n - 1);
  T = (tList[[-1]] - tList[[1]]);
  phasor = Exp[-I * omega * tList];
  (dt/T) * Total[sList * phasor]
];

(* ===========================
   Step 5: Steady-state window protocol (FROZEN defaults)
   =========================== *)

nDiscardPeriodsDefault = 5;
nMeasurePeriodsDefault = 10;

WindowIndicesByPeriods[tList_List, omega_, nDiscard_Integer : nDiscardPeriodsDefault, nMeasure_Integer : nMeasurePeriodsDefault] := Module[
  {t0, Tper, tStart, tEnd, idx},
  If[Length[tList] < 8, Return[$Failed]];
  Tper = 2 Pi/omega;
  t0 = tList[[1]];
  tStart = t0 + nDiscard*Tper;
  tEnd   = tStart + nMeasure*Tper;

  idx = Flatten @ Position[tList, _?(# >= tStart && # <= tEnd &)];
  If[idx === {} || Length[idx] < 8, Return[$Failed]];
  idx
];

ComplexAmplitudeLockInWindow[tList_List, sList_List, omega_, nDiscard_Integer : nDiscardPeriodsDefault, nMeasure_Integer : nMeasurePeriodsDefault] := Module[
  {idx, tt, ss},
  idx = WindowIndicesByPeriods[tList, omega, nDiscard, nMeasure];
  If[idx === $Failed, Return[$Failed]];
  tt = tList[[idx]];
  ss = sList[[idx]];
  ComplexAmplitudeLockIn[tt, ss, omega]
];

(* ===========================
   Step 5: Assemble uMat/jMat (complex amplitudes) from time series
   =========================== *)

BuildPortAmplitudeMatrices[runs_List, omega_, nDiscard_Integer : nDiscardPeriodsDefault, nMeasure_Integer : nMeasurePeriodsDefault] := Module[
  {uCols, jCols},
  uCols = Table[
    Table[
      ComplexAmplitudeLockInWindow[runs[[k, "t"]], runs[[k, "uPorts"]][[i]], omega, nDiscard, nMeasure],
      {i, Length[runs[[k, "uPorts"]]]}
    ],
    {k, Length[runs]}
  ];
  jCols = Table[
    Table[
      ComplexAmplitudeLockInWindow[runs[[k, "t"]], runs[[k, "jPorts"]][[i]], omega, nDiscard, nMeasure],
      {i, Length[runs[[k, "jPorts"]]]}
    ],
    {k, Length[runs]}
  ];
  {Transpose[uCols], Transpose[jCols]}
];

(* Estimate Z^eff from matrices of complex amplitudes.
   uMat and jMat: each column corresponds to one driven port. *)
EstimateZeff[uMat_, jMat_] := Module[{Z},
  Z = LinearSolve[Transpose[uMat], Transpose[jMat]];
  Transpose[Z]
];

(* ===========================
   Step 5: Robust Zeff solve policy (FROZEN default behavior)
   =========================== *)

Options[EstimateZeffRobust] = {
  "CondThresh" -> 10^10,
  "RidgeLambda" -> 0.0,
  "SVTol" -> 10^-12,
  "Verbose" -> True
};

EstimateZeffRobust[uMat_, jMat_, opts : OptionsPattern[]] := Module[
  {condThresh, lambda, svtol, verbose, sv, cond, pinv, Z, resid},
  condThresh = OptionValue["CondThresh"];
  lambda = OptionValue["RidgeLambda"];
  svtol = OptionValue["SVTol"];
  verbose = OptionValue["Verbose"];

  sv = SingularValueList[uMat];
  cond = If[Min[sv] == 0, Infinity, Max[sv]/Min[sv]];

  If[verbose,
    Print["[Zeff] singular values(u): ", sv];
    Print["[Zeff] cond(u) ~ ", N[cond]];
  ];

  If[cond < condThresh && lambda == 0,
    Z = EstimateZeff[uMat, jMat],
    If[lambda == 0,
      pinv = PseudoInverse[uMat, Tolerance -> svtol],
      pinv = Inverse[ConjugateTranspose[uMat].uMat + lambda IdentityMatrix[Length[uMat[[1]]]]].ConjugateTranspose[uMat]
    ];
    Z = jMat . pinv;
  ];

  resid = Norm[jMat - Z.uMat]/Max[1.*^-30, Norm[jMat]];
  If[verbose, Print["[Zeff] relative residual = ", N[resid]]];
  <|"Z" -> Z, "cond" -> cond, "resid" -> resid|>
];

(* ---------------------------
   Config bundle (human-readable single object)
   --------------------------- *)
Paper7Config[] := <|
  "Version" -> "definitions.wl v1.0",
  "PrimaryVconfFamily" -> "F1_modulated_brane_trap",
  "SmoothStep" -> "S(u)=(1+Tanh[u])/2",
  "Gate" -> "G=GateR(R3;a,dr)*GateW(w;L,dw) with SmoothAbs",
  "GeometryDOFs" -> "{a,L} only (baseline)",
  "AxisConventions" -> "Throat axis = w, brane plane = w=0",
  "GeometryMeasures" -> "Tube=B^3(a)×[0,L], V=(4π/3)a^3L, A=(4πa^2)L+2(4π/3)a^3",
  "Constraint" -> "FixedNormN",
  "ProjectionWeight" -> "W(w)=|χ0(w)|^2 (Ω_out ground state)",
  "ProjectionWindow" -> "Numerics: |w| <= WprojCut = wProjNSig * ellW",
  "ProjectionRenorm" -> "Use WwProj(w)=W(w)/WprojNorm on |w|<=WprojCut; tail mass = 1-WprojNorm",
  "Gamma" -> "R3=rPort, w=0",
  "GammaMeasure" -> "dμ = r_port^2 sinθ dθ dφ",
  "GammaNormal" -> "outward radial on Γ",
  "GammaQuadrature" -> "θ: Gauss-Legendre (nTheta=" <> ToString[nThetaDefault] <> "), φ: uniform (nPhi=" <> ToString[nPhiDefault] <> ")",
  "MeasurementRegion" -> "Γ: sphere R3=rPort at w=0",
  "FluxPrimary" -> "leakage, outward-oriented: Jw_out = sign(w) Jw at w=±Wcut over R3<Rmeasure",
  "LeakageConvention" -> "j_out+(t)=∫ Jw(w=+Wcut)d^3x, j_out-(t)=∫[-Jw(w=-Wcut)]d^3x, j_leak=j_out+ + j_out-",
  "Ports" -> "Y_lm complex basis, ordering = PortsList[lmax], lmax baseline = 2",
  "Drive" -> "Vconf -> Vconf + eps cos(ωt) exp(-w^2/wDrive^2) exp(-(R3-rPort)^2/rDrive^2) Y_lm(theta(x,y,z),phi(x,y))",
  "nWDefault" -> nWDefault,
  "nDiscardPeriods" -> nDiscardPeriodsDefault,
  "nMeasurePeriods" -> nMeasurePeriodsDefault,
  "ZeffCondThresh" -> ("CondThresh" /. Options[EstimateZeffRobust]),
  "ZeffRidgeLambda" -> ("RidgeLambda" /. Options[EstimateZeffRobust]),
  "ZeffSVTol" -> ("SVTol" /. Options[EstimateZeffRobust]),
  "ResponseOperator" -> "j_mouth(ω)=Zeff(ω) u(ω); leakage j_leak reported separately",
  "ResponseExtraction" -> StringTemplate[
     "Lock-in: discard `` periods, measure `` periods; nw=``; Zeff via EstimateZeffRobust (CondThresh=`` RidgeLambda=`` SVTol=``)"
   ][nDiscardPeriodsDefault, nMeasurePeriodsDefault, nWDefault,
     ("CondThresh" /. Options[EstimateZeffRobust]),
     ("RidgeLambda" /. Options[EstimateZeffRobust]),
     ("SVTol" /. Options[EstimateZeffRobust])
   ]
|>;

Paper7FreezeSignature[] := Hash[Paper7Config[], "SHA256"];

Paper7FreezeSummary[] := Module[{cfg = Paper7Config[]},
  Print["=== Paper 7 Operational Definitions (Executable) ==="];
  Print["Version: ", cfg["Version"]];
  Print["Primary Vconf: ", cfg["PrimaryVconfFamily"]];
  Print["Geometry DOFs: ", cfg["GeometryDOFs"]];
  Print["Axis conventions: ", cfg["AxisConventions"]];
  Print["Constraint: ", cfg["Constraint"]];
  Print["Projection weight: ", cfg["ProjectionWeight"]];
  Print["Projection window: ", cfg["ProjectionWindow"]];
  Print["Projection renorm: ", cfg["ProjectionRenorm"]];
  Print["Gamma: ", cfg["Gamma"]];
  Print["Gamma measure: ", cfg["GammaMeasure"]];
  Print["Gamma normal: ", cfg["GammaNormal"]];
  Print["Gamma quadrature: ", cfg["GammaQuadrature"]];
  Print["Measurement region: ", cfg["MeasurementRegion"]];
  Print["Flux primary: ", cfg["FluxPrimary"]];
  Print["Ports: ", cfg["Ports"]];
  Print["Drive: ", cfg["Drive"]];
  Print["Response operator: ", cfg["ResponseOperator"]];
  Print["Response extraction: ", cfg["ResponseExtraction"]];
  Print["Zeff cond thresh: ", cfg["ZeffCondThresh"]];
  Print["Zeff ridge lambda: ", cfg["ZeffRidgeLambda"]];
  Print["Zeff sv tol: ", cfg["ZeffSVTol"]];
  Print["Freeze Signature (SHA256): ", Paper7FreezeSignature[]];
];

(* ============================================================
   Paper 7: Numeric config exporter (always works)
   -cfg exports *frozen numeric scan params* from Paper7NumericParams[].
   If the corresponding symbol (e.g. dr) is already numeric at export time,
   that runtime numeric value overrides the frozen value.

   Usage:
     math -script mathematica/inner_throat/master_checks.wl -- -cfg paper7_defs.json
     math -script mathematica/inner_throat/master_checks.wl -- -cfg-only -cfg paper7_defs.json
   ============================================================ *)

ClearAll[
  Paper7NumericParams,
  Paper7CfgParseFile, Paper7CfgHasFlag,
  Paper7CfgPairs, Paper7CfgNumValue,
  Paper7CfgBuildNumericAssociation, Paper7CfgExportNumericJSON,
  Paper7CfgMaybeExportNumericJSON
];

(* ---- 1) ONE PLACE TO PUT FROZEN NUMBERS (fill these in once) ---- *)
Paper7NumericParams[] := <|
  "hbar" -> 1.0,
  "m"    -> 1.0,
  "K"    -> 1.0,

  "dr"   -> 0.25,
  "dw"   -> 0.25,
  "epsW" -> 1.*^-3,
  "V0"   -> 50.0,
  "p"    -> 4.0,

  "OmOut" -> 2.0,
  "OmIn"  -> 0.5,

  "Pvac" -> 0.1,
  "Sig"  -> 0.1,

  "wProjNSig" -> 3.0,
  "Wcut"      -> 3.0,
  "Rmeasure"  -> 2.5,
  "rPort"     -> 2.0,
  "epsAng"    -> 1.*^-6,

  "wDrive" -> 0.5,
  "rDrive" -> 0.5
|>;

(* ---- 2) CLI helpers ---- *)
Paper7CfgHasFlag[flag_String] := Module[{cmd = If[ListQ[$ScriptCommandLine], $ScriptCommandLine, {}]},
  MemberQ[cmd, flag] || AnyTrue[cmd, StringStartsQ[#, flag <> "="] &]
];

Paper7CfgParseFile[] := Module[{cmd = If[ListQ[$ScriptCommandLine], $ScriptCommandLine, {}], tok, pos, out},
  tok = SelectFirst[cmd, StringStartsQ[#, "-cfg="] &, Missing["NotFound"]];
  If[tok =!= Missing["NotFound"],
    out = StringDrop[tok, StringLength["-cfg="]];
    Return[If[StringLength[out] > 0, out, "paper7_defs.json"]];
  ];
  pos = FirstPosition[cmd, "-cfg", Missing["NotFound"]];
  If[pos === Missing["NotFound"], Return[Missing["NotFound"]]];
  out = If[
    pos[[1]] < Length[cmd] &&
      StringQ[cmd[[pos[[1]] + 1]]] &&
      !StringStartsQ[cmd[[pos[[1]] + 1]], "-"],
    cmd[[pos[[1]] + 1]],
    "paper7_defs.json"
  ];
  out
];

(* ---- 3) Key<->Symbol mapping (runtime override if numeric) ---- *)
Paper7CfgPairs[] := {
  "hbar" -> hbar, "m" -> m, "K" -> K,
  "dr" -> dr, "dw" -> dw, "epsW" -> epsW, "V0" -> V0, "p" -> p,
  "OmOut" -> OmOut, "OmIn" -> OmIn,
  "Pvac" -> Pvac, "Sig" -> Sig,
  "wProjNSig" -> wProjNSig, "Wcut" -> Wcut, "Rmeasure" -> Rmeasure,
  "rPort" -> rPort, "epsAng" -> epsAng,
  "wDrive" -> wDrive, "rDrive" -> rDrive
};

Paper7CfgNumValue[key_String, sym_Symbol, frozen_Association] := Module[{v},
  v = Quiet @ Check[Evaluate[sym], $Failed];
  Which[
    NumericQ[v], N[v],
    KeyExistsQ[frozen, key] && NumericQ[frozen[key]], N[frozen[key]],
    True, Null
  ]
];

Paper7CfgBuildNumericAssociation[] := Module[{frozen, pairs, params, missing},
  frozen = Paper7NumericParams[];
  pairs = Paper7CfgPairs[];

  params = Association @ Map[
    Function[{kv}, First[kv] -> Paper7CfgNumValue[First[kv], Last[kv], frozen]],
    pairs
  ];

  missing = Keys @ Select[params, # === Null &];

  <|
    "freezeSHA256" -> ToString[Paper7FreezeSignature[]],
    "operational"  -> Paper7Config[],
    "params"       -> params,
    "missing_numeric" -> missing,
    "ready"        -> (missing === {}),
    "cwd"          -> Directory[]
  |>
];

Paper7CfgExportNumericJSON[file_String] := Module[{payload},
  payload = Paper7CfgBuildNumericAssociation[];
  Export[file, payload, "JSON"]
];

Paper7CfgMaybeExportNumericJSON[] := Module[{file, only},
  file = Paper7CfgParseFile[];
  If[file === Missing["NotFound"], Return[]];

  only = Paper7CfgHasFlag["-cfg-only"];

  Print["[Paper7] -cfg: exporting numeric scan config to: ", file];
  Paper7CfgExportNumericJSON[file];

  If[only,
    Print["[Paper7] -cfg-only: exiting after export."];
    Exit[0];
  ];
];

(* Trigger on load *)
Paper7CfgMaybeExportNumericJSON[];


(* End of file *)
