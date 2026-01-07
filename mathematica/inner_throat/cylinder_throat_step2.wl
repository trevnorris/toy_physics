(* ============================================================
   cylinder_throat_step2_v7.wl  — Paper 7 Step 2
   Patch over v6:
     - Fix resonance exclusion for λ>0 by using ω-spacing derived from κ:
         Δω_mode ≈ (π/L) (c^2 |κ|/|ω|)
       (reduces over-exclusion near cutoff)
     - Cap exclusion width to a fraction of band width to prevent NaN wipeout.
   ============================================================ *)

ClearAll["Global`*"];

(* ----------------------------
   Robust CLI detection
   ---------------------------- *)
scriptArgs = Quiet@Check[$ScriptCommandLine, {}];
allArgs = DeleteDuplicates@Flatten@{scriptArgs, $CommandLine};

cliMode = ($FrontEnd === Null) || MemberQ[allArgs, "-cli"];
doSpectrum = MemberQ[allArgs, "-spectrum"];
doCriteria = MemberQ[allArgs, "-criteria"];
doLocality = MemberQ[allArgs, "-locality"];
runScan    = MemberQ[allArgs, "-scan"];

If[cliMode, Print["[CLI mode enabled]"];];

(* ----------------------------
   User parameters
   ---------------------------- *)
c      = 1.0;
a      = 1.0;

Lovera0 = 1.85;
L0      = Lovera0*a;

wallBC   = "Neumann";
bottomBC = "Neumann";

(* Damping policy (scale-invariant recommended) *)
gammaMode  = "RelativeOmega";  (* "RelativeOmega" | "Absolute" *)
deltaRel   = 0.02;
gammaAbs   = 0.01;
gammaFloor = 1.*10^-6;

outerModel = "AutoScaled";

modeList = {
  {0, If[wallBC === "Dirichlet", 1, 0]},
  {0, If[wallBC === "Dirichlet", 2, 1]}
};

scanLovera = Range[1.2, 2.6, 0.02];

etaBandOmega0 = {0.35, 1.35};
etaBandLambda = {1.05, 1.35};

nΩ = 700;

excludeMultGamma   = 6.0;
excludeFracSpacing = 0.03;

(* NEW: cap exclusion to avoid erasing the whole ω band *)
excludeCapFracBand = 0.35;  (* max Δω = 0.35*(ωhi-ωlo) *)

denTol = 5.*10^-4;
ZCap   = 10^6;

etaLocalMax = 0.25;
seriesOrder = 8;
tolLocal    = 0.10;

epsMismatch = 10^-12;
epsEnergy   = 10^-12;
denomFloorMismatch = 10^-8;

(* ----------------------------
   Formatting helpers
   ---------------------------- *)
NumStr[x_, prec_:12] := Module[{nf},
  nf = NumberForm[
    N[x, prec],
    {Infinity, prec},
    NumberPadding -> {"", ""},
    ExponentFunction -> (Null &)
  ];
  ToString[nf, OutputForm]
];

PrintTable[header_List, rows_List] := Module[{},
  Print[StringRiffle[header, "\t"]];
  Scan[(Print[StringRiffle[(ToString[#, OutputForm] & /@ #), "\t"]]) &, rows];
];

FiniteQ[z_] := (z =!= Indeterminate) && z =!= $Failed &&
               NumericQ[z] && z =!= ComplexInfinity && z =!= Infinity;

MedianSafe[list_] := If[Length[list] > 0, Median[list], Indeterminate];

(* ----------------------------
   Damped complex frequency w(ω)
   ---------------------------- *)
Wfreq[ω_?NumericQ] := Which[
  gammaMode === "Absolute",
    ω + I*gammaAbs,
  gammaMode === "RelativeOmega",
    ω*(1 + I*deltaRel) + I*gammaFloor,
  True,
    ω + I*gammaAbs
];

GammaEff[ω_?NumericQ] := Im[Wfreq[ω]];

(* ----------------------------
   Bessel derivative
   ---------------------------- *)
BesselJPrime[m_Integer, x_?NumericQ] := (BesselJ[m - 1, x] - BesselJ[m + 1, x])/2;

(* ----------------------------
   Robust radial roots
   ---------------------------- *)
RadialRoot::bc  = "Unknown boundary condition `1`.";
RadialRoot::dir = "Dirichlet roots use n>=1; got n=`1`.";

RadialRoot[bc_, m_Integer, n_Integer] := Module[{left, right, δ, sol},
  δ = 10^-10;
  Which[
    bc === "Dirichlet",
      If[n < 1, Message[RadialRoot::dir, n]; Return[$Failed]];
      BesselJZero[m, n],

    bc === "Neumann",
      Which[
        m == 0,
          If[n == 0, 0, BesselJZero[1, n]],

        m > 0,
          left  = If[n == 0, 0, BesselJZero[m, n]];
          right = BesselJZero[m, n + 1];
          sol = FindRoot[
            BesselJPrime[m, x] == 0,
            {x, left + If[left == 0, 10^-4, δ], right - δ},
            WorkingPrecision -> 60, AccuracyGoal -> 20, PrecisionGoal -> 20,
            MaxIterations -> 200
          ];
          x /. sol,

        True, $Failed
      ],

    True,
      Message[RadialRoot::bc, bc];
      $Failed
  ]
];

Lambda[m_Integer, n_Integer] := RadialRoot[wallBC, m, n]/a;

(* ----------------------------
   κ branch control
   ---------------------------- *)
KappaBranch[κ_?NumericQ] := Module[{k = κ},
  If[Re[k] < 0, k = -k];
  If[Re[k] == 0 && Im[k] < 0, k = -k];
  k
];

Kappa[ω_?NumericQ, m_Integer, n_Integer] := Module[{w, k, λ, κ},
  w = Wfreq[ω];
  k = w/c;
  λ = Lambda[m, n];
  κ = Sqrt[k^2 - λ^2];
  KappaBranch[κ]
];

PoleDenom[κ_?NumericQ, L_?NumericQ] := Which[
  bottomBC === "Dirichlet", Sin[κ*L],
  bottomBC === "Neumann",  Cos[κ*L],
  True, 1
];

Zmn[ω_?NumericQ, m_Integer, n_Integer, L_?NumericQ] := Module[{κ, den, Z},
  κ = Kappa[ω, m, n];
  den = PoleDenom[κ, L];
  If[Abs[den] < denTol, Return[Indeterminate]];
  Z = Which[
    bottomBC === "Dirichlet", κ*Cot[κ*L],
    bottomBC === "Neumann",  -κ*Tan[κ*L],
    True, $Failed
  ];
  If[Z === $Failed, Return[$Failed]];
  If[Abs[Z] > ZCap, Sign[Re[Z]]*ZCap + I*Sign[Im[Z]]*ZCap, Z]
];

qFund := If[bottomBC === "Dirichlet", 1, 0];

ResonanceOmega[m_Integer, n_Integer, q_Integer, L_?NumericQ] := Module[{λ, κres},
  λ = Lambda[m, n];
  κres = Which[
    bottomBC === "Dirichlet", (q*Pi)/L,
    bottomBC === "Neumann",  ((q + 1/2)*Pi)/L,
    True, $Failed
  ];
  c*Sqrt[λ^2 + κres^2]
];

(* ----------------------------
   Mode-dependent band selection
   ---------------------------- *)
OmegaBand[m_Integer, n_Integer, L_?NumericQ] := Module[{λ, ω0, ωlo, ωhi},
  λ = Lambda[m, n];
  If[Abs[λ] < 10^-12,
    ω0  = ResonanceOmega[m, n, qFund, L];
    ωlo = etaBandOmega0[[1]]*ω0;
    ωhi = etaBandOmega0[[2]]*ω0;,
    ωlo = etaBandLambda[[1]]*λ;
    ωhi = etaBandLambda[[2]]*λ;
  ];
  If[ωlo >= ωhi, {Indeterminate, Indeterminate}, {ωlo, ωhi}]
];

ResonancesInBand[m_Integer, n_Integer, L_?NumericQ, ωlo_?NumericQ, ωhi_?NumericQ] := Module[{qs, res},
  qs = Range[0, 200];
  res = ResonanceOmega[m, n, #, L] & /@ qs;
  Select[res, ωlo <= # <= ωhi &]
];

(* NEW: mode-correct ω spacing from κ spacing π/L:
   Δω ≈ (π/L) * (dω/dκ) = (π/L) * (c^2 |κ|/|ω|) *)
ResSpacingOmega[m_Integer, n_Integer, L_?NumericQ, ωref_?NumericQ] := Module[{κref, ωabs},
  ωabs = Max[Abs[ωref], 10^-12];
  κref = Kappa[ωref, m, n];
  κref = Max[Abs[Re[κref]], 10^-12];
  (Pi/L) * (c^2*κref/ωabs)
];

MakeOmegaGrid[m_Integer, n_Integer, L_?NumericQ] := Module[
  {band, ωlo, ωhi, ωgrid, res, ωref, γref, Δγ, Δs, Δ, cap, keepQ},
  band = OmegaBand[m, n, L];
  {ωlo, ωhi} = band;
  If[!FiniteQ[ωlo] || !FiniteQ[ωhi], Return[{}]];

  ωgrid = Subdivide[ωlo, ωhi, nΩ];
  res = ResonancesInBand[m, n, L, ωlo, ωhi];

  ωref = (ωlo + ωhi)/2;
  γref = GammaEff[ωref];
  Δγ   = excludeMultGamma*γref;
  Δs   = excludeFracSpacing*ResSpacingOmega[m, n, L, ωref];

  Δ = Max[Δγ, Δs];

  (* cap to avoid erasing whole band *)
  cap = excludeCapFracBand*(ωhi - ωlo);
  Δ = Min[Δ, cap];

  keepQ[ω_] := If[Length[res] == 0, True, Min[Abs[ω - #] & /@ res] > Δ];
  Select[ωgrid, keepQ]
];

(* ----------------------------
   Outer DtN model
   ---------------------------- *)
ZOuterMode[ω_?NumericQ, m_Integer, n_Integer] := Module[
  {λ = Lambda[m, n], w, k},
  w = Wfreq[ω];
  k = w/c;
  Which[
    outerModel === "SphereMonopole", I*k - 1/a,
    outerModel === "SphereScaled",  I*k,
    outerModel === "CylContinuation", I*Kappa[ω, m, n],
    outerModel === "Auto",
      If[Abs[λ] < 10^-12, I*k - 1/a, I*Kappa[ω, m, n]],
    outerModel === "AutoScaled",
      If[Abs[λ] < 10^-12, I*k, I*Kappa[ω, m, n]],
    True,
      If[Abs[λ] < 10^-12, I*k, I*Kappa[ω, m, n]]
  ]
];

(* ----------------------------
   Criteria
   ---------------------------- *)
CriterionDtNMismatch[m_Integer, n_Integer, LoveraVal_?NumericQ] := Module[
  {L, ωgrid, Zin, Zout, vals, keep, denom2},
  L = LoveraVal*a;
  ωgrid = MakeOmegaGrid[m, n, L];
  If[Length[ωgrid] < 30, Return[Indeterminate]];

  Zin = Zmn[#, m, n, L] & /@ ωgrid;
  keep = Map[FiniteQ, Zin];
  If[Count[keep, True] < 30, Return[Indeterminate]];

  Zout = ZOuterMode[#, m, n] & /@ ωgrid;

  vals = Table[
    denom2 = Max[Abs[Zout[[i]]]^2, denomFloorMismatch] + epsMismatch;
    Abs[Zin[[i]] - Zout[[i]]]^2/denom2
    ,
    {i, Length[ωgrid]}
  ];

  Mean[Pick[vals, keep]]
];

EStoredFast[ω_?NumericQ, m_Integer, n_Integer, L_?NumericQ] := Module[
  {κ, den, β, factor},
  κ = Kappa[ω, m, n];
  den = PoleDenom[κ, L];
  If[Abs[den] < denTol, Return[Indeterminate]];
  β = Im[κ];
  factor = If[Abs[β] < 10^-10, L, Sinh[2*β*L]/(2*β)];
  (Abs[κ]^2/Abs[den]^2) * factor
];

CriterionEnergyReduced[m_Integer, n_Integer, LoveraVal_?NumericQ] := Module[
  {L, ωgrid, vals},
  L = LoveraVal*a;
  ωgrid = MakeOmegaGrid[m, n, L];
  If[Length[ωgrid] < 30, Return[Indeterminate]];

  vals = Table[
    Module[{ω = ωgrid[[i]], κ, E},
      κ = Kappa[ω, m, n];
      E = EStoredFast[ω, m, n, L];
      If[FiniteQ[E], E/(Abs[κ]^2*L + epsEnergy), Indeterminate]
    ],
    {i, Length[ωgrid]}
  ];

  vals = Select[vals, FiniteQ];
  If[Length[vals] < 30, Indeterminate, Mean[vals]]
];

(* ----------------------------
   Mode table
   ---------------------------- *)
ModeTableRows[] := Module[{rows, mMax, nMax, nStart},
  mMax = Max[modeList[[All, 1]]];
  nMax = Max[modeList[[All, 2]]];
  nStart = If[wallBC === "Dirichlet", 1, 0];
  rows = Flatten[
    Table[
      Table[
        Module[{lam = Lambda[m, n]},
          {ToString[m], ToString[n], NumStr[lam, 12], NumStr[c*lam, 12]}
        ],
        {n, nStart, nMax}
      ],
      {m, 0, mMax}
    ],
    1
  ];
  rows
];

If[cliMode,
  Print["=== Cylinder DtN Step 2 v7 Summary ==="];
  Print["wallBC=", wallBC, "  bottomBC=", bottomBC,
        "  a=", NumStr[a, 6], "  L0=", NumStr[L0, 6], "  (L0/a)=", NumStr[L0/a, 6]];
  Print["gammaMode=", gammaMode,
        "  (deltaRel=", NumStr[deltaRel, 6], ", gammaAbs=", NumStr[gammaAbs, 6], ", gammaFloor=", NumStr[gammaFloor, 6], ")"];
  Print["nΩ=", nΩ, "  etaBandOmega0=", etaBandOmega0, "  etaBandLambda=", etaBandLambda];
  Print["exclusion: multGamma=", excludeMultGamma, "  fracSpacing=", excludeFracSpacing,
        "  capFracBand=", excludeCapFracBand];
  Print["outerModel=", outerModel];
  Print["modes: ", modeList];
  Print["--- Mode table (subset) ---"];
  PrintTable[{"m","n","lambda","omega_cutoff"}, ModeTableRows[]];
];

(* ----------------------------
   Criteria scan
   ---------------------------- *)
RunCriteriaScan[] := Module[
  {dataMis, dataERed, goodMis, goodERed, bestMis, bestERed,
   m, n, ωrefList, γOverωList, medγOverω},

  Print["--- Step 2B v7: Criteria scans (resonance-excluded) ---"];
  Print["Bands: λ=0 uses ω0(L); λ>0 uses λ (cutoff-safe)."];
  Print["Damping: ", gammaMode];
  Print["Criteria: mismatch(min), energyRed(min)"];

  Do[
    m = mn[[1]]; n = mn[[2]];
    Print[""];
    Print["### Mode (m,n)=(", m, ",", n, ") ###"];

    ωrefList = Table[
      Module[{L = Lov*a, band, ωlo, ωhi},
        band = OmegaBand[m, n, L];
        {ωlo, ωhi} = band;
        If[FiniteQ[ωlo] && FiniteQ[ωhi], (ωlo + ωhi)/2, Indeterminate]
      ],
      {Lov, scanLovera}
    ];
    γOverωList = Select[
      Table[
        If[FiniteQ[ωrefList[[i]]] && ωrefList[[i]] > 0,
          GammaEff[ωrefList[[i]]]/ωrefList[[i]],
          Indeterminate
        ],
        {i, Length[ωrefList]}
      ],
      FiniteQ
    ];
    medγOverω = MedianSafe[γOverωList];

    dataMis  = Table[{Lov, CriterionDtNMismatch[m, n, Lov]},   {Lov, scanLovera}];
    dataERed = Table[{Lov, CriterionEnergyReduced[m, n, Lov]}, {Lov, scanLovera}];

    goodMis  = Select[dataMis,  FiniteQ[#[[2]]] &];
    goodERed = Select[dataERed, FiniteQ[#[[2]]] &];

    bestMis  = If[Length[goodMis]  > 0, First@MinimalBy[goodMis,  Last], {Indeterminate, Indeterminate}];
    bestERed = If[Length[goodERed] > 0, First@MinimalBy[goodERed, Last], {Indeterminate, Indeterminate}];

    Print["median(γ_eff/ω) over scan (ref mid-band) = ", NumStr[medγOverω, 10]];
    Print["criterion\tbest_L/a\tvalue"];
    Print["mismatch(min)\t", NumStr[bestMis[[1]], 6],  "\t", NumStr[bestMis[[2]], 12]];
    Print["energyRed(min)\t", NumStr[bestERed[[1]], 6], "\t", NumStr[bestERed[[2]], 12]];

    If[runScan,
      Print["--- Scan table (TSV): L/a  mismatch  energyRed ---"];
      Print["L/a\tmismatch\tenergyRed"];
      Do[
        Module[{Lov = scanLovera[[i]], vm, ver},
          vm  = dataMis[[i, 2]];
          ver = dataERed[[i, 2]];
          Print[
            NumStr[Lov, 6], "\t",
            If[FiniteQ[vm],  NumStr[vm,  10], "NaN"], "\t",
            If[FiniteQ[ver], NumStr[ver, 10], "NaN"]
          ];
        ],
        {i, Length[scanLovera]}
      ];
    ];
    ,
    {mn, modeList}
  ];
];

If[cliMode && doCriteria, RunCriteriaScan[]];

(* ----------------------------
   Locality check (unchanged, γ=0 series)
   ---------------------------- *)
LocalityCheck[] := Module[
  {m, n, λ, L, ω0, ωmax, ωgrid, ωs, ZexactExpr, ZserPoly,
   Zexact, Zseries, relerr, relerrList, maxErr, goodPts, ωgoodMax},

  {m, n} = {0, If[wallBC === "Dirichlet", 1, 0]};
  λ = Lambda[m, n];
  If[Abs[λ] > 10^-12,
    Print["--- Locality check skipped: selected mode has λ!=0 (non-analytic at ω=0)."];
    Return[$Failed];
  ];

  L = L0;
  ω0 = ResonanceOmega[m, n, If[bottomBC === "Dirichlet", 1, 0], L];
  ωmax = etaLocalMax*ω0;

  ωgrid = Subdivide[ωmax/250., ωmax, 250];
  ωs = Unique["ω"];

  ZexactExpr[ωsym_] := Module[{k},
    k = (ωsym/c);
    Which[
      bottomBC === "Neumann",  -k*Tan[k*L],
      bottomBC === "Dirichlet", k*Cot[k*L],
      True, Indeterminate
    ]
  ];

  ZserPoly = Normal@Series[ZexactExpr[ωs], {ωs, 0, seriesOrder}];

  Zexact[ω_?NumericQ]  := N[ZexactExpr[ω], 50];
  Zseries[ω_?NumericQ] := N[ZserPoly /. ωs -> ω, 50];

  relerr[ω_?NumericQ] := Module[{ze = Zexact[ω], zs = Zseries[ω], denom},
    denom = Max[Abs[ze], 10^-12];
    Abs[zs - ze]/denom
  ];

  relerrList = Table[{ω, relerr[ω]}, {ω, ωgrid}];
  maxErr = Max[relerrList[[All, 2]]];

  goodPts = Select[relerrList, #[[2]] <= tolLocal &];
  ωgoodMax = If[Length[goodPts] > 0, Max[goodPts[[All, 1]]], 0.0];

  Print["--- Locality check (series about ω=0, λ=0 mode) ---"];
  Print["mode (m,n)=(", m, ",", n, ")  L/a=", NumStr[L/a, 6]];
  Print["ω0=", NumStr[ω0, 12], "  testing ω∈[", NumStr[First[ωgrid], 12], ", ", NumStr[ωmax, 12], "] (etaLocalMax=", etaLocalMax, ")"];
  Print["seriesOrder=", seriesOrder, "  tolLocal=", tolLocal];
  Print["max relerr over test range = ", NumStr[maxErr, 10]];
  Print["largest ω with relerr <= tolLocal = ", NumStr[ωgoodMax, 12],
        "  (fraction of ω0 = ", NumStr[If[ω0 > 0, ωgoodMax/ω0, 0], 8], ")"];
  Print["Z_series(ω) = ", ToString[(ZserPoly /. ωs -> ω), InputForm]];
];

If[cliMode && doLocality, LocalityCheck[]];

If[cliMode && !(doSpectrum || doCriteria || doLocality),
  Print["(no -spectrum / -criteria / -locality flag provided; nothing else to run)"];
];
