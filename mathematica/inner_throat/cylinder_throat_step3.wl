(* ============================================================
   cylinder_throat_step3_2.wl  — Paper 7 Step 3.1b
   Resonance-resolved near–far matching scan using ka-anchored overlap bands.

   Purpose:
     Implement the Paper-7 themed anchor robustly for lambda>0 modes:
       Choose L/a that best matches the inner throat DtN eigenvalue
       Z_in(ω;L) to an outer reference Z_out(ω) over an overlap band
       defined in ka = ω a / c.

   Step 3.1 upgrade:
     * For lambda>0 channels, replace band-averaged mismatch with a
       resonance-resolved metric:
         - predict resonance frequencies ω_q(L)
         - sample at ω_q(1±ε) (controlled offsets from poles)
         - aggregate mismatch with a trimmed mean
       This dramatically reduces sensitivity to small changes in the
       overlap-band endpoints.

   CLI:
     math -script mathematica/inner_throat/cylinder_throat_step3_1.wl -- -cli -scan

   Optional flags:
     -spectrum     Print cutoff + a few resonance estimates per mode.
     -plot         (Notebook only) plot mismatch + coverage.

   Output (CLI):
     For each mode:
       - TSV table: L/a, mismatch, coverage, Nω, Nres, medAbsZin, medAbsZout
       - Best L/a and robustness under ±bandJitterFrac scaling.
   ============================================================ *)

ClearAll["Global`*"];

(* ----------------------------
   Robust CLI detection
   ---------------------------- *)
scriptArgs = Quiet@Check[$ScriptCommandLine, {}];
allArgs = DeleteDuplicates@Flatten@{scriptArgs, $CommandLine};

cliMode    = ($FrontEnd === Null) || MemberQ[allArgs, "-cli"];
runScan    = MemberQ[allArgs, "-scan"];
doSpectrum = MemberQ[allArgs, "-spectrum"];
doPlotFlag = MemberQ[allArgs, "-plot"];

plotMode = (!cliMode) || doPlotFlag;

If[cliMode, Print["[CLI mode enabled]"];];

(* ----------------------------
   User parameters
   ---------------------------- *)
c = 1.0;
a = 1.0;

(* Reference L/a (used when not scanning) *)
Lovera0 = 1.85;
L0 = Lovera0*a;

wallBC   = "Neumann";   (* "Neumann" or "Dirichlet" at r=a *)
bottomBC = "Neumann";   (* "Neumann" or "Dirichlet" at z=-L *)

(* Damping policy: keep scans scale-invariant *)
gammaMode  = "RelativeOmega";  (* "RelativeOmega" | "Absolute" *)
deltaRel   = 0.02;             (* ω -> ω(1+i δ) + i γ_floor *)
gammaAbs   = 0.01;             (* used only for gammaMode="Absolute" *)
gammaFloor = 1.*10^-6;

(* Inner modes to scan *)
modeList = {
  {0, If[wallBC === "Dirichlet", 1, 0]},
  {0, If[wallBC === "Dirichlet", 2, 1]}
};

(* L/a scan range *)
scanLovera = Range[1.2, 2.6, 0.02];

(* Overlap band definition in ka = ω a / c *)
kaBandMonopole = {0.10, 0.60};   (* λ=0 channel *)
etaCutBand     = {1.01, 1.17};   (* λ>0: ka in [η1 x_cut, η2 x_cut] where x_cut = λ a *)

bandJitterFrac = 0.05;           (* robustness: scale band endpoints by 1±bandJitterFrac *)

(* --- Metric mode ---
   "Grid"       : band-averaged mismatch on ω-grid (Step 3 behavior)
   "Resonance"  : resonance-resolved mismatch (Step 3.1 behavior)
   "Hybrid"     : λ=0 uses Grid; λ>0 uses Resonance (recommended)
*)
metricMode = "Hybrid";

(* ----------------------------
   Step 3.2 additions (resonance-labeled fixed-q sampling)
   ---------------------------- *)
resSamplingMode = "FixedQ";         (* "FixedQ" recommended; "AutoBand" mimics 3.1b *)
qFixedList = {0,1};              (* axial indices used when resSamplingMode=FixedQ *)
applyKaFilterHighMode = False;       (* if True, filter ω samples to etaCutBand*ωcut *)
maxSamplePoints = 60;                (* cap on ω samples used in resonance metric *)
qRobustSets = {{0}, {0, 1}, {0, 1, 2}, {0, 1, 2, 3}};
detuneScaleRobust = {0.8, 1.0, 1.2};

(* ω grid resolution (Grid metric only) *)
nΩ = 900;

(* Resonance-resolved sampling controls (Resonance metric)

   IMPORTANT: the ka overlap band for λ>0 modes can be narrow (e.g. 1.04–1.17 × cutoff),
   which may include only ~1 pole for typical L/a. If we only sample at ω_res(1±ε),
   we can end up with too few in-band samples and the scan returns NaN.

   Fix (Step 3.1b): sample a *fan* of detunings around each resonance:
     ω = ω_res * (1 ± δ_j), for several δ_j values.
   This yields O(10–50) stable sample points even if only one resonance lies near the band.
*)

(* Detuning fan (relative): ω = ω_res (1 ± δ).  Include enough δ values that even a narrow band
   can capture O(20–40) in-band samples from the nearest few resonances. *)
detuneList = {0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20};
resScaleFactors = Join[1 - detuneList, 1 + detuneList];
maxDetune = Max[detuneList];

minResonances = 1;               (* require at least this many resonances contributing samples *)
maxResonances = 24;              (* cap; if more, subsample evenly *)
minSamplePoints = 12;            (* min ω sample points after safety filtering *)

(* Pole guard: require denominator not too small *)
poleGuard = 3.0*deltaRel;        (* scales with loss fraction *)
denTol    = 5.*10^-4;            (* hard veto used inside Z_in *)

(* Data hygiene / robustness *)
minOmegaPoints = 80;             (* minimum ω points after pole guard (Grid metric) *)
coverageMin    = 0.30;           (* minimum keep fraction after finiteness checks *)
trimFrac       = 0.10;           (* trimmed mean fraction on each side *)

ZCap = 10^8;                     (* cap to prevent numerical blowups *)
denomFloorMismatch = 10^-8;      (* mismatch denom floor *)

(* Outer DtN model: reference continuation in the overlap region *)
outerModel = "AutoScaled";
(*  options:
      "SphereMonopole" : Z_out = i k - 1/a  (λ=0)
      "SphereScaled"   : Z_out = i k        (λ=0)
      "CylContinuation": Z_out = i κ(ω)     (all modes)
      "Auto"           : λ=0 -> i k - 1/a, λ>0 -> i κ
      "AutoScaled"     : λ=0 -> i k,       λ>0 -> i κ
*)

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

FiniteQ[z_] := NumericQ[z] && FreeQ[z, Indeterminate | ComplexInfinity | Infinity | DirectedInfinity | $Failed];

MedianSafe[list_] := If[Length[list] > 0, Median[list], Indeterminate];

TrimmedMeanLocal[list_, frac_:0.10] := Module[{s, k},
  If[Length[list] < 5, Return[Mean[list]]];
  s = Sort[list];
  k = Max[1, Floor[frac*Length[s]]];
  If[2*k >= Length[s], Mean[s], Mean[Take[s, {k + 1, -k}]]]
];

SubsampleEven[list_List, max_Integer] := Module[{n = Length[list], idx},
  If[n <= max, Return[list]];
  idx = Round[Subdivide[1, n, max]];
  list[[DeleteDuplicates[idx]]]
];

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

(* ----------------------------
   Bessel derivative (stable recurrence)
   ---------------------------- *)
BesselJPrime[m_Integer, x_?NumericQ] := (BesselJ[m - 1, x] - BesselJ[m + 1, x])/2;

(* ----------------------------
   Radial roots (memoized)
   ---------------------------- *)
RadialRoot::bc  = "Unknown boundary condition `1`.";
RadialRoot::dir = "Dirichlet roots use n>=1; got n=`1`.";

RadialRoot[bc_, m_Integer, n_Integer] := RadialRoot[bc, m, n] = Module[
  {left, right, δ, sol, guess},
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
          If[n == 0,
            left = 0;
            right = BesselJZero[m, 1];
            guess = 0.5*right;,
            left = BesselJZero[m, n];
            right = BesselJZero[m, n + 1];
            guess = (left + right)/2;
          ];

          sol = FindRoot[
            BesselJPrime[m, x] == 0,
            {x, guess, left + If[left == 0, 10^-6, δ], right - δ},
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

Lambda[m_Integer, n_Integer] := Lambda[m, n] = RadialRoot[wallBC, m, n]/a;

(* ----------------------------
   κ(ω) branch control
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

(* ----------------------------
   DtN eigenvalue Z_in(ω;L) with pole guarding
   ---------------------------- *)
PoleDenom[κ_?NumericQ, L_?NumericQ] := Which[
  bottomBC === "Dirichlet", Sin[κ*L],
  bottomBC === "Neumann",  Cos[κ*L],
  True, 1
];

Zin[ω_?NumericQ, m_Integer, n_Integer, L_?NumericQ] := Module[{κ, den, Z},
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

PoleSafeQ[ω_?NumericQ, m_Integer, n_Integer, L_?NumericQ] := Module[{κ, den},
  κ = Kappa[ω, m, n];
  den = PoleDenom[κ, L];
  Abs[den] > poleGuard
];

(* ----------------------------
   Outer DtN reference Z_out(ω)
   ---------------------------- *)
Zout[ω_?NumericQ, m_Integer, n_Integer] := Module[
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
   ka-anchored overlap band
   ---------------------------- *)
KaBand[m_Integer, n_Integer, bandScale_:1.0] := Module[{λ, xcut, xlo, xhi},
  λ = Lambda[m, n];
  If[Abs[λ] < 10^-12,
    xlo = bandScale*kaBandMonopole[[1]];
    xhi = bandScale*kaBandMonopole[[2]];,
    xcut = λ*a;
    xlo = bandScale*etaCutBand[[1]]*xcut;
    xhi = bandScale*etaCutBand[[2]]*xcut;
  ];
  {xlo, xhi}
];

OmegaBandFromKa[m_Integer, n_Integer, bandScale_:1.0] := Module[{xlo, xhi},
  {xlo, xhi} = KaBand[m, n, bandScale];
  {(c/a)*xlo, (c/a)*xhi}
];

(* ----------------------------
   Grid ω generation (Grid metric)
   ---------------------------- *)
MakeOmegaGridKa[m_Integer, n_Integer, L_?NumericQ, bandScale_:1.0] := Module[
  {ωlo, ωhi, ωgrid},
  {ωlo, ωhi} = OmegaBandFromKa[m, n, bandScale];
  If[!(FiniteQ[ωlo] && FiniteQ[ωhi]) || ωlo <= 0 || ωhi <= ωlo, Return[{}]];
  ωgrid = Subdivide[ωlo, ωhi, nΩ];
  Select[ωgrid, PoleSafeQ[#, m, n, L] &]
];

(* ----------------------------
   Resonance list + samples (Resonance metric)
   ---------------------------- *)

ResonanceOmega[m_Integer, n_Integer, q_Integer, L_?NumericQ] := Module[{λ, κres},
  λ = Lambda[m, n];
  κres = Which[
    bottomBC === "Dirichlet", (q*Pi)/L,
    bottomBC === "Neumann",  ((q + 1/2)*Pi)/L,
    True, $Failed
  ];
  c*Sqrt[λ^2 + κres^2]
];

ResonanceQRange[m_Integer, n_Integer, L_?NumericQ, bandScale_:1.0] := Module[
  {ωlo, ωhi, ωloEff, ωhiEff, λ, klo2, khi2, κlo, κhi, qMin, qMax},

  {ωlo, ωhi} = OmegaBandFromKa[m, n, bandScale];
  λ = Lambda[m, n];

  (*
    Use an *expanded* ω window so that resonances just outside the band
    can still contribute samples when we detune by up to maxDetune.
      ω = ω_res (1 ± δ),  δ ∈ detuneList
    A resonance can contribute in-band samples if:
      ω_res (1 + maxDetune) >= ω_lo   and   ω_res (1 - maxDetune) <= ω_hi
    => ω_res ∈ [ ω_lo/(1+maxDetune),  ω_hi/(1-maxDetune) ]
  *)
  ωloEff = ωlo/(1 + maxDetune);
  ωhiEff = ωhi/(1 - maxDetune);

  (* If ωhiEff is below cutoff, there are no propagating resonances. *)
  If[ωhiEff <= c*Abs[λ], Return[{}]];

  klo2 = Max[0, (ωloEff/c)^2 - λ^2];
  khi2 = Max[0, (ωhiEff/c)^2 - λ^2];
  κlo = Sqrt[klo2];
  κhi = Sqrt[khi2];

  Which[
    bottomBC === "Neumann",
      qMin = Ceiling[(κlo*L/Pi) - 1/2];
      qMax = Floor[(κhi*L/Pi) - 1/2];
      qMin = Max[qMin, 0];
      qMax = Max[qMax, -1];,

    bottomBC === "Dirichlet",
      qMin = Ceiling[(κlo*L/Pi)];
      qMax = Floor[(κhi*L/Pi)];
      qMin = Max[qMin, 1];
      qMax = Max[qMax, 0];,

    True,
      Return[{}]
  ];

  If[qMax < qMin, {}, Range[qMin, qMax]]
];


MakeOmegaSamplesResonance[m_Integer, n_Integer, L_?NumericQ, bandScale_:1.0] := Module[
  {lam, omegaCut, qList, omegasByRes, omegaAll, omegaSafe, omegaLo, omegaHi, resScales, nRes},

  lam = Lambda[wallBC, m, n];
  omegaCut = c*lam;

  (* Detune factors: avoid exact poles by not including 0 detune. *)
  resScales = Flatten@Table[{1 - d, 1 + d}, {d, detuneList}];
  resScales = Select[DeleteDuplicates[resScales], # > 0 &];

  (* Choose the resonance list. *)
  qList = Which[
    resSamplingMode === "FixedQ", qFixedList,
    True, (* fallback to 3.1b behavior: use any resonances that land in etaCutBand (optionally scaled) *)
      Module[{etaLo, etaHi, qMaxTry, qs, okQ},
        etaLo = (etaCutBand[[1]]*bandScale);
        etaHi = (etaCutBand[[2]]*bandScale);
        qMaxTry = 200;
        qs = Range[0, qMaxTry];
        okQ[q_] := Module[{w = ResonanceOmega[m, n, q, L]}, (etaLo*omegaCut <= w <= etaHi*omegaCut)];
        Select[qs, okQ]
      ]
  ];

  (* Optionally filter ω samples to the overlap band relative to cutoff. *)
  If[applyKaFilterHighMode,
    omegaLo = (etaCutBand[[1]]*bandScale)*omegaCut;
    omegaHi = (etaCutBand[[2]]*bandScale)*omegaCut;
  ,
    omegaLo = -Infinity; omegaHi = Infinity;
  ];

  omegasByRes = Table[
    Module[{w0 = ResonanceOmega[m, n, q, L], ws},
      ws = w0*resScales;
      ws = Select[ws, omegaLo <= # <= omegaHi &];
      ws
    ],
    {q, qList}
  ];

  omegaAll = Flatten[omegasByRes];

  (* Pole safety and basic finiteness checks. *)
  omegaSafe = Select[omegaAll, PoleSafeQ[#, m, n, L] &];
  omegaSafe = Sort[DeleteDuplicates[omegaSafe]];

  (* Keep runtime bounded: cap the number of ω points used. *)
  omegaSafe = SubsampleEven[omegaSafe, maxSamplePoints];

  (* Count how many resonances contributed at least one ω sample. *)
  nRes = Count[omegasByRes, lst_ /; (Length@Select[lst, MemberQ[omegaSafe, #] &] > 0)];

  {omegaSafe, nRes}
];

(* ----------------------------
   Matching stats (Grid and Resonance)
   Return: {mismatch, coverage, Nω_keep, medAbsZin, medAbsZout, Nres}
   ---------------------------- *)

MatchStatsGrid[m_Integer, n_Integer, LoveraVal_?NumericQ, bandScale_:1.0] := Module[
  {L, ωgrid, ZinVals, ZoutVals, keep, zinKeep, zoutKeep,
   mismatchVals, cov, mm, medIn, medOut},

  L = LoveraVal*a;
  ωgrid = MakeOmegaGridKa[m, n, L, bandScale];
  If[Length[ωgrid] < minOmegaPoints, Return[Indeterminate]];

  ZinVals  = Zin[#, m, n, L] & /@ ωgrid;
  ZoutVals = Zout[#, m, n] & /@ ωgrid;

  keep = MapThread[FiniteQ[#1] && FiniteQ[#2] &, {ZinVals, ZoutVals}];
  cov  = N[Count[keep, True]]/Length[ωgrid];
  If[cov < coverageMin, Return[Indeterminate]];

  zinKeep = Pick[ZinVals, keep];
  zoutKeep= Pick[ZoutVals, keep];

  mismatchVals = Table[
    Abs[zinKeep[[i]] - zoutKeep[[i]]]^2 / (Abs[zoutKeep[[i]]]^2 + denomFloorMismatch),
    {i, Length[zinKeep]}
  ];

  mm = TrimmedMeanLocal[mismatchVals, trimFrac];
  medIn  = MedianSafe[Abs /@ zinKeep];
  medOut = MedianSafe[Abs /@ zoutKeep];

  {mm, cov, Length[zinKeep], medIn, medOut, 0}
];

MatchStatsResonance[m_Integer, n_Integer, LoveraVal_?NumericQ, bandScale_:1.0] := Module[
  {L, ωsamples, nRes, ZinVals, ZoutVals, keep, zinKeep, zoutKeep,
   mismatchVals, cov, mm, medIn, medOut, planned},

  L = LoveraVal*a;
  {ωsamples, nRes} = MakeOmegaSamplesResonance[m, n, L, bandScale];
  planned = Length[ωsamples];
  If[planned < minSamplePoints, Return[Indeterminate]];

  ZinVals  = Zin[#, m, n, L] & /@ ωsamples;
  ZoutVals = Zout[#, m, n] & /@ ωsamples;

  keep = MapThread[FiniteQ[#1] && FiniteQ[#2] &, {ZinVals, ZoutVals}];
  cov  = N[Count[keep, True]]/planned;
  If[cov < coverageMin, Return[Indeterminate]];

  zinKeep = Pick[ZinVals, keep];
  zoutKeep= Pick[ZoutVals, keep];

  mismatchVals = Table[
    Abs[zinKeep[[i]] - zoutKeep[[i]]]^2 / (Abs[zoutKeep[[i]]]^2 + denomFloorMismatch),
    {i, Length[zinKeep]}
  ];

  mm = TrimmedMeanLocal[mismatchVals, trimFrac];
  medIn  = MedianSafe[Abs /@ zinKeep];
  medOut = MedianSafe[Abs /@ zoutKeep];

  {mm, cov, Length[zinKeep], medIn, medOut, nRes}
];

MatchStats[m_Integer, n_Integer, LoveraVal_?NumericQ, bandScale_:1.0] := Module[{λ},
  λ = Lambda[m, n];
  Which[
    metricMode === "Grid",      MatchStatsGrid[m, n, LoveraVal, bandScale],
    metricMode === "Resonance", MatchStatsResonance[m, n, LoveraVal, bandScale],
    metricMode === "Hybrid",
      If[Abs[λ] < 10^-12,
        MatchStatsGrid[m, n, LoveraVal, bandScale],
        MatchStatsResonance[m, n, LoveraVal, bandScale]
      ],
    True,
      MatchStatsResonance[m, n, LoveraVal, bandScale]
  ]
];

(* ----------------------------
   Scan utilities
   ---------------------------- *)
BestFromScan[data_] := Module[{finite, idx},
  If[!ListQ[data], Return[Indeterminate]];
  finite = Select[data, FiniteQ[#[[2]]] &];
  If[Length[finite] == 0, Return[Indeterminate]];
  idx = First@Ordering[finite[[All, 2]], 1];
  finite[[idx]]
];

RunModeScan[m_Integer, n_Integer, bandScale_:1.0] := Module[{rows, stats},
  rows = Table[
    stats = MatchStats[m, n, Lov, bandScale];
    If[stats === Indeterminate,
      {Lov, Indeterminate, Indeterminate, 0, 0, Indeterminate, Indeterminate},
      {Lov, stats[[1]], stats[[2]], stats[[3]], stats[[6]], stats[[4]], stats[[5]]}
    ],
    {Lov, scanLovera}
  ];
  rows
];

(* ----------------------------
   Quick spectrum/cutoff summary
   ---------------------------- *)
CutoffOmega[m_Integer, n_Integer] := c*Lambda[m, n];

(* ----------------------------
   Header
   ---------------------------- *)
If[cliMode,
  Print["=== Cylinder DtN Step 3.2 Summary ==="]; 
  Print["wallBC=", wallBC, "  bottomBC=", bottomBC, "  a=", NumStr[a, 6], "  c=", NumStr[c, 6]];
  Print["L0=", NumStr[L0, 6], "  (L0/a)=", NumStr[Lovera0, 6]];
  Print["gammaMode=", gammaMode, "  deltaRel=", NumStr[deltaRel, 6], "  gammaAbs=", NumStr[gammaAbs, 6], "  gammaFloor=", NumStr[gammaFloor, 6]];
  Print["outerModel=", outerModel, "  metricMode=", metricMode, "  nΩ=", nΩ];
  Print["resSamplingMode=", resSamplingMode, "  qFixedList=", qFixedList, "  applyKaFilterHighMode=", applyKaFilterHighMode, "  maxSamplePoints=", maxSamplePoints];
  Print["kaBandMonopole=", kaBandMonopole, "  etaCutBand=", etaCutBand, "  bandJitterFrac=", bandJitterFrac];
  Print["detuneList=", detuneList, "  (#scales=", Length[resScaleFactors], ")  maxDetune=", NumStr[maxDetune, 6]];
  Print["minRes=", minResonances, "  maxRes=", maxResonances, "  minSamplePoints=", minSamplePoints];
  Print["poleGuard=", NumStr[poleGuard, 6], "  denTol=", NumStr[denTol, 6], "  minOmegaPoints=", minOmegaPoints, "  coverageMin=", NumStr[coverageMin, 4]];
  Print["modes: ", modeList];
];

If[cliMode && doSpectrum,
  Print["--- Spectrum/cutoff summary at L0 ---"];
  Do[
    Module[{m = mn[[1]], n = mn[[2]], ωcut, kaLoHi, ωLoHi, res},
      ωcut = CutoffOmega[m, n];
      kaLoHi = KaBand[m, n, 1.0];
      ωLoHi  = OmegaBandFromKa[m, n, 1.0];
      Print["Mode (", m, ",", n, ") : cutoff ω_cut=", NumStr[ωcut, 12],
            " ; band ka=", kaLoHi, " ; band ω={", NumStr[ωLoHi[[1]], 10], ", ", NumStr[ωLoHi[[2]], 10], "}"];
      res = Table[ResonanceOmega[m, n, q, L0], {q, 0, 4}];
      Print["  ω_res(q=0..4) @ L0: ", StringRiffle[NumStr[#, 10] & /@ res, ", "]];
    ],
    {mn, modeList}
  ];
];

If[cliMode,
  If[runScan,
    Print["--- Step 3.2: hybrid matching scans (Grid λ=0; Fixed-q resonance λ>0) ---"];,
    Print["--- Step 3.1b: single-point evaluation at L0 (no -scan) ---"]; 
  ];
];

resultsByMode = Association[];

Do[
  Module[{m = mn[[1]], n = mn[[2]], rows, best, stats0, bandScales, bests, drift},

    If[cliMode, Print[""]; Print["### Mode (m,n)=(", m, ",", n, ") ###"];];

    If[runScan,
      rows = RunModeScan[m, n, 1.0];

      If[cliMode,
        Print["L/a\tmismatch\tcoverage\tNω\tNres\tmed|Zin|\tmed|Zout|"];
        Scan[
          (Print[
            NumStr[#[[1]], 6], "\t",
            If[FiniteQ[#[[2]]], NumStr[#[[2]], 12], "NaN"], "\t",
            If[FiniteQ[#[[3]]], NumStr[#[[3]], 6], "NaN"], "\t",
            ToString[#[[4]], OutputForm], "\t",
            ToString[#[[5]], OutputForm], "\t",
            If[FiniteQ[#[[6]]], NumStr[#[[6]], 12], "NaN"], "\t",
            If[FiniteQ[#[[7]]], NumStr[#[[7]], 12], "NaN"]
          ]) &,
          rows
        ];
      ];

      best = BestFromScan[rows];

      If[cliMode,
        If[best === Indeterminate,
          Print["Best L/a: Indeterminate (insufficient finite points)"];,
          Print["Best L/a=", NumStr[best[[1]], 6],
                "  mismatch=", NumStr[best[[2]], 12],
                "  coverage=", NumStr[best[[3]], 6],
                "  Nω=", best[[4]],
                "  Nres=", best[[5]]
          ];
        ];
      ];

      bandScales = {1 - bandJitterFrac, 1.0, 1 + bandJitterFrac};
      bests = Table[
        Module[{rowsB, b},
          rowsB = RunModeScan[m, n, s];
          b = BestFromScan[rowsB];
          {s, b}
        ],
        {s, bandScales}
      ];

      If[cliMode,
        Print["--- Robustness: band scaling (ka endpoints × s) ---"];
        Scan[
          (If[#[[2]] === Indeterminate,
             Print["s=", NumStr[#[[1]], 4], "  best_L/a=Indeterminate"];,
             Print["s=", NumStr[#[[1]], 4],
                   "  best_L/a=", NumStr[#[[2, 1]], 6],
                   "  mismatch=", NumStr[#[[2, 2]], 12],
                   "  coverage=", NumStr[#[[2, 3]], 6],
                   "  Nres=", ToString[#[[2, 5]], OutputForm]
             ];
          ]) &,
          bests
        ];
      ];

      drift = If[AllTrue[bests[[All, 2]], # =!= Indeterminate &],
                 Max[bests[[All, 2, 1]]] - Min[bests[[All, 2, 1]]],
                 Indeterminate];

      If[cliMode,
        Print["Band-jitter drift Δ(L/a)=", If[drift === Indeterminate, "Indeterminate", NumStr[drift, 6]]];
      ];

      resultsByMode[{m, n}] = <|"rows" -> rows, "best" -> best, "bests" -> bests, "drift" -> drift|>;

      ,
      stats0 = MatchStats[m, n, Lovera0, 1.0];
      If[cliMode,
        If[stats0 === Indeterminate,
          Print["At L0: Indeterminate (coverage/points failed)"];,
          Print["At L0: mismatch=", NumStr[stats0[[1]], 12],
                "  coverage=", NumStr[stats0[[2]], 6],
                "  Nω=", stats0[[3]],
                "  Nres=", stats0[[6]],
                "  med|Zin|=", NumStr[stats0[[4]], 12],
                "  med|Zout|=", NumStr[stats0[[5]], 12]
          ];
        ];
      ];
      resultsByMode[{m, n}] = <|"stats0" -> stats0|>;
    ];
  ],
  {mn, modeList}
];

(* ----------------------------
   Plotting (Notebook only)
   ---------------------------- *)
If[plotMode && runScan,
  Do[
    Module[{m = mn[[1]], n = mn[[2]], rows, good, p1, p2},
      If[!KeyExistsQ[resultsByMode, {m, n}], Continue[]];
      rows = resultsByMode[{m, n}]["rows"];
      good = Select[rows, FiniteQ[#[[2]]] && FiniteQ[#[[3]]] &];
      If[Length[good] == 0, Continue[]];

      p1 = ListLinePlot[
        good[[All, {1, 2}]],
        PlotRange -> All,
        AxesLabel -> {"L/a", "mismatch"},
        PlotLabel -> Row[{"Mode (", m, ",", n, ") mismatch vs L/a"}]
      ];

      p2 = ListLinePlot[
        good[[All, {1, 3}]],
        PlotRange -> {0, 1},
        AxesLabel -> {"L/a", "coverage"},
        PlotLabel -> Row[{"Mode (", m, ",", n, ") coverage vs L/a"}]
      ];

      Print[p1];
      Print[p2];
    ],
    {mn, modeList}
  ];
];
