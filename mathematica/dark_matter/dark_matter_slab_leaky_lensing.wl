(* ===================================================================== *)
(* dark_matter_slab_test_revised_v2.wl                                    *)
(* ===================================================================== *)
(* Paper-grade validation harness for the slab (finite bulk thickness ±H): *)
(*   - Two equivalent representations (images + KK/modes)                 *)
(*   - Correct far-field normalization: F ~ GM/(H r)                       *)
(*   - Convergence + log-slope diagnostics                                *)
(*   - Tail-corrected image force sum (fixes large-r mismatch)            *)
(*   - Optional far-field disk sanity check + (slow) full disk mode build *)
(* ===================================================================== *)

ClearAll["Global`*"];
$Assumptions = {r > 0, H > 0};

SetOptions[Plot, ImageSize -> 500, PlotRange -> All, PlotTheme -> "Scientific"];
SetOptions[ListLinePlot, ImageSize -> 500, PlotRange -> All, PlotTheme -> "Scientific"];

(* ------------------------------- *)
(* Run toggles (set False for speed) *)
(* ------------------------------- *)
RunConvergenceStudy = True;
RunCrossCheckTables = True;
RunMainPlots = False;
RunKneeEstimate = True;
RunOptionalDiskFar = True;
RunOptionalDiskFull = False; (* slow *)
RunLeakyStudy = False;

(* Cross-check plot tuning.
   After tail-correction, images and modes agree to ~machine precision, so the
   raw relative-error curve is dominated by round-off and Plot's adaptive sampling.
   These controls make that plot more interpretable (fixed grid + log-smoothing). *)
MakeRelDiffPlot = False;
RelDiffRMin = 0.3;          (* avoid the r≈0 singular region for the discrepancy plot *)
RelDiffGridN = 400;         (* fixed log-spaced samples for the discrepancy plot *)
RelDiffSmoothWindow = 11;   (* odd integer; moving-average smoothing in log10 space *)
RelDiffFloor = 1.*^-16;     (* display floor to avoid log(0) artifacts *)
RelDiffUseHighPrecision = False; (* set True for a cleaner curve (slower) *)
RelDiffWorkingPrecision = 50;     (* used only if RelDiffUseHighPrecision=True *)

(* ------------------------------- *)
(* Units & conventions             *)
(* ------------------------------- *)
(* Use dimensionless units: Phi_N = 1/r, F_N = 1/r^2. *)
G = 1.0;

fmt[x_?NumericQ, spec_: {10, 8}] := ToString[NumberForm[x, spec]];

Print["\n=============================================================="];
Print[" TEST 1: SLAB VACUUM (FINITE BULK THICKNESS)"];
Print["=============================================================="];
Print["Model: Neumann boundaries at z=±H implemented by image charges at z=2 n H."];
Print["Goal: verify crossover F(r) ~ 1/r^2 (r<<H)  ->  F(r) ~ 1/(H r) (r>>H)."];

(* ===================================================================== *)
(* 1A) Method-of-images force on the brane (z=0)                           *)
(* ===================================================================== *)

PotentialSlabImages[r_?NumericQ, H_?NumericQ, nMax_Integer?NonNegative] := Module[
  {term0, rest},
  term0 = 1.0/r;
  rest = If[nMax == 0, 0.0, 2.0*Sum[1.0/Sqrt[r^2 + (2.0*n*H)^2], {n, 1, nMax}]];
  term0 + rest
];

ForceSlabImages[r_?NumericQ, H_?NumericQ, nMax_Integer?NonNegative] := Module[
  {term0, rest},
  term0 = 1.0/r^2;
  rest = If[nMax == 0, 0.0, 2.0*Sum[r/(r^2 + (2.0*n*H)^2)^(3/2), {n, 1, nMax}]];
  term0 + rest
];

(* Tail correction for the images force sum.
   For spacing Δu = 2H and tail starting at u0 ≈ 2H (nMax + 1/2):
     tail ≈ (1/H) ∫_{u0}^{∞} r/(r^2 + u^2)^{3/2} du
          = (1/H) [ 1/r - u0/(r Sqrt[r^2+u0^2]) ].
*)
TailForceImages[r_?NumericQ, H_?NumericQ, nMax_Integer?NonNegative] := Module[
  {u0},
  u0 = 2.0*H*(nMax + 0.5);
  (1.0/H)*(1.0/r - u0/(r*Sqrt[r^2 + u0^2]))
];

ForceSlabImagesTail[r_?NumericQ, H_?NumericQ, nMax_Integer?NonNegative] :=
  ForceSlabImages[r, H, nMax] + TailForceImages[r, H, nMax];

(* ===================================================================== *)
(* 1B) KK / mode representation (fast at large r)                          *)
(* ===================================================================== *)

(* For Neumann walls at z=±H and a brane source at z=0:
     F(r) = (1/H) [ 1/r + 2 Σ_{m>=1} k_m K1(k_m r) ],  k_m = m π / H.
*)
ForceSlabMode[r_?NumericQ, H_?NumericQ, mMax_Integer?Positive] := Module[
  {sum, k},
  sum = 1.0/r;
  sum += 2.0*Sum[(k = (Pi*m)/H; k*BesselK[1, k*r]), {m, 1, mMax}];
  sum/H
];

(* Potential differences (avoids log gauge issue):
     Phi(r)-Phi(r0) = (1/H)[ -Log(r/r0) + 2 Σ (K0(k r)-K0(k r0)) ].
*)
PotentialSlabModeDelta[r_?NumericQ, r0_?NumericQ, H_?NumericQ, mMax_Integer?Positive] := Module[
  {sum, k},
  sum = -Log[r/r0];
  sum += 2.0*Sum[(k = (Pi*m)/H; (BesselK[0, k*r] - BesselK[0, k*r0])), {m, 1, mMax}];
  sum/H
];

(* ===================================================================== *)
(* Diagnostics                                                            *)
(* ===================================================================== *)

LogSlope[f_, r_?NumericQ, eps_: 1.0*^-4] := Module[
  {rp, rm},
  rp = r*(1.0 + eps);
  rm = r*(1.0 - eps);
  (Log[f[rp]] - Log[f[rm]])/(Log[rp] - Log[rm])
];

RelativeError[x_?NumericQ, xRef_?NumericQ] := Abs[(x - xRef)/xRef];

(* ===================================================================== *)
(* Parameters                                                             *)
(* ===================================================================== *)

HVal = 10.0;
MVal = 1.0;

rMin = 0.1;
rMax = 300.0;

nMaxMain = 5000;   (* images truncation *)
mMaxMain = 300;    (* mode truncation *)

Print["\n--- Parameters ---"];
Print["H = ", HVal, "   |   M = ", MVal, "   |   nMax(images) = ", nMaxMain, "   |   mMax(modes) = ", mMaxMain];
Print["Plot range r ∈ [", rMin, ", ", rMax, "]"];

(* Force models (use tail-corrected images) *)
FImages[r_?NumericQ] := G*MVal*ForceSlabImagesTail[r, HVal, nMaxMain];
FMode[r_?NumericQ] := G*MVal*ForceSlabMode[r, HVal, mMaxMain];

vFromForce[f_] := Function[{rr}, Sqrt[rr*f[rr]]];

vImages = vFromForce[FImages];
vMode = vFromForce[FMode];

vNewton[r_?NumericQ] := Sqrt[G*MVal/r];
vFlatAsym[r_?NumericQ] := Sqrt[G*MVal/HVal];

(* ===================================================================== *)
(* PART 1C: Leaky-slab parameterizations (boundary sensitivity)           *)
(* ===================================================================== *)

LeakyModel = "weightedImages";   (* "weightedImages" or "liftedZero" *)

(* --- Model A: weighted images (partial reflection) --- *)
ForceSlabImagesWeighted[r_?NumericQ, H_?NumericQ, nMax_Integer?NonNegative, w_?NumericQ] := Module[
  {term0, rest},
  term0 = 1.0/r^2;
  rest = If[nMax == 0, 0.0, 2.0*Sum[(w^n)*r/(r^2 + (2.0*n*H)^2)^(3/2), {n, 1, nMax}]];
  term0 + rest
];

(* For w<1 the weighted sum converges quickly; tail correction is usually unnecessary.
   If you run w extremely close to 1 (e.g. 0.9999), increase nMaxWeighted. *)
nMaxWeighted = 10000;

FWeighted[w_?NumericQ] := Function[{rr}, G*MVal*ForceSlabImagesWeighted[rr, HVal, nMaxWeighted, w]];

(* --- Model B: lifted zero mode (screened / Robin-like) --- *)
(* Replace the Neumann zero-mode term 1/r by k0 K1(k0 r), which matches 1/r for k0 r << 1
   and screens exponentially for k0 r >> 1. This is a convenient phenomenological proxy
   for Robin boundary conditions that lift the k=0 mode. *)
ForceSlabModeLiftedZero[r_?NumericQ, H_?NumericQ, mMax_Integer?Positive, k0_?NumericQ] := Module[
  {sum, k},
  sum = If[k0 == 0.0, 1.0/r, k0*BesselK[1, k0*r]];
  sum += 2.0*Sum[(k = (Pi*m)/H; k*BesselK[1, k*r]), {m, 1, mMax}];
  sum/H
];

FLiftedZero[k0_?NumericQ] := Function[{rr}, G*MVal*ForceSlabModeLiftedZero[rr, HVal, mMaxMain, k0]];

(* ===================================================================== *)
(* PART 1D: Lensing / deflection integral from a given force law          *)
(* ===================================================================== *)

cVal = 1.0; (* restore physical c later; matches Paper II convention when c0->c *)

DeflectionFromForce[f_, b_?NumericQ] := Module[
  {integrand},
  integrand[u_?NumericQ] := f[b*Sec[u]]*Sec[u];
  (2.0*b/cVal^2)*NIntegrate[integrand[u], {u, -Pi/2, Pi/2},
    Method -> "GlobalAdaptive", MaxRecursion -> 12, AccuracyGoal -> 8, PrecisionGoal -> 8
  ]
];

(* Convenience wrappers *)
DeflMode[b_?NumericQ] := DeflectionFromForce[FMode, b];
DeflImages[b_?NumericQ] := DeflectionFromForce[FImages, b];

(* Reference asymptotes *)
DeflNewton[b_?NumericQ] := 4.0*G*MVal/(b*cVal^2);
DeflPlateau[b_?NumericQ] := 2.0*Pi*G*MVal/(HVal*cVal^2);



(* ===================================================================== *)
(* Critical checks                                                        *)
(* ===================================================================== *)

Print["\n--- Critical Checks (point mass) ---"];

rSmall = 0.05*HVal;
rLarge = 100.0*HVal;

nearCheckImages = FImages[rSmall]*rSmall^2;
nearCheckMode = FMode[rSmall]*rSmall^2;

farCheckImages = FImages[rLarge]*HVal*rLarge;
farCheckMode = FMode[rLarge]*HVal*rLarge;

Print["Near-field check (r << H):  r^2 F(r) -> 1"];
Print["  Images: r=", rSmall, "  r^2 F=", fmt[nearCheckImages, {12, 9}]];
Print["  Modes : r=", rSmall, "  r^2 F=", fmt[nearCheckMode, {12, 9}]];

Print["Far-field check (r >> H):  H r F(r) -> 1"];
Print["  Images: r=", rLarge, "  H r F=", fmt[farCheckImages, {12, 9}]];
Print["  Modes : r=", rLarge, "  H r F=", fmt[farCheckMode, {12, 9}]];

Print["(If H r F -> 1, the far-field normalization is GM/(H r).)"];

(* ===================================================================== *)
(* Convergence study                                                      *)
(* ===================================================================== *)

Print["\n--- Convergence Study ---"];
If[RunConvergenceStudy,

  rTestPts = {0.1*HVal, 0.3*HVal, 1.0*HVal, 3.0*HVal, 10.0*HVal, 30.0*HVal, 100.0*HVal, 300.0*HVal, 100.0*HVal*10};

  nRef = 20000;
  mRef = 800;

  FImagesRef[r_?NumericQ] := G*MVal*ForceSlabImagesTail[r, HVal, nRef];
  FModeRef[r_?NumericQ] := G*MVal*ForceSlabMode[r, HVal, mRef];

  nList = {200, 500, 1000, 2000, 5000};
  mList = {30, 50, 100, 200, 300};

  convImages = Table[
    With[{fr = FImagesRef[rt]},
      Prepend[Table[RelativeError[G*MVal*ForceSlabImagesTail[rt, HVal, nn], fr], {nn, nList}], rt]
    ],
    {rt, rTestPts}
  ];

  convModes = Table[
    With[{fr = FModeRef[rt]},
      Prepend[Table[RelativeError[G*MVal*ForceSlabMode[rt, HVal, mm], fr], {mm, mList}], rt]
    ],
    {rt, rTestPts}
  ];

  Print["Relative error vs reference (Images+tail). Columns: r, then nMax = ", nList];
  Print@TableForm[convImages, TableHeadings -> {None, Prepend[nList, "r"]}];

  Print["\nRelative error vs reference (Modes). Columns: r, then mMax = ", mList];
  Print@TableForm[convModes, TableHeadings -> {None, Prepend[mList, "r"]}];
];

(* ===================================================================== *)
(* Cross-check: images vs modes                                           *)
(* ===================================================================== *)

Print["\n--- Cross-check: Images vs Modes (main truncations) ---"];
If[RunCrossCheckTables,

  crossPts = {1., 3., 10., 30., 100., 300., 1000.};

  cross = Table[
    {rt, FImages[rt], FMode[rt], RelativeError[FMode[rt], FImages[rt]]},
    {rt, crossPts}
  ];

  Print@TableForm[cross, TableHeadings -> {None, {"r", "F_images", "F_modes", "rel err (modes vs images)"}}];
];

(* ===================================================================== *)
(* Plots                                                                  *)
(* ===================================================================== *)

Print["\n--- Plots ---"];
If[RunMainPlots,

  plotForce = LogLogPlot[
    {FImages[r], FMode[r], 1/r^2, 1/(HVal*r)},
    {r, rMin, rMax},
    PlotLegends -> {"F(r) images (tail)", "F(r) modes", "reference ~ r^-2", "reference ~ (H r)^-1"},
    AxesLabel -> {"r", "F(r)"},
    PlotLabel -> "Slab Force: crossover from 1/r^2 to 1/(H r)"
  ];

  plotV = Plot[
    {vImages[r], vMode[r], vNewton[r], vFlatAsym[r]},
    {r, rMin, rMax},
    PlotLegends -> {"v(r) images (tail)", "v(r) modes", "Kepler (Newton)", "Flat asymptote sqrt(GM/H)"},
    AxesLabel -> {"r", "v(r)"},
    PlotLabel -> "Rotation Curve: Keplerian -> Flat crossover"
  ];

  plotSlope = Plot[
    {LogSlope[FImages, r], LogSlope[FMode, r]},
    {r, rMin, rMax},
    PlotRange -> {-2.2, -0.8},
    PlotLegends -> {"d ln F / d ln r (images)", "d ln F / d ln r (modes)"},
    AxesLabel -> {"r", "slope"},
    PlotLabel -> "Log-slope diagnostic (should go -2 -> -1)"
  ];

  plotRelDiff = LogLogPlot[
    Evaluate[
      If[!MakeRelDiffPlot,
        (* Fallback: simple expression (may look jagged near the noise floor) *)
        RelativeError[FMode[r], FImages[r]],
        (* Preferred: fixed log-grid + log-space smoothing (more interpretable) *)
        Module[{r0, rGrid, errRaw, errFloor, win, rSmooth, errSmooth, dataRaw, dataSmooth},
          r0 = Max[rMin, RelDiffRMin];
          rGrid = Exp[Subdivide[Log[r0], Log[rMax], RelDiffGridN]];
          errRaw = RelativeError[FMode[#], FImages[#]] & /@ rGrid;
          errFloor = Max[#, RelDiffFloor] & /@ errRaw;
          win = Max[3, If[OddQ[RelDiffSmoothWindow], RelDiffSmoothWindow, RelDiffSmoothWindow + 1]];
          rSmooth = MovingAverage[rGrid, win];
          errSmooth = 10^(MovingAverage[Log10[errFloor], win]);
          dataRaw = Transpose[{rGrid, errFloor}];
          dataSmooth = Transpose[{rSmooth, errSmooth}];
          (* Return a graphics-ready expression: we'll plot the two datasets via Interpolation *)
          {Interpolation[dataRaw, InterpolationOrder -> 1][r], Interpolation[dataSmooth, InterpolationOrder -> 1][r]}
        ]
      ]
    ],
    {r, Max[rMin, RelDiffRMin], rMax},
    PlotRange -> {RelDiffFloor, 1.*^-10},
    PlotLegends -> If[MakeRelDiffPlot,
      {"raw rel err (floored)", "smoothed (log MA)"},
      {"|F_mode - F_img| / F_img"}
    ],
    AxesLabel -> {"r", "relative error"},
    PlotLabel -> "Mode vs Image discrepancy (noise-floor cross-check)"
  ];

  Print[plotForce];
  Print[plotV];
  Print[plotSlope];
  Print[plotRelDiff];
];

(* ===================================================================== *)
(* Leaky boundary studies + lensing plots                                 *)
(* ===================================================================== *)

If[RunLeakyStudy,

  Print["\n=============================================================="];
  Print[" LEAKY SLAB STUDY (BOUNDARY SENSITIVITY)"];
  Print["=============================================================="];
  Print["Model = ", LeakyModel];

  rPlotMin = rMin; rPlotMax = rMax;
  bPlotMin = 0.1*HVal; bPlotMax = 100.0*HVal;

  If[LeakyModel == "weightedImages",
    wList = {1.0, 0.99, 0.95, 0.90};
    fList = FWeighted /@ wList;
    vList = vFromForce /@ fList;
    dList = (Function[{bb}, DeflectionFromForce[#, bb]] &) /@ fList;

    plotLeakySlope = Plot[
      Evaluate[LogSlope[#, r] & /@ fList],
      {r, rPlotMin, rPlotMax},
      PlotRange -> {-2.2, -0.8},
      PlotLegends -> (("w=" <> ToString[#]) & /@ wList),
      AxesLabel -> {"r", "d ln F / d ln r"},
      PlotLabel -> "Leaky slab (weighted images): log-slope"
    ];

    plotLeakyV = Plot[
      Evaluate[#[r] & /@ vList],
      {r, rPlotMin, rPlotMax},
      PlotLegends -> (("w=" <> ToString[#]) & /@ wList),
      AxesLabel -> {"r", "v(r)"},
      PlotLabel -> "Leaky slab (weighted images): rotation curves"
    ];

    plotLeakyDefl = LogLogPlot[
      Evaluate[#[b] & /@ dList],
      {b, bPlotMin, bPlotMax},
      PlotLegends -> (("w=" <> ToString[#]) & /@ wList),
      AxesLabel -> {"b", "Δθ(b)"},
      PlotLabel -> "Leaky slab (weighted images): deflection vs impact parameter"
    ];
  ];

  If[LeakyModel == "liftedZero",
    k0List = {0.0, 0.03/HVal, 0.1/HVal, 0.3/HVal}; (* k0=0 -> Neumann *)
    fList = FLiftedZero /@ k0List;
    vList = vFromForce /@ fList;
    dList = (Function[{bb}, DeflectionFromForce[#, bb]] &) /@ fList;

    plotLeakySlope = Plot[
      Evaluate[LogSlope[#, r] & /@ fList],
      {r, rPlotMin, rPlotMax},
      PlotRange -> {-2.2, -0.8},
      PlotLegends -> (("k0=" <> ToString[#]) & /@ k0List),
      AxesLabel -> {"r", "d ln F / d ln r"},
      PlotLabel -> "Leaky slab (lifted zero mode): log-slope"
    ];

    plotLeakyV = Plot[
      Evaluate[#[r] & /@ vList],
      {r, rPlotMin, rPlotMax},
      PlotLegends -> (("k0=" <> ToString[#]) & /@ k0List),
      AxesLabel -> {"r", "v(r)"},
      PlotLabel -> "Leaky slab (lifted zero mode): rotation curves"
    ];

    plotLeakyDefl = LogLogPlot[
      Evaluate[#[b] & /@ dList],
      {b, bPlotMin, bPlotMax},
      PlotLegends -> (("k0=" <> ToString[#]) & /@ k0List),
      AxesLabel -> {"b", "Δθ(b)"},
      PlotLabel -> "Leaky slab (lifted zero mode): deflection vs impact parameter"
    ];
  ];

  Print[plotLeakySlope];
  Print[plotLeakyV];
  Print[plotLeakyDefl];

  (* Baseline deflection curve with asymptotes (Neumann slab) *)
  plotDeflection = LogLogPlot[
    {DeflMode[b], DeflNewton[b], DeflPlateau[b]},
    {b, bPlotMin, bPlotMax},
    PlotLegends -> {"Δθ(b) from F_mode", "Newtonian 4GM/(b c^2)", "Slab plateau 2πGM/(H c^2)"},
    AxesLabel -> {"b", "Δθ(b)"},
    PlotLabel -> "Deflection: Newton 1/b → slab plateau"
  ];
  Print[plotDeflection];
];


(* ===================================================================== *)
(* Knee estimate                                                          *)
(* ===================================================================== *)

Print["\n--- Knee Scale Estimate ---"];
If[RunKneeEstimate,

  slopeTarget = -1.5;
  rScan = Exp[Subdivide[Log[rMin], Log[rMax], 250]];
  slopeVals = LogSlope[FMode, #] & /@ rScan;
  idx = First@Ordering[Abs[slopeVals - slopeTarget], 1];
  rKneeGuess = rScan[[idx]];

  rKnee = Quiet@Check[
    r /. FindRoot[LogSlope[FMode, r] - slopeTarget == 0, {r, rKneeGuess}],
    rKneeGuess
  ];

  Print["Estimated knee where slope ~ -3/2: r_knee ≈ ", fmt[rKnee, {10, 6}], "  (compare to H = ", HVal, ")"];
];

(* ===================================================================== *)
(* PART 2: Wake force (separate toy mechanism)                             *)
(* ===================================================================== *)

Print["\n=============================================================="];
Print[" TEST 2: WAKE FORCE (VELOCITY-DEPENDENT ADDITIVE TERM)"];
Print["=============================================================="];
Print["Toy equation (circular orbit):  v^2/r = GM/r^2 + C*(v/r)  => v^2 - C v - GM/r = 0"];
Print["As r -> ∞: v -> C (flat set by coupling). Distinct from slab geometry."];

vWakeModel[r_?NumericQ, mass_?NumericQ, coupling_?NumericQ] :=
  (coupling + Sqrt[coupling^2 + 4.0*G*mass/r])/2.0;

massWake = 100.0;
couplingWake = 2.0;

plotWake = Plot[
  {vWakeModel[r, massWake, couplingWake], Sqrt[G*massWake/r], couplingWake},
  {r, 1, rMax},
  PlotLegends -> {"Wake model", "Pure Newton", "asymptote v=C"},
  AxesLabel -> {"r", "v(r)"},
  PlotLabel -> "Wake-force rotation curve (separate mechanism)"
];

Print[plotWake];

Print["\n--- Wake critical check ---"];
vInf = vWakeModel[10^8, massWake, couplingWake];
Print["v(r=1e8) ≈ ", fmt[vInf, {12, 9}], "   (expected -> C = ", couplingWake, ")"];

(* ===================================================================== *)
(* OPTIONAL: Extended source far-field sanity check (exponential disk)     *)
(* ===================================================================== *)

If[RunOptionalDiskFar,

  Print["\n=============================================================="];
  Print[" OPTIONAL: EXTENDED SOURCE FAR-FIELD CHECKS"];
  Print["=============================================================="];
  Print["Far-field slab regime (axisymmetric): a(r) ≈ G M_enc(r)/(H r), so v^2(r) ≈ G M_enc(r)/H."];

  SigmaExp[R_?NumericQ, Mtot_?NumericQ, Rd_?NumericQ] := (Mtot/(2.0*Pi*Rd^2))*Exp[-R/Rd];

  (* Closed-form enclosed mass for exponential disk *)
  MencExpClosed[r_?NumericQ, Mtot_?NumericQ, Rd_?NumericQ] := Mtot*(1.0 - Exp[-r/Rd]*(1.0 + r/Rd));

  aSlabFar[r_?NumericQ, H_?NumericQ, Mtot_?NumericQ, Rd_?NumericQ] := G*MencExpClosed[r, Mtot, Rd]/(H*r);
  vSlabFar[r_?NumericQ, H_?NumericQ, Mtot_?NumericQ, Rd_?NumericQ] := Sqrt[r*aSlabFar[r, H, Mtot, Rd]];

  MtotDisk = 50.0;
  RdDisk = 8.0;
  vInfDisk = Sqrt[G*MtotDisk/HVal];

  plotDiskFar = Plot[
    {vSlabFar[r, HVal, MtotDisk, RdDisk], vInfDisk},
    {r, 0.5, rMax},
    PlotLegends -> Placed[{"v(r) far-field slab law", "v_inf = sqrt(G M_tot / H)"}, Right],
    AxesLabel -> {"r", "v(r)"},
    PlotLabel -> "Exponential disk: far-field slab prediction"
  ];

  Print[plotDiskFar];

  Print["\nDone. (For full near-field disk curve, enable RunOptionalDiskFull=True.)"];
];

(* ===================================================================== *)
(* OPTIONAL (SLOW): Full near-field exponential disk via mode integrals    *)
(* ===================================================================== *)
(* Uses identity: ∫_0^{2π} K0(k |r-r'|) dθ = 2π I0(k min) K0(k max).        *)
(* For k_m = m π / H:
     A_m(r)=∫_0^r dR R Σ(R) I0(k_m R)
     B_m(r)=∫_r^∞ dR R Σ(R) K0(k_m R)
     a_m(r)=(4π k_m/H)[ K1(k_m r) A_m(r) - I1(k_m r) B_m(r) ]
   Total: a(r)=G M_enc/(H r) + G Σ_{m>=1} a_m(r).
*)

If[RunOptionalDiskFull,

  Print["\n=============================================================="];
  Print[" OPTIONAL (SLOW): FULL DISK CURVE USING MODE INTEGRALS"];
  Print["=============================================================="];

  mMaxDiskFull = 60;
  rMaxDiskFull = rMax;
  rDiskGrid = Exp[Subdivide[Log[0.5], Log[rMaxDiskFull], 60]];

  ClearAll[Aint, Bint];

  Aint[k_?NumericQ, rr_?NumericQ] := Aint[k, rr] =
    NIntegrate[
      R*SigmaExp[R, MtotDisk, RdDisk]*BesselI[0, k*R],
      {R, 0, rr},
      Method -> "GlobalAdaptive",
      MaxRecursion -> 12,
      AccuracyGoal -> 8,
      PrecisionGoal -> 8,
      WorkingPrecision -> 40
    ];

  Bint[k_?NumericQ, rr_?NumericQ] := Bint[k, rr] =
    NIntegrate[
      R*SigmaExp[R, MtotDisk, RdDisk]*BesselK[0, k*R],
      {R, rr, Infinity},
      Method -> "GlobalAdaptive",
      MaxRecursion -> 12,
      AccuracyGoal -> 8,
      PrecisionGoal -> 8,
      WorkingPrecision -> 40
    ];

  aDiskFull[rr_?NumericQ] := Module[
    {a0, sum, k},
    a0 = G*MencExpClosed[rr, MtotDisk, RdDisk]/(HVal*rr);
    sum = Sum[
      (k = (Pi*m)/HVal;
       (4.0*Pi*k/HVal)*(BesselK[1, k*rr]*Aint[k, rr] - BesselI[1, k*rr]*Bint[k, rr])
      ),
      {m, 1, mMaxDiskFull}
    ];
    a0 + G*sum
  ];

  vDiskFull[rr_?NumericQ] := Sqrt[rr*aDiskFull[rr]];

  diskFullData = Table[{rr, vDiskFull[rr]}, {rr, rDiskGrid}];

  plotDiskFull = ListLinePlot[
    diskFullData,
    PlotLegends -> Placed[{"v(r) full disk (0-mode + m>0)"}, Right],
    AxesLabel -> {"r", "v(r)"},
    PlotLabel -> "Exponential disk: full slab prediction (slow)"
  ];

  Print[plotDiskFull];
];


(*"
Ouptput:

==============================================================
 TEST 1: SLAB VACUUM (FINITE BULK THICKNESS)
==============================================================
Model: Neumann boundaries at z=±H implemented by image charges at z=2 n H.
Goal: verify crossover F(r) ~ 1/r^2 (r<<H)  ->  F(r) ~ 1/(H r) (r>>H).

--- Parameters ---
H = 10.   |   M = 1.   |   nMax(images) = 5000   |   mMax(modes) = 300
Plot range r ∈ [0.1, 300.]

--- Critical Checks (point mass) ---
Near-field check (r << H):  r^2 F(r) -> 1
  Images: r=0.5  r^2 F=1.000037534
  Modes : r=0.5  r^2 F=1.000037534
Far-field check (r >> H):  H r F(r) -> 1
  Images: r=1000.  H r F=1.000000000
  Modes : r=1000.  H r F=1.000000000
(If H r F -> 1, the far-field normalization is GM/(H r).)

--- Convergence Study ---
Relative error vs reference (Images+tail). Columns: r, then nMax = {200, 500, 1000, 2000, 5000}
TableForm[{{1., 1.909011767326864*^-14, 4.4395622495973576*^-16, 0., 0., 0.}, {3., 5.178745771995719*^-13, 1.33836933088188*^-14, 7.435385171566001*^-16, 0., 0.}, {10., 1.577012023950812*^-11, 4.061758444736008*^-13, 2.588999635620653*^-14, 2.1221308488693875*^-15, 9.90327729472381*^-16}, {30., 1.738940842495338*^-10, 4.4780323526552605*^-12, 2.798038855699347*^-13, 1.729271225873667*^-14, 1.0401631433826568*^-15}, {100., 1.930684128195123*^-9, 4.978829848399423*^-11, 3.11903280980529*^-12, 1.9493955061283061*^-13, 5.4210108624257685*^-15}, {300., 1.7162034651480634*^-8, 4.4719517802368814*^-10, 2.8048689672960224*^-11, 1.7500649367174679*^-12, 4.0332320816460557*^-14}, {1000., 1.662964697315885*^-7, 4.857925314186856*^-9, 3.0993708034082745*^-10, 1.9480944635219544*^-11, 4.988685246148927*^-13}, {3000., 5.728385171402695*^-7, 3.6148311313666727*^-8, 2.6551374624318658*^-9, 1.731447830163158*^-10, 4.4723339615026355*^-12}, {10000., 1.3810829226526847*^-7, 8.825560657178261*^-8, 1.7861526907315924*^-8, 1.676691893507362*^-9, 4.855599429475355*^-11}}, TableHeadings -> {None, {r, 200, 500, 1000, 2000, 5000}}]

Relative error vs reference (Modes). Columns: r, then mMax = {30, 50, 100, 200, 300}
TableForm[{{1., 0.00018426484198499413, 4.292783116216137*^-7, 8.856926687946726*^-14, 0., 0.}, {3., 1.383849103514626*^-12, 0., 0., 0., 0.}, {10., 0., 0., 0., 0., 0.}, {30., 0., 0., 0., 0., 0.}, {100., 0., 0., 0., 0., 0.}, {300., 0., 0., 0., 0., 0.}, {1000., 0., 0., 0., 0., 0.}, {3000., 0., 0., 0., 0., 0.}, {10000., 0., 0., 0., 0., 0.}}, TableHeadings -> {None, {r, 30, 50, 100, 200, 300}}]

--- Cross-check: Images vs Modes (main truncations) ---
TableForm[{{1., 1.0002995450516297, 1.0002995450516299, 2.2197811247986788*^-16}, {3., 0.11198710614926968, 0.11198710614926964, 3.7176925857830005*^-16}, {10., 0.012261662448154562, 0.012261662448154548, 1.1318031193970056*^-15}, {30., 0.003335483451827395, 0.003335483451827389, 1.8202855009196472*^-15}, {100., 0.001000000000000329, 0.001000000000000323, 6.0715321659168274*^-15}, {300., 0.0003333333333333485, 0.0003333333333333334, 4.53738609185163*^-14}, {1000., 0.00010000000000004988, 0.0001, 4.987329993430832*^-13}}, TableHeadings -> {None, {r, F_images, F_modes, rel err (modes vs images)}}]

--- Plots ---

--- Knee Scale Estimate ---
Estimated knee where slope ~ -3/2: r_knee ≈ 10.496773  (compare to H = 10.)

==============================================================
 TEST 2: WAKE FORCE (VELOCITY-DEPENDENT ADDITIVE TERM)
==============================================================
Toy equation (circular orbit):  v^2/r = GM/r^2 + C*(v/r)  => v^2 - C v - GM/r = 0
As r -> ∞: v -> C (flat set by coupling). Distinct from slab geometry.

--- Wake critical check ---
v(r=1e8) ≈ 2.000000500   (expected -> C = 2.)

==============================================================
 OPTIONAL: EXTENDED SOURCE FAR-FIELD CHECKS
==============================================================
Far-field slab regime (axisymmetric): a(r) ≈ G M_enc(r)/(H r), so v^2(r) ≈ G M_enc(r)/H.

Done. (For full near-field disk curve, enable RunOptionalDiskFull=True.)
"*)
