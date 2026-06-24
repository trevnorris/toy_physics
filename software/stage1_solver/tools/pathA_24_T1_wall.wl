(* pathA_24 T1 Mathematica lead script for the frozen polar-arrow wall.

   The script transcribes only the T0 freeze-action block, checks its full
   SHA-256 hash, derives the static wall functional with rho(w), and runs
   spectrum fixtures that can return stable, unstable, or connected-vacuum
   unprotected outcomes.
*)

ClearAll[
  scriptPath, stage1Root, projectRoot, scratchDir, freezeReport,
  expectedHash, expectedLPolLines, freezeText, freezeBlock, blockLines,
  hashInput, freezeHash, linePositions, forbiddenFragments,
  forbiddenPresent, fail, exprString, dimTuple, dimString, assertDim,
  dimsReport, Mdim, Ldim, Tdim, dimless, dm, drho, dkConst, dhbar, da,
  dwDim, drhop, dgradP, dcs2, dEnergyDensity, dSigma, qpDim, uDim,
  opGradDim, opPotDim, kpDim, bDim, rho, rho0, kConst, mass, alen,
  hbar, lHalf, pNorm2, pPrime2, rhoPrime, cs2, cs20, u, staticDensity,
  dU, d2U, stationaryPoints, kP, bDepth, sigmaSpread, dSigmaDL,
  sigmaRadial, radialNegativeOmega2, confinementGapLimit,
  finiteDifferenceSpectrum, countNegative, spectrumTol, nGrid, xMax,
  phi4Vals, unstableVals, lGeo, qGeo, connectedVals, radialVals,
  phi4Control, unstableControl, connectedControl, radialBaseline,
  controlValues, allControlsPass, imposedTermScan, sourceFiles,
  termNeedles, scanNeedles, labels, tauThreshold, report, outPath,
  exportResult
];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
exprString[x_] := ToString[InputForm[FullSimplify[x]]];

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_24_T1_wall.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
projectRoot = ParentDirectory[ParentDirectory[stage1Root]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
freezeReport = FileNameJoin[{stage1Root, "reports", "pathA_24_T0_freeze.md"}];

expectedHash = "8fa41ac51e88a1464a4a5b22c6fe64fc218cf36ba2e3583d26b97c994e5da064";
expectedLPolLines = {
  "  L_pol =",
  "      (1/2) m rho a^2 (D_t^v P^i)(D_t^v P^i)",
  "    - (1/2) m rho c_s^2(rho) a^2 (partial_j P^i)(partial_j P^i)",
  "    - (1/4) m rho c_s^2(rho) (P^i P^i - 1)^2."
};

freezeText = Import[freezeReport, "Text"];
freezeBlock = StringCases[
  freezeText,
  "```freeze-action\n" ~~ Shortest[block__] ~~ "\n```" :> block
];
If[Length[freezeBlock] =!= 1, fail["freeze-action block not found uniquely"]];
freezeBlock = First[freezeBlock];
hashInput = freezeBlock <> "\n";
freezeHash = ToLowerCase[IntegerString[Hash[hashInput, "SHA256"], 16, 64]];
If[freezeHash =!= expectedHash,
  fail["freeze-action SHA-256 mismatch: expected " <> expectedHash <> ", got " <> freezeHash]
];
blockLines = StringSplit[freezeBlock, "\n"];
linePositions = Flatten[FirstPosition[blockLines, #, Missing["notfound"]] & /@ expectedLPolLines];
If[MemberQ[linePositions, Missing["notfound"]],
  fail["frozen L_pol line missing or changed"]
];
If[linePositions =!= Sort[linePositions], fail["frozen L_pol lines out of order"]];
forbiddenFragments = {"gamma_w", "w_hat", "P dot w", "delta(w)", "Z(w)|nabla P|", "A_i P^i"};
forbiddenPresent = Select[forbiddenFragments, StringContainsQ[freezeBlock, #] &];
If[forbiddenPresent =!= {}, fail["forbidden frozen-action fragments present: " <> ToString[forbiddenPresent]]];
If[! StringContainsQ[freezeBlock, "no explicit easy axis"],
  fail["baseline no-explicit-easy-axis line missing from freeze-action"]
];

dimTuple[v_] := Normal[v];
dimString[v_] := Module[{labels = {"M", "L", "T"}, powers = Normal[v], parts},
  parts = MapThread[
    If[#2 == 0, Nothing, If[#2 == 1, #1, #1 <> "^" <> ToString[#2]]] &,
    {labels, powers}
  ];
  If[parts === {}, "1", StringRiffle[parts, " "]]
];
assertDim[name_, actual_, expected_] := If[FullSimplify[actual - expected] =!= {0, 0, 0},
  fail[name <> ": expected " <> ToString[dimTuple[expected]] <> ", got " <> ToString[dimTuple[actual]]]
];

Mdim = {1, 0, 0}; Ldim = {0, 1, 0}; Tdim = {0, 0, 1}; dimless = {0, 0, 0};
dm = Mdim;
drho = -4 Ldim;
dkConst = Mdim + 18 Ldim - 2 Tdim;
dhbar = Mdim + 2 Ldim - Tdim;
da = Ldim;
dwDim = Ldim;
drhop = drho - Ldim;
dgradP = -Ldim;
dcs2 = dkConst + 4 drho - dm;
dEnergyDensity = Mdim - 2 Ldim - 2 Tdim;
dSigma = dEnergyDensity + dwDim;
qpDim = 2 dhbar - dm - drho + 2 drhop;
uDim = dkConst + 5 drho;
opGradDim = dm + drho + dcs2 + 2 da + 2 dgradP;
opPotDim = dm + drho + dcs2;
kpDim = dm + drho + dcs2 + 2 da;
bDim = dm + drho + dcs2;

assertDim["c_s^2(rho)", dcs2, 2 Ldim - 2 Tdim];
assertDim["GNLS quantum pressure", qpDim, dEnergyDensity];
assertDim["GNLS U(rho)", uDim, dEnergyDensity];
assertDim["OP gradient", opGradDim, dEnergyDensity];
assertDim["OP potential", opPotDim, dEnergyDensity];
assertDim["surface tension", dSigma, Mdim - Ldim - 2 Tdim];
assertDim["spread sigma K_P/L", kpDim - Ldim, Mdim - Ldim - 2 Tdim];
assertDim["radial saddle sigma B*a", bDim + da, Mdim - Ldim - 2 Tdim];
assertDim["omega^2", dcs2 - 2 da, -2 Tdim];

dimsReport = <|
  "energy_density" -> <|"triple_MLT" -> dimTuple[dEnergyDensity], "string" -> dimString[dEnergyDensity]|>,
  "sigma_brane" -> <|"triple_MLT" -> dimTuple[dSigma], "string" -> dimString[dSigma]|>,
  "terms" -> <|
    "gnls_quantum_pressure" -> <|"triple_MLT" -> dimTuple[qpDim], "string" -> dimString[qpDim]|>,
    "gnls_U" -> <|"triple_MLT" -> dimTuple[uDim], "string" -> dimString[uDim]|>,
    "op_gradient" -> <|"triple_MLT" -> dimTuple[opGradDim], "string" -> dimString[opGradDim]|>,
    "op_potential" -> <|"triple_MLT" -> dimTuple[opPotDim], "string" -> dimString[opPotDim]|>
  |>
|>;

$Assumptions = rho >= 0 && rho0 > 0 && kConst > 0 && mass > 0 && alen > 0 && hbar > 0 && lHalf > 0;
cs2 = 5 kConst rho^4/mass;
cs20 = 5 kConst rho0^4/mass;
u = kConst rho^5/4;
staticDensity = FullSimplify[
  hbar^2 rhoPrime^2/(8 mass rho)
  + u
  + 1/2 mass rho cs2 alen^2 pPrime2
  + 1/4 mass rho cs2 (pNorm2 - 1)^2
];
dU = D[u, rho];
d2U = D[u, {rho, 2}];
stationaryPoints = Solve[dU == 0, rho, Reals][[All, 1, 2]];
If[stationaryPoints =!= {0}, fail["unexpected U stationary points: " <> ToString[stationaryPoints]]];
If[FullSimplify[d2U == 5 kConst rho^3] =!= True, fail["unexpected U''"]];

kP = FullSimplify[mass rho0 cs20 alen^2];
bDepth = FullSimplify[mass rho0 cs20];
sigmaSpread = FullSimplify[kP Pi^2/(4 lHalf)];
dSigmaDL = FullSimplify[D[sigmaSpread, lHalf]];
If[Block[{$Assumptions = rho0 > 0 && kConst > 0 && mass > 0 && alen > 0},
    Limit[sigmaSpread, lHalf -> Infinity]
  ] =!= 0,
  fail["spread sigma limit did not vanish"]
];
sigmaRadial = FullSimplify[(2 Sqrt[2]/3) bDepth alen];
radialNegativeOmega2 = FullSimplify[-cs20/(2 alen^2)];
confinementGapLimit = Block[
  {$Assumptions = rho0 > 0 && kConst > 0 && mass > 0 && alen > 0},
  Limit[cs20 Pi^2/(4 lHalf^2), lHalf -> Infinity]
];
If[confinementGapLimit =!= 0, fail["confinement gap spread limit did not vanish"]];

spectrumTol = 10^-3;
finiteDifferenceSpectrum[pot_, xmin_, xmax_, n_Integer, count_Integer] := Module[
  {h, xs, diag, off, mat, vals},
  h = N[(xmax - xmin)/(n + 1)];
  xs = Table[N[xmin + h i], {i, 1, n}];
  diag = 2/h^2 + (pot /@ xs);
  off = ConstantArray[-1/h^2, n - 1];
  mat = SparseArray[
    {Band[{1, 1}] -> diag, Band[{1, 2}] -> off, Band[{2, 1}] -> off},
    {n, n}
  ];
  vals = Sort[Eigenvalues[N[Normal[mat]]]];
  Take[vals, count]
];
countNegative[vals_] := Count[vals, x_ /; x < -spectrumTol];

nGrid = 401;
xMax = 14.;
phi4Vals = finiteDifferenceSpectrum[
  Function[x, 2 - 3 Sech[x/Sqrt[2]]^2],
  -xMax, xMax, nGrid, 4
];
unstableVals = finiteDifferenceSpectrum[
  Function[x, -1],
  -xMax, xMax, nGrid, 4
];
lGeo = 10.;
qGeo = Pi/(2 lGeo);
connectedVals = finiteDifferenceSpectrum[
  Function[x, -N[qGeo]^2],
  -lGeo, lGeo, nGrid, 4
];
radialVals = finiteDifferenceSpectrum[
  Function[y, -2 Sech[y]^2],
  -xMax, xMax, nGrid, 4
];

phi4Control = <|
  "operator" -> "-d^2/dx^2 + 2 - 3 sech^2(x/sqrt(2))",
  "pi0" -> "Z2",
  "lowest_eigenvalues" -> phi4Vals,
  "negative_count" -> countNegative[phi4Vals],
  "expected" -> "no negative modes; translation near-zero; protected",
  "pass" -> (countNegative[phi4Vals] == 0 && First[phi4Vals] > -spectrumTol)
|>;
unstableControl = <|
  "operator" -> "-d^2/dx^2 - 1",
  "lowest_eigenvalues" -> unstableVals,
  "negative_count" -> countNegative[unstableVals],
  "expected" -> "negative mode found",
  "pass" -> (countNegative[unstableVals] >= 1)
|>;
connectedControl = <|
  "operator" -> "-d^2/dw^2 - (pi/(2L))^2, Dirichlet endpoints",
  "L" -> lGeo,
  "pi0" -> 0,
  "lowest_eigenvalues" -> connectedVals,
  "negative_count_below_tol" -> countNegative[connectedVals],
  "expected" -> "clean clamped spectrum but topological diagnosis flags unwinding",
  "pass" -> (countNegative[connectedVals] == 0 && Abs[First[connectedVals]] < 5.0*^-4)
|>;
radialBaseline = <|
  "operator_dimensionless" -> "-d^2/dy^2 - 2 sech^2(y)",
  "analytic_lowest_eigenvalue_dimensionless" -> "-1",
  "lowest_eigenvalues_dimensionless" -> radialVals,
  "negative_count" -> countNegative[radialVals],
  "expected" -> "negative transverse mode for the amplitude-collapse saddle",
  "pass" -> (countNegative[radialVals] >= 1 && Abs[First[radialVals] + 1] < 5.0*^-3)
|>;
controlValues = {phi4Control, unstableControl, connectedControl};
allControlsPass = And @@ Lookup[controlValues, "pass"] && radialBaseline["pass"];
If[! TrueQ[allControlsPass], fail["spectrum fixture failure"]];

sourceFiles = {
  FileNameJoin[{projectRoot, "research", "pde", "paper", "pde.tex"}],
  FileNameJoin[{projectRoot, "docs", "conceptual_foundation.md"}],
  FileNameJoin[{stage1Root, "reports", "pathA_23_stage0_action_and_contracts.md"}]
};
termNeedles = <|
  "V_conf" -> {"V_conf", "V_{\\mathrm{conf}}"},
  "Z(w)" -> {"Z(w)"},
  "W(w)" -> {"W(w)"},
  "B_ell" -> {"B_\\ell", "B_ℓ"},
  "k_w" -> {"k_w"}
|>;
imposedTermScan = AssociationMap[
  Function[key,
    Module[{hits = {}, textLines, rel},
      scanNeedles = termNeedles[key];
      Do[
        If[FileExistsQ[src],
          textLines = StringSplit[Import[src, "Text"], "\n"];
          rel = FileNameDrop[src, FileNameDepth[projectRoot]];
          Do[
            If[AnyTrue[scanNeedles, StringContainsQ[textLines[[i]], #] &],
              AppendTo[hits, rel <> ":" <> ToString[i]]; Break[]
            ],
            {i, Length[textLines]}
          ]
        ],
        {src, sourceFiles}
      ];
      hits
    ]
  ],
  Keys[termNeedles]
];

labels = <|
  "T1a" -> "BRANE_IMPOSED_NOT_DERIVED",
  "T1b" -> "FAIL_NO_WALL_PROFILE",
  "T1c" -> "FAIL_WALL_UNWINDS_SPHERE_VACUA",
  "T1d" -> "W_EMERGENT_BUT_UNSTABLE",
  "T1e" -> "FAIL_NO_BOUND_ZERO_MODES",
  "rollup" -> "T1_FAIL_NO_STABLE_WALL"
|>;
tauThreshold = "1000 * tau_Hubble";

report = <|
  "schema" -> "stage1_pathA_24_T1_wall_mathematica/v1",
  "engine" -> "mathematica",
  "pass" -> True,
  "freeze_fidelity_guard" -> <|
    "path" -> freezeReport,
    "sha256" -> freezeHash,
    "checked_lines" -> expectedLPolLines,
    "forbidden_fragments_absent" -> forbiddenFragments
  |>,
  "labels" -> labels,
  "dimensions" -> dimsReport,
  "derivations" -> <|
    "static_energy_density" -> exprString[staticDensity],
    "U" -> <|
      "expr" -> exprString[u],
      "dU" -> exprString[dU],
      "d2U" -> exprString[d2U],
      "stationary_points_nonnegative" -> (exprString /@ stationaryPoints),
      "single_nonnegative_minimum" -> "rho=0; no second scalar vacuum"
    |>,
    "coefficients" -> <|
      "c_s0_squared" -> exprString[cs20],
      "K_P0" -> exprString[kP],
      "B0" -> exprString[bDepth]
    |>,
    "fixed_background_minimizing_sequence" -> <|
      "domain" -> "w in [-Lhalf, Lhalf]",
      "profile" -> "P_L(w)=(sin(pi*(w+Lhalf)/(2*Lhalf)),0,0,-cos(pi*(w+Lhalf)/(2*Lhalf)))",
      "rho" -> "rho0",
      "core_at_w0" -> "P_L(0)=(1,0,0,0), flat/in-plane for this finite-box representative",
      "sigma_L" -> exprString[sigmaSpread],
      "d_sigma_d_Lhalf" -> exprString[dSigmaDL],
      "limit_Lhalf_to_infinity" -> "0",
      "width_status" -> "arbitrary_box_size_dependent_no_natural_localized_width"
    |>,
    "radial_amplitude_kink_saddle" -> <|
      "profile" -> "P=(0,0,0,tanh(w/(sqrt(2)*a)))",
      "sigma_saddle" -> exprString[sigmaRadial],
      "core_texture" -> "amplitude_collapse_P=0_not_flat_orientation",
      "transverse_negative_omega_squared" -> exprString[radialNegativeOmega2],
      "not_minimizer_reason" -> "finite fixed sigma but beaten by sigma_L -> 0 spread sequence"
    |>,
    "topology" -> <|
      "vacuum_manifold" -> "{P in R^4 | P.P = 1} = S^3",
      "sphere_dimension" -> 3,
      "pi0" -> 0,
      "connected" -> True,
      "topologically_protected" -> False
    |>,
    "unwinding" -> <|
      "deltaE_unwind" -> "0",
      "tau_unwind" -> "not_computed_no_metastable_minimum",
      "tau_threshold_preregistered_before_scan" -> tauThreshold,
      "local_minimum_against_delocalization" -> False
    |>,
    "confinement" -> <|
      "gap_spread_limit_omega_squared" -> exprString[confinementGapLimit],
      "bound_zero_modes_from_wall" -> False
    |>
  |>,
  "spectra" -> <|
    "tolerance_negative_threshold" -> spectrumTol,
    "controls" -> <|
      "phi4_kink_topologically_stable" -> phi4Control,
      "known_unstable_saddle_phi_at_zero" -> unstableControl,
      "connected_vacuum_clamped_antipodal_geodesic" -> connectedControl
    |>,
    "baseline" -> <|
      "finite_box_geodesic_sensitivity" -> connectedControl,
      "radial_soft_spin_saddle_transverse" -> radialBaseline
    |>
  |>,
  "imposed_existing_brane_terms_source_scan" -> imposedTermScan
|>;

If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir, CreateIntermediateDirectories -> True]];
outPath = FileNameJoin[{scratchDir, "pathA_24_T1_wall_mathematica.json"}];
exportResult = Export[outPath, report, "RawJSON"];
If[exportResult =!= outPath || ! FileExistsQ[outPath] || FileByteCount[outPath] <= 0,
  fail["RawJSON export failed for " <> outPath]
];

Print["PATHA_24_T1_WALL_MATHEMATICA: PASS"];
Print["freeze_sha256: ", freezeHash];
Print["T1_rollup: ", labels["rollup"]];
Print["T1a: ", labels["T1a"]];
Print["T1b: ", labels["T1b"]];
Print["T1c: ", labels["T1c"]];
Print["T1d: ", labels["T1d"]];
Print["T1e: ", labels["T1e"]];
Print["sigma_L: ", report["derivations", "fixed_background_minimizing_sequence", "sigma_L"]];
Print["sigma_L_limit: ", report["derivations", "fixed_background_minimizing_sequence", "limit_Lhalf_to_infinity"]];
Print["radial_saddle_sigma: ", report["derivations", "radial_amplitude_kink_saddle", "sigma_saddle"]];
Print["pi0_vacuum: ", report["derivations", "topology", "pi0"]];
Print["deltaE_unwind: ", report["derivations", "unwinding", "deltaE_unwind"]];
Print["tau_threshold: ", tauThreshold];
Print["confinement_gap_limit_omega2: ", report["derivations", "confinement", "gap_spread_limit_omega_squared"]];
Print["controls_pass: ", And @@ Lookup[controlValues, "pass"]];
Print["wrote: ", outPath];
Exit[0]
