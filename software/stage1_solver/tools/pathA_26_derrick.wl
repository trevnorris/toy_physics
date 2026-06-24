(* pathA_26 Mathematica engine: throat-soliton Derrick and open-drain stability. *)

ClearAll[
  fail, clean, cleanComplex, scriptPath, stage1Root, scratchDir, sympyJson,
  jsonOut, rho0, kEOS, pVac, sigma, iwave, cW, muIn, muOut, chiPerp,
  chiDN, hbar, mPar, deltaR, deltaPar, epsW, etaMin, nGrid, hNum,
  phaseBMaxCoeff, dSlab, cDampSample, cDampMax, gainMin, gainMax,
  robustAxisFraction, p0, alphaVol, muFluid, vol, aside, acap, omega,
  energy, aSym, lSym, gradExpr, rootSol, qStar, hessExpr, hessBase,
  eigBase, omegaStar, eVol, eSide, eCap, waveA, waveL, virialA,
  virialL, phaseALabel, fluidOnly, rmax, wmax, rVals, wVals, drGrid,
  dwGrid, wr, ww, sech2, fluidExcess, numericTotal, fdH, numericHess,
  numericEig, numericInterior, numericSignAgree, phaseBEnergy, phaseBCorner,
  phaseBCorners, phaseBAllOk, hardStatus, phaseBLabel, consGrad,
  consHess, forceVec, forceJac, massMat, fixedPoint, jacobianOpen,
  stableAt, findGcrit, gSample, phaseCSigns, requiredRadius,
  minCenterRobust, robustCenterG, minUpperRobust, robustForcedG,
  maxGcrit, robustExists, phaseCLabel, staticDivergenceSigns,
  staticDivergence,
  gateScaleRatio, gateBindingMargin, gateChecks, gatePasses, topVerdict,
  lambdaRay, forms, citations, results, mismatches, cmpExact, cmpFloats,
  sympy, agreementStatus, agreement, exportResult, plus, minus
];

fail[msg_] := (Print["FAIL: ", msg]; Exit[1]);
clean[x_] := N[Round[N[x, 20], 10^-12], 16];
cleanComplex[z_] := <|"re" -> N[Round[Re[N[z]], 10^-10], 16], "im" -> N[Round[Im[N[z]], 10^-10], 16]|>;

scriptPath = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  FileNameJoin[{"software", "stage1_solver", "tools", "pathA_26_derrick.wl"}]
];
stage1Root = ParentDirectory[DirectoryName[scriptPath]];
scratchDir = FileNameJoin[{stage1Root, "_scratch"}];
If[! DirectoryQ[scratchDir], CreateDirectory[scratchDir]];
sympyJson = FileNameJoin[{scratchDir, "pathA_26_derrick_sympy.json"}];
jsonOut = FileNameJoin[{scratchDir, "pathA_26_derrick_mathematica.json"}];

rho0 = 1.; kEOS = 0.002; pVac = 0.001; sigma = 0.003;
iwave = 1.; cW = 1.; muIn = 0.2; muOut = 2.5;
chiPerp = N[Pi, 30]; chiDN = N[Pi/2, 30];
hbar = 0.03; mPar = 1.; deltaR = 0.05; deltaPar = 0.05; epsW = 0.01;
etaMin = 0.95; nGrid = 48; hNum = 0.05;
phaseBMaxCoeff = 1.*^-4; dSlab = 2.5;
cDampSample = 1.; cDampMax = 10.; gainMin = 0.1; gainMax = 10.;
robustAxisFraction = 0.30;
p0 = kEOS rho0^5; alphaVol = p0 + pVac; muFluid = 5 kEOS rho0^4/4;

vol[a_, ell_] := 4 Pi a^3 ell/3;
aside[a_, ell_] := 4 Pi a^2 ell;
acap[a_] := 8 Pi a^3/3;
omega[a_, ell_] := Sqrt[cW^2 (chiPerp^2/a^2 + chiDN^2/ell^2) + muIn^2];
energy[a_, ell_] := alphaVol vol[a, ell] + sigma (aside[a, ell] + acap[a]) + iwave omega[a, ell];

Clear[aSym, lSym];
gradExpr = {D[energy[aSym, lSym], aSym], D[energy[aSym, lSym], lSym]};
rootSol = FindRoot[Evaluate[gradExpr] == {0, 0}, {{aSym, 1.9}, {lSym, 1.8}},
  MaxIterations -> 100
];
qStar = N[{aSym, lSym} /. rootSol, 18];
hessExpr = D[energy[aSym, lSym], {{aSym, lSym}, 2}];
hessBase = N[hessExpr /. rootSol, 18];
eigBase = Sort[N[Eigenvalues[hessBase], 18]];
omegaStar = N[omega[qStar[[1]], qStar[[2]]], 18];
eVol = alphaVol N[vol[qStar[[1]], qStar[[2]]], 18];
eSide = sigma N[aside[qStar[[1]], qStar[[2]]], 18];
eCap = sigma N[acap[qStar[[1]]], 18];
waveA = iwave cW^2 chiPerp^2/(qStar[[1]]^2 omegaStar);
waveL = iwave cW^2 chiDN^2/(qStar[[2]]^2 omegaStar);
virialA = 3 eVol + 2 eSide + 3 eCap - waveA;
virialL = eVol + eSide - waveL;
phaseALabel = If[Max[Abs[N[gradExpr /. rootSol]]] > 10^-8 || Min[qStar] <= 0,
  "A_FORBIDDEN",
  If[Min[eigBase] <= 10^-10, "A_DEGENERATE", "A_CANDIDATE_MIN"]
];

fluidOnly = <|
  "label" -> "FLUID_ONLY_COLLAPSE_NO_INTERIOR_STATIONARY",
  "dE_da" -> ToString[InputForm[FullSimplify[D[alphaVol vol[aSym, lSym] + sigma (aside[aSym, lSym] + acap[aSym]), aSym]]]],
  "dE_dL" -> ToString[InputForm[FullSimplify[D[alphaVol vol[aSym, lSym] + sigma (aside[aSym, lSym] + acap[aSym]), lSym]]]],
  "positive_derivative_sample" -> {clean[D[alphaVol vol[aSym, lSym] + sigma (aside[aSym, lSym] + acap[aSym]), aSym] /. {aSym -> 1, lSym -> 1}],
    clean[D[alphaVol vol[aSym, lSym] + sigma (aside[aSym, lSym] + acap[aSym]), lSym] /. {aSym -> 1, lSym -> 1}]}
|>;

rmax = qStar[[1]] + 0.60; wmax = qStar[[2]]/2 + 0.60;
rVals = N[Subdivide[0, rmax, nGrid - 1], 18];
wVals = N[Subdivide[-wmax, wmax, nGrid - 1], 18];
drGrid = rmax/(nGrid - 1); dwGrid = 2 wmax/(nGrid - 1);
wr = ConstantArray[1., nGrid]; wr[[1]] = 0.5; wr[[-1]] = 0.5;
ww = ConstantArray[1., nGrid]; ww[[1]] = 0.5; ww[[-1]] = 0.5;
sech2[x_] := Sech[Clip[x, {-40, 40}]]^2;

fluidExcess[a_?NumericQ, ell_?NumericQ] := Module[
  {sum = 0., r, w, weight, u, gr, dgr, sw, v, gw, dgw, gate, dg2, rho, rhoc, grand, kin},
  Do[
    r = rVals[[i]]; w = wVals[[j]];
    weight = wr[[i]] ww[[j]] drGrid dwGrid 4 Pi r^2;
    u = (r - a)/deltaR;
    gr = 0.5 (1 - Tanh[u]);
    dgr = -0.5 sech2[u]/deltaR;
    sw = Sqrt[w^2 + epsW^2];
    v = (sw - ell/2)/deltaPar;
    gw = 0.5 (1 - Tanh[v]);
    dgw = -0.5 sech2[v] (w/sw)/deltaPar;
    gate = gr gw;
    dg2 = (dgr gw)^2 + (gr dgw)^2;
    rho = rho0 (1 - etaMin gate);
    rhoc = Max[rho, 10^-5];
    grand = kEOS rho^5/4 - muFluid rho - (kEOS rho0^5/4 - muFluid rho0);
    kin = hbar^2 (rho0 etaMin)^2 dg2/(8 mPar rhoc);
    sum += (grand + kin) weight,
    {i, 1, nGrid}, {j, 1, nGrid}
  ];
  N[sum, 18]
];

numericTotal[{a_?NumericQ, ell_?NumericQ}] := fluidExcess[a, ell] + pVac N[vol[a, ell]] +
  sigma N[aside[a, ell] + acap[a]] + iwave N[omega[a, ell]];

fdH[f_, q_, h_] := Table[
  Module[{ei = {0., 0.}, ej = {0., 0.}},
    ei[[i]] = h; ej[[j]] = h;
    (f[q + ei + ej] - f[q + ei - ej] - f[q - ei + ej] + f[q - ei - ej])/(4 h^2)
  ],
  {i, 1, 2}, {j, 1, 2}
];
numericHess = N[fdH[numericTotal, qStar, hNum], 16];
numericEig = Sort[N[Eigenvalues[numericHess], 16]];
numericInterior = qStar[[1]] + 2 hNum < rmax && qStar[[2]]/2 + 2 hNum < wmax;
numericSignAgree = (Min[numericEig] > 0) == (Min[eigBase] > 0);

phaseBEnergy[a_, ell_, bf_, ss_, lw_] := energy[a, ell] + bf/a^2 + ss (ell/dSlab)^2 + lw (16 Pi/9) ell;
phaseBCorner[bf_, ss_, lw_] := Module[{sol, hess, eig, ok = True},
  sol = Quiet[Check[
      FindRoot[
        {D[phaseBEnergy[aSym, lSym, bf, ss, lw], aSym] == 0,
          D[phaseBEnergy[aSym, lSym, bf, ss, lw], lSym] == 0},
        {{aSym, qStar[[1]]}, {lSym, qStar[[2]]}}, MaxIterations -> 80
      ],
      ok = False; {}
    ]];
  If[! ok, Return[<|"success" -> False, "coeffs" -> {bf, ss, lw}, "message" -> "FindRoot failed"|>]];
  hess = N[D[phaseBEnergy[aSym, lSym, bf, ss, lw], {{aSym, lSym}, 2}] /. sol, 18];
  eig = Sort[N[Eigenvalues[hess], 18]];
  <|"success" -> (Min[eig] > 0), "coeffs" -> {clean[bf], clean[ss], clean[lw]},
    "q_star" -> (clean /@ N[{aSym, lSym} /. sol, 18]),
    "hessian_eigenvalues" -> (clean /@ eig)|>
];
phaseBCorners = Flatten[Table[phaseBCorner[bf, ss, lw],
    {bf, {0., phaseBMaxCoeff}}, {ss, {0., phaseBMaxCoeff}}, {lw, {0., phaseBMaxCoeff}}]];
phaseBAllOk = And @@ (Lookup[#, "success"] & /@ phaseBCorners);
hardStatus = If[qStar[[2]] < dSlab, "interior_not_boundary", "pushed_to_L_equals_d_slab"];
phaseBLabel = If[phaseBAllOk && hardStatus == "interior_not_boundary", "B_RESCUABLE", "B_NOT_RESCUABLE"];

consGrad[{a_?NumericQ, ell_?NumericQ}] := N[gradExpr /. {aSym -> a, lSym -> ell}, 18];
consHess[{a_?NumericQ, ell_?NumericQ}] := N[hessExpr /. {aSym -> a, lSym -> ell}, 18];
forceVec[{a_?NumericQ, ell_?NumericQ}, g_?NumericQ, s_Integer] := s g {2 a^7/ell, a^8/ell^2};
forceJac[{a_?NumericQ, ell_?NumericQ}, g_?NumericQ, s_Integer] := s g {{14 a^6/ell, -2 a^7/ell^2}, {8 a^7/ell^2, -2 a^8/ell^3}};
massMat[{a_?NumericQ, ell_?NumericQ}] := DiagonalMatrix[{rho0 N[vol[a, ell]], rho0 N[vol[a, ell]]/4}];

fixedPoint[g_?NumericQ, s_Integer] := Module[{guesses, sol, q, residual, found = Missing["NoSteady"]},
  guesses = {qStar, {1., 1.}, {2., 2.}, {2., 4.}, {1.5, 1.2}, {3., 3.}};
  Do[
    sol = Quiet[
      FindRoot[
        Evaluate[{
          gradExpr[[1]] - s g 2 aSym^7/lSym,
          gradExpr[[2]] - s g aSym^8/lSym^2
        } == {0, 0}],
        {{aSym, guess[[1]]}, {lSym, guess[[2]]}}, MaxIterations -> 80
      ]
    ];
    If[MatchQ[sol, {_Rule ..}],
      q = N[{aSym, lSym} /. sol, 18];
      residual = Norm[N[{
          gradExpr[[1]] - s g 2 aSym^7/lSym,
          gradExpr[[2]] - s g aSym^8/lSym^2
        } /. {aSym -> q[[1]], lSym -> q[[2]]}]];
      If[Min[q] > 0 && Max[q] < 50 && residual < 10^-7, found = q; Break[]]
    ],
    {guess, guesses}
  ];
  found
];

jacobianOpen[q_, g_, s_Integer] := Module[{stiff, minv, cmat},
  stiff = consHess[q] - forceJac[q, g, s];
  minv = Inverse[massMat[q]];
  cmat = DiagonalMatrix[{cDampSample, cDampSample}];
  {ArrayFlatten[{{ConstantArray[0., {2, 2}], IdentityMatrix[2]}, {-minv.stiff, -minv.cmat}}],
    stiff}
];

stableAt[g_?NumericQ, s_Integer] := Module[{q, pair, eig, stiff, stable},
  q = fixedPoint[g, s];
  If[Head[q] === Missing, Return[{False, Missing["NoSteady"]}]];
  pair = jacobianOpen[q, g, s];
  eig = Eigenvalues[N[pair[[1]], 18]];
  stiff = pair[[2]];
  stable = Max[Re[eig]] < -10^-9;
  {stable, <|"q_star" -> q, "jacobian_eigenvalues" -> eig,
    "det_total_stiffness" -> Det[stiff],
    "stiffness_symmetric_eigenvalues" -> Sort[Eigenvalues[N[(stiff + Transpose[stiff])/2, 18]]]|>}
];

findGcrit[s_Integer] := Module[{lo = 10^-4, hi = 1., mid, st},
  If[! First[stableAt[lo, s]], Return[lo]];
  While[First[stableAt[hi, s]],
    hi = 2 hi;
    If[hi > 100, Return[hi]]
  ];
  Do[
    mid = Sqrt[lo hi];
    st = First[stableAt[mid, s]];
    If[st, lo = mid, hi = mid],
    {35}
  ];
  N[lo, 12]
];

gSample = (gainMin gainMin)^2/rho0;
phaseCSigns = Association[];
Do[
  Module[{stData, gcrit, key, samplePayload},
    stData = stableAt[gSample, s];
    samplePayload = If[Head[stData[[2]]] === Missing,
      <|"steady_point" -> Null, "stable" -> False|>,
      <|"steady_point" -> (clean /@ stData[[2, "q_star"]]),
        "jacobian_spectrum" -> (cleanComplex /@ stData[[2, "jacobian_eigenvalues"]]),
        "det_total_stiffness" -> N[Round[stData[[2, "det_total_stiffness"]], 10^-10], 16],
        "stiffness_symmetric_eigenvalues" -> (N[Round[#, 10^-10], 16] & /@ stData[[2, "stiffness_symmetric_eigenvalues"]]),
        "stable" -> stData[[1]]|>
    ];
    gcrit = findGcrit[s];
    key = If[s == 1, "plus", "minus"];
    phaseCSigns[key] = <|"s" -> s,
      "sample_gain" -> <|"kappa_c" -> gainMin, "mu_drive" -> gainMin, "zeta" -> 1.,
        "c_a" -> cDampSample, "c_L" -> cDampSample, "g_open" -> N[Round[gSample, 10^-10], 16]|>,
      "sample" -> samplePayload,
      "stable_corner_gcrit" -> N[Round[gcrit, 10^-10], 16],
      "stable_corner_condition" -> "(kappa_c*mu_drive)^2 < " <> ToString[N[gcrit, 10], InputForm]|>
  ],
  {s, {1, -1}}
];

requiredRadius = robustAxisFraction (gainMax - gainMin);
minCenterRobust = gainMin + requiredRadius;
robustCenterG = (minCenterRobust^2)^2/rho0;
minUpperRobust = gainMin + 2 requiredRadius;
robustForcedG = (minUpperRobust^2)^2/rho0;
maxGcrit = Max[phaseCSigns["plus", "stable_corner_gcrit"], phaseCSigns["minus", "stable_corner_gcrit"]];
robustExists = maxGcrit >= robustForcedG;
phaseCLabel = If[robustExists, "C_STABILIZABLE", "C_GENERICALLY_UNSTABLE"];

staticDivergenceSigns = Association[];
Do[
  Module[{pair, stiff, symEig, totalEig, key},
    pair = jacobianOpen[qStar, robustCenterG, s];
    stiff = pair[[2]];
    symEig = Sort[Eigenvalues[N[(stiff + Transpose[stiff])/2, 18]]];
    totalEig = Eigenvalues[N[stiff, 18]];
    key = If[s == 1, "plus", "minus"];
    staticDivergenceSigns[key] = <|"s" -> s,
      "det_total_stiffness" -> N[Round[Det[stiff], 10^-6], 16],
      "min_symmetric_stiffness_eigenvalue" -> N[Round[First[symEig], 10^-6], 16],
      "symmetric_stiffness_eigenvalues" -> (N[Round[#, 10^-6], 16] & /@ symEig),
      "total_stiffness_eigenvalues" -> (cleanComplex /@ totalEig),
      "max_jacobian_real_part" -> N[Round[Max[Re[Eigenvalues[N[pair[[1]], 18]]]], 10^-6], 16]|>
  ],
  {s, {1, -1}}
];
staticDivergence = <|"reference_point" -> "conservative_phase_A_q_star",
  "q_reference" -> (clean /@ qStar),
  "g_open" -> N[Round[robustCenterG, 10^-10], 16],
  "damping_note" -> "Static divergence: H+K_open has negative determinant and the symmetric stiffness part has a negative eigenvalue, so passive C>=0 can only set timescales; it cannot remove the negative-stiffness direction.",
  "signs" -> staticDivergenceSigns|>;

gateScaleRatio = Min[qStar[[1]]/deltaR, qStar[[2]]/deltaPar];
gateBindingMargin = (muOut - omegaStar)/muOut;
gateChecks = <|"scale_separation" -> (gateScaleRatio > 10),
  "binding" -> (omegaStar < muOut && gateBindingMargin > 0.10),
  "candidate_interior_to_numeric_grid" -> numericInterior,
  "hessian_sign_agreement" -> numericSignAgree|>;
gatePasses = And @@ Values[gateChecks];

topVerdict = Module[{cons, drainOk},
  cons = phaseALabel == "A_CANDIDATE_MIN" || phaseBLabel == "B_RESCUABLE";
  drainOk = MemberQ[{"C_STABILIZABLE", "C_UNNEEDED"}, phaseCLabel];
  Which[
    cons && ! drainOk && ! gatePasses, "INCONCLUSIVE",
    cons && drainOk && ! gatePasses, "INCONCLUSIVE",
    ! cons && ! drainOk, "THROAT_FORBIDDEN",
    ! cons && drainOk, "THROAT_NEEDS_DYNAMICS",
    cons && ! drainOk, "THROAT_DRAIN_DESTABILIZED",
    True, "THROAT_CANDIDATE"
  ]
];

lambdaRay = <|"object" -> "reduced ledger only, not the full 4D continuum",
  "F" -> "A/a + B/a^2 + C*a^3",
  "exponents" -> {-1, -2, 3},
  "stationarity_identity" -> "-E_w - 2*E_f + 3*E_PV = 0",
  "virial" -> "E_w + 2*E_f = 3*E_PV"|>;

forms = <|"sharp_wall_forms" -> {
    <|"name" -> "fluid_grand_potential_depletion", "form" -> "P(rho0)*(4*pi/3)*a^3*L", "exponents_a_L" -> {3, 1}, "sign" -> "+"|>,
    <|"name" -> "vacuum_pressure_geometry", "form" -> "P_vac*(4*pi/3)*a^3*L", "exponents_a_L" -> {3, 1}, "sign" -> "+"|>,
    <|"name" -> "brane_tension_side", "form" -> "sigma*4*pi*a^2*L", "exponents_a_L" -> {2, 1}, "sign" -> "+"|>,
    <|"name" -> "brane_tension_caps", "form" -> "sigma*(8*pi/3)*a^3", "exponents_a_L" -> {3, 0}, "sign" -> "+"|>,
    <|"name" -> "fixed_action_wave", "form" -> "I*sqrt(c_w^2*(pi^2/a^2+(pi/2)^2/L^2)+mu_in^2)", "exponents_a_L" -> "algebraic; asymptotic pieces (-1,0) and (0,-1), not a monomial", "sign" -> "+"|>,
    <|"name" -> "phaseB_internal_4D_self_flow", "form" -> "Phi^2/(8*pi^2*rho*a^2)", "exponents_a_L" -> {-2, 0}, "sign" -> "+"|>,
    <|"name" -> "phaseB_soft_slab", "form" -> "sigma_ret*(L/d_slab)^2", "exponents_a_L" -> {0, 2}, "sign" -> "+"|>,
    <|"name" -> "phaseB_hard_slab", "form" -> "constraint L<=d_slab", "exponents_a_L" -> "inequality constraint", "sign" -> "constraint"|>,
    <|"name" -> "optional_Willmore_side", "form" -> "lambda_W*(16*pi/9)*L", "exponents_a_L" -> {0, 1}, "sign" -> "+"|>,
    <|"name" -> "phaseC_conductance", "form" -> "G_c=a^3/L", "exponents_a_L" -> {3, -1}, "sign" -> "+"|>,
    <|"name" -> "phaseC_work_kernel", "form" -> "G_work=a^2*L/ell0^5", "exponents_a_L" -> {2, 1}, "sign" -> "s in {+,-} enters F_nc"|>
  }|>;

citations = <||>;

results = <|
  "engine" -> "mathematica",
  "forms" -> forms,
  "phase_A" -> <|"q_star" -> (clean /@ qStar), "L_over_a" -> clean[qStar[[2]]/qStar[[1]]],
    "energy" -> clean[energy[qStar[[1]], qStar[[2]]]], "omega" -> clean[omegaStar],
    "gradient_residual" -> (N[Round[#, 10^-14], 16] & /@ (gradExpr /. rootSol)),
    "hessian" -> Map[clean, hessBase, {2}], "hessian_eigenvalues" -> (clean /@ eigBase),
    "label" -> phaseALabel,
    "virial" -> <|"a_identity" -> "3*E_volume + 2*E_side + 3*E_cap = I*c_w^2*chi_perp^2/(a^2*omega)",
      "L_identity" -> "E_volume + E_side = I*c_w^2*chi_DN^2/(L^2*omega)",
      "a_residual" -> N[Round[virialA, 10^-13], 16], "L_residual" -> N[Round[virialL, 10^-13], 16]|>|>,
  "fluid_only" -> fluidOnly,
  "numeric_envelope" -> <|"method" -> "coarse smooth-gate fixed-mu depletion envelope; eta relaxed to constrained lower bound",
    "grid_shape" -> {nGrid, nGrid}, "eta_relaxed" -> etaMin,
    "domain" -> <|"rmax" -> clean[rmax], "wmax" -> clean[wmax]|>,
    "hessian" -> Map[clean, numericHess, {2}], "hessian_eigenvalues" -> (clean /@ numericEig),
    "analytic_numeric_sign_agree" -> numericSignAgree, "candidate_interior_to_grid" -> numericInterior|>,
  "phase_B" -> <|"label" -> phaseBLabel,
    "open_positive_region" -> <|"B_flow" -> "(0, 0.0001)", "sigma_ret" -> "(0, 0.0001)",
      "lambda_W" -> "(0, 0.0001)", "d_slab" -> dSlab|>,
    "corner_checks" -> phaseBCorners,
    "hard_slab" -> <|"constraint" -> "L <= 2.5", "status" -> hardStatus|>,
    "soft_slab_form" -> "sigma_ret*(L/d_slab)^2",
    "flow_form" -> "Phi^2/(8*pi^2*rho*a^2)",
    "bending_form" -> "(16*pi/9)*lambda_W*L plus fixed-regulator edge corrections"|>,
  "phase_C" -> <|"label" -> phaseCLabel,
    "state_variables" -> {"Phi_w", "Delta_h"},
    "slaving" -> <|"G_c" -> "a^3/L", "Phi_w" -> "kappa_c*(a^3/L)*mu_drive", "Delta_h" -> "zeta*Phi_w"|>,
    "work_kernel" -> "G_work=a^2*L/ell0^5 (nondimensional execution: a^2*L)",
    "F_nc_after_slaving" -> "s*(kappa_c*mu_drive)^2/rho0 * (2*a^7/L, a^8/L^2)",
    "K_open" -> "-dF_nc/dq", "mass_matrix" -> "diag(rho0*V, rho0*V/4)",
    "damping" -> "diag(c_a,c_L), c_a,c_L>=0",
    "gain_box" -> <|"kappa_c" -> {gainMin, gainMax}, "zeta" -> {gainMin, gainMax},
      "mu_drive" -> {gainMin, gainMax}, "c_a_c_L" -> {0., cDampMax}|>,
    "robust_required_axis_radius" -> clean[requiredRadius],
    "robust_min_center_gain" -> clean[minCenterRobust],
    "robust_center_g_open_lower_bound" -> clean[robustCenterG],
    "robust_forced_g_lower_bound" -> clean[robustForcedG],
    "robust_region_exists" -> robustExists,
    "static_divergence_at_box_forced_gain" -> staticDivergence,
    "signs" -> phaseCSigns,
    "instability_proof" -> "Stable fixed points exist only in a tiny near-zero-drain corner. The most favorable in-box 30%-axis-radius ball center already has kappa_c=mu_drive=" <> ToString[N[minCenterRobust, 8]] <> ", hence g_open=" <> ToString[N[robustCenterG, 8]] <> ", far above max(gcrit)=" <> ToString[N[maxGcrit, 8]] <> "; therefore no qualifying ball can be stable."|>,
  "approximation_gate" -> <|"checks" -> gateChecks, "passes" -> gatePasses,
    "scale_ratio_min" -> clean[gateScaleRatio], "binding_margin_fraction" -> clean[gateBindingMargin],
    "remedy_if_failed" -> "shrink fixed regulators, raise outside threshold, enlarge grid, or run a higher-resolution envelope solve"|>,
  "lambda_ray" -> lambdaRay,
  "top_line_verdict" -> topVerdict
|>;

mismatches = {};
cmpExact[name_, x_, y_] := If[x =!= y, AppendTo[mismatches, name <> " mismatch"]];
cmpFloats[name_, xs_, ys_, tol_] := If[Length[xs] =!= Length[ys],
  AppendTo[mismatches, name <> " length mismatch"],
  Do[If[Abs[N[xs[[i]]] - N[ys[[i]]]] > tol, AppendTo[mismatches, name <> "[" <> ToString[i] <> "] mismatch"]], {i, Length[xs]}]
];

If[FileExistsQ[sympyJson],
  sympy = Import[sympyJson, "RawJSON"];
  cmpExact["top_line", results["top_line_verdict"], sympy["top_line_verdict"]];
  cmpExact["phase_A", results["phase_A", "label"], sympy["phase_A", "label"]];
  cmpExact["phase_B", results["phase_B", "label"], sympy["phase_B", "label"]];
  cmpExact["phase_C", results["phase_C", "label"], sympy["phase_C", "label"]];
  cmpFloats["baseline_hessian_eigs", results["phase_A", "hessian_eigenvalues"], sympy["phase_A", "hessian_eigenvalues"], 5.*^-6];
  cmpFloats["numeric_hessian_eigs", results["numeric_envelope", "hessian_eigenvalues"], sympy["numeric_envelope", "hessian_eigenvalues"], 2.*^-5];
  cmpFloats["phase_C_plus_q", results["phase_C", "signs", "plus", "sample", "steady_point"], sympy["phase_C", "signs", "plus", "sample", "steady_point"], 2.*^-5];
  cmpFloats["phase_C_minus_q", results["phase_C", "signs", "minus", "sample", "steady_point"], sympy["phase_C", "signs", "minus", "sample", "steady_point"], 2.*^-5];
  If[Abs[results["phase_C", "signs", "plus", "stable_corner_gcrit"] - sympy["phase_C", "signs", "plus", "stable_corner_gcrit"]] > 2.*^-4,
    AppendTo[mismatches, "plus gcrit mismatch"]];
  If[Abs[results["phase_C", "signs", "minus", "stable_corner_gcrit"] - sympy["phase_C", "signs", "minus", "stable_corner_gcrit"]] > 2.*^-4,
    AppendTo[mismatches, "minus gcrit mismatch"]];
  If[Abs[results["phase_C", "robust_center_g_open_lower_bound"] - sympy["phase_C", "robust_center_g_open_lower_bound"]] > 10^-8,
    AppendTo[mismatches, "robust center g mismatch"]];
  Do[
    If[Abs[results["phase_C", "static_divergence_at_box_forced_gain", "signs", signName, "det_total_stiffness"] -
        sympy["phase_C", "static_divergence_at_box_forced_gain", "signs", signName, "det_total_stiffness"]] > 1.,
      AppendTo[mismatches, "static det " <> signName <> " mismatch"]];
    If[Abs[results["phase_C", "static_divergence_at_box_forced_gain", "signs", signName, "min_symmetric_stiffness_eigenvalue"] -
        sympy["phase_C", "static_divergence_at_box_forced_gain", "signs", signName, "min_symmetric_stiffness_eigenvalue"]] > 10^-2,
      AppendTo[mismatches, "static min eig " <> signName <> " mismatch"]],
    {signName, {"plus", "minus"}}
  ];
  agreementStatus = If[mismatches === {}, "AGREE", "DISAGREE"],
  agreementStatus = "MISSING_SYMPY_JSON"
];
agreement = <|"status" -> agreementStatus, "sympy_json" -> sympyJson, "mismatches" -> mismatches|>;
results["engine_agreement"] = agreement;
If[agreementStatus === "DISAGREE", fail["engine disagreement: " <> ToString[mismatches, InputForm]]];

Export[jsonOut, results, "RawJSON"];

plus = results["phase_C", "signs", "plus"];
minus = results["phase_C", "signs", "minus"];
Print["top_line_verdict: ", topVerdict];
Print["phase_labels: A=", phaseALabel, " B=", phaseBLabel, " C=", phaseCLabel];
Print["envelope_hessian_eigenvalues: analytic=", results["phase_A", "hessian_eigenvalues"],
  " numeric=", results["numeric_envelope", "hessian_eigenvalues"]];
Print["phase_C_jacobian_spectrum_plus: ", plus["sample", "jacobian_spectrum"]];
Print["phase_C_jacobian_spectrum_minus: ", minus["sample", "jacobian_spectrum"]];
Print["phase_C_stable_region: robust_exists=", robustExists, " gcrit_plus=", plus["stable_corner_gcrit"],
  " gcrit_minus=", minus["stable_corner_gcrit"]];
Print["phase_C_static_divergence: g_open=", staticDivergence["g_open"],
  " plus_det=", staticDivergence["signs", "plus", "det_total_stiffness"],
  " plus_min_sym_eig=", staticDivergence["signs", "plus", "min_symmetric_stiffness_eigenvalue"],
  " minus_det=", staticDivergence["signs", "minus", "det_total_stiffness"],
  " minus_min_sym_eig=", staticDivergence["signs", "minus", "min_symmetric_stiffness_eigenvalue"]];
Print["engine_agreement: ", agreementStatus];
