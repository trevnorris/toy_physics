(* master_checks.wl
   Run: math -script master_checks.wl
*)

ClearAll["Global`*"];
Get[FileNameJoin[{DirectoryName[$InputFileName], "definitions.wl"}]];
Paper7FreezeSummary[];

(* ============================================================
   H2 / Step-2 Support: Quick normalization checks for W(w)
   ============================================================ *)

ClearAll[CheckWwNormalization];

CheckWwNormalization[] := Module[
  {W, normInf, normCut, normDef, dInf, dCut, asmInf, asmCut},

  If[DownValues[Ww] === {}, Print["[CheckWwNormalization] ERROR: Ww[w_] is not defined."]; Return[$Failed]];
  If[DownValues[WprojCut] === {}, Print["[CheckWwNormalization] ERROR: WprojCut[] is not defined."]; Return[$Failed]];
  If[DownValues[WprojNorm] === {}, Print["[CheckWwNormalization] ERROR: WprojNorm[] is not defined."]; Return[$Failed]];

  W = Ww[w];

  asmInf = {ellW[] > 0};
  asmCut = {ellW[] > 0, wProjNSig > 0};

  Print["--- W(w) normalization checks ---"];

  normInf = Quiet @ FullSimplify[Integrate[W, {w, -Infinity, Infinity}], Assumptions -> asmInf];
  dInf    = Quiet @ FullSimplify[normInf - 1, Assumptions -> asmInf];

  Print[" ∫_{-∞}^{∞} Ww(w) dw = ", normInf];
  If[dInf === 0,
    Print[" PASS: infinite-range normalization"],
    Print[" FAIL: infinite-range normalization (Δ = ", dInf, ")"]
  ];

  normCut = Quiet @ FullSimplify[Integrate[W, {w, -WprojCut[], WprojCut[]}], Assumptions -> asmCut];
  normDef = Quiet @ FullSimplify[WprojNorm[], Assumptions -> asmCut];
  dCut    = Quiet @ FullSimplify[normCut - normDef, Assumptions -> asmCut];

  Print[" ∫_{|w|≤WprojCut} Ww(w) dw = ", normCut];
  Print[" WprojNorm[]                 = ", normDef];

  If[dCut === 0,
    Print[" PASS: truncated integral matches WprojNorm[]"],
    Print[" FAIL: truncated integral mismatch (Δ = ", dCut, ")"]
  ];

  If[DownValues[WprojTailMass] =!= {},
    Print[" Tail mass (WprojTailMass[]) = ", Quiet @ FullSimplify[WprojTailMass[], Assumptions -> asmCut]]
  ];

  Print["--- end W(w) checks ---"];
  <|"normInf" -> normInf, "normCut" -> normCut, "normDef" -> normDef, "dInf" -> dInf, "dCut" -> dCut|>
];

CheckWwNormalization[];

Print["=== Paper 7 Master Checks: Hard-Mode 4D ==="];

(* ---- 1. Setup & Coordinates ---- *)
space = {x, y, z, w};
hbar = Symbol["hbar"]; m = Symbol["m"];
K = Symbol["K"];
c = Symbol["c"];
a = Symbol["a"]; L = Symbol["L"];
(* Dynamic symbols for moving wall checks *)
aObs = Symbol["aObs"]; LObs = Symbol["LObs"];

VconfF = Vconf[x, y, z, w, a, L];
mu2F   = mu2A[x, y, z, w, a, L];

(* ---- 2. Helpers ---- *)
computeAtomicEL[Latomic_, S_, St_, Sgrad_List, map_] := Module[
  {term0, termT, termSpace, dLdS, dLdSt, dLdSgrad},
  dLdS     = D[Latomic, S];
  dLdSt    = D[Latomic, St];
  dLdSgrad = Table[D[Latomic, Sgrad[[i]]], {i, Length[Sgrad]}];

  term0 = dLdS /. map;
  termT = D[ (dLdSt /. map), t ];
  termSpace = Sum[D[(dLdSgrad[[i]] /. map), space[[i]]], {i, 1, 4}];

  Expand[term0 - termT - termSpace]
];

assertZero[name_, expr_] := Module[{r},
  r = Expand[expr];
  If[r === 0,
    Print["[PASS] ", name]; True,
    r = TimeConstrained[Simplify[r], 2.0, r];
    If[r === 0,
       Print["[PASS] ", name, " (via Simplify)"];
       True,
       Print["[FAIL] ", name];
       Print["  Residual (short): ", Short[r, 3]];
       False
    ]
  ]
];

assertSmall[name_, val_, tol_] := Module[{v = N[val]},
  If[NumericQ[v] && Abs[v] < tol,
    Print["[PASS] ", name, " (", v, ")"]; True,
    Print["[FAIL] ", name, " (", v, ")"]; False
  ]
];

(* ============================================================
   CHECK 0: Gamma Quadrature Sanity (Step 3)
   ============================================================ *)
Print["-- Section 0: Gamma Quadrature Sanity --"];
Block[{rPort = 1, lmaxCheck = 2, Anum, Aexact, absErr, relErr, ports, m, diag, normMat, maxOff},
  Anum = GammaIntegrate[Function[{th, ph}, 1]];
  Aexact = 4 Pi rPort^2;
  absErr = N[Anum - Aexact];
  relErr = N[Anum/Aexact - 1];
  Print["Gamma area: Anum = ", Anum, ", Aexact = ", Aexact];
  Print["abs error = ", absErr, ", rel error = ", relErr];
  assertSmall["Gamma area relative error", relErr, 1*^-8];

  ports = PortsList[lmaxCheck];
  m = Table[
    GammaIntegrate[
      Function[{th, ph}, PortP[keyi][th, ph] * Conjugate[PortP[keyj][th, ph]]]
    ],
    {keyi, ports}, {keyj, ports}
  ];
  diag = Diagonal[m];
  normMat = Table[
    If[i == j, 0, Abs[m[[i, j]]]/Sqrt[Abs[diag[[i]]] * Abs[diag[[j]]]]],
    {i, Length[ports]}, {j, Length[ports]}
  ];
  maxOff = Max[N[Flatten[Abs[normMat]]]];
  Print["Port orthonormal check: max off-diagonal (normalized) = ", maxOff];
  assertSmall["Port orthonormal off-diagonal", maxOff, 1*^-6];
];

(* ============================================================
   CHECK 0.5: Drive Envelope Sanity (Step 4)
   ============================================================ *)
Print["-- Section 0.5: Drive Envelope Sanity --"];

assertZero["DriveEnvelope is even in w",
  Simplify[DriveEnvelope[x,y,z,w] - DriveEnvelope[x,y,z,-w]]
];

Block[{rPort = 1.0, wDrive = 0.4, rDrive = 0.2, epsAng = 10^-9},
  Module[{center, offW, offR},
    center = N @ DriveEnvelope[rPort, 0, 0, 0];
    assertSmall["DriveEnvelope center equals 1 (|1-center|)",
      Abs[center - 1], 10^-12
    ];

    offW = N @ DriveEnvelope[rPort, 0, 0, 5 wDrive];
    assertSmall["DriveEnvelope decays for |w| >> wDrive", offW, 10^-8];

    offR = N @ DriveEnvelope[rPort + 5 rDrive, 0, 0, 0];
    assertSmall["DriveEnvelope decays for |R3-rPort| >> rDrive", offR, 10^-8];
  ]
];

(* ============================================================
   CHECK 0.75: Zeff extraction sanity (Step 5)
   ============================================================ *)
Print["-- Section 0.75: Zeff Extraction Sanity --"];

Block[
  {nPorts = 5, omega = 3.0, nDiscard = 2, nMeasure = 6,
   Tper, t0, t1, nSamp, tList, Ztrue, runs, uTS, jTS, k, res, Zhat, uMat, jMat},

  Tper = 2 Pi/omega;
  t0 = 0.0;
  t1 = (nDiscard + nMeasure) * Tper;
  nSamp = 4000;
  tList = N @ Subdivide[t0, t1, nSamp - 1];

  Ztrue = Table[If[i == j, 1 + 0.2 I, 0.05 (i - j)], {i, nPorts}, {j, nPorts}];

  runs = Table[
    uTS = ConstantArray[0.0, {nPorts, Length[tList]}];
    uTS[[k]] = Cos[omega * tList];
    jTS = Ztrue . uTS;
    <|"t" -> tList, "uPorts" -> uTS, "jPorts" -> jTS|>,
    {k, 1, nPorts}
  ];

  {uMat, jMat} = BuildPortAmplitudeMatrices[runs, omega, nDiscard, nMeasure];

  res = EstimateZeffRobust[uMat, jMat, "Verbose" -> False];
  Zhat = res["Z"];

  Print["Zeff residual (rel) = ", N[res["resid"]]];
  assertSmall["Zeff residual small", res["resid"], 1*^-10];
  assertSmall["Zeff matches Ztrue (Frobenius rel err)",
    Norm[Zhat - Ztrue]/Norm[Ztrue], 1*^-10
  ];
];

(* ============================================================
   CHECK 1: GNLS EOM
   ============================================================ *)
Print["-- Section 1: GNLS EOM --"];
psiF  = psi[x,y,z,w,t];
psibF = psib[x,y,z,w,t];

atomMap = {
  P   -> psiF,  Pt  -> D[psiF, t],  Px  -> D[psiF, x],  Py  -> D[psiF, y],  Pz  -> D[psiF, z],  Pw  -> D[psiF, w],
  Pb  -> psibF, Pbt -> D[psibF, t], Pbx -> D[psibF, x], Pby -> D[psibF, y], Pbz -> D[psibF, z], Pbw -> D[psibF, w],
  Rho -> psibF * psiF,
  Vc  -> VconfF
};

LgnlsAtomic =
  (I hbar/2) (Pb * Pt - P * Pbt) -
  (hbar^2/(2 m)) (Px*Pbx + Py*Pby + Pz*Pbz + Pw*Pbw) -
  Vc * (Pb * P) -
  (K/4) * (Pb * P)^5;

eomPsi = computeAtomicEL[LgnlsAtomic, Pb, Pbt, {Pbx, Pby, Pbz, Pbw}, atomMap];
eomPsib = computeAtomicEL[LgnlsAtomic, P, Pt, {Px, Py, Pz, Pw}, atomMap];

lapPsi  = D[psiF, {x,2}] + D[psiF, {y,2}] + D[psiF, {z,2}] + D[psiF, {w,2}];
lapPsib = D[psibF, {x,2}] + D[psibF, {y,2}] + D[psibF, {z,2}] + D[psibF, {w,2}];
rhoReal = psibF * psiF;

expectPsi  = I hbar D[psiF, t] + (hbar^2/(2 m)) lapPsi - VconfF psiF - (5 K/4) rhoReal^4 psiF;
expectPsib = -I hbar D[psibF, t] + (hbar^2/(2 m)) lapPsib - VconfF psibF - (5 K/4) rhoReal^4 psibF;

assertZero["GNLS EOM for psi", eomPsi - expectPsi];
assertZero["GNLS EOM for psib", eomPsib - expectPsib];

(* ============================================================
   CHECK 2: Continuity
   ============================================================ *)
Print["-- Section 2: Continuity --"];
rhsPsi  = - (hbar^2/(2 m)) lapPsi + VconfF psiF + (5 K/4) rhoReal^4 psiF;
rhsPsib =   (hbar^2/(2 m)) lapPsib - VconfF psibF - (5 K/4) rhoReal^4 psibF;

dtRho = ( (1/(I hbar))*rhsPsi * psibF ) + ( psiF * (1/(I hbar))*rhsPsib );
divJ = (hbar/(2 I m)) * ( psibF * lapPsi - psiF * lapPsib );

assertZero["Continuity Equation", dtRho + divJ];

(* ============================================================
   CHECK 2.5: Dynamic Continuity (NEW)
   ============================================================ *)
Print["-- Section 2.5: Dynamic Continuity (Time-Dependent V) --"];
(* Ensure continuity holds even if V depends on t *)
VconfDyn25 = Vconf[x, y, z, w, aObs[t], LObs[t]];

rhsPsiDyn25  = - (hbar^2/(2 m)) lapPsi + VconfDyn25 * psiF + (5 K/4) rhoReal^4 * psiF;
rhsPsibDyn25 =   (hbar^2/(2 m)) lapPsib - VconfDyn25 * psibF - (5 K/4) rhoReal^4 * psibF;

dtRhoDyn = ( (1/(I hbar))*rhsPsiDyn25 * psibF ) + ( psiF * (1/(I hbar))*rhsPsibDyn25 );
(* divJ is same kinematic form *)
assertZero["Dynamic Continuity (dtRho + divJ with V(t))", dtRhoDyn + divJ];


(* ============================================================
   CHECK 3: Breathing Forces (Static)
   ============================================================ *)
Print["-- Section 3: Breathing Forces (Static) --"];
HfluidAtomic = (hbar^2/(2 m)) (Px*Pbx + Py*Pby + Pz*Pbz + Pw*Pbw) + Vc*(Pb*P) + (K/4)*(Pb*P)^5;

mapA = {
   Vc -> VconfF,
   P -> psiF, Pb -> psibF,
   Px->D[psiF,x], Pbx->D[psibF,x], Py->D[psiF,y], Pby->D[psibF,y],
   Pz->D[psiF,z], Pbz->D[psibF,z], Pw->D[psiF,w], Pbw->D[psibF,w]
};

negdHda = -D[ (HfluidAtomic /. mapA), a ];
expectedForceA = - (psibF * psiF) * D[VconfF, a];

negdHdL = -D[ (HfluidAtomic /. mapA), L ];
expectedForceL = - (psibF * psiF) * D[VconfF, L];

assertZero["Fluid Force (a)", negdHda - expectedForceA];
assertZero["Fluid Force (L)", negdHdL - expectedForceL];

(* ============================================================
   CHECK 4: Wave EOM (Strict Sign Fixed)
   ============================================================ *)
Print["-- Section 4: Wave EOM (Strict) --"];
AF = A[x,y,z,w,t];
waveMap = {
  ASym -> AF, At -> D[AF, t],
  Ax -> D[AF, x], Ay -> D[AF, y], Az -> D[AF, z], Aw -> D[AF, w],
  Mu2 -> mu2F
};

LwaveAtomic = 1/2 * (At^2 - c^2*(Ax^2 + Ay^2 + Az^2 + Aw^2) - Mu2 * ASym^2);
eomWave = computeAtomicEL[LwaveAtomic, ASym, At, {Ax, Ay, Az, Aw}, waveMap];

lapA = D[AF, {x,2}] + D[AF, {y,2}] + D[AF, {z,2}] + D[AF, {w,2}];
standardWaveEq = D[AF, {t,2}] - c^2 lapA + mu2F AF;

(* EL returns -(Eq), so eom + Eq == 0 *)
assertZero["Wave EOM matches standard form", eomWave + standardWaveEq];

(* ============================================================
   CHECK 5: Wave Forces
   ============================================================ *)
Print["-- Section 5: Wave Forces --"];
HwaveAtomic = 1/2 * (At^2 + c^2*(Ax^2 + Ay^2 + Az^2 + Aw^2) + Mu2 * ASym^2);

negdHdaWave = -D[ (HwaveAtomic /. waveMap), a ];
expectForceWaveA = - (AF^2 / 2) * D[mu2F, a];

negdHdLWave = -D[ (HwaveAtomic /. waveMap), L ];
expectForceWaveL = - (AF^2 / 2) * D[mu2F, L];

assertZero["Wave Force (a)", negdHdaWave - expectForceWaveA];
assertZero["Wave Force (L)", negdHdLWave - expectForceWaveL];

(* ============================================================
   CHECK 5.5: Combined Forces (Robust Regression Guard)
   ============================================================ *)
Print["-- Section 5.5: Combined Total Force (Robust) --"];
(* Calculate H_tot first, THEN differentiate. This ensures no summation errors. *)

HtotAtomic = HfluidAtomic + HwaveAtomic;
fullMap = Join[mapA, waveMap];

(* Total Forces from Hamiltonian Derivative *)
ForceTotA = -D[ HtotAtomic /. fullMap, a ];
ForceTotL = -D[ HtotAtomic /. fullMap, L ];

(* Expected Sums *)
ExpectedTotA = expectedForceA + expectForceWaveA;
ExpectedTotL = expectedForceL + expectForceWaveL;

assertZero["Total Force A (Calculated from H_tot)", ForceTotA - ExpectedTotA];
assertZero["Total Force L (Calculated from H_tot)", ForceTotL - ExpectedTotL];


(* ============================================================
   CHECK 6: Thermo & Sound Speed (Robust)
   ============================================================ *)
Print["-- Section 6: Thermodynamics & Sound Speed --"];
rhoSym = Symbol["rho"];
Udens[r_] = (K/4) * r^5;
enthalpy[r_] = D[Udens[r], r];
Pressure[r_] = Simplify[r * enthalpy[r] - Udens[r]];

assertZero["Pressure Consistency", Pressure[rhoSym] - K * rhoSym^5];

csSqThermo = D[Pressure[rhoSym], rhoSym];
csSqHydro  = rhoSym * D[enthalpy[rhoSym], rhoSym];

assertZero["Sound Speed Consistency (dP/drho == rho h')", csSqThermo - csSqHydro];

(* ============================================================
   CHECK 7: Moving Wall Energy Balance (True Dynamic Check)
   ============================================================ *)
Print["-- Section 7: Fluid Moving Wall Work Balance --"];

(* 1. Define Fluxes (Symbolic form) *)
FluxX = -(hbar^2/(2 m)) (Pbt * Px + Pt * Pbx);
FluxY = -(hbar^2/(2 m)) (Pbt * Py + Pt * Pby);
FluxZ = -(hbar^2/(2 m)) (Pbt * Pz + Pt * Pbz);
FluxW = -(hbar^2/(2 m)) (Pbt * Pw + Pt * Pbw);

(* 2. Define Time-Dependent Geometry *)
VconfDyn = Vconf[x, y, z, w, aObs[t], LObs[t]];

(* 3. Update EOMs to see Dynamic Potential *)
rhsPsiDyn  = - (hbar^2/(2 m)) lapPsi  + VconfDyn * psiF  + (5 K/4) rhoReal^4 * psiF;
rhsPsibDyn =   (hbar^2/(2 m)) lapPsib - VconfDyn * psibF - (5 K/4) rhoReal^4 * psibF;

(* Base temporal rules *)
baseRulesEOM = {
  D[psiF, t]  -> (1/(I hbar)) rhsPsiDyn,
  D[psibF, t] -> (1/(I hbar)) rhsPsibDyn
};

(* Hardened rules: Explicitly add mixed derivatives *)
(* This prevents failures if Mathematica expands D[H,t] into D[psi,x,t] *)
mixedRules = Flatten@Table[{
   D[psiF, t, var]  -> D[(1/(I hbar)) rhsPsiDyn, var],
   D[psiF, var, t]  -> D[(1/(I hbar)) rhsPsiDyn, var],
   D[psibF, t, var] -> D[(1/(I hbar)) rhsPsibDyn, var],
   D[psibF, var, t] -> D[(1/(I hbar)) rhsPsibDyn, var]
}, {var, space}];

rulesEOMDyn = Join[baseRulesEOM, mixedRules];

(* 4. Inject Dynamic V into Hamiltonian and Fluxes *)
HrealDyn = HfluidAtomic /. atomMap /. VconfF -> VconfDyn;
Sx = FluxX /. atomMap;
Sy = FluxY /. atomMap;
Sz = FluxZ /. atomMap;
Sw = FluxW /. atomMap;

divS = D[Sx, x] + D[Sy, y] + D[Sz, z] + D[Sw, w];
WorkTermDyn = rhoReal * D[VconfDyn, t];

balanceDyn = D[HrealDyn, t] + divS - WorkTermDyn;
finalBalance = balanceDyn /. rulesEOMDyn /. rulesEOMDyn;

assertZero["Moving Wall Balance (dH/dt + divS = rho dV/dt)", finalBalance];

(* ============================================================
   CHECK 7.5: Wave Moving Wall Work Balance
   ============================================================ *)
Print["-- Section 7.5: Wave Moving Wall Work Balance --"];
(* Canonical Wave Fluxes: S_i = - c^2 * A_t * A_i *)
SwaveX = -c^2 * At * Ax;
SwaveY = -c^2 * At * Ay;
SwaveZ = -c^2 * At * Az;
SwaveW = -c^2 * At * Aw;

(* Dynamic Geometry for Wave *)
mu2Dyn = mu2A[x, y, z, w, aObs[t], LObs[t]];

(* Update Wave EOM for dynamic mu2 *)
rhsWaveDyn = c^2 * lapA - mu2Dyn * AF;
baseWaveRule = { D[AF, {t, 2}] -> rhsWaveDyn };

(* Hardened rules for wave mixed derivatives *)
mixedWaveRules = Flatten@Table[{
   D[AF, {t, 2}, var] -> D[rhsWaveDyn, var],
   D[AF, var, {t, 2}] -> D[rhsWaveDyn, var]
}, {var, space}];
rulesWaveEOM = Join[baseWaveRule, mixedWaveRules];

(* Map to dynamic H *)
HwaveDyn = HwaveAtomic /. waveMap /. mu2F -> mu2Dyn;

(* Fluxes *)
SxW = SwaveX /. waveMap;
SyW = SwaveY /. waveMap;
SzW = SwaveZ /. waveMap;
SwW = SwaveW /. waveMap;

divSW = D[SxW, x] + D[SyW, y] + D[SzW, z] + D[SwW, w];

(* Expected Work Term: 1/2 A^2 * d(mu2)/dt *)
WorkTermWave = (1/2) * AF^2 * D[mu2Dyn, t];

(* Balance *)
balanceWave = D[HwaveDyn, t] + divSW - WorkTermWave;
finalBalanceWave = balanceWave /. rulesWaveEOM /. rulesWaveEOM;

assertZero["Wave Moving Wall Balance (dH/dt + divS = 1/2 A^2 dmu^2/dt)", finalBalanceWave];


(* ============================================================
   CHECK 8: Madelung / Hydro Consistency
   ============================================================ *)
Print["-- Section 8: Madelung / Hydro Consistency --"];
NonlinearTerm = (5 K/4) * rhoReal^4 * psiF;
EffectivePotential = NonlinearTerm / psiF;
ExpectedEnthalpy = enthalpy[rhoSym] /. rhoSym -> rhoReal;

assertZero["Madelung Potential Matches Enthalpy", EffectivePotential - ExpectedEnthalpy];

(* ============================================================
   CHECK 9: Dynamic Geometry Chain Rule (Explicit)
   ============================================================ *)
Print["-- Section 9: Dynamic Geometry Chain Rule (Explicit Check) --"];

aSym = Symbol["aSym"];
LSym = Symbol["LSym"];
Vabstract = Vconf[x, y, z, w, aSym, LSym];

dVda = D[Vabstract, aSym];
dVdL = D[Vabstract, LSym];
subs = {aSym -> aObs[t], LSym -> LObs[t]};

rhs = (dVda /. subs) * D[aObs[t], t] + (dVdL /. subs) * D[LObs[t], t];
lhs = D[Vconf[x, y, z, w, aObs[t], LObs[t]], t];

assertZero["Chain Rule Identity", lhs - rhs];

(* ============================================================
   CHECK 10: Explicit Potential Localization (Refined)
   ============================================================ *)
Print["-- Section 10: Explicit Potential Localization --"];

R = Symbol["R"];
V0 = Symbol["V0"];
delta = Symbol["delta"];
uSym = Symbol["uSym"];

VTestRadial = V0 * (1/2)*(1 + Tanh[(R - a)/delta]);
ForceDensityRadial = -D[VTestRadial, a];

(* 1. Check Far Field (R -> Infinity) *)
assump = {V0 > 0, delta > 0, a > 0};
limFar = Simplify[Limit[ForceDensityRadial, R -> Infinity], Assumptions -> assump];
assertZero["Force density vanishes at infinity", limFar];

(* 2. Check Deep Interior using substitution (R-a)/delta -> uSym *)
ForceDensityU = ForceDensityRadial /. (R - a)/delta -> uSym;
limInsideU = Limit[ForceDensityU, uSym -> -Infinity];
assertZero["Force density vanishes deep inside (u -> -inf)", limInsideU];

(* 3. Check Peak Location *)
valPeak = Simplify[ForceDensityRadial /. R -> a];
assertZero["Force density peaks at wall (R=a)", valPeak - (V0/(2 delta))];

Print[""];
Print["=== Master Checks Complete ==="];

(*"
Output:

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
--- W(w) normalization checks ---
 ∫_{-∞}^{∞} Ww(w) dw = 1
 PASS: infinite-range normalization
 ∫_{|w|≤WprojCut} Ww(w) dw = Erf[wProjNSig]
 WprojNorm[]                 = Erf[wProjNSig]
 PASS: truncated integral matches WprojNorm[]
 Tail mass (WprojTailMass[]) = Erfc[wProjNSig]
--- end W(w) checks ---
=== Paper 7 Master Checks: Hard-Mode 4D ===
-- Section 0: Gamma Quadrature Sanity --
Gamma area: Anum = 12.566370551478308, Aexact = 4*Pi
abs error = -6.28808649594248*^-8, rel error = -5.003900183098153*^-9
[PASS] Gamma area relative error (-5.003900183098153*^-9)
Port orthonormal check: max off-diagonal (normalized) = 1.097638764081165*^-8
[PASS] Port orthonormal off-diagonal (1.097638764081165*^-8)
-- Section 0.5: Drive Envelope Sanity --
[PASS] DriveEnvelope is even in w
[PASS] DriveEnvelope center equals 1 (|1-center|) (0.)
[PASS] DriveEnvelope decays for |w| >> wDrive (1.388794386496407*^-11)
[PASS] DriveEnvelope decays for |R3-rPort| >> rDrive (1.388794386496407*^-11)
-- Section 0.75: Zeff Extraction Sanity --
Zeff residual (rel) = 0.
[PASS] Zeff residual small (0.)
[PASS] Zeff matches Ztrue (Frobenius rel err) (1.214250135728164*^-16)
-- Section 1: GNLS EOM --
[PASS] GNLS EOM for psi
[PASS] GNLS EOM for psib
-- Section 2: Continuity --
[PASS] Continuity Equation
-- Section 2.5: Dynamic Continuity (Time-Dependent V) --
[PASS] Dynamic Continuity (dtRho + divJ with V(t))
-- Section 3: Breathing Forces (Static) --
[PASS] Fluid Force (a)
[PASS] Fluid Force (L)
-- Section 4: Wave EOM (Strict) --
[PASS] Wave EOM matches standard form
-- Section 5: Wave Forces --
[PASS] Wave Force (a)
[PASS] Wave Force (L)
-- Section 5.5: Combined Total Force (Robust) --
[PASS] Total Force A (Calculated from H_tot)
[PASS] Total Force L (Calculated from H_tot)
-- Section 6: Thermodynamics & Sound Speed --
[PASS] Pressure Consistency
[PASS] Sound Speed Consistency (dP/drho == rho h')
-- Section 7: Fluid Moving Wall Work Balance --
[PASS] Moving Wall Balance (dH/dt + divS = rho dV/dt)
-- Section 7.5: Wave Moving Wall Work Balance --
[PASS] Wave Moving Wall Balance (dH/dt + divS = 1/2 A^2 dmu^2/dt)
-- Section 8: Madelung / Hydro Consistency --
[PASS] Madelung Potential Matches Enthalpy
-- Section 9: Dynamic Geometry Chain Rule (Explicit Check) --
[PASS] Chain Rule Identity
-- Section 10: Explicit Potential Localization --
[PASS] Force density vanishes at infinity
[PASS] Force density vanishes deep inside (u -> -inf)
[PASS] Force density peaks at wall (R=a)

=== Master Checks Complete ===
"*)
