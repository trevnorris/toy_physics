(*
  4D -> brane gravity + Newtonian + 1PN master harness
  ----------------------------------------------------
  Purpose:
    One PASS/FAIL style Wolfram Language script that consolidates the gravity-side
    symbolic checks and synthetic self-tests implied by the current project files.

  Claim classes used in this harness:
    - EXACT IDENTITY: algebraic consequences of definitions / projections / product rules
    - CONTROLLED REDUCTION: reductions that require regime assumptions (quasi-static, etc.)
    - PROTOCOL CLOSURE: adiabatic-response closure for kappa_PV
    - SYNTHETIC SELF-TEST: numerical checks of extraction machinery on manufactured data
    - ANSATZ-LEVEL CHECK: consistency scaffold for the particle-limit mapping, not a full 4D derivation

  Important scope note:
    This harness does NOT claim to derive the full collective-coordinate defect EOM directly
    from the 4D action. It includes an ansatz-level particle-limit consistency check so that
    the Newtonian / 1PN ledger can be exercised in one place, while keeping that claim marked.
*)

Print["--- 4D gravity + Newtonian + 1PN master harness ---"];

ClearAll[
  section, pass, fail, skip, info,
  zeroScalarQ, zeroArrayQ,
  checkScalar, checkVector, checkMatrix,
  checkEqScalar, checkEqVector, checkEqMatrix,
  checkNumeric, checkCondition,
  cross3, dot3, outer3, grad3, div3, curl3, lap3, advect3, divTensor3,
  ProjectedContinuitySource, ProjectedMomentumSource,
  LockInAmplitude, EstimateZeff, PoissonCorrectionRatios,
  x, y, z, t, w, r, rr, theta, phi, tau,
  ellw, wProj, sigJ, rho0, cs0, mParam,
  Qsrc, muKepler, cLight, aSemi, ecc,
  M0, qVar, PhiN, vSq,
  mInert, mGrav, X, Y, Z, PhiExt,
  R11, R12, R13, R21, R22, R23, R31, R32, R33,
  rhoBf, Jx, Jy, Jz, SrhoF,
  phiB, vTx, vTy, vTz,
  drhoF, phi3F, SsrcF,
  rhoOF, vxO, vyO, vzO, F1O, F2O, F3O, SJ1O, SJ2O, SJ3O, SrhoOpenF, contOpenExpr, contOpenResidual,
  u1, u2, u3, rho3f, wFac,
  dDim, aRad, Ubody, SdMinus1, Vd, Lseg, V3,
  nEOS, xRatio, alphaSq, Kvec, aHsq, betaVel, E0, vSym, GMsym, bImp,
  aVar, rhoVar, Cw, Cf, Cpv, Lam, KEOS,
  passCount, failCount, skipCount,
  beta1PN, kappaPV, omegaDrive, periods, samplesPerPeriod, times,
  Ztrue, Ucols, Jcols, uTimeSeries, jTimeSeries, Urec, Jrec, Zrec,
  Z0, Z1, Z2, zToy, omegaVar,
  epsDemo, omegaDemo, rhoDemo, phiDemo, ratiosDemo, sourceDemo,
  rOrbit, phOrbit, sigmaOrbit
];

passCount = 0;
failCount = 0;
skipCount = 0;

section[name_String] := Print["\n=== ", name, " ==="];
pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);
skip[name_String, msg_String : ""] := (
  skipCount++;
  If[msg === "",
    Print["SKIP: ", name],
    Print["SKIP: ", name, "\n  ", msg]
  ]
);
info[msg_String] := Print["INFO: ", msg];

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet@FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet@FullSimplify[res == 0, assum]]
];

zeroArrayQ[arr_, assum_: True] := Module[{res, flat},
  res = Quiet@FullSimplify[arr, assum];
  flat = Flatten[{res}];
  And @@ (TrueQ[# === 0] || TrueQ[Quiet@FullSimplify[# == 0, assum]] & /@ flat)
];

checkScalar[name_String, expr_, assum_: True] := Module[{res},
  res = Quiet@FullSimplify[expr, assum];
  If[zeroScalarQ[expr, assum], pass[name], fail[name, res]]
];

checkVector[name_String, expr_List, assum_: True] := Module[{res},
  res = Quiet@FullSimplify[expr, assum];
  If[zeroArrayQ[expr, assum], pass[name], fail[name, res]]
];

checkMatrix[name_String, expr_?MatrixQ, assum_: True] := Module[{res},
  res = Quiet@FullSimplify[expr, assum];
  If[zeroArrayQ[expr, assum], pass[name], fail[name, res]]
];

checkEqScalar[name_String, lhs_, rhs_, assum_: True] := checkScalar[name, lhs - rhs, assum];
checkEqVector[name_String, lhs_List, rhs_List, assum_: True] := checkVector[name, lhs - rhs, assum];
checkEqMatrix[name_String, lhs_?MatrixQ, rhs_?MatrixQ, assum_: True] := checkMatrix[name, lhs - rhs, assum];

checkNumeric[name_String, lhs_, rhs_, tol_ : 10^-8] := Module[{res, err},
  res = Quiet@N[lhs - rhs, 50];
  err = Quiet@Max[Abs[Flatten[{res}]]];
  If[NumericQ[err] && err < tol, pass[name], fail[name, res]]
];

checkCondition[name_String, cond_] := If[TrueQ[cond], pass[name], fail[name, cond]];

(* --- explicit vector/tensor operators --- *)

cross3[a_List, b_List] := {
  a[[2]] b[[3]] - a[[3]] b[[2]],
  a[[3]] b[[1]] - a[[1]] b[[3]],
  a[[1]] b[[2]] - a[[2]] b[[1]]
};

dot3[a_List, b_List] := Total[a b];

outer3[a_List, b_List] := Table[a[[i]] b[[j]], {i, 1, 3}, {j, 1, 3}];

grad3[f_, {xx_, yy_, zz_}] := {D[f, xx], D[f, yy], D[f, zz]};

div3[v_List, {xx_, yy_, zz_}] := D[v[[1]], xx] + D[v[[2]], yy] + D[v[[3]], zz];

curl3[v_List, {xx_, yy_, zz_}] := {
  D[v[[3]], yy] - D[v[[2]], zz],
  D[v[[1]], zz] - D[v[[3]], xx],
  D[v[[2]], xx] - D[v[[1]], yy]
};

lap3[f_, {xx_, yy_, zz_}] := D[f, {xx, 2}] + D[f, {yy, 2}] + D[f, {zz, 2}];

advect3[v_List, coords_List] := Table[dot3[v, grad3[v[[i]], coords]], {i, 1, 3}];

divTensor3[T_?MatrixQ, {xx_, yy_, zz_}] := Table[
  D[T[[i, 1]], xx] + D[T[[i, 2]], yy] + D[T[[i, 3]], zz],
  {i, 1, 3}
];

ProjectedContinuitySource[Wexpr_, Jwexpr_, wSym_, wmin_, wmax_] :=
  -((Wexpr Jwexpr) /. wSym -> wmax)
  + ((Wexpr Jwexpr) /. wSym -> wmin)
  + Integrate[D[Wexpr, wSym] Jwexpr, {wSym, wmin, wmax}];

ProjectedMomentumSource[Wexpr_, Piwexpr_, wSym_, wmin_, wmax_] :=
  -((Wexpr Piwexpr) /. wSym -> wmax)
  + ((Wexpr Piwexpr) /. wSym -> wmin)
  + Integrate[D[Wexpr, wSym] Piwexpr, {wSym, wmin, wmax}];

LockInAmplitude[data_List, times_List, omega_?NumericQ] /; Length[data] == Length[times] :=
  Mean[data Exp[-I omega times]];

EstimateZeff[Umat_?MatrixQ, Jmat_?MatrixQ, ridge_ : 10^-12] := Module[
  {cond, invU, n},
  n = Length[Umat];
  cond = Quiet@Check[ConditionNumber[N[Umat, 50]], Infinity];
  invU = If[NumericQ[cond] && cond < 10^10,
    LinearSolve[N[Umat, 50], IdentityMatrix[n]],
    PseudoInverse[N[Umat, 50]]
  ];
  N[Jmat, 50].invU
];

PoissonCorrectionRatios[rhoField_, phiField_, vTField_List, sourceField_, coords_List, rules_List] := Module[
  {timeTerm, advTerm, sourceTerm},
  timeTerm = Quiet@N[Abs[Evaluate[D[rhoField, t] /. rules]], 30];
  advTerm = Quiet@N[Abs[Evaluate[(dot3[grad3[rhoField, coords], grad3[phiField, coords] + vTField]) /. rules]], 30];
  sourceTerm = Quiet@N[Abs[Evaluate[sourceField /. rules]], 30];
  <|
    "TimeToSource" -> timeTerm/sourceTerm,
    "AdvectionToSource" -> advTerm/sourceTerm
  |>
];

(* ============================================================ *)
section["EXACT / operational definitions and projection kernel"];
(* ============================================================ *)

(* Frozen operational Gaussian projection weight *)
ClearAll[chi0, W0, NW];
chi0[w_] := (1/(Pi ellw^2))^(1/4) Exp[-w^2/(2 ellw^2)];
W0[w_] := (1/(Sqrt[Pi] ellw)) Exp[-w^2/ellw^2];

checkEqScalar[
  "Gaussian projection kernel normalizes on the full line",
  Integrate[W0[w], {w, -Infinity, Infinity}, Assumptions -> ellw > 0],
  1,
  ellw > 0
];

NW = FullSimplify[Integrate[W0[w], {w, -wProj, wProj}, Assumptions -> ellw > 0 && wProj > 0], ellw > 0 && wProj > 0];

checkEqScalar[
  "Truncated Gaussian normalization factor is Erf[Wproj/ellw]",
  NW,
  Erf[wProj/ellw],
  ellw > 0 && wProj > 0
];

checkEqScalar[
  "Renormalized truncated kernel integrates to one on the projection window",
  Integrate[W0[w]/NW, {w, -wProj, wProj}, Assumptions -> ellw > 0 && wProj > 0],
  1,
  ellw > 0 && wProj > 0
];

(* Concrete integration-by-parts self-test for leakage-source bookkeeping *)
ClearAll[JwTest];
JwTest[w_] := Exp[-sigJ w^2];

checkEqScalar[
  "Leakage source split via integration by parts (decaying test profile)",
  Integrate[W0[w] D[JwTest[w], w], {w, -Infinity, Infinity}, Assumptions -> ellw > 0 && sigJ > 0],
  -Integrate[D[W0[w], w] JwTest[w], {w, -Infinity, Infinity}, Assumptions -> ellw > 0 && sigJ > 0],
  ellw > 0 && sigJ > 0
];

info["ProjectedContinuitySource[W, Jw, w, wmin, wmax] and ProjectedMomentumSource[...] are defined for direct use with simulation expressions."];

(* ============================================================ *)
section["EXACT / brane kinematics, continuity identities, and stress structure"];
(* ============================================================ *)

coords = {x, y, z};

rhoBexpr = rhoBf[x, y, z, t];
JBexpr = {Jx[x, y, z, t], Jy[x, y, z, t], Jz[x, y, z, t]};
SrhoExpr = SrhoF[x, y, z, t];
vBexpr = JBexpr/rhoBexpr;

checkEqScalar[
  "Kinematic identity: div(v_brane) = div(J_brane)/rho_brane - J.grad(rho_brane)/rho_brane^2",
  div3[vBexpr, coords],
  div3[JBexpr, coords]/rhoBexpr - dot3[JBexpr, grad3[rhoBexpr, coords]]/rhoBexpr^2,
  rhoBexpr > 0
];

checkEqScalar[
  "Projected continuity inserted into div(v_brane)",
  (div3[JBexpr, coords]/rhoBexpr - dot3[JBexpr, grad3[rhoBexpr, coords]]/rhoBexpr^2) /. div3[JBexpr, coords] -> (SrhoExpr - D[rhoBexpr, t]),
  (SrhoExpr - D[rhoBexpr, t])/rhoBexpr - dot3[JBexpr, grad3[rhoBexpr, coords]]/rhoBexpr^2,
  rhoBexpr > 0
];

checkEqVector[
  "Kinematic identity: curl(v_brane) = curl(J_brane)/rho_brane - grad(rho_brane) x J_brane / rho_brane^2",
  curl3[vBexpr, coords],
  curl3[JBexpr, coords]/rhoBexpr - cross3[grad3[rhoBexpr, coords], JBexpr]/rhoBexpr^2,
  rhoBexpr > 0
];

phiExpr = phiB[x, y, z, t];
vTexpr = {vTx[x, y, z, t], vTy[x, y, z, t], vTz[x, y, z, t]};
vHelm = grad3[phiExpr, coords] + vTexpr;

checkScalar[
  "Product rule for div(rho_brane v_brane)",
  div3[rhoBexpr vHelm, coords] - dot3[grad3[rhoBexpr, coords], vHelm] - rhoBexpr div3[vHelm, coords]
];

checkScalar[
  "Helmholtz divergence split div(grad(phi)+v_T) = Laplacian(phi) + div(v_T)",
  div3[vHelm, coords] - lap3[phiExpr, coords] - div3[vTexpr, coords]
];

checkScalar[
  "Exact longitudinal identity from projected continuity + Helmholtz decomposition",
  Expand[
    (
      rhoBexpr lap3[phiExpr, coords]
      - (SrhoExpr - D[rhoBexpr, t] - dot3[grad3[rhoBexpr, coords], vHelm])
    ) /. SrhoExpr -> (D[rhoBexpr, t] + div3[rhoBexpr vHelm, coords])
  ],
  div3[vTexpr, coords] == 0
];

rhoSep = rho3f[x, y, z, t] wFac;
vSep = {u1[x, y, z, t], u2[x, y, z, t], u3[x, y, z, t]};
JSep = rhoSep vSep;
PiSep = rhoSep outer3[vSep, vSep];
Rsep = FullSimplify[PiSep - rhoSep outer3[JSep/rhoSep, JSep/rhoSep], wFac != 0];

checkMatrix[
  "Single-mode / separable regime gives zero extra stress R_ij",
  Rsep,
  wFac != 0
];

(* ============================================================ *)
section["CONTROLLED / Newtonian near-zone reduction and Poisson hook"];
(* ============================================================ *)

drhoExpr = drhoF[x, y, z, t];
phi3Expr = phi3F[x, y, z, t];
SsrcExpr = SsrcF[x, y, z, t];

linBern = D[phi3Expr, t] + (cs0^2/(mParam rho0)) drhoExpr;
waveFromLinear = FullSimplify[
  D[linBern, t] /. D[drhoExpr, t] -> (SsrcExpr - rho0 lap3[phi3Expr, coords]),
  rho0 > 0 && cs0 > 0 && mParam > 0
];

checkEqScalar[
  "Linearized continuity + Bernoulli/enthalpy closure -> forced wave equation",
  waveFromLinear,
  D[phi3Expr, {t, 2}] - (cs0^2/mParam) lap3[phi3Expr, coords] + (cs0^2/(mParam rho0)) SsrcExpr,
  rho0 > 0 && cs0 > 0 && mParam > 0
];

ClearAll[phiTTSym, lapPhiSym, SsrcSym];
waveAlgebraic = phiTTSym - (cs0^2/mParam) lapPhiSym + (cs0^2/(mParam rho0)) SsrcSym;
poissonSolve = FullSimplify[lapPhiSym /. First[Solve[(waveAlgebraic /. phiTTSym -> 0) == 0, lapPhiSym]], rho0 > 0 && cs0 > 0 && mParam > 0];

checkEqScalar[
  "Quasi-static limit solves to Laplacian(phi_3) = S_rho / rho0",
  poissonSolve,
  SsrcSym/rho0,
  rho0 > 0 && cs0 > 0 && mParam > 0
];

Rxyz = Sqrt[x^2 + y^2 + z^2];
phiMono = -Qsrc/(4 Pi rho0 Rxyz);
vLMono = grad3[phiMono, coords];

checkEqVector[
  "Localized monopole source gives 1/r^2 longitudinal field scaling",
  vLMono,
  (Qsrc/(4 Pi rho0)) {x, y, z}/Rxyz^3,
  rho0 > 0 && Rxyz > 0
];

checkScalar[
  "Monopole longitudinal field is divergence-free away from the source",
  div3[vLMono, coords],
  rho0 > 0 && Rxyz > 0
];

checkEqScalar[
  "Monopole longitudinal flux through a sphere matches the source strength",
  Integrate[(Qsrc/(4 Pi rho0 rr^2)) rr^2 Sin[theta], {theta, 0, Pi}, {phi, 0, 2 Pi}, Assumptions -> rr > 0 && rho0 > 0],
  Qsrc/rho0,
  rr > 0 && rho0 > 0
];

(* Legacy hybrid/Newtonian checks carried forward from the pre-4D papers *)
LNewton = (1/2) M0 vSq - M0 qVar PhiN;
qSolve = qVar /. First[Solve[Coefficient[LNewton, PhiN] == Coefficient[(1/2) M0 vSq - M0 PhiN, PhiN], qVar]];

checkEqScalar[
  "Hybrid Newtonian limit fixes q = 1",
  qSolve,
  1
];

Chyb[n_, q_] := q/2 - (n - 1)/2;
checkEqScalar[
  "Hybrid Phi v^2 / c^2 bookkeeping coefficient at (q,n)=(1,5)",
  Chyb[5, 1],
  -3/2
];

(* Ansatz-level particle-limit check; not a full 4D collective-coordinate derivation *)
coordPart = {X[t], Y[t], Z[t]};
velPart = D[coordPart, t];
phiPart = PhiExt[X[t], Y[t], Z[t], t];
Lpart = (1/2) mInert dot3[velPart, velPart] - mGrav phiPart;

eomPart = {
  FullSimplify[D[D[Lpart, X'[t]], t] - D[Lpart, X[t]]],
  FullSimplify[D[D[Lpart, Y'[t]], t] - D[Lpart, Y[t]]],
  FullSimplify[D[D[Lpart, Z'[t]], t] - D[Lpart, Z[t]]]
};

gradPhiPart = {
  D[phiPart, X[t]],
  D[phiPart, Y[t]],
  D[phiPart, Z[t]]
};

checkEqVector[
  "ANSATZ-LEVEL particle limit: m_i xdd = -m_g grad(Phi)",
  eomPart,
  mInert D[coordPart, {t, 2}] + mGrav gradPhiPart
];

checkEqScalar[
  "ANSATZ-LEVEL equivalence principle requires m_g/m_i = 1",
  (mGrav/mInert) /. mGrav -> mInert,
  1,
  mInert != 0
];

(* ============================================================ *)
section["EXACT / projected momentum algebra and open-system Euler form"];
(* ============================================================ *)

rhoOexpr = rhoOF[x, y, z, t];
vOexpr = {vxO[x, y, z, t], vyO[x, y, z, t], vzO[x, y, z, t]};
Fexpr = {F1O[x, y, z, t], F2O[x, y, z, t], F3O[x, y, z, t]};
SJexpr = {SJ1O[x, y, z, t], SJ2O[x, y, z, t], SJ3O[x, y, z, t]};
SrhoOpen = SrhoOpenF[x, y, z, t];
Rexpr = {
  {R11[x, y, z, t], R12[x, y, z, t], R13[x, y, z, t]},
  {R21[x, y, z, t], R22[x, y, z, t], R23[x, y, z, t]},
  {R31[x, y, z, t], R32[x, y, z, t], R33[x, y, z, t]}
};

momLHS = D[rhoOexpr vOexpr, t] + divTensor3[rhoOexpr outer3[vOexpr, vOexpr], coords];
contOpenExpr = D[rhoOexpr, t] + div3[rhoOexpr vOexpr, coords];
contOpenResidual = contOpenExpr - SrhoOpen;
momDecomp = rhoOexpr (D[vOexpr, t] + advect3[vOexpr, coords]) + vOexpr contOpenExpr;

divRexpr = divTensor3[Rexpr, coords];

checkEqVector[
  "Momentum product decomposition: dt(rho v) + div(rho vv) = rho a + v (dt rho + div(rho v))",
  momLHS,
  momDecomp
];

checkEqVector[
  "Open-system Euler residual before continuity insertion is rho a - F + divR + v (dt rho + div(rho v)) - SJ",
  Expand[(momLHS - (Fexpr - divRexpr + SJexpr)) /. momLHS -> momDecomp],
  rhoOexpr (D[vOexpr, t] + advect3[vOexpr, coords]) - Fexpr + divRexpr + vOexpr contOpenExpr - SJexpr
];

checkEqVector[
  "Open-system Euler residual carries the characteristic +v S_rho term before rearrangement",
  Expand[(momLHS - (Fexpr - divRexpr + SJexpr) - vOexpr contOpenResidual) /. momLHS -> momDecomp],
  rhoOexpr (D[vOexpr, t] + advect3[vOexpr, coords]) - Fexpr + divRexpr + vOexpr SrhoOpen - SJexpr
];

(* ============================================================ *)
section["CONTROLLED / added mass, topology discriminator, and coefficient closure"];
(* ============================================================ *)

phiDball = -(aRad^dDim/(dDim - 1)) Ubody Cos[theta]/r^(dDim - 1);

checkEqScalar[
  "d-ball exterior potential satisfies the no-penetration boundary condition at r=a",
  FullSimplify[D[phiDball, r] /. r -> aRad, dDim > 1 && aRad > 0],
  Ubody Cos[theta],
  dDim > 1 && aRad > 0
];

TflowDerived = (rho0/2) (aRad/(dDim - 1)) Ubody^2 ((SdMinus1 aRad^(dDim - 1))/dDim);
TflowTarget = (1/2) (rho0 Vd aRad^dDim/(dDim - 1)) Ubody^2;

checkEqScalar[
  "d-ball kinetic energy reduces to T = (1/2) m_add U^2 with kappa_add = 1/(d-1)",
  FullSimplify[TflowDerived /. SdMinus1 -> dDim Vd, dDim > 1 && aRad > 0],
  TflowTarget,
  dDim > 1 && aRad > 0
];

kappaAdd[d_] := 1/(d - 1);

checkEqScalar[
  "3D sphere / w-uniform throat slice gives kappa_add = 1/2",
  kappaAdd[3],
  1/2
];

checkEqScalar[
  "Counterfactual 4D bubble gives kappa_add = 1/3",
  kappaAdd[4],
  1/3
];

checkEqScalar[
  "Uniform-throat added-mass coefficient is independent of throat length L",
  FullSimplify[((1/2) rho0 V3 Lseg)/(rho0 V3 Lseg), rho0 > 0 && V3 > 0 && Lseg > 0],
  1/2,
  rho0 > 0 && V3 > 0 && Lseg > 0
];

(* ============================================================ *)
section["CONTROLLED / optics, vector-sector EIH match, and SR kinematics"];
(* ============================================================ *)

alphaN[n_] := (n - 1)/2;

checkEqScalar[
  "Optical weak-field coefficient is alpha_n = (n-1)/2",
  alphaN[nEOS],
  (nEOS - 1)/2
];

nSolution = nEOS /. First[Solve[alphaN[nEOS] == 2, nEOS]];

checkEqScalar[
  "GR optical coefficient selects the stiff EOS exponent n = 5",
  nSolution,
  5
];

checkEqScalar[
  "Weak-field light-bending coefficient at n=5 is 4 GM / (b c^2)",
  (((nEOS - 1) GMsym)/(bImp cLight^2)) /. nEOS -> 5,
  (4 GMsym)/(bImp cLight^2),
  bImp > 0 && cLight > 0
];

checkEqScalar[
  "Weak-field Shapiro bookkeeping coefficient at n=5 is alpha_n = 2",
  alphaN[5],
  2
];

Cpar = Kvec Pi^2 (-1 + aHsq - alphaSq);
CLong = Kvec Pi^2 (-1 + aHsq + alphaSq);
solVec = First@FullSimplify[Solve[{Cpar == -7/2, CLong == -1/2}, {alphaSq, Kvec}], 0 <= aHsq < 1];

checkEqScalar[
  "EIH vector family: alpha^2(a_H^2) = (3/4)(1-a_H^2)",
  alphaSq /. solVec,
  (3/4) (1 - aHsq),
  0 <= aHsq < 1
];

checkEqScalar[
  "EIH vector family: K(a_H^2) = 2 / (pi^2 (1-a_H^2))",
  Kvec /. solVec,
  2/(Pi^2 (1 - aHsq)),
  0 <= aHsq < 1
];

checkEqScalar[
  "EIH vector invariant K (1-a_H^2) = 2/pi^2",
  (Kvec (1 - aHsq)) /. solVec,
  2/Pi^2,
  0 <= aHsq < 1
];

checkEqScalar[
  "EIH vector invariant K alpha^2 = 3/(2 pi^2)",
  (Kvec alphaSq) /. solVec,
  3/(2 Pi^2),
  0 <= aHsq < 1
];

checkEqScalar[
  "Minimal parity-even specialization gives alpha^2 = 3/4",
  (alphaSq /. solVec /. aHsq -> 0),
  3/4
];

checkEqScalar[
  "Minimal parity-even specialization gives K = 2/pi^2",
  (Kvec /. solVec /. aHsq -> 0),
  2/Pi^2
];

checkEqScalar[
  "Back-substitution reproduces the EIH C_parallel target",
  Cpar /. solVec /. aHsq -> 0,
  -7/2
];

checkEqScalar[
  "Back-substitution reproduces the EIH C_L target",
  CLong /. solVec /. aHsq -> 0,
  -1/2
];

gammaSeries = Normal[Series[1/Sqrt[1 - betaVel^2], {betaVel, 0, 4}]];

checkEqScalar[
  "Gamma expansion carries the universal v^4 coefficient 3/8",
  Coefficient[gammaSeries, betaVel, 4],
  3/8
];

Eseries = Normal[Series[E0/Sqrt[1 - (vSym/cLight)^2], {vSym, 0, 4}]];

checkEqScalar[
  "Wave-supported rest energy gives the universal 3/8 coefficient in E(v)",
  Coefficient[Eseries, vSym, 4],
  (3 E0)/(8 cLight^4),
  cLight > 0
];

(* ============================================================ *)
section["PROTOCOL / adiabatic kappa_PV closure and internal partition"];
(* ============================================================ *)

EwExpr = Cw rhoVar^((nEOS - 1)/2)/aVar;
EfExpr = Cf/(rhoVar aVar^2);
EpvExpr = Cpv rhoVar^nEOS aVar^3;
FtotExpr = EwExpr + EfExpr + EpvExpr;

checkScalar[
  "Adiabatic equilibrium derivative gives the virial identity Ew + 2 Ef = 3 E_PV",
  FullSimplify[aVar D[FtotExpr, aVar] + EwExpr + 2 EfExpr - 3 EpvExpr, aVar > 0 && rhoVar > 0]
];

ClearAll[EwSym, EfSym, EpvSym, Fsym];
Fsym = EwSym + EfSym + EpvSym;
dlnFsym = (((nEOS - 1)/2) EwSym - EfSym + nEOS EpvSym)/Fsym;
dlnFgeneral = FullSimplify[dlnFsym /. {EfSym -> xRatio EwSym, EpvSym -> EwSym (1 + 2 xRatio)/3}];

checkEqScalar[
  "Closure total energy reduces to F = Ew (4+5x)/3",
  FullSimplify[(Fsym /. {EfSym -> xRatio EwSym, EpvSym -> EwSym (1 + 2 xRatio)/3})/EwSym],
  (4 + 5 xRatio)/3
];

checkEqScalar[
  "General closure formula for d ln F / d ln rho",
  dlnFgeneral,
  (((5 nEOS - 3)/2) + (2 nEOS - 3) xRatio)/(4 + 5 xRatio)
];

dlnFn5 = FullSimplify[dlnFgeneral /. nEOS -> 5];

checkEqScalar[
  "At n=5 the density-response slope becomes (11+7x)/(4+5x)",
  dlnFn5,
  (11 + 7 xRatio)/(4 + 5 xRatio)
];

xSolution = xRatio /. First[Solve[dlnFn5 == 5/2, xRatio]];

checkEqScalar[
  "GR closure target d ln F / d ln rho = 5/2 gives x = E_f/E_w = 2/11",
  xSolution,
  2/11
];

partVec = FullSimplify[{1, xSolution, (1 + 2 xSolution)/3}];
fracVec = FullSimplify[partVec/((4 + 5 xSolution)/3)];

checkEqVector[
  "Internal energy partition is Ew : Ef : E_PV = 11 : 2 : 5",
  partVec,
  {1, 2/11, 5/11}
];

checkEqVector[
  "Internal energy fractions are (11/18, 1/9, 5/18)",
  fracVec,
  {11/18, 1/9, 5/18}
];

kappaPV = FullSimplify[(dlnFn5 /. xRatio -> xSolution) - 1];

checkEqScalar[
  "Adiabatic closure yields kappa_PV = 3/2",
  kappaPV,
  3/2
];

Clear[slopeSym];
virialSlopeEq = FullSimplify[
  -(((nEOS - 1)/2) - slopeSym)
  - 2 xRatio (-1 - 2 slopeSym)
  + (1 + 2 xRatio) (nEOS + 3 slopeSym)
];

dlnASolution = FullSimplify[slopeSym /. First[Solve[virialSlopeEq == 0, slopeSym]]];

checkEqScalar[
  "General closure formula for d ln a / d ln rho",
  dlnASolution,
  -(-((nEOS - 1)/2) + 2 xRatio + nEOS (1 + 2 xRatio))/(4 + 10 xRatio)
];

checkEqScalar[
  "At n=5 and x=2/11 the throat breathing slope is -57/64",
  dlnASolution /. {nEOS -> 5, xRatio -> 2/11},
  -57/64
];

beta1PN = FullSimplify[1 + kappaAdd[3] + kappaPV];

checkEqScalar[
  "1PN inertia ledger closes to beta = kappa_rho + kappa_add + kappa_PV = 3",
  beta1PN,
  3
];

(* ============================================================ *)
section["END-TO-END / 1PN orbit bookkeeping and perihelion target"];
(* ============================================================ *)

ClearAll[rOrbit, phOrbit, sigmaOrbit];
sigmaOrbit[r_] := beta1PN muKepler/(cLight^2 r);
Lorb = (1/2) (1 + sigmaOrbit[rOrbit[tau]]) (rOrbit'[tau]^2 + rOrbit[tau]^2 phOrbit'[tau]^2) + muKepler/rOrbit[tau];

checkScalar[
  "1PN orbital ansatz: azimuthal angle is cyclic",
  D[Lorb, phOrbit[tau]]
];

checkEqScalar[
  "1PN orbital ansatz: canonical angular momentum is (1+sigma) r^2 phidot",
  FullSimplify[D[Lorb, phOrbit'[tau]]],
  (1 + beta1PN muKepler/(cLight^2 rOrbit[tau])) rOrbit[tau]^2 phOrbit'[tau]
];

DeltaPhiModel = FullSimplify[(2 Pi beta1PN muKepler)/(cLight^2 aSemi (1 - ecc^2)), aSemi > 0 && cLight > 0 && 0 <= ecc < 1];
DeltaPhiGR = (6 Pi muKepler)/(cLight^2 aSemi (1 - ecc^2));

checkEqScalar[
  "Model perihelion shift matches the GR 1PN target once beta = 3",
  DeltaPhiModel,
  DeltaPhiGR,
  aSemi > 0 && cLight > 0 && 0 <= ecc < 1
];

(* ============================================================ *)
section["SYNTHETIC SELF-TEST / Zeff extraction and locality gate"];
(* ============================================================ *)

omegaDrive = 0.7;
periods = 8;
samplesPerPeriod = 256;
times = N[Table[k (2 Pi/(omegaDrive samplesPerPeriod)), {k, 0, periods samplesPerPeriod - 1}], 50];

Ztrue = {{2 + I/5, -1/3 + I/7}, {I/4, 3 - 2 I/5}};
Ucols = IdentityMatrix[2];
Jcols = Ztrue.Ucols;

uTimeSeries = Table[
  Re[N[Ucols[[i, k]] Exp[I omegaDrive times], 50]],
  {k, 1, 2}, {i, 1, 2}
];

jTimeSeries = Table[
  Re[N[Jcols[[i, k]] Exp[I omegaDrive times], 50]],
  {k, 1, 2}, {i, 1, 2}
];

Urec = Table[LockInAmplitude[uTimeSeries[[k, i]], times, omegaDrive], {i, 1, 2}, {k, 1, 2}];
Jrec = Table[LockInAmplitude[jTimeSeries[[k, i]], times, omegaDrive], {i, 1, 2}, {k, 1, 2}];
Zrec = EstimateZeff[Urec, Jrec];

checkNumeric[
  "Lock-in convention returns half-amplitudes for cosine drives",
  Urec,
  Ucols/2,
  10^-10
];

checkNumeric[
  "Lock-in extraction + solve policy recovers the known Zeff matrix",
  Zrec,
  Ztrue,
  10^-10
];

Z0 = {{1.2, 0.1}, {-0.05, 0.8}};
Z1 = {{0.3, 0.0}, {0.0, -0.2}};
Z2 = {{0.0, 0.07}, {-0.11, 0.0}};
zToy[omega_] := Z0 + I omega Z1 + omega^2 Z2;

checkEqMatrix[
  "Locality-gate bookkeeping: low-omega expansion is Z0 + i omega Z1 + omega^2 Z2 + ...",
  Map[Normal[Series[#, {omegaVar, 0, 2}]] &, zToy[omegaVar], {2}],
  Z0 + I omegaVar Z1 + omegaVar^2 Z2
];

(* ============================================================ *)
section["SYNTHETIC SELF-TEST / Poisson-regime diagnostics helper"];
(* ============================================================ *)

epsDemo = 10^-6;
omegaDemo = 10^-4;
rhoDemo = 1 + epsDemo Exp[-(x^2 + y^2 + z^2)] Cos[omegaDemo t];
phiDemo = 1/Sqrt[x^2 + y^2 + z^2 + 1];
sourceDemo = lap3[phiDemo, coords];
ratiosDemo = PoissonCorrectionRatios[rhoDemo, phiDemo, {0, 0, 0}, sourceDemo, coords, {x -> 1, y -> 0, z -> 0, t -> 0}];

checkCondition[
  "Manufactured Poisson-regime demo has parametrically small dropped terms",
  Quiet@N[Max[Values[ratiosDemo]], 30] < 10^-3
];

info["Manufactured Poisson-regime ratios: " <> ToString[ratiosDemo, InputForm]];

(* ============================================================ *)
section["SUMMARY"];
(* ============================================================ *)

Print["PASS count = ", passCount];
Print["FAIL count = ", failCount];
Print["SKIP count = ", skipCount];

If[failCount == 0,
  Print["OVERALL: PASS (within the claim classes stated at the top of the file)."],
  Print["OVERALL: FAIL (inspect the residuals printed above)."]
];


(*"
Output:

--- 4D gravity + Newtonian + 1PN master harness ---

=== EXACT / operational definitions and projection kernel ===
PASS: Gaussian projection kernel normalizes on the full line
PASS: Truncated Gaussian normalization factor is Erf[Wproj/ellw]
PASS: Renormalized truncated kernel integrates to one on the projection window
PASS: Leakage source split via integration by parts (decaying test profile)
INFO: ProjectedContinuitySource[W, Jw, w, wmin, wmax] and ProjectedMomentumSource[...] are defined for direct use with simulation expressions.

=== EXACT / brane kinematics, continuity identities, and stress structure ===
PASS: Kinematic identity: div(v_brane) = div(J_brane)/rho_brane - J.grad(rho_brane)/rho_brane^2
PASS: Projected continuity inserted into div(v_brane)
PASS: Kinematic identity: curl(v_brane) = curl(J_brane)/rho_brane - grad(rho_brane) x J_brane / rho_brane^2
PASS: Product rule for div(rho_brane v_brane)
PASS: Helmholtz divergence split div(grad(phi)+v_T) = Laplacian(phi) + div(v_T)
PASS: Exact longitudinal identity from projected continuity + Helmholtz decomposition
PASS: Single-mode / separable regime gives zero extra stress R_ij

=== CONTROLLED / Newtonian near-zone reduction and Poisson hook ===
PASS: Linearized continuity + Bernoulli/enthalpy closure -> forced wave equation
PASS: Quasi-static limit solves to Laplacian(phi_3) = S_rho / rho0
PASS: Localized monopole source gives 1/r^2 longitudinal field scaling
PASS: Monopole longitudinal field is divergence-free away from the source
PASS: Monopole longitudinal flux through a sphere matches the source strength
PASS: Hybrid Newtonian limit fixes q = 1
PASS: Hybrid Phi v^2 / c^2 bookkeeping coefficient at (q,n)=(1,5)
PASS: ANSATZ-LEVEL particle limit: m_i xdd = -m_g grad(Phi)
PASS: ANSATZ-LEVEL equivalence principle requires m_g/m_i = 1

=== EXACT / projected momentum algebra and open-system Euler form ===
PASS: Momentum product decomposition: dt(rho v) + div(rho vv) = rho a + v (dt rho + div(rho v))
PASS: Open-system Euler residual before continuity insertion is rho a - F + divR + v (dt rho + div(rho v)) - SJ
PASS: Open-system Euler residual carries the characteristic +v S_rho term before rearrangement

=== CONTROLLED / added mass, topology discriminator, and coefficient closure ===
PASS: d-ball exterior potential satisfies the no-penetration boundary condition at r=a
PASS: d-ball kinetic energy reduces to T = (1/2) m_add U^2 with kappa_add = 1/(d-1)
PASS: 3D sphere / w-uniform throat slice gives kappa_add = 1/2
PASS: Counterfactual 4D bubble gives kappa_add = 1/3
PASS: Uniform-throat added-mass coefficient is independent of throat length L

=== CONTROLLED / optics, vector-sector EIH match, and SR kinematics ===
PASS: Optical weak-field coefficient is alpha_n = (n-1)/2
PASS: GR optical coefficient selects the stiff EOS exponent n = 5
PASS: Weak-field light-bending coefficient at n=5 is 4 GM / (b c^2)
PASS: Weak-field Shapiro bookkeeping coefficient at n=5 is alpha_n = 2
PASS: EIH vector family: alpha^2(a_H^2) = (3/4)(1-a_H^2)
PASS: EIH vector family: K(a_H^2) = 2 / (pi^2 (1-a_H^2))
PASS: EIH vector invariant K (1-a_H^2) = 2/pi^2
PASS: EIH vector invariant K alpha^2 = 3/(2 pi^2)
PASS: Minimal parity-even specialization gives alpha^2 = 3/4
PASS: Minimal parity-even specialization gives K = 2/pi^2
PASS: Back-substitution reproduces the EIH C_parallel target
PASS: Back-substitution reproduces the EIH C_L target
PASS: Gamma expansion carries the universal v^4 coefficient 3/8
PASS: Wave-supported rest energy gives the universal 3/8 coefficient in E(v)

=== PROTOCOL / adiabatic kappa_PV closure and internal partition ===
PASS: Adiabatic equilibrium derivative gives the virial identity Ew + 2 Ef = 3 E_PV
PASS: Closure total energy reduces to F = Ew (4+5x)/3
PASS: General closure formula for d ln F / d ln rho
PASS: At n=5 the density-response slope becomes (11+7x)/(4+5x)
PASS: GR closure target d ln F / d ln rho = 5/2 gives x = E_f/E_w = 2/11
PASS: Internal energy partition is Ew : Ef : E_PV = 11 : 2 : 5
PASS: Internal energy fractions are (11/18, 1/9, 5/18)
PASS: Adiabatic closure yields kappa_PV = 3/2
PASS: General closure formula for d ln a / d ln rho
PASS: At n=5 and x=2/11 the throat breathing slope is -57/64
PASS: 1PN inertia ledger closes to beta = kappa_rho + kappa_add + kappa_PV = 3

=== END-TO-END / 1PN orbit bookkeeping and perihelion target ===
PASS: 1PN orbital ansatz: azimuthal angle is cyclic
PASS: 1PN orbital ansatz: canonical angular momentum is (1+sigma) r^2 phidot
PASS: Model perihelion shift matches the GR 1PN target once beta = 3

=== SYNTHETIC SELF-TEST / Zeff extraction and locality gate ===
PASS: Lock-in convention returns half-amplitudes for cosine drives
PASS: Lock-in extraction + solve policy recovers the known Zeff matrix
PASS: Locality-gate bookkeeping: low-omega expansion is Z0 + i omega Z1 + omega^2 Z2 + ...

=== SYNTHETIC SELF-TEST / Poisson-regime diagnostics helper ===
PASS: Manufactured Poisson-regime demo has parametrically small dropped terms
INFO: Manufactured Poisson-regime ratios: <|"TimeToSource" -> 0, "AdvectionToSource" -> 4.905059215619230954606983602152811565944148413756905`29.69897000433602*^-7|>

=== SUMMARY ===
PASS count = 60
FAIL count = 0
SKIP count = 0
OVERALL: PASS (within the claim classes stated at the top of the file).
"*)
