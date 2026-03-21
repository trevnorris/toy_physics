
ClearAll["Global`*"];
Print["--- 2PN Family-1 TF inertia scale from parent PDE ---"];

passCount = 0;
failCount = 0;
section[name_String] := Print["\n=== ", name, " ==="];
pass[name_String] := (passCount++; Print["PASS: ", name]);
fail[name_String, res_] := (failCount++; Print["FAIL: ", name, "\n  residual -> ", res]);

zeroScalarQ[expr_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[expr, assum];
  TrueQ[res === 0] || TrueQ[Quiet @ FullSimplify[res == 0, assum]]
];

checkEqScalar[name_String, lhs_, rhs_, assum_: True] := Module[{res},
  res = Quiet @ FullSimplify[lhs - rhs, assum];
  If[zeroScalarQ[lhs - rhs, assum], pass[name], fail[name, res]]
];

assum = nEOS > 1 && a0 > 0 && Lam > 0 && Sigma > 0 && rhoGeom > 0 && betaGeom > 0 &&
  mPsi > 0 && KEOS > 0 && OmegaIn > 0;

section["1) Family-1 TF bulk inertia from the parent n-polytrope"];

alpha = 1/(nEOS - 1);

c0 = Assuming[assum,
  FullSimplify[(1/2) Beta[1/2, alpha + 1]]
];
c2 = Assuming[assum,
  FullSimplify[(1/2) Beta[3/2, alpha + 1]]
];

c2OverC0Closed = Assuming[assum,
  FullSimplify[c2/c0]
];
mLLClosed = Assuming[assum,
  FullSimplify[c2/(4 c0)]
];

rhoCenter = Assuming[assum,
  FullSimplify[((nEOS - 1) mPsi OmegaIn^2 (Lam a0)^2/(8 nEOS KEOS))^alpha]
];
rhoEffTF = Assuming[assum,
  FullSimplify[c0 rhoCenter]
];

mHatGeneral = {{3/5, 0}, {0, mLLClosed}};

Print["alpha = ", alpha];
Print["c0(n) = ", c0];
Print["c2(n) = ", c2];
Print["c2/c0 = ", c2OverC0Closed];
Print["mLL(n) = ", mLLClosed];
Print["rhoCenter = ", rhoCenter];
Print["rhoEffTF = ", rhoEffTF];
Print["mHatGeneral = ", MatrixForm[mHatGeneral]];

checkEqScalar[
  "Exact TF ratio c2/c0 = (n-1)/(3 n - 1)",
  c2OverC0Closed,
  (nEOS - 1)/(3 nEOS - 1),
  assum
];

checkEqScalar[
  "Exact TF axial inertia coefficient mLL = (n-1)/(4(3 n - 1))",
  mLLClosed,
  (nEOS - 1)/(4 (3 nEOS - 1)),
  assum
];

section["2) Frozen n = 5 specialization"];

c0n5 = Assuming[KEOS > 0 && mPsi > 0 && OmegaIn > 0 && a0 > 0 && Lam > 0,
  FullSimplify[c0 /. nEOS -> 5]
];
rhoCenterN5 = Assuming[KEOS > 0 && mPsi > 0 && OmegaIn > 0 && a0 > 0 && Lam > 0,
  FullSimplify[rhoCenter /. nEOS -> 5]
];
rhoEffN5 = Assuming[KEOS > 0 && mPsi > 0 && OmegaIn > 0 && a0 > 0 && Lam > 0,
  FullSimplify[rhoEffTF /. nEOS -> 5]
];
mHatN5 = Assuming[True, FullSimplify[mHatGeneral /. nEOS -> 5]];

Print["c0(5) = ", c0n5];
Print["N[c0(5),20] = ", N[c0n5, 20]];
Print["rhoCenter(5) = ", rhoCenterN5];
Print["rhoEffTF(5) = ", rhoEffN5];
Print["mHatTF(5) = ", MatrixForm[mHatN5]];

checkEqScalar[
  "TF n=5 axial inertia coefficient is 1/14",
  mHatN5[[2, 2]],
  1/14
];

section["3) Geometry breathing response with TF inertia"];

a =.; L =.;
a = Symbol["a"]; L = Symbol["L"];

V = (4 Pi/3) a^3 L;
A = 4 Pi a^2 L + (8 Pi/3) a^3;

sigma = Sigma/a0^3;
Pvac = rhoGeom Sigma/a0^4;
kappab = betaGeom Sigma/a0;

Egeom = Expand[Pvac V + sigma A + kappab a^2/L];
H = D[Egeom, {{a, L}, 2}];
g = {D[V, a], D[V, L]};

subs0 = {a -> a0, L -> Lam a0};
H0 = Assuming[assum, FullSimplify[H /. subs0]];
V0 = Assuming[assum, FullSimplify[V /. subs0]];
g0 = Assuming[assum, FullSimplify[g /. subs0]];

hBar = Assuming[assum, FullSimplify[(a0^2/Sigma) H0]];
gBar = Assuming[assum, FullSimplify[(a0/V0) g0]];

Delta0 = Assuming[assum,
  FullSimplify[(gBar . LinearSolve[hBar, gBar])/Sigma]
];
Delta2TF = Assuming[assum,
  FullSimplify[(gBar . LinearSolve[hBar, mHatN5 . LinearSolve[hBar, gBar]])/Sigma]
];
lamEffTF = Assuming[assum,
  FullSimplify[Delta0/Delta2TF]
];
Ytf = Assuming[assum,
  FullSimplify[(gBar . LinearSolve[hBar - s mHatN5, gBar])/Sigma]
];

Print["hBar = ", MatrixForm[hBar]];
Print["gBar = ", gBar];
Print["Delta0 = ", Delta0];
Print["Delta2TF = ", Delta2TF];
Print["lamEffTF = ", lamEffTF];
Print["Ytf(s) = ", Ytf];

section["4) EM-worked point with TF inertia"];

x01 = SetPrecision[2.40482555769577276862163187933, 50];
LamEM = N[Sqrt[2] Pi/x01, 50];
rhoEx = 1/10;
betaEx = 12;
target = 109/280;

SigmaStar = N[(Delta0 /. {a0 -> 1, Lam -> LamEM, rhoGeom -> rhoEx, betaGeom -> betaEx, Sigma -> 1})/target, 60];

hNum = N[hBar /. {a0 -> 1, Lam -> LamEM, rhoGeom -> rhoEx, betaGeom -> betaEx}, 60];
mNum = N[mHatN5, 60];
gNum = N[gBar /. {a0 -> 1, Lam -> LamEM}, 60];

{valsRaw, vecsRaw} = Eigensystem[LinearSolve[mNum, hNum]];
ord = Ordering[valsRaw];
vals = valsRaw[[ord]];
vecs = vecsRaw[[ord]];

vecsMass = Table[
  vecs[[i]]/Sqrt[vecs[[i]].mNum.vecs[[i]]],
  {i, Length[vecs]}
];

residues = Table[
  N[(gNum.vecsMass[[i]])^2/(SigmaStar vals[[i]]), 40],
  {i, Length[vals]}
];
residueFractions = N[residues/Total[residues], 40];

lamEffTFNum = N[lamEffTF /. {a0 -> 1, Lam -> LamEM, rhoGeom -> rhoEx, betaGeom -> betaEx}, 40];

sGrid = Table[x, {x, 0, 0.1 vals[[1]], 0.1 vals[[1]]/399}];
exactGrid = Table[
  Total@Table[residues[[i]]/(1 - sg/vals[[i]]), {i, Length[vals]}],
  {sg, sGrid}
];
padeGrid = Table[target/(1 - sg/lamEffTFNum), {sg, sGrid}];
relErr = Max[Abs[(padeGrid - exactGrid)/exactGrid]];

rhoEffWorked = Assuming[KEOS > 0 && OmegaIn > 0 && mPsi > 0,
  FullSimplify[rhoEffN5 /. {Lam -> LamEM, a0 -> 1}]
];
omegaScaleWorked = Assuming[KEOS > 0 && OmegaIn > 0 && mPsi > 0 && Sigma > 0,
  FullSimplify[Sigma/(rhoEffWorked ((4 Pi/3) LamEM))]
];
omegaMinusSq = N[vals[[1]] omegaScaleWorked, 40];
omegaPlusSq = N[vals[[2]] omegaScaleWorked, 40];
omegaEffSq = N[lamEffTFNum omegaScaleWorked, 40];

Print["LamEM = ", LamEM];
Print["SigmaStar = ", SigmaStar];
Print["mHatTF(worked point) = ", MatrixForm[mNum]];
Print["Dimensionless poles lambda_i = ", N[vals, 20]];
Print["Residues = ", residues];
Print["Residue fractions = ", residueFractions];
Print["lamEffTFNum = ", lamEffTFNum];
Print["Max relative error on 0 <= s <= 0.1 lambda_- = ", relErr];
Print["rhoEffWorked = ", rhoEffWorked];
Print["Omega_-^2 = ", omegaMinusSq];
Print["Omega_+^2 = ", omegaPlusSq];
Print["Omega_eff^2 = ", omegaEffSq];

checkEqScalar[
  "Static geometry channel still sums to 109/280 at the EM worked point",
  N[Total[residues], 30],
  N[target, 30]
];

If[And @@ Thread[residues > 0], pass["TF dynamic geometry residues stay positive"], fail["TF dynamic geometry residues stay positive", residues]];
If[relErr < 10^-3, pass["TF one-pole Padé error remains below 1e-3 on the low-frequency band"], fail["TF one-pole Padé error remains below 1e-3 on the low-frequency band", relErr]];

section["Summary"];
Print["PASS count = ", passCount];
Print["FAIL count = ", failCount];
