
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

(*"
Output:

--- 2PN Family-1 TF inertia scale from parent PDE ---

=== 1) Family-1 TF bulk inertia from the parent n-polytrope ===
alpha = (-1 + nEOS)^(-1)
c0(n) = Beta[1/2, nEOS/(-1 + nEOS)]/2
c2(n) = Beta[3/2, nEOS/(-1 + nEOS)]/2
c2/c0 = (-1 + nEOS)/(-1 + 3*nEOS)
mLL(n) = (-1 + nEOS)/(-4 + 12*nEOS)
rhoCenter = ((a0^2*Lam^2*mPsi*(-1 + nEOS)*OmegaIn^2)/(KEOS*nEOS))^(-1 + nEOS)^(-1)/8^(-1 + nEOS)^(-1)
rhoEffTF = 2^((2 + nEOS)/(1 - nEOS))*((a0^2*Lam^2*mPsi*(-1 + nEOS)*OmegaIn^2)/(KEOS*nEOS))^(-1 + nEOS)^(-1)*Beta[1/2, nEOS/(-1 + nEOS)]
mHatGeneral = MatrixForm[{{3/5, 0}, {0, (-1 + nEOS)/(-4 + 12*nEOS)}}]
PASS: Exact TF ratio c2/c0 = (n-1)/(3 n - 1)
PASS: Exact TF axial inertia coefficient mLL = (n-1)/(4(3 n - 1))

=== 2) Frozen n = 5 specialization ===
c0(5) = (Sqrt[Pi]*Gamma[5/4])/(2*Gamma[7/4])
N[c0(5),20] = 0.87401918476403993682161319663185145628`20.
rhoCenter(5) = ((mPsi/KEOS)^(1/4)*Sqrt[a0*Lam*OmegaIn])/10^(1/4)
rhoEffTF(5) = ((mPsi/KEOS)^(1/4)*Sqrt[a0*Lam*OmegaIn]*Sqrt[Pi]*Gamma[5/4])/(2*10^(1/4)*Gamma[7/4])
mHatTF(5) = MatrixForm[{{3/5, 0}, {0, 1/14}}]
PASS: TF n=5 axial inertia coefficient is 1/14

=== 3) Geometry breathing response with TF inertia ===
hBar = MatrixForm[{{(2*betaGeom)/Lam + 8*Pi*(2 + Lam + Lam*rhoGeom), (-2*betaGeom)/Lam^2 + 4*Pi*(2 + rhoGeom)}, {(-2*betaGeom)/Lam^2 + 4*Pi*(2 + rhoGeom), (2*betaGeom)/Lam^3}}]
gBar = {3, Lam^(-1)}
Delta0 = (-4*betaGeom + Lam*Pi*(-2 + Lam*(5 + 2*rhoGeom)))/(2*Pi*(Lam^3*Pi*(2 + rhoGeom)^2 - betaGeom*(2 + Lam*(3 + 2*rhoGeom)))*Sigma)
Delta2TF = (4*betaGeom^2*(42 + 5*Lam^2) + Lam^4*Pi^2*(248 - 40*Lam*(4 + rhoGeom) + 42*rhoGeom*(4 + rhoGeom) + 5*Lam^2*(4 + rhoGeom)^2) - 4*betaGeom*Lam^2*Pi*(42*(2 + rhoGeom) + 5*Lam*(-4 + Lam*(4 + rhoGeom))))/(1120*Pi^2*(Lam^3*Pi*(2 + rhoGeom)^2 - betaGeom*(2 + Lam*(3 + 2*rhoGeom)))^2*Sigma)
lamEffTF = (560*Pi*(Lam^3*Pi*(2 + rhoGeom)^2 - betaGeom*(2 + Lam*(3 + 2*rhoGeom)))*(-4*betaGeom + Lam*Pi*(-2 + Lam*(5 + 2*rhoGeom))))/(4*betaGeom^2*(42 + 5*Lam^2) + Lam^4*Pi^2*(248 - 40*Lam*(4 + rhoGeom) + 42*rhoGeom*(4 + rhoGeom) + 5*Lam^2*(4 + rhoGeom)^2) - 4*betaGeom*Lam^2*Pi*(42*(2 + rhoGeom) + 5*Lam*(-4 + Lam*(4 + rhoGeom))))
Ytf(s) = -((560*(4*betaGeom + Lam*Pi*(2 - Lam*(5 + 2*rhoGeom))) - 3*Lam*(14 + 15*Lam^2)*s)/((2*betaGeom*(-560*Pi*(2 + 3*Lam + 2*Lam*rhoGeom) + (42 + 5*Lam^2)*s) + Lam^3*(1120*Pi^2*(2 + rhoGeom)^2 + 40*Pi*(2 + Lam + Lam*rhoGeom)*s - 3*s^2))*Sigma))

=== 4) EM-worked point with TF inertia ===
LamEM = 1.8474865771201280510433744839396122428696573796896114673321817811568650327843`50.
SigmaStar = 0.20761432918354888540403250898444970518124118990882318258536279134129844964869`48.15554576568318
mHatTF(worked point) = MatrixForm[{{0.6`60., 0}, {0, 0.07142857142857142857142857142857142857142857142857142857142857142857142857143`60.}}]
Dimensionless poles lambda_i = {5.92556257692685858382223855943745125857`20., 237.91117494303324065318581486681708803893`20.}
Residues = {0.002628002865702182717908131119550467474353909167098464977`40., 0.3866577114200121029963775831661638182399318051186158207373`40.}
Residue fractions = {0.006750833049510194137745657921781017365312794190711653152`39.69897000433602, 0.993249166950489805862254342078218982634687205809288346848`39.69897000433602}
lamEffTFNum = 188.1769589801712634979284218699305818526037636010761944631566`40.
Max relative error on 0 <= s <= 0.1 lambda_- = 0.00007100969970122965
rhoEffWorked = 0.66805406566754947933377993028331689290644414273125724893935683781830504186007`50.30102999566398*(mPsi/KEOS)^(1/4)*Sqrt[OmegaIn]
Omega_-^2 = (1.1461674649753605162726260465006270842628882262896971236478`40.*(KEOS/mPsi)^(1/4)*Sigma)/Sqrt[OmegaIn]
Omega_+^2 = (46.018592282791748388736666001371046906470351267877636869586`40.*(KEOS/mPsi)^(1/4)*Sigma)/Sqrt[OmegaIn]
Omega_eff^2 = (36.3986212686211154220866044111710401757808298741717865596416`40.*(KEOS/mPsi)^(1/4)*Sigma)/Sqrt[OmegaIn]
PASS: Static geometry channel still sums to 109/280 at the EM worked point
PASS: TF dynamic geometry residues stay positive
PASS: TF one-pole Padé error remains below 1e-3 on the low-frequency band

=== Summary ===
PASS count = 6
FAIL count = 0
"*)
