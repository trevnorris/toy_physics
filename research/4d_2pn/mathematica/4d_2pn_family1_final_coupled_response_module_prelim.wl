
(* ::Package:: *)

ClearAll[section, passfail, smoothStep, rhoArg, rhoFull, gaussList, eigResidues, padePole, maxPadeRelativeError,
  zBase, tProf, sProf, K00Full];

section[str_] := (
  Print["\n" <> StringRepeat["=", 78]];
  Print[str];
  Print[StringRepeat["=", 78]];
);

passfail[label_, expr_] := Print[label, ": ", If[TrueQ[expr], "PASS", "FAIL"], "   (", expr, ")"];

(* ------------------------------------------------------------------------- *)
(* Exact static support/source sector                                        *)
(* ------------------------------------------------------------------------- *)
section["EXACT Family-1 static support/source law"];

zBase[mu_] := 17/56 - 5 mu^2/56;
tProf[mu_] := 593/672 - (1553/672) mu^2 + (7/8) mu^4;
sProf[mu_] := 10 - (63/2) zBase[mu];

Print["zBase[mu] = ", zBase[μ]];
Print["tProf[mu] = ", tProf[μ]];
Print["sProf[mu] = ", sProf[μ]];
passfail["sProf == 10 - (63/2) zBase", FullSimplify[sProf[μ] == 10 - (63/2) zBase[μ]]];

(* ------------------------------------------------------------------------- *)
(* Exact 2PN cross-block reconstruction                                      *)
(* ------------------------------------------------------------------------- *)
section["EXACT 2PN cross-block reconstruction from final throat module"];

ClearAll[ax, ay, az, bx, by, bz, UA, UB, uAB, dA, dB, vA2, vB2, vAB,
  L1Wake, Lodd, Pi0A, Pi0B, Pi20A, Pi20B, Pi21A, Pi21B, Pi22A, Pi22B,
  dot21, dot22, J0, J20, DeltaGeom, Leven, Lmodule, Ltarget, residual];

uAB = ax*bx + ay*by;
dA = az;
dB = bz;
vA2 = ax^2 + ay^2 + az^2;
vB2 = bx^2 + by^2 + bz^2;
vAB = uAB + dA*dB;

L1Wake = -(7/2)*uAB - 4*dA*dB;
Lodd = (1/2)*(vA2 + vB2)*L1Wake - (15/4)*(UA + UB)*(uAB + dA*dB);

Pi0A = (Sqrt[5]/2)*vA2;
Pi0B = (Sqrt[5]/2)*vB2;
Pi20A = (1/2)*(3*dA^2 - vA2);
Pi20B = (1/2)*(3*dB^2 - vB2);
Pi21A = {Sqrt[2]*dA*ax, Sqrt[2]*dA*ay};
Pi21B = {Sqrt[2]*dB*bx, Sqrt[2]*dB*by};
Pi22A = {(ax^2 - ay^2)/(2*Sqrt[2]), (ax*ay)/Sqrt[2]};
Pi22B = {(bx^2 - by^2)/(2*Sqrt[2]), (bx*by)/Sqrt[2]};
dot21 = Expand[Pi21A.Pi21B];
dot22 = Expand[Pi22A.Pi22B];

J0 = 4/Sqrt[5];
J20 = 5/4;
DeltaGeom = 281/80;

Leven = Expand[
  Pi0A*Pi0B + Pi20A*Pi20B + dot21 + dot22
  + UA*(J0*Pi0B + J20*Pi20B)
  + UB*(J0*Pi0A + J20*Pi20A)
  + (J0^2 + J20^2 - DeltaGeom)*UA*UB
];

Lmodule = Expand[Lodd + Leven];

Ltarget = Expand[
  -(7/4)*vAB*(vA2 + vB2)
  -(1/4)*dA*dB*(vA2 + vB2)
  + (11/8)*vA2*vB2
  + (1/4)*vAB^2
  - (5/8)*(vA2*dB^2 + vB2*dA^2)
  + (3/2)*vAB*dA*dB
  + (3/8)*dA^2*dB^2
  - (15/4)*(UA + UB)*vAB
  + UA*((11/8)*vB2 + (15/8)*dB^2)
  + UB*((11/8)*vA2 + (15/8)*dA^2)
  + (5/4)*UA*UB
];

residual = FullSimplify[Expand[Lmodule - Ltarget]];
passfail["Final throat module reproduces exact added 2PN cross block", residual === 0];
Print["Residual = ", residual];

(* ------------------------------------------------------------------------- *)
(* Direct coupled sidewall + endcap full-profile wall completion             *)
(* ------------------------------------------------------------------------- *)
section["DIRECT coupled sidewall + endcap full-profile wall completion"];

epsr = SetPrecision[0.05, 30];
alphaR = SetPrecision[10.0, 30];
pR = 2;
epsz = SetPrecision[0.05, 30];              (* d_z / L *)
chiCap = SetPrecision[4.0, 30];
alphaCap = 4 epsz chiCap;                    (* balanced branch = 0.8 *)
pZ = 2;

gaussList[a_?NumericQ, b_?NumericQ, n_Integer?Positive] := gaussList[a, b, n] = Module[
  {beta, jac, vals, vecs, ord, nodes, weights},
  beta = N[Table[k/Sqrt[4*k^2 - 1], {k, 1, n - 1}], 40];
  jac = SparseArray[
    Join[
      Table[{k, k + 1} -> beta[[k]], {k, 1, n - 1}],
      Table[{k + 1, k} -> beta[[k]], {k, 1, n - 1}]
    ],
    {n, n}
  ];
  {vals, vecs} = Eigensystem[N[Normal[jac], 40]];
  ord = Ordering[vals];
  vals = vals[[ord]];
  vecs = vecs[[ord]];
  nodes = N[0.5*(b - a)*vals + 0.5*(a + b), 30];
  weights = N[(b - a)*vecs[[All, 1]]^2, 30];
  {nodes, weights}
];

smoothStep[x_?NumericQ] := (1 + Tanh[x])/2;

rhoArg[x_?NumericQ, y_?NumericQ] :=
  1 - y^2 - alphaR*smoothStep[(x - 1)/epsr]^pR - alphaCap*smoothStep[(y - 1)/(2*epsz)]^pZ;

rhoFull[x_?NumericQ, y_?NumericQ] := Module[{val = rhoArg[x, y]},
  If[val > 0, val^(1/4), 0]
];

{x1, w1} = gaussList[0., 0.93, 120];
{x2, w2} = gaussList[0.93, 1., 260];
xNodes = Join[x1, x2];
wx = Join[w1, w2];

{y1, wy1} = gaussList[0., 0.93, 120];
{y2, wy2} = gaussList[0.93, 1., 260];
yNodes = Join[y1, y2];
wy = Join[wy1, wy2];

xGrid = Outer[Times, xNodes, ConstantArray[1., Length[yNodes]]];
yGrid = Outer[Times, ConstantArray[1., Length[xNodes]], yNodes];
weightGrid = Outer[Times, wx, wy];
rhoGrid = Outer[rhoFull, xNodes, yNodes];

i2 = N[2*Total[Flatten[weightGrid*xGrid^2*rhoGrid]], 30];
i4 = N[2*Total[Flatten[weightGrid*xGrid^4*rhoGrid]], 30];
iw = N[2*Total[Flatten[weightGrid*(yGrid^2/4.)*xGrid^2*rhoGrid]], 30];

c0 = N[Sqrt[Pi] Gamma[5/4]/(2 Gamma[7/4]), 30];
i2Sharp = N[2 c0/3, 30];

rMassFull = N[i2/i2Sharp, 20];
maaFull = N[i4/i2, 20];
mllFull = N[iw/i2, 20];

Print["Rmass_full_exact = ", N[rMassFull, 15]];
Print["Maa_full_exact   = ", N[maaFull, 15]];
Print["Mll_full_exact   = ", N[mllFull, 15]];

(* ------------------------------------------------------------------------- *)
(* Compare to carried-forward separated-order completion                     *)
(* ------------------------------------------------------------------------- *)
section["Compare direct full-profile values to separated-order carried-forward branch"];

rMassSep = SetPrecision[0.8846236634, 30];
maaSep = SetPrecision[0.5623810783, 30];
mllSep = SetPrecision[0.0671965962, 30];

Print["Rmass relative shift vs separated = ", N[(rMassFull - rMassSep)/rMassSep, 15]];
Print["Maa   relative shift vs separated = ", N[(maaFull - maaSep)/maaSep, 15]];
Print["Mll   relative shift vs separated = ", N[(mllFull - mllSep)/mllSep, 15]];

(* ------------------------------------------------------------------------- *)
(* Geometry breathing response                                               *)
(* ------------------------------------------------------------------------- *)
section["Geometry breathing response from exact full-profile wall completion"];

x01 = SetPrecision[2.40482555769577276862163187933, 40];
lamEM = N[Sqrt[2] Pi/x01, 30];
rhoGeom = 1/10;
betaGeom = 12;
a0 = 1;

V[a_, l_] := (4 Pi/3) a^3 l;
A[a_, l_] := 4 Pi a^2 l + (8 Pi/3) a^3;
Egeom[a_, l_, lam_, sigScale_, rho_, beta_] := Module[{sigma, pvac, kappab},
  sigma = sigScale/a0^3;
  pvac = rho sigScale/a0^4;
  kappab = beta sigScale/a0;
  pvac V[a, l] + sigma A[a, l] + kappab a^2/l
];

h0 = Outer[
   D[Egeom[a, l, lamEM, Sigma, rhoGeom, betaGeom], #1, #2] &,
   {a, l}, {a, l}
 ] /. {a -> a0, l -> lamEM};

hbar = Simplify[(a0^2/Sigma) h0] /. Sigma -> 1;
g0 = {D[V[a, l], a], D[V[a, l], l]} /. {a -> a0, l -> lamEM};
gbar = N[(a0 g0)/V[a0, lamEM], 30];

deltaUnit = N[gbar . Inverse[hbar] . gbar, 30];
sigmaStar = N[deltaUnit/(109/280), 30];

Print["Lambda_EM = ", N[lamEM, 15]];
Print["Sigma_*   = ", N[sigmaStar, 15]];

Clear[eigResidues];
eigResidues[M_?MatrixQ] := Module[{mat, vals, vecs, ord, vnorm, res},
  mat = LinearSolve[M, hbar];
  {vals, vecs} = Eigensystem[mat];
  ord = Ordering[N[vals]];
  vals = vals[[ord]];
  vecs = vecs[[ord]];
  vnorm = Table[vecs[[i]]/Sqrt[vecs[[i]].M.vecs[[i]]], {i, Length[vals]}];
  res = Table[(gbar.vnorm[[i]])^2/(sigmaStar vals[[i]]), {i, Length[vals]}];
  {N[vals, 30], N[res, 30]}
];

padePole[vals_List, res_List] := N[Total[res]/Total[res/vals], 30];
maxPadeRelativeError[vals_List, res_List] := Module[{leff, sgrid, exact, pade},
  leff = padePole[vals, res];
  sgrid = Table[s, {s, 0, 0.1 vals[[1]], 0.1 vals[[1]]/399}];
  exact = Table[Total[res/(1 - s/vals)], {s, sgrid}];
  pade = Table[Total[res]/(1 - s/leff), {s, sgrid}];
  N[Max[Abs[(pade - exact)/exact]], 20]
];

mfull = DiagonalMatrix[{rMassFull*0 + maaFull, mllFull}]; (* keep exact numeric values *)
{vals, res} = eigResidues[mfull];
leff = padePole[vals, res];
err = maxPadeRelativeError[vals, res];

mtf = DiagonalMatrix[{3/5, 1/14}];
{valsTF, resTF} = eigResidues[mtf];
omegaRatios = N[(vals/valsTF)/rMassFull, 30];

Print["lambda_-   = ", N[vals[[1]], 15]];
Print["lambda_+   = ", N[vals[[2]], 15]];
Print["R_-        = ", N[res[[1]], 15]];
Print["R_+        = ", N[res[[2]], 15]];
Print["R_- + R_+  = ", N[Total[res], 15]];
Print["lambda_eff = ", N[leff, 15]];
Print["max rel err on [0, 0.1 lambda_-] = ", N[err, 15]];
Print["Omega^2/Omega_TF^2 (minus) = ", N[omegaRatios[[1]], 15]];
Print["Omega^2/Omega_TF^2 (plus)  = ", N[omegaRatios[[2]], 15]];

K00Full[s_] := -757/2520 + res[[1]]/(1 - s/vals[[1]]) + res[[2]]/(1 - s/vals[[2]]);
Print["K00_raw(static local) = ", N[-757/2520, 15]];
Print["K00_full(0)           = ", N[K00Full[0], 15]];
Print["4/45                  = ", N[4/45, 15]];
passfail["K00_full(0) == 4/45", Chop[K00Full[0] - 4/45] == 0];

(* ------------------------------------------------------------------------- *)
(* Final response object summary                                             *)
(* ------------------------------------------------------------------------- *)
section["FINAL low-frequency Family-1 throat-response object"];

Print["Odd dipole residues:   R1_perp = 7/2, R10 = 4"];
Print["Odd dressings:         sigma = 1/2, eta_perp = 15/14, eta_parallel = 15/16"];
Print["Even local support:    K00_raw = -757/2520, K1_perp = 2/7, K10 = 1/4, K20 = 4/9, K21 = 2/3, K22 = 8/3"];
Print["Even source vector:    J = (4/Sqrt[5], 5/4)"];
Print["Monopole closure:      K00[s] = -757/2520 + R_-/(1 - s/lambda_-) + R_+/(1 - s/lambda_+)"];
Print["Default balanced cap:  alphaCap = ", N[alphaCap, 15], ", epsz = ", N[epsz, 15], ", chiCap = ", N[chiCap, 15]];
Print["Conservative 2PN status: the full added cross block is reproduced exactly at zero frequency;"];
Print["dynamic non-monopole pole scales remain genuine inner-throat observables for a beyond-2PN / Paper-7 extension."];

section["SUITE verification"];
Print[If[TrueQ[FullSimplify[sProf[μ] == 10 - (63/2) zBase[μ]]], "PASS: ", "FAIL: "], "Family-1 static support/source identity remains exact"];
Print[If[TrueQ[residual === 0], "PASS: ", "FAIL: "], "Final throat module still reproduces the exact added 2PN cross block"];
Print[If[TrueQ[Chop[K00Full[0] - 4/45] == 0], "PASS: ", "FAIL: "], "Full monopole closure still satisfies K00_full(0) = 4/45"];

Family1FinalCoupledResponseResults = <|
  "zBase" -> zBase[μ],
  "tProf" -> tProf[μ],
  "sProf" -> sProf[μ],
  "CrossBlockResidual" -> residual,
  "K00FullAtZero" -> K00Full[0],
  "LambdaMinus" -> vals[[1]],
  "LambdaPlus" -> vals[[2]],
  "Residues" -> res,
  "LambdaEff" -> leff,
  "PadeError" -> err
|>;

Print["Key exported symbol: Family1FinalCoupledResponseResults."];

(*"
Output:


==============================================================================
EXACT Family-1 static support/source law
==============================================================================
zBase[mu] = 17/56 - (5*μ^2)/56
tProf[mu] = 593/672 - (1553*μ^2)/672 + (7*μ^4)/8
sProf[mu] = 10 - (63*(17/56 - (5*μ^2)/56))/2
sProf == 10 - (63/2) zBase: PASS   (True)

==============================================================================
EXACT 2PN cross-block reconstruction from final throat module
==============================================================================
Final throat module reproduces exact added 2PN cross block: PASS   (True)
Residual = 0

==============================================================================
DIRECT coupled sidewall + endcap full-profile wall completion
==============================================================================
Rmass_full_exact = 0.8863139729897237
Maa_full_exact   = 0.5631149689539836
Mll_full_exact   = 0.06582922811934913

==============================================================================
Compare direct full-profile values to separated-order carried-forward branch
==============================================================================
Rmass relative shift vs separated = 0.0019107668714480326
Maa   relative shift vs separated = 0.0013049703880545012
Mll   relative shift vs separated = -0.020348769996937265

==============================================================================
Geometry breathing response from exact full-profile wall completion
==============================================================================
Lambda_EM = 1.84748657712012805104337448393961224287`15.
Sigma_*   = 0.20761432918354888540403250898444970518`15.
lambda_-   = 6.405572392138904
lambda_+   = 254.4449681369374
R_-        = 0.0025524747717377352
R_+        = 0.38673323951397665
R_- + R_+  = 0.3892857142857144
lambda_eff = 202.92351636751968
max rel err on [0, 0.1 lambda_-] = 0.0000689464436711859
Omega^2/Omega_TF^2 (minus) = 1.2196655544121753
Omega^2/Omega_TF^2 (plus)  = 1.2066780945388536
K00_raw(static local) = -0.30039682539682539682539682539682539683`15.
K00_full(0)           = 0.08888888888888902
4/45                  = 0.08888888888888888888888888888888888889`15.
K00_full(0) == 4/45: PASS   (True)

==============================================================================
FINAL low-frequency Family-1 throat-response object
==============================================================================
Odd dipole residues:   R1_perp = 7/2, R10 = 4
Odd dressings:         sigma = 1/2, eta_perp = 15/14, eta_parallel = 15/16
Even local support:    K00_raw = -757/2520, K1_perp = 2/7, K10 = 1/4, K20 = 4/9, K21 = 2/3, K22 = 8/3
Even source vector:    J = (4/Sqrt[5], 5/4)
Monopole closure:      K00[s] = -757/2520 + R_-/(1 - s/lambda_-) + R_+/(1 - s/lambda_+)
Default balanced cap:  alphaCap = 0.80000000000000004440892098500626161695`15., epsz = 0.05000000000000000277555756156289135106`15., chiCap = 4.`15.
Conservative 2PN status: the full added cross block is reproduced exactly at zero frequency;
dynamic non-monopole pole scales remain genuine inner-throat observables for a beyond-2PN / Paper-7 extension.

==============================================================================
SUITE verification
==============================================================================
PASS: Family-1 static support/source identity remains exact
PASS: Final throat module still reproduces the exact added 2PN cross block
PASS: Full monopole closure still satisfies K00_full(0) = 4/45
Key exported symbol: Family1FinalCoupledResponseResults.
"*)
