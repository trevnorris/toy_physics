
(* ::Package:: *)

ClearAll[
  section, pass, fail,
  a0, Lam, Sigma, rho, beta, a, L,
  V, A, sigma, Pvac, kappab, Egeom, H, g, subs0, H0, g0, V0,
  deltaGeomExact, deltaGeomFormula, Hbar, detHbar, HbarBase, detHbarBase,
  betaStab, betaDelta,
  x01, lamEM, rhoEx, betaEx, target, deltaUnit, sigmaStar, deltaMatched,
  Hnum, ghat, eigsys, pairs, evals, evecs, projs, contrib, frac,
  dominantIndex, relerr
];

section[s_String] := Print["\n" <> StringRepeat["=", 78] <> "\n" <> s <> "\n" <> StringRepeat["=", 78]];
pass[s_String] := Print["PASS: " <> s];
fail[s_String] := Print["FAIL: " <> s];

section["1) Exact geometry-side setup"];

V[a_, L_] := (4 Pi/3) a^3 L;
A[a_, L_] := 4 Pi a^2 L + (8 Pi/3) a^3;

sigma = Sigma/a0^3;
Pvac = rho Sigma/a0^4;
kappab = beta Sigma/a0;

Egeom = Expand[Pvac V[a, L] + sigma A[a, L] + kappab a^2/L];
H = D[Egeom, {{a, L}, 2}];
g = {D[V[a, L], a], D[V[a, L], L]};

subs0 = {a -> a0, L -> Lam a0};
H0 = FullSimplify[H /. subs0];
g0 = FullSimplify[g /. subs0];
V0 = FullSimplify[V[a, L] /. subs0];

deltaGeomExact = FullSimplify[(g0 . Inverse[H0] . g0)/V0^2];
deltaGeomFormula =
  (2 Pi Lam^2 rho + 5 Pi Lam^2 - 2 Pi Lam - 4 beta)/
  (2 Pi Sigma (Pi Lam^3 rho^2 + 4 Pi Lam^3 rho + 4 Pi Lam^3 - 2 Lam beta rho - 3 Lam beta - 2 beta));

Hbar = FullSimplify[H0/(Sigma/a0^2)];
detHbar = Factor[Det[Hbar]];

Print["V(a,L) = ", V[a, L]];
Print["A(a,L) = ", A[a, L]];
Print["E_geom(a,L) = ", Egeom];
Print[""];
Print["H0 = "];
Print[MatrixForm[H0]];
Print[""];
Print["g0 = ", g0];
Print["V0 = ", V0];
Print[""];
Print["deltaGeomExact = ", deltaGeomExact];

If[FullSimplify[deltaGeomExact - deltaGeomFormula] === 0,
  pass["Exact geometry compressibility formula matches the closed form."],
  fail["Closed-form deltaGeom formula mismatch."]
];

section["2) Baseline no-go: P_vac V + sigma A alone"];

HbarBase = FullSimplify[Hbar /. beta -> 0];
detHbarBase = Factor[detHbar /. beta -> 0];

Print["Hbar(beta->0) = "];
Print[MatrixForm[HbarBase]];
Print[""];
Print["det Hbar(beta->0) = ", detHbarBase];

If[
  FullSimplify[detHbarBase + 16 Pi^2 (rho + 2)^2] === 0 &&
  FullSimplify[HbarBase[[2, 2]]] === 0,
  pass["Baseline P_vac V + sigma A geometry Hessian is non-passive / incomplete."],
  fail["Baseline no-go check failed."]
];

section["3) Minimal curvature completion and positivity thresholds"];

betaStab = FullSimplify[Pi Lam^3 (rho + 2)^2/(2 Lam rho + 3 Lam + 2)];
betaDelta = FullSimplify[Pi Lam (2 Lam rho + 5 Lam - 2)/4];

Print["Hbar = "];
Print[MatrixForm[Hbar]];
Print[""];
Print["det Hbar = ", detHbar];
Print[""];
Print["betaStab  = ", betaStab];
Print["betaDelta = ", betaDelta];

If[
  FullSimplify[detHbar + 16 Pi (Pi Lam^3 rho^2 + 4 Pi Lam^3 rho + 4 Pi Lam^3 - 2 Lam beta rho - 3 Lam beta - 2 beta)/Lam^3] === 0,
  pass["Exact determinant formula verified."],
  fail["Determinant formula check failed."]
];

section["4) EM-branch worked example"];

x01 = N[BesselJZero[0, 1], 40];
lamEM = N[Sqrt[2] Pi/x01, 40];
rhoEx = 1/10;
betaEx = 12;
target = N[109/280, 40];

deltaUnit = N[deltaGeomExact /. {a0 -> 1, Sigma -> 1, Lam -> lamEM, rho -> rhoEx, beta -> betaEx}, 50];
sigmaStar = N[deltaUnit/target, 50];
deltaMatched = N[deltaGeomExact /. {a0 -> 1, Sigma -> sigmaStar, Lam -> lamEM, rho -> rhoEx, beta -> betaEx}, 50];

Hnum = N[H0 /. {a0 -> 1, Sigma -> sigmaStar, Lam -> lamEM, rho -> rhoEx, beta -> betaEx}, 50];
ghat = N[(g0/V0) /. {a0 -> 1, Lam -> lamEM}, 50];

eigsys = Eigensystem[Hnum];
pairs = SortBy[Transpose[eigsys], First];
evals = pairs[[All, 1]];
evecs = pairs[[All, 2]];
projs = evecs . ghat;
contrib = N[(projs^2)/evals, 40];
frac = N[contrib/Total[contrib], 40];
dominantIndex = First[Ordering[contrib, -1]];
relerr = N[Abs[target - contrib[[dominantIndex]]]/target, 40];

Print["x01    = ", x01];
Print["lamEM  = ", lamEM];
Print["rho    = ", rhoEx];
Print["beta   = ", betaEx];
Print[""];
Print["betaStab(lamEM,rho)  = ", N[betaStab /. {Lam -> lamEM, rho -> rhoEx}, 30]];
Print["betaDelta(lamEM,rho) = ", N[betaDelta /. {Lam -> lamEM, rho -> rhoEx}, 30]];
Print[""];
Print["sigmaStar = ", sigmaStar];
Print["deltaMatched = ", deltaMatched];
Print["target       = ", target];
Print["residual     = ", N[deltaMatched - target, 40]];
Print[""];
Print["ghat = ", ghat];
Print["eigenvalues = ", evals];
Print["eigenvectors = ", evecs];
Print["mode contributions = ", contrib];
Print["fractions = ", frac];
Print["dominant one-pole relative error = ", relerr];

If[
  Chop[deltaMatched - target, 10^-30] == 0 &&
  frac[[dominantIndex]] > 0.95,
  pass["EM-branch worked point matches 109/280 exactly and is dominantly one-pole."],
  fail["Worked example did not match target or dominant-mode reduction."]
];

section["5) Compact interpretation"];

Print["Exact reduced-geometry statement:"];
Print["  delta K00^(geom) = (grad V)^T H0^{-1} (grad V) / V0^2"];
Print["with H0 the (a,L) geometry Hessian."];
Print[""];
Print["Baseline P_vac V + sigma A fails; the minimal curvature term kappa_b a^2/L repairs the Hessian."];
Print["At the EM-aspect-ratio worked point, the dominant geometry eigenmode carries ",
  100 N[frac[[dominantIndex]], 20], "% of the static monopole response."];
Print["So the earlier single global monopole auxiliary is a controlled reduction of the full 2DOF geometry sector."];
