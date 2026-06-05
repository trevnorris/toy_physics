(* moving_throat_pde_stage009_projected_maxwell_near_throat_mathematica_audit.wl *)
(* Independent Mathematica verification of unit 009. *)

ClearAll["Global`*"];
$HistoryLength = 0;

fmt[expr_] := ToString[InputForm[expr]];

assertZero[label_String, expr_] := Module[{res},
  res = FullSimplify[expr, Assumptions -> $Assumptions];
  If[res =!= 0,
    Print["FAIL: ", label, " :: ", fmt[res]]; Exit[1],
    Print["OK: ", label, " residual = ", fmt[res]]
  ]
];

(* --- M1: half-line integration-by-parts recombination --- *)
ClearAll[w, u, ell, q0, q1, q2, q3, q4];
$Assumptions = ell > 0 && Element[{q0, q1, q2, q3, q4}, Reals];
Qw = q0 + q1 w + q2 w^2/2 + q3 w^3/6 + q4 w^4/24;
Wel = Exp[-w/ell]/ell;
avgQ = Assuming[$Assumptions, Integrate[Wel Qw, {w, 0, Infinity}]];
avgDQ = Assuming[$Assumptions, Integrate[Wel D[Qw, w], {w, 0, Infinity}]];
bdry = Assuming[$Assumptions, Limit[Wel Qw, w -> Infinity] - (Wel Qw /. w -> 0)];
minusWp = Assuming[$Assumptions, -Integrate[D[Wel, w] Qw, {w, 0, Infinity}]];
assertZero["M1a half-line IBP recombination", bdry + minusWp - avgDQ];
assertZero["M1b half-line dQ expansion",
  avgDQ - (q1 + ell q2 + ell^2 q3 + ell^3 q4)];
assertZero["M1c half-line Q expansion",
  avgQ - (q0 + ell q1 + ell^2 q2 + ell^3 q3 + ell^4 q4)];
unitGammaKernel = u Exp[-u];
gammaNorm = Assuming[$Assumptions, Integrate[unitGammaKernel, {u, 0, Infinity}]];
gammaMu1 = Assuming[$Assumptions, Integrate[u unitGammaKernel, {u, 0, Infinity}]];
gammaMouthKernel = (w/ell) Exp[-w/ell]/ell;
avgQGamma = Assuming[$Assumptions,
  Integrate[gammaMouthKernel Qw, {w, 0, Infinity}]];
avgQGammaLeading = Assuming[$Assumptions,
  Normal[Series[avgQGamma, {ell, 0, 1}]]];
assertZero["M1d Gamma one-sided kernel normalization", gammaNorm - 1];
assertZero["M1e Gamma one-sided first moment", gammaMu1 - 2];
assertZero["M1f Gamma mouth first-moment leading term",
  avgQGammaLeading - (q0 + ell gammaMu1 q1)];

(* --- M2: even-kernel Taylor moments (unit Gaussian) --- *)
ClearAll[u, sigma, q0, q1, q2, q3, q4, q5];
$Assumptions = Element[{q0, q1, q2, q3, q4, q5}, Reals];
polyInScaledU =
  q0 + q1 sigma u + q2 (sigma u)^2/2 + q3 (sigma u)^3/6 +
    q4 (sigma u)^4/24 + q5 (sigma u)^5/120;
unitGaussian = Exp[-u^2/2]/Sqrt[2 Pi];
avgQgauss = Assuming[$Assumptions,
  Integrate[unitGaussian (polyInScaledU /. sigma -> 1), {u, -Infinity, Infinity}]];
avgDQgauss = Assuming[$Assumptions,
  Integrate[unitGaussian D[polyInScaledU /. sigma -> 1, u], {u, -Infinity, Infinity}]];
assertZero["M2a Gaussian even-kernel Q",
  avgQgauss - (q0 + q2/2 + q4/8)];
assertZero["M2b Gaussian even-kernel dQ",
  avgDQgauss - (q1 + q3/2 + q5/8)];

(* --- M3: zero-mode effective-parameter series, mouth case --- *)
ClearAll[ell, mu0, xi, s0, s1, s2, z0, z1, z2, h0, h1, h2];
$Assumptions =
  ell > 0 && z0 != 0 && h0 != 0 &&
    Element[{mu0, xi, s0, s1, s2, z0, z1, z2, h0, h1, h2}, Reals];
mouthSource = s0 + ell s1 + ell^2 s2;
mouthGaugeWeight = h0 + ell h1 + ell^2 h2;
mouthFieldWeight = z0 + ell z1 + ell^2 z2;
muHalfSer = Assuming[$Assumptions,
  Normal[Series[mu0 mouthSource/mouthFieldWeight, {ell, 0, 1}]]];
xiHalfSer = Assuming[$Assumptions,
  Normal[Series[xi mouthFieldWeight/mouthGaugeWeight, {ell, 0, 1}]]];
assertZero["M3a mouth mu_eff to O(ell)",
  muHalfSer - (mu0 s0/z0 + ell (mu0 s1/z0 - mu0 s0 z1/z0^2))];
assertZero["M3b mouth xi_eff to O(ell)",
  xiHalfSer - (xi z0/h0 + ell (xi z1/h0 - xi z0 h1/h0^2))];

(* --- M4: zero-mode effective-parameter series, symmetric case --- *)
ClearAll[sigma, m2, mu0, xi, s0, s2, z0, z2, h0, h2];
$Assumptions =
  sigma > 0 && z0 != 0 && h0 != 0 &&
    Element[{m2, mu0, xi, s0, s2, z0, z2, h0, h2}, Reals];
symSource = s0 + m2 sigma^2 s2/2;
symFieldWeight = z0 + m2 sigma^2 z2/2;
symGaugeWeight = h0 + m2 sigma^2 h2/2;
muSymSer = Assuming[$Assumptions,
  Normal[Series[mu0 symSource/symFieldWeight, {sigma, 0, 2}]]];
xiSymSer = Assuming[$Assumptions,
  Normal[Series[xi symFieldWeight/symGaugeWeight, {sigma, 0, 2}]]];
assertZero["M4a symmetric mu_eff to O(sigma^2)",
  muSymSer - (mu0 s0/z0 + (m2 sigma^2/2) (mu0 s2/z0 - mu0 s0 z2/z0^2))];
assertZero["M4b symmetric xi_eff to O(sigma^2)",
  xiSymSer - (xi z0/h0 + (m2 sigma^2/2) (xi z2/h0 - xi z0 h2/h0^2))];

(* --- M5: Gaussian-localizer asymptotic fingerprints --- *)
ClearAll[w, sigma, ell, lam, r];
$Assumptions = sigma > 0 && ell > 0 && lam > 0 && r > 0;
observerGaussian = Exp[-w^2/(2 sigma^2)]/(Sqrt[2 Pi] sigma);
observerMouth = Exp[-w/ell]/ell;
localizer = Exp[-w^2/lam^2];
ISym = Assuming[$Assumptions,
  FullSimplify[Integrate[observerGaussian localizer, {w, -Infinity, Infinity}]]];
ISymSer = Assuming[$Assumptions, Normal[Series[ISym, {sigma, 0, 4}]]];
assertZero["M5a symmetric Gaussian asymptotic",
  ISymSer - (1 - sigma^2/lam^2 + 3 sigma^4/(2 lam^4))];
IMou = Assuming[$Assumptions,
  FullSimplify[Integrate[observerMouth localizer, {w, 0, Infinity}]]];
IMouR = Assuming[$Assumptions, FullSimplify[IMou /. ell -> 1/r]];
IMouSer = Assuming[$Assumptions,
  FullSimplify[Normal[Series[IMouR, {r, Infinity, 7}]] /. r -> 1/ell]];
assertZero["M5b mouth Gaussian asymptotic",
  IMouSer - (1 - 2 ell^2/lam^2 + 12 ell^4/lam^4 - 120 ell^6/lam^6)];

Print["STATUS: PASS"];
Exit[0];
