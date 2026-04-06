(* ::Package:: *)

ClearAll[
  section, pass, fail,
  mu, P2, xi, chi, DeltaSChi, tau0, tau1, pi0, pi1, pi2,
  zBase, tProfile, zTarget, tTarget, sourceTarget, sourceFromZ,
  s0, s1, s2,
  K1perp, K10sym, a0, a2, zFromDipole, pIso, qAx, J0, J20,
  densityLM, Kchannel, K00raw, K11, K10, K20, K21, K22, deltaK00
];

section[s_String] := Print["\n" <> StringRepeat["=", 78] <> "\n" <> s <> "\n" <> StringRepeat["=", 78]];
pass[s_String] := Print["PASS: " <> s];
fail[s_String] := Print["FAIL: " <> s];

section["1) Carried-forward local Family-1 support/source data"];

P2 = LegendreP[2, mu];

xi = 1176/1553;
chi = Expand[1 - xi mu^2 + 1/2 xi^2 mu^4];
DeltaSChi = Expand[D[(1 - mu^2) D[chi, mu], mu]];

tau0 = -1714441/790272;
tau1 = 2411809/790272;

pi0 = 17671223/241790304;
pi1 = 60295225/241790304;
pi2 = 12059045/967161216;

zBase = Expand[pi0 + pi1 chi + pi2 DeltaSChi];
tProfile = Expand[tau0 + tau1 chi];

zTarget = Expand[17/56 - 5 mu^2/56];
tTarget = Expand[593/672 - 1553 mu^2/672 + 7 mu^4/8];

Print["chi(mu)       = ", chi];
Print["Delta chi(mu) = ", DeltaSChi];
Print[""];
Print["zBase(mu)     = ", zBase];
Print["t(mu)         = ", tProfile];
Print[""];
Print["z residual    = ", FullSimplify[zBase - zTarget]];
Print["t residual    = ", FullSimplify[tProfile - tTarget]];

If[FullSimplify[zBase - zTarget] === 0 && FullSimplify[tProfile - tTarget] === 0,
  pass["Carried-forward Family-1 support profiles reproduced exactly."],
  fail["Support profile reconstruction failed."]
];

section["2) Exact source-from-support closure"];

sourceTarget = Expand[7/16 + 45 mu^2/16];
sourceFromZ = Expand[10 - 63 zBase/2];

s0 = 177262811/23027648;
s1 = -180885675/23027648;
s2 = -36177135/92110592;

Print["S_target(mu)  = ", sourceTarget];
Print["S_from_z(mu)  = ", sourceFromZ];
Print[""];
Print["source residual = ", FullSimplify[sourceFromZ - sourceTarget]];
Print[""];
Print["Exact flare-basis coefficient relation:"];
Print["  s0 - [10 - 63 pi0/2] = ", FullSimplify[s0 - (10 - 63 pi0/2)]];
Print["  s1 + 63 pi1/2        = ", FullSimplify[s1 + 63 pi1/2]];
Print["  s2 + 63 pi2/2        = ", FullSimplify[s2 + 63 pi2/2]];

If[
  FullSimplify[sourceFromZ - sourceTarget] === 0 &&
  FullSimplify[s0 - (10 - 63 pi0/2)] === 0 &&
  FullSimplify[s1 + 63 pi1/2] === 0 &&
  FullSimplify[s2 + 63 pi2/2] === 0,
  pass["Source sector is fixed exactly by the base support profile."],
  fail["Source-from-support closure failed."]
];

section["3) Dipole support residues determine the scalar source sector"];

a0 = FullSimplify[(K10sym + 2 K1perp)/3];
a2 = FullSimplify[5 (K10sym - K1perp)/3];
zFromDipole = Expand[a0 + a2 P2];

pIso = FullSimplify[10 - 63 a0/2];
qAx = FullSimplify[-63 a2/2];

Print["zBase(mu) reconstructed from dipole support residues:"];
Print["  a0 = ", a0];
Print["  a2 = ", a2];
Print["  zFromDipole(mu) = ", zFromDipole];
Print[""];
Print["Wall/source coefficients determined by dipole support:"];
Print["  pIso = ", pIso];
Print["  qAx  = ", qAx];
Print[""];
Print["At the passed dipole support values K1perp=2/7, K10=1/4:"];
Print["  z residual = ", FullSimplify[(zFromDipole /. {K1perp -> 2/7, K10sym -> 1/4}) - zTarget]];
Print["  pIso residual = ", FullSimplify[(pIso /. {K1perp -> 2/7, K10sym -> 1/4}) - 11/8]];
Print["  qAx residual  = ", FullSimplify[(qAx /. {K1perp -> 2/7, K10sym -> 1/4}) - 15/8]];

If[
  FullSimplify[(zFromDipole /. {K1perp -> 2/7, K10sym -> 1/4}) - zTarget] === 0 &&
  FullSimplify[(pIso /. {K1perp -> 2/7, K10sym -> 1/4}) - 11/8] === 0 &&
  FullSimplify[(qAx /. {K1perp -> 2/7, K10sym -> 1/4}) - 15/8] === 0,
  pass["Dipole support residues reconstruct the exact base and source profiles."],
  fail["Dipole-to-source reconstruction failed."]
];

section["4) Canonical scalar source vector from the dipole support data"];

J0 = FullSimplify[(2/Sqrt[5]) (pIso + qAx/3)];
J20 = FullSimplify[(2/3) qAx];

Print["Canonical mapping:"];
Print["  J0  = (2/Sqrt[5]) (pIso + qAx/3)"];
Print["  J20 = (2/3) qAx"];
Print[""];
Print["J0(K1perp,K10)  = ", J0];
Print["J20(K1perp,K10) = ", J20];
Print[""];
Print["At the passed dipole support values:"];
Print["  J0 residual  = ", FullSimplify[(J0 /. {K1perp -> 2/7, K10sym -> 1/4}) - 4/Sqrt[5]]];
Print["  J20 residual = ", FullSimplify[(J20 /. {K1perp -> 2/7, K10sym -> 1/4}) - 5/4]];

If[
  FullSimplify[(J0 /. {K1perp -> 2/7, K10sym -> 1/4}) - 4/Sqrt[5]] === 0 &&
  FullSimplify[(J20 /. {K1perp -> 2/7, K10sym -> 1/4}) - 5/4] === 0,
  pass["Canonical scalar source vector is fixed by the dipole support residues."],
  fail["Canonical source-vector reconstruction failed."]
];

section["5) Full static local wall operator + monopole auxiliary completion"];

densityLM[ell_Integer, m_Integer] := Module[{P, norm},
  P = LegendreP[ell, m, mu];
  norm = Rationalize[(2 ell + 1)/2 * Factorial[ell - m]/Factorial[ell + m], 0];
  FullSimplify[norm P^2]
];

Kchannel[ell_Integer, m_Integer] := Module[{w},
  w = densityLM[ell, m];
  FullSimplify[
    Integrate[zBase w, {mu, -1, 1}] +
    (ell (ell + 1) - 2) Integrate[tProfile w, {mu, -1, 1}]
  ]
];

K00raw = Kchannel[0, 0];
K11 = Kchannel[1, 1];
K10 = Kchannel[1, 0];
K20 = Kchannel[2, 0];
K21 = Kchannel[2, 1];
K22 = Kchannel[2, 2];

deltaK00 = FullSimplify[4/45 - K00raw];

Print["Local Family-1 wall operator:"];
Print["  O_loc = zBase(mu) + 1/2 { -Delta_S - 2, t(mu) }"];
Print[""];
Print["Channel values from O_loc:"];
Print["  K00 raw = ", K00raw];
Print["  K11     = ", K11];
Print["  K10     = ", K10];
Print["  K20     = ", K20];
Print["  K21     = ", K21];
Print["  K22     = ", K22];
Print[""];
Print["Monopole add-on needed to reach the carried-forward target 4/45:"];
Print["  deltaK00 = ", deltaK00];
Print[""];
Print["Minimal full static even-wall completion:"];
Print["  O_full = O_loc + deltaK00 * P_00"];
Print["with P_00 the monopole projector."];
Print[""];
Print["Equivalent auxiliary breathing-channel realizations:"];
Print["  residue-normalized pole:   Y_mono(omega) = ", deltaK00, "/(1 - omega^2/Omega_mono^2)"];
Print["  stiffness-normalized aux:  E[b,Psi00] = 1/2 b^2 - Sqrt[", deltaK00, "] b Psi00"];

If[
  K00raw === -757/2520 &&
  K11 === 2/7 &&
  K10 === 1/4 &&
  K20 === 4/9 &&
  K21 === 2/3 &&
  K22 === 8/3 &&
  deltaK00 === 109/280,
  pass["Local wall operator and monopole auxiliary completion are exact."],
  fail["Static local wall operator or monopole completion failed."]
];

section["6) Compact closure statement"];

Print["Exact static Family-1 even-sector closure:"];
Print["  1) local support operator"];
Print["       O_loc = zBase + 1/2{ -Delta_S - 2, t }"];
Print["  2) source fixed by the base support profile"];
Print["       S = 10 - (63/2) zBase"];
Print["  3) one global monopole breathing add-on"];
Print["       deltaK00 = ", deltaK00];
Print[""];
Print["So the entire static even wall sector is:"];
Print["  [local Family-1 support/source constitutive law] + [one monopole auxiliary channel]"];
