ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

pass[name_String] := Print["PASS: ", name];
fmt[expr_] := ToString[InputForm[expr]];

fail[name_String, detail_: Missing["NotAvailable"]] := (
  Print["FAIL: ", name];
  If[!MissingQ[detail], Print["  residual -> ", fmt[detail]]];
  Exit[1];
);

expectZero[name_String, expr_] := Module[{res},
  res = FullSimplify[Together[Expand[expr]], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["I. TANH-WALL SHELL MOMENTS"];

Clear[t];
$Assumptions = Element[t, Reals];

iF = FullSimplify[Integrate[(1/4)*(1 - t^2), {t, -1, 1}], Assumptions -> $Assumptions];
iG = FullSimplify[Integrate[t^2*(1 - t^2), {t, -1, 1}], Assumptions -> $Assumptions];
Print["I_f = ", fmt[iF]];
Print["I_g = ", fmt[iG]];
expectZero["I_f - 1/3", iF - 1/3];
expectZero["I_g - 4/15", iG - 4/15];

banner["II. D/N HALF-WAVE NORMALIZATION DATA"];

Clear[zTube, lW];
$Assumptions = Element[{zTube, lW}, Reals] && lW > 0;

chi = Sqrt[2/lW]*Sin[Pi*zTube/(2*lW)];
chiNorm = FullSimplify[Integrate[chi^2, {zTube, 0, lW}], Assumptions -> $Assumptions];
chiGrad = FullSimplify[Integrate[D[chi, zTube]^2, {zTube, 0, lW}], Assumptions -> $Assumptions];
chiInt = FullSimplify[Integrate[chi, {zTube, 0, lW}], Assumptions -> $Assumptions];
chiPrime0 = FullSimplify[(D[chi, zTube] /. zTube -> 0), Assumptions -> $Assumptions];

Print["Integral chi^2 dz = ", fmt[chiNorm]];
Print["Integral (chi')^2 dz = ", fmt[chiGrad]];
Print["Integral chi dz = ", fmt[chiInt]];
Print["chi'(0) = ", fmt[chiPrime0]];
expectZero["D/N norm check", chiNorm - 1];
expectZero["D/N stiffness check", chiGrad - Pi^2/(4*lW^2)];

banner["III. SHELL STIFFNESS K_s"];

Clear[a, ell, hWall, hbar, mPsi, rhoW, cSw];
$Assumptions =
  Element[{a, ell, hWall, hbar, mPsi, rhoW, cSw}, Reals] &&
  a > 0 && ell > 0 && hWall > 0 && hbar > 0 && mPsi > 0 && rhoW > 0 && cSw > 0;

kS = FullSimplify[4*Pi*a^2*(hWall*ell*iF + (hbar^2/(4*mPsi*rhoW))*(iG/ell)), Assumptions -> $Assumptions];
kSExpected = 4*Pi*a^2*(hWall*ell/3 + hbar^2/(15*mPsi*rhoW*ell));
Print["K_s = ", fmt[kS]];
expectZero["K_s closed form", kS - kSExpected];

healingSub = {hWall -> mPsi*cSw^2/rhoW, ell -> hbar/(2*mPsi*cSw)};
kSHeal = FullSimplify[kS /. healingSub, Assumptions -> $Assumptions];
Print["K_s on healing lock = ", fmt[kSHeal]];
expectZero["healing-lock K_s", kSHeal - 3*Pi*a^2*hbar^2/(5*mPsi*rhoW*(hbar/(2*mPsi*cSw)))];

banner["IV. BILINEAR GNLS SHELL/MIXED HYBRIDIZATION"];

Clear[rho0, varrhoS, v0, qStar, aQ, sAmp, qAmp, mMass];
$Assumptions = Element[{rho0, varrhoS, v0, qStar, aQ, sAmp, qAmp, mMass}, Reals] && mMass != 0;

expr = Expand[(1/2)*mMass*(rho0 + sAmp*varrhoS)*(v0 - (qStar/mMass)*qAmp*aQ)^2];
sqCoeff = Coefficient[Coefficient[expr, sAmp, 1], qAmp, 1];
Print["sq coefficient = ", fmt[sqCoeff]];
expectZero["bilinear sq coefficient", sqCoeff + qStar*varrhoS*v0*aQ];

banner["V. EFFECTIVE MIXED STIFFNESS AND MOUTH COUPLINGS"];

Clear[mu0, zQ, cSound, tM];
$Assumptions =
  Element[{a, ell, lW, mu0, zQ, cSound, tM, qStar, v0}, Reals] &&
  a > 0 && ell > 0 && lW > 0 && mu0 > 0 && zQ > 0 && cSound > 0 && tM > 0;

kQ = FullSimplify[(zQ/mu0)*(Pi^2*cSound^2/(4*lW^2)), Assumptions -> $Assumptions];
gQ = FullSimplify[(zQ/mu0)*chiPrime0, Assumptions -> $Assumptions];
jS = FullSimplify[4*Pi*a^2*ell*iF, Assumptions -> $Assumptions];
gS = FullSimplify[tM*jS, Assumptions -> $Assumptions];
iQ = FullSimplify[chiInt, Assumptions -> $Assumptions];
iSqUniform = FullSimplify[jS*iQ, Assumptions -> $Assumptions];
lamUniform = FullSimplify[qStar*v0*iSqUniform, Assumptions -> $Assumptions];

Print["K_q = ", fmt[kQ]];
Print["g_q = ", fmt[gQ]];
Print["J_s = ", fmt[jS]];
Print["g_s = ", fmt[gS]];
Print["I_sq (uniform-core closure) = ", fmt[iSqUniform]];
Print["lambda (uniform-core closure) = ", fmt[lamUniform]];

expectZero["g_q closed form", gQ - (zQ/mu0)*Pi/(Sqrt[2]*lW^(3/2))];
expectZero["J_s closed form", jS - 4*Pi*a^2*ell/3];
expectZero["I_q closed form", iQ - 2*Sqrt[2*lW]/Pi];
expectZero["lambda uniform closure", lamUniform - (8*Sqrt[2]*qStar*v0*a^2*ell*Sqrt[lW])/3];

Print[""];
Print["Stage 118 Mathematica audit passed."];

Exit[0];
