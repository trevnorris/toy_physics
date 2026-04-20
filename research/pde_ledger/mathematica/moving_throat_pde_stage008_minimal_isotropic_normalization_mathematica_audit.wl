ClearAll["Global`*"];
$HistoryLength = 0;

banner[title_String] := (
  Print[""];
  Print[StringRepeat["=", 88]];
  Print[title];
  Print[StringRepeat["=", 88]];
);

subbanner[title_String] := (
  Print[""];
  Print[StringRepeat["-", 88]];
  Print[title];
  Print[StringRepeat["-", 88]];
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

Clear[k, varpi, cCoupling, omegaU, omegaW, gU, gW, r, gConst, cSpeed, cs, a, mhat, x];
$Assumptions = Element[{k, varpi, cCoupling, omegaU, omegaW, gU, gW, r, gConst, cSpeed, cs, a, mhat, x}, Reals] &&
  k > 0 && varpi > 0 && omegaU > 0 && omegaW > 0 && gConst > 0 && cSpeed > 0 && cs > 0 && a > 0 && mhat > 0 && x >= 0;

zeroFrequencyCoefficients[] := Module[{delta, q, p, b0, z0, n0, d0},
  banner["SECTION I — MINIMAL ISOTROPIC ZERO-FREQUENCY COEFFICIENTS"];
  delta = FullSimplify[omegaU^2*omegaW^2 - r^2, Assumptions -> $Assumptions];
  q = FullSimplify[gU^2*omegaW^2 + 2*gU*gW*r + gW^2*omegaU^2, Assumptions -> $Assumptions];
  p = FullSimplify[omegaU^2*gW + r*gU, Assumptions -> $Assumptions];
  b0 = FullSimplify[cCoupling^2/varpi^2, Assumptions -> $Assumptions];
  z0 = FullSimplify[q/delta, Assumptions -> $Assumptions];
  n0 = FullSimplify[p^2/delta^2, Assumptions -> $Assumptions];
  d0 = FullSimplify[k - b0 - z0, Assumptions -> $Assumptions];
  Print["Delta = ", fmt[delta]];
  Print["Q     = ", fmt[q]];
  Print["P     = ", fmt[p]];
  Print["B0    = ", fmt[b0]];
  Print["Z0    = ", fmt[z0]];
  Print["N0    = ", fmt[n0]];
  Print["D0    = ", fmt[d0]];
  {delta, q, p, b0, z0, n0, d0}
];

normalizationFormula[] := Module[{delta, q, p, b0, z0, n0, d0, p0, p0Compact, target, residual},
  banner["SECTION II — EXACT MINIMAL ISOTROPIC NORMALIZATION FORMULA"];
  {delta, q, p, b0, z0, n0, d0} = zeroFrequencyCoefficients[];
  p0 = FullSimplify[n0/d0, Assumptions -> $Assumptions];
  p0Compact = FullSimplify[p^2/(delta*(k*delta - delta*cCoupling^2/varpi^2 - q)), Assumptions -> $Assumptions];

  subbanner["II.1 — P0 in raw and compact form"];
  Print["P0 raw    = ", fmt[p0]];
  Print["P0 compact= ", fmt[p0Compact]];
  expectZero["P0 - P0_compact", p0 - p0Compact];

  target = FullSimplify[54*gConst*cs^5/(5*a^5*cSpeed^5), Assumptions -> $Assumptions];
  residual = FullSimplify[mhat^2*p0Compact - target, Assumptions -> $Assumptions];

  subbanner["II.2 — Exact target equation"];
  Print["Target residual = ", fmt[residual]];
];

stabilityAndPositivity[] := Module[{delta, q, p, b0, z0, n0, d0, compactDenom},
  banner["SECTION III — STABILITY AND POSITIVITY STRUCTURE"];
  {delta, q, p, b0, z0, n0, d0} = zeroFrequencyCoefficients[];
  compactDenom = FullSimplify[delta*d0, Assumptions -> $Assumptions];
  expectZero["Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", compactDenom - (k*delta - delta*cCoupling^2/varpi^2 - q)];
  expectZero["N0 - P^2/Delta^2", n0 - p^2/delta^2];
  Print["If Delta > 0 and D0 > 0, then P0 > 0 whenever P != 0."];
];

monotonicDerivatives[] := Module[{delta, q, p, b0, z0, n0, d0, p0, dP0dK, dP0dX},
  banner["SECTION IV — EXACT MONOTONIC DERIVATIVES"];
  {delta, q, p, b0, z0, n0, d0} = zeroFrequencyCoefficients[];
  p0 = FullSimplify[n0/(k - x - q/delta), Assumptions -> $Assumptions];
  dP0dK = FullSimplify[D[p0, k], Assumptions -> $Assumptions];
  dP0dX = FullSimplify[D[p0, x], Assumptions -> $Assumptions];

  subbanner["IV.1 — Derivatives with respect to K and X = C^2/varpi^2"];
  Print["dP0/dK = ", fmt[dP0dK]];
  Print["dP0/dX = ", fmt[dP0dX]];

  expectZero["dP0/dK + N0/(K - X - Q/Delta)^2", dP0dK + n0/(k - x - q/delta)^2];
  expectZero["dP0/dX - N0/(K - X - Q/Delta)^2", dP0dX - n0/(k - x - q/delta)^2];
  expectZero["dP0/dX + dP0/dK", dP0dX + dP0dK];
];

normalizationFormula[];
stabilityAndPositivity[];
monotonicDerivatives[];

Print[""];
Print["Stage 8 Mathematica audit passed."];

Exit[0];
