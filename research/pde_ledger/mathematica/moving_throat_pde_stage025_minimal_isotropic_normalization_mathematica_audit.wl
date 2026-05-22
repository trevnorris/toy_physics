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
samplePoint = {k -> 2, varpi -> 1, cCoupling -> 1, omegaU -> 2, omegaW -> 2,
               r -> 1, gU -> 1, gW -> 1};

zeroFrequencyCoefficients[] := Module[{delta, q, p, b0, z0, n0, d0},
  banner["Section I: Zero-Frequency Coefficients"];
  delta = Factor[omegaU^2*omegaW^2 - r^2];
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

normalizationFormula[] := Module[{delta, q, p, b0, z0, n0, d0, p0, p0Combined, p0Compact, p0Raw, p0Compactv, target, residual, mhatSq, mhatSqSample, solvability},
  banner["Section II: Minimal Isotropic Normalization Formula"];
  {delta, q, p, b0, z0, n0, d0} = zeroFrequencyCoefficients[];
  p0 = FullSimplify[n0/d0, Assumptions -> $Assumptions];
  p0Combined = Together[n0/d0];
  p0Compact = Apart[p0Combined, k];

  subbanner["II.1: raw P0 versus Mathematica-canonical form"];
  Print["P0 raw    = ", fmt[p0]];
  Print["P0 compact= ", fmt[p0Compact]];
  expectZero["P0 - P0_compact", p0 - p0Compact];

  p0Raw = FullSimplify[p0 /. samplePoint];
  p0Compactv = FullSimplify[p0Compact /. samplePoint];
  Print["P0 raw at sample     = ", fmt[p0Raw]];
  Print["P0 compact at sample = ", fmt[p0Compactv]];
  If[p0Raw =!= 1/3, fail["P0 raw at sample point != 1/3", p0Raw]];
  If[p0Compactv =!= 1/3, fail["P0 compact at sample point != 1/3", p0Compactv]];
  pass["P0 numerical at sample point"];

  (*
    target coefficient 54/5 is carried forward from
    scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:321-342,
    where Gamma5_port = a^5/(27*cs^5) and gamma_GR = 2*gConst/(5*cSpeed^5)
    imply mhat^2*P0 = gamma_GR/Gamma5_port = 54*gConst*cs^5/(5*a^5*cSpeed^5).
  *)
  target = FullSimplify[54*gConst*cs^5/(5*a^5*cSpeed^5), Assumptions -> $Assumptions];
  residual = FullSimplify[mhat^2*p0Compact - target, Assumptions -> $Assumptions];

  subbanner["II.2: target equation and positive solution"];
  Print["Target residual = ", fmt[residual]];
  mhatSq = FullSimplify[target / p0Compact, Assumptions -> $Assumptions];
  mhatSqSample = mhatSq /. Join[samplePoint, {gConst -> 1, cs -> 1, a -> 1, cSpeed -> 1}];
  Print["mhat^2 on sample = ", fmt[mhatSqSample]];
  If[!TrueQ[mhatSqSample > 0], fail["mhat^2 on sample is not positive", mhatSqSample]];
  pass["II.2 target equation solvability at sample point"];
  solvability = Reduce[
    mhat^2 == target/p0Compact && mhat > 0 && delta > 0 && d0 > 0,
    mhat, Reals
  ];
  Print["II.2 solvability (Reduce form) = ", fmt[solvability]];
  If[solvability === False, fail["II.2 target equation unsolvable for mhat > 0"]];
  pass["II.2 target equation solvable for mhat > 0"];
];

stabilityAndPositivity[] := Module[{delta, q, p, b0, z0, n0, d0, compactDenom, deltaVal, d0Val, p0PosVal},
  banner["Section III: Stability Sample Positivity"];
  {delta, q, p, b0, z0, n0, d0} = zeroFrequencyCoefficients[];
  compactDenom = FullSimplify[delta*d0, Assumptions -> $Assumptions];
  expectZero["Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)", compactDenom - (k*delta - delta*cCoupling^2/varpi^2 - q)];
  expectZero["N0 - P^2/Delta^2", n0 - p^2/delta^2];
  deltaVal = FullSimplify[delta /. samplePoint];
  d0Val = FullSimplify[d0 /. samplePoint];
  p0PosVal = FullSimplify[(n0/d0) /. samplePoint];
  Print["Delta on sample = ", fmt[deltaVal]];
  Print["D0    on sample = ", fmt[d0Val]];
  Print["P0    on sample = ", fmt[p0PosVal]];
  If[!TrueQ[deltaVal > 0], fail["Delta on sample is not positive", deltaVal]];
  If[!TrueQ[d0Val > 0], fail["D0 on sample is not positive", d0Val]];
  If[!TrueQ[p0PosVal > 0], fail["P0 on sample is not positive", p0PosVal]];
  pass["Stability positivity on sample point"];
  Print["If Delta > 0 and D0 > 0, then P0 > 0 whenever P != 0."];
];

monotonicDerivatives[] := Module[{delta, q, p, b0, z0, n0, d0, p0, h, dP0dK, dP0dX, dP0dKDirect, dP0dXDirect, sampleIV, dP0dKVal, dP0dXVal},
  banner["Section IV: K and X Slope Signs"];
  {delta, q, p, b0, z0, n0, d0} = zeroFrequencyCoefficients[];
  p0 = FullSimplify[n0/(k - x - q/delta), Assumptions -> $Assumptions];
  dP0dK = FullSimplify[Limit[((p0 /. k -> k + h) - p0)/h, h -> 0], Assumptions -> $Assumptions];
  dP0dX = FullSimplify[Limit[((p0 /. x -> x + h) - p0)/h, h -> 0], Assumptions -> $Assumptions];
  dP0dKDirect = FullSimplify[D[p0, k], Assumptions -> $Assumptions];
  dP0dXDirect = FullSimplify[D[p0, x], Assumptions -> $Assumptions];

  subbanner["IV.1: quotient-limit derivatives in K and X"];
  Print["dP0/dK = ", fmt[dP0dK]];
  Print["dP0/dX = ", fmt[dP0dX]];

  expectZero["Limit dP0/dK - D[p0,k]", dP0dK - dP0dKDirect];
  expectZero["Limit dP0/dX - D[p0,x]", dP0dX - dP0dXDirect];
  expectZero["dP0/dK + N0/(K - X - Q/Delta)^2", dP0dK + n0/(k - x - q/delta)^2];
  expectZero["dP0/dX - N0/(K - X - Q/Delta)^2", dP0dX - n0/(k - x - q/delta)^2];
  expectZero["dP0/dX + dP0/dK", dP0dX + dP0dK];
  sampleIV = Join[samplePoint, {x -> 1}];
  dP0dKVal = FullSimplify[dP0dK /. sampleIV];
  dP0dXVal = FullSimplify[dP0dX /. sampleIV];
  Print["dP0/dK on sample = ", fmt[dP0dKVal]];
  Print["dP0/dX on sample = ", fmt[dP0dXVal]];
  If[!TrueQ[dP0dKVal < 0], fail["dP0/dK on sample is not negative", dP0dKVal]];
  If[!TrueQ[dP0dXVal > 0], fail["dP0/dX on sample is not positive", dP0dXVal]];
  If[dP0dKVal + dP0dXVal =!= 0, fail["dP0/dK + dP0/dX on sample is not zero", dP0dKVal + dP0dXVal]];
  pass["Monotonic sign checks on sample point"];
];

normalizationFormula[];
stabilityAndPositivity[];
monotonicDerivatives[];

Print[""];
Print["Stage 8 Mathematica audit passed."];

Exit[0];
