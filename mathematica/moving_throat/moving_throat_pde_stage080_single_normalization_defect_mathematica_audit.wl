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

banner["STAGE 080 — SINGLE NORMALIZATION DEFECT"];

Clear[gConst, cLight, cSound, aRad, k0, omegaQ, nQ];
$Assumptions =
  Element[{gConst, cLight, cSound, aRad, k0, omegaQ, nQ}, Reals] &&
  gConst > 0 && cLight > 0 && cSound > 0 && aRad > 0 && k0 > 0 && omegaQ > 0 && nQ > 0;

k2 = FullSimplify[k0/(4*omegaQ^2), Assumptions -> $Assumptions];
k4 = FullSimplify[k0/(4*omegaQ^4), Assumptions -> $Assumptions];
gamma5 = FullSimplify[9*k2^(5/2)/k0^(3/2), Assumptions -> $Assumptions];
gamma5Expected = FullSimplify[9*k0/(32*omegaQ^5), Assumptions -> $Assumptions];

Print["K2 = ", fmt[k2]];
Print["K4 = ", fmt[k4]];
Print["Gamma5 = ", fmt[gamma5Expected]];
expectZero["Gamma5 - 9 K0/(32 Omega^5)", gamma5 - gamma5Expected];

k0Target = FullSimplify[64*gConst*omegaQ^5/(45*cLight^5), Assumptions -> $Assumptions];
omegaGeom = FullSimplify[3*cSound/(2*aRad), Assumptions -> $Assumptions];
k0TargetGeom = FullSimplify[k0Target /. omegaQ -> omegaGeom, Assumptions -> $Assumptions];

Print["K0_target = ", fmt[k0Target]];
Print["K0_target (Omega=3 c_s/(2a)) = ", fmt[k0TargetGeom]];
expectZero["geometric target reduction", k0TargetGeom - 54*gConst*cSound^5/(5*aRad^5*cLight^5)];

k2Target = FullSimplify[k0Target/(4*omegaQ^2), Assumptions -> $Assumptions];
k4Target = FullSimplify[k0Target/(4*omegaQ^4), Assumptions -> $Assumptions];
gamma5Target = FullSimplify[9*k2Target^(5/2)/k0Target^(3/2), Assumptions -> $Assumptions];
expectZero["Gamma5_target - 2G/(5c^5)", gamma5Target - 2*gConst/(5*cLight^5)];

r0 = FullSimplify[(nQ*k0Target)/k0Target - 1, Assumptions -> $Assumptions];
r2 = FullSimplify[((k2 /. k0 -> nQ*k0Target)/k2Target) - 1, Assumptions -> $Assumptions];
r4 = FullSimplify[((k4 /. k0 -> nQ*k0Target)/k4Target) - 1, Assumptions -> $Assumptions];
r5 = FullSimplify[((gamma5 /. k0 -> nQ*k0Target)/gamma5Target) - 1, Assumptions -> $Assumptions];

Print["R0 = ", fmt[Factor[r0]]];
Print["R2 = ", fmt[Factor[r2]]];
Print["R4 = ", fmt[Factor[r4]]];
Print["R5 = ", fmt[Factor[r5]]];

expectZero["R0 - (N_Q - 1)", r0 - (nQ - 1)];
expectZero["R2 - (N_Q - 1)", r2 - (nQ - 1)];
expectZero["R4 - (N_Q - 1)", r4 - (nQ - 1)];
expectZero["R5 - (N_Q - 1)", r5 - (nQ - 1)];

Print[""];
Print["Stage 080 Mathematica audit passed."];

Exit[0];
