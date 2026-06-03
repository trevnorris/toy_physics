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

banner["INDEPENDENT SHELL-INTEGRAL DERIVATION (concrete Gaussian profile)"];

(* Concrete symmetric profile: f(u) = Exp[-u^2]; constant compressibility h' = 1. *)
fProf[u_] := Exp[-u^2];
fpProf[u_] := D[fProf[u], u];   (* = -2 u Exp[-u^2] *)
hConst = 1;

(* Direct moment integrals of (f')^2 / h'. *)
j1Num = Integrate[(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];
j2Num = Integrate[xi*(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];
j3Num = Integrate[xi^2*(fpProf[xi])^2/hConst, {xi, -Infinity, Infinity}];

Print["J1_num = ", fmt[j1Num]];
Print["J2_num = ", fmt[j2Num]];
Print["J3_num = ", fmt[j3Num]];

(* Claim (3): J2 = 0 for a centred symmetric layer (parity). *)
expectZero["independent: J2 vanishes for symmetric profile", j2Num];

(* Claim (2): expanding the full shell weight gives the (1, 2, 1) coefficient pattern. *)
i1Direct = Integrate[(fpProf[xi])^2/hConst*(aSym + ellSym*xi)^2, {xi, -Infinity, Infinity}];
i1Direct = FullSimplify[i1Direct, Assumptions -> aSym > 0 && ellSym > 0];
i1Poly   = aSym^2*j1Num + 2*aSym*ellSym*j2Num + ellSym^2*j3Num;
expectZero["independent: I1 polynomial expansion matches direct integral",
           i1Direct - i1Poly];

(* Claim (6): for constant h' = H_w, J1 = I_f / H_w with I_f = integral of (f')^2. *)
ifMomDirect = Integrate[(fpProf[xi])^2, {xi, -Infinity, Infinity}];
hwSym = 1;  (* hConst above; documented numerically here *)
expectZero["independent: J1 = I_f / H_w under constant compressibility",
           j1Num - ifMomDirect/hwSym];

banner["STAGE 065 — THIN-WALL CONFINEMENT BRANCH"];

Clear[v0, ell, a, kx, j1, j2, j3, peReq, delta0, deltaInf, tx, len, kappa, hw, ifMom];
$Assumptions =
  Element[{v0, ell, a, kx, j1, j2, j3, peReq, delta0, deltaInf, tx, len, kappa, hw, ifMom}, Reals] &&
  v0 > 0 && ell > 0 && a > 0 && kx > 0 && j1 > 0 && j3 > 0 && peReq > 0 &&
  delta0 > 0 && deltaInf > 0 && tx > 0 && len > 0 && kappa > 0 && hw > 0 && ifMom > 0;

gPhi = FullSimplify[v0/ell, Assumptions -> $Assumptions];
i1 = FullSimplify[4*Pi*ell*(a^2*j1 + 2*a*ell*j2 + ell^2*j3), Assumptions -> $Assumptions];
i1Sym = FullSimplify[i1 /. j2 -> 0, Assumptions -> $Assumptions];

Print["g_phi = ", fmt[gPhi]];
Print["I_1 = ", fmt[i1]];
Print["I_1 | J_2=0 = ", fmt[i1Sym]];

gEq = FullSimplify[gPhi^2*i1/kx, Assumptions -> $Assumptions];
gEqSym = FullSimplify[gPhi^2*i1Sym/kx, Assumptions -> $Assumptions];
gEqTw = FullSimplify[4*Pi*a^2*j1*v0^2/(kx*ell), Assumptions -> $Assumptions];

Print["G_eq = ", fmt[gEq]];
Print["G_eq | J_2=0 = ", fmt[gEqSym]];
Print["G_eq^(tw) = ", fmt[gEqTw]];
expectZero[
  "thin-wall remainder after dropping O(ell/a) correction",
  gEqSym - gEqTw - 4*Pi*v0^2*ell*j3/kx
];

banner["PARENT WALL THRESHOLDS FOR V0"];

gFail = FullSimplify[peReq/(kappa*deltaInf), Assumptions -> $Assumptions];
gSuff = FullSimplify[peReq/(kappa*delta0), Assumptions -> $Assumptions];
gEqCoeff = FullSimplify[gEqTw/v0^2, Assumptions -> $Assumptions];

v0FailSq = FullSimplify[gFail/gEqCoeff, Assumptions -> $Assumptions];
v0SuffSq = FullSimplify[gSuff/gEqCoeff, Assumptions -> $Assumptions];

Print["V0_fail^2 = ", fmt[v0FailSq]];
Print["V0_suff^2 = ", fmt[v0SuffSq]];

subsKappa = {kappa -> kx*len^2/tx};
v0FailSqGeom = FullSimplify[v0FailSq /. subsKappa, Assumptions -> $Assumptions];
v0SuffSqGeom = FullSimplify[v0SuffSq /. subsKappa, Assumptions -> $Assumptions];

Print["V0_fail^2 with kappa inserted = ", fmt[v0FailSqGeom]];
Print["V0_suff^2 with kappa inserted = ", fmt[v0SuffSqGeom]];

expectZero[
  "K_X cancellation in V0_fail^2",
  v0FailSqGeom - tx*ell*peReq/(4*Pi*a^2*len^2*j1*deltaInf)
];
expectZero[
  "K_X cancellation in V0_suff^2",
  v0SuffSqGeom - tx*ell*peReq/(4*Pi*a^2*len^2*j1*delta0)
];

banner["CONSTANT-COMPRESSIBILITY WALL LAYER"];

v0FailConst = FullSimplify[v0FailSqGeom /. j1 -> ifMom/hw, Assumptions -> $Assumptions];
v0SuffConst = FullSimplify[v0SuffSqGeom /. j1 -> ifMom/hw, Assumptions -> $Assumptions];

Print["V0_fail^2 | H~const = ", fmt[v0FailConst]];
Print["V0_suff^2 | H~const = ", fmt[v0SuffConst]];
expectZero[
  "constant-H fail threshold",
  v0FailConst - hw*tx*ell*peReq/(4*Pi*a^2*len^2*ifMom*deltaInf)
];
expectZero[
  "constant-H suff threshold",
  v0SuffConst - hw*tx*ell*peReq/(4*Pi*a^2*len^2*ifMom*delta0)
];

Print[""];
Print["Stage 065 Mathematica audit passed."];

Exit[0];
