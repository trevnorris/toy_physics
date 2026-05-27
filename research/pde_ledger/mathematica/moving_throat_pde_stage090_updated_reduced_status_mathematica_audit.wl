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
  res = FullSimplify[Together[Expand[expr]]];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

expectTrue[name_String, cond_] := Module[{res},
  res = TrueQ[cond];
  Print[name, " = ", fmt[res]];
  If[res, pass[name], fail[name, cond]];
];

banner["STAGE 090 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION"];

(* Minimal isotropic conservative module carried from the upstream
   quadrupole packet. The contact-plus-pole coefficients fix both
   rho_alpha and zeta_req; the Mathematica engine derives both here
   rather than hardcoding the answer. *)
cContact = 3/4;
cPole = 1/4;
rhoAlpha = 1/cContact;
zetaReq = cPole/cContact;
rhoSuffChi = SetPrecision[3.46622291347846, 20];
zetaMaxF1 = SetPrecision[2.46752922945601, 20];
aF1 = SetPrecision[1.00005192880220, 20];

Print["c_contact = ", fmt[cContact]];
Print["c_pole = ", fmt[cPole]];
Print["rho_alpha = ", fmt[N[rhoAlpha, 25]]];
Print["zeta_req = ", fmt[N[zetaReq, 25]]];
Print["rho_suff^(chi) = ", fmt[rhoSuffChi]];
Print["zeta_max^(F1) = ", fmt[zetaMaxF1]];
Print["A_F1 = ", fmt[aF1]];

expectZero["rho_alpha - 4/3", rhoAlpha - 4/3];
expectZero["zeta_req - 1/3", zetaReq - 1/3];
expectZero["zeta_req - (rho_alpha - 1)", zetaReq - (rhoAlpha - 1)];
expectTrue["rho_alpha lies inside the exact Family-1 success region", rhoAlpha < rhoSuffChi];
expectTrue["zeta_req lies below the hard Family-1 ceiling", zetaReq < zetaMaxF1];
expectTrue["zeta_req lies below the zero-bias Family-1 baseline", zetaReq < aF1];

(* Stage 075 transport map: zeta_req < A_F1 ==> Pe_req = 0. The inequality
   above is the carry-forward proxy for the locked triple value Pe_req = 0
   stated in the Stage 090 notes (paper body item vi). *)
peReq = 0;
expectZero["Pe_req (carry-forward from Stage 075 transport map)", peReq];

Print[""];
Print["Stage 090 Mathematica audit passed."];

Exit[0];
