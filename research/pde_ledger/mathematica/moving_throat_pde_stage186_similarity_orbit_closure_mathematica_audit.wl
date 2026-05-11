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

banner["STAGE 169 — EXACT MICROSCOPIC SIMILARITY ORBIT"];

Clear[chi, delta, eParam, fParam, lam1, c1, gam1, kU, kEta, kW, mu1, tau1];
$Assumptions = Element[{chi, delta, eParam, fParam, lam1, c1, gam1, kU, kEta, kW, mu1, tau1}, Reals] &&
  chi > 0 && delta > 0;

dlogCtr = (1 + delta)*(gam1 + c1 - kU) + (1 + chi)*(tau1 - kU);
dlogCnt = (2 + 2*eParam)*lam1 + 2*eParam*gam1 + mu1 + (fParam - eParam)*kU - kEta - (2 + eParam)*kW - fParam*tau1;
dlogEta = 2*c1 - kU - kEta;

Print["delta log C_tr = ", fmt[Expand[dlogCtr]]];
Print["delta log C_nt = ", fmt[Expand[dlogCnt]]];
Print["delta log eps_eta = ", fmt[Expand[dlogEta]]];

mMat = {
  {0, 1 + delta, 1 + delta, -(2 + chi + delta), 0, 0, 0, 1 + chi},
  {2 + 2*eParam, 0, 2*eParam, fParam - eParam, -1, -(2 + eParam), 1, -fParam},
  {0, 2, 0, -1, -1, 0, 0, 0}
};
Print[""];
Print["M_* = ", fmt[MatrixForm[mMat]]];

minor = mMat[[All, {8, 5, 7}]];
detMinor = FullSimplify[Det[minor], Assumptions -> $Assumptions];
Print[""];
Print["minor det(tau1,kEta,mu1) = ", fmt[detMinor]];
If[!TrueQ[detMinor === 1 + chi], fail["Unexpected rank-3 minor determinant", detMinor]];

banner["Linear compatibility solve"];
sol = First[Solve[{dlogCtr == 0, dlogEta == 0, dlogCnt == 0}, {tau1, kEta, mu1}]];
tauFormula = FullSimplify[tau1 /. sol, Assumptions -> $Assumptions];
kEtaFormula = FullSimplify[kEta /. sol, Assumptions -> $Assumptions];
muFormula = FullSimplify[mu1 /. sol, Assumptions -> $Assumptions];

Print["tau1 = ", fmt[tauFormula]];
Print["kEta = ", fmt[kEtaFormula]];
Print["mu1  = ", fmt[muFormula]];

expectZero[
  "tau1 - [kU - ((1+delta)/(1+chi)) (gam1+c1-kU)]",
  tauFormula - (kU - (1 + delta)*(gam1 + c1 - kU)/(1 + chi))
];
expectZero["kEta - (2 c1 - kU)", kEtaFormula - (2*c1 - kU)];
expectZero[
  "mu1 - Stage 185 form",
  muFormula - (
    (2*c1 - kU) + 2*kW - 2*lam1 -
    eParam*(2*gam1 + 2*lam1 - kU - kW) -
    fParam*(1 + delta)*(gam1 + c1 - kU)/(1 + chi)
  )
];

banner["Finite five-parameter similarity orbit"];
Clear[lamSym, cSym, gamSym, uSym, wSym];
$Assumptions = Element[{chi, delta, eParam, fParam, lamSym, cSym, gamSym, uSym, wSym}, Reals] &&
  chi > 0 && delta > 0;

etaExp = 2*cSym - uSym;
tauExp = uSym - (1 + delta)*(gamSym + cSym - uSym)/(1 + chi);
muExp = etaExp + 2*wSym - 2*lamSym - eParam*(2*gamSym + 2*lamSym - uSym - wSym) + fParam*(tauExp - uSym);

Print["K_eta exponent = ", fmt[etaExp]];
Print["T_U exponent   = ", fmt[tauExp]];
Print["mu_W exponent  = ", fmt[FullSimplify[muExp, Assumptions -> $Assumptions]]];

ctrOrbit = (1 + delta)*(gamSym + cSym - uSym) + (1 + chi)*(tauExp - uSym);
cntOrbit = (2 + 2*eParam)*lamSym + 2*eParam*gamSym + muExp + (fParam - eParam)*uSym - etaExp - (2 + eParam)*wSym - fParam*tauExp;
etaOrbit = 2*cSym - uSym - etaExp;

expectZero["finite orbit preserves C_tr", ctrOrbit];
expectZero["finite orbit preserves C_nt", cntOrbit];
expectZero["finite orbit preserves eps_eta", etaOrbit];

banner["Linearization reproduces compatibility ledger"];
basisSubs = {lamSym -> lam1, cSym -> c1, gamSym -> gam1, uSym -> kU, wSym -> kW};
expectZero["linearized tau formula", (tauExp /. basisSubs) - tauFormula];
expectZero["linearized kEta formula", (etaExp /. basisSubs) - kEtaFormula];
expectZero["linearized mu formula", (muExp /. basisSubs) - muFormula];

Print[""];
Print["Conclusion:"];
Print["  The three Stage 185 compatibility equations are exactly the tangent-space"];
Print["  equations of a finite five-parameter multiplicative similarity orbit."];
Print["  The coherent weak-axisymmetric zero-defect branch is therefore codimension 3,"];
Print["  not an isolated fine-tuning point."];

Print[""];
Print["Stage 186 Mathematica audit passed."];

Exit[0];
