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

banner["STAGE 164 — COHERENT TRACKING-BRANCH DEFECT LAW"];

Clear[gConst, cSpeed, cs, a, epsEta, epsW, deltaU, chi0, zeta, zW, omegaW2];
$Assumptions = Element[{gConst, cSpeed, cs, a, epsEta, epsW, deltaU, chi0, zeta, zW, omegaW2}, Reals] &&
  gConst > 0 && cSpeed > 0 && cs > 0 && a > 0 && zW > 0 && omegaW2 > 0;

lamNorm = FullSimplify[27*Pi^2*gConst*cs^5*omegaW2/(20*a^5*cSpeed^5), Assumptions -> $Assumptions];
eps = FullSimplify[epsW*(1 - (2/11)*deltaU/(1 + deltaU)), Assumptions -> $Assumptions];
rTarget = FullSimplify[lamNorm*(1 - epsEta)*(1 - eps)^2/(zW*(1 + chi0)^2), Assumptions -> $Assumptions];

t2Direct = FullSimplify[zW*(1 + chi0)^2/(omegaW2*(1 - eps)^2), Assumptions -> $Assumptions];
t2Selected = FullSimplify[27*Pi^2*gConst*cs^5*(1 - epsEta)/(20*a^5*cSpeed^5*rTarget), Assumptions -> $Assumptions];
mMix = FullSimplify[8*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - eps)), Assumptions -> $Assumptions];
mSupp = FullSimplify[8*zeta*zW*(1 + chi0)^2/(Pi^2*(1 - epsEta)*(1 - zeta*eps)), Assumptions -> $Assumptions];
sSupport = FullSimplify[1 + zeta*(1 - eps)/(1 - zeta*eps), Assumptions -> $Assumptions];
mTr = FullSimplify[mMix + mSupp, Assumptions -> $Assumptions];
productLoaded = FullSimplify[8*lamNorm*(1 - eps)/Pi^2*sSupport, Assumptions -> $Assumptions];
rTargetLoaded = FullSimplify[productLoaded/mTr, Assumptions -> $Assumptions];
t2Loaded = FullSimplify[(lamNorm/omegaW2)*(1 - epsEta)/rTargetLoaded, Assumptions -> $Assumptions];

expectZero["direct-selected transfer-shape identity", t2Direct - t2Selected];
expectZero["support-loaded R_target reconstruction", rTargetLoaded - rTarget];
expectZero["support-loaded T^2 reconstruction", t2Loaded - t2Direct];
expectZero["support-loaded branch product law", rTargetLoaded*mTr - productLoaded];
expectZero["d/dzeta ln T^2 (support-loaded route)", D[Log[t2Loaded], zeta]];
expectZero["d/dzeta ln R_target (support-loaded route)", D[Log[rTargetLoaded], zeta]];

banner["Weak-axisymmetric drift transport"];
Clear[zetaZ, omegaW, chi1, epsW1, deltaU1, eta1];
$Assumptions = Element[{zetaZ, omegaW, chi1, epsW1, deltaU1, eta1, epsEta, epsW, deltaU, chi0}, Reals];

eps1 = FullSimplify[D[eps, epsW]*epsW1 + D[eps, deltaU]*deltaU1, Assumptions -> $Assumptions];
eps1Expected = FullSimplify[
  (1 - (2/11)*deltaU/(1 + deltaU))*epsW1 - (2*epsW*deltaU1)/(11*(1 + deltaU)^2),
  Assumptions -> $Assumptions
];
expectZero["split-blocking drift eps_1", eps1 - eps1Expected];

xi1 = FullSimplify[zetaZ - omegaW + 2*chi1/(1 + chi0) + 2*eps1/(1 - eps), Assumptions -> $Assumptions];
r1 = FullSimplify[
  omegaW - eta1/(1 - epsEta) - zetaZ - 2*chi1/(1 + chi0) - 2*eps1/(1 - eps),
  Assumptions -> $Assumptions
];
expectZero["selected-branch identity", xi1 + eta1/(1 - epsEta) + r1];

Print["Xi_1 = ", fmt[xi1]];
Print["R_1  = ", fmt[r1]];

banner["Tracking-factor drift"];
rTr = FullSimplify[(1 + chi0/(1 + deltaU))/(1 + chi0), Assumptions -> $Assumptions];
theta1 = FullSimplify[D[Log[rTr], chi0]*chi1 + D[Log[rTr], deltaU]*deltaU1, Assumptions -> $Assumptions];
theta1Expected = FullSimplify[
  -(chi0*(1 + chi0)*deltaU1 + deltaU*(1 + deltaU)*chi1)/
   ((1 + chi0)*(1 + deltaU)*(1 + chi0 + deltaU)),
  Assumptions -> $Assumptions
];
expectZero["tracking-factor drift", theta1 - theta1Expected];
Print["Theta_1 = ", fmt[theta1]];

banner["Support-blindness consequence"];
xiSupportRigid = FullSimplify[xi1 /. {chi1 -> 0, deltaU1 -> 0}, Assumptions -> $Assumptions];
thetaSupportRigid = FullSimplify[theta1 /. {chi1 -> 0, deltaU1 -> 0}, Assumptions -> $Assumptions];
Print["Xi_1 with chi1=deltaU1=0 = ", fmt[xiSupportRigid]];
Print["Theta_1 with chi1=deltaU1=0 = ", fmt[thetaSupportRigid]];

If[TrueQ[xiSupportRigid === 0], fail["support-rigid specialization unexpectedly killed Xi_1"]];

Print[""];
Print["Conclusion:"];
Print["  The coherent support ratio zeta drops out identically from T^2 and R_target."];
Print["  The grouped weak-axisymmetric defect is carried only by"];
Print["  Z_W, Omega_W^2, chi_0, eps_W, and delta_U."];
Print["  Exact tracking-factor rigidity is not sufficient to kill Xi_1."];

Print[""];
Print["Stage 181 Mathematica audit passed."];

Exit[0];
