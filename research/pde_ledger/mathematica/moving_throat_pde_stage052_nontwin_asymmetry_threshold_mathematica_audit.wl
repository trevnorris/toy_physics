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
  (* Strip ConditionalExpression wrapper: under $Assumptions, a result of
     the form ConditionalExpression[0, cond] is identically zero on the
     declared domain.  Solve[]/Reduce[] often introduce these wrappers when
     auxiliary inequalities are nontrivial. *)
  res = res /. ConditionalExpression[e_, _] :> e;
  res = FullSimplify[res, Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 035 — NON-TWIN ASYMMETRY THRESHOLD"];

Clear[piTr, cMix, eps, kW, kPhi0, omega0];
$Assumptions =
  Element[{piTr, cMix, eps, kW, kPhi0, omega0}, Reals] &&
  piTr > 0 && cMix > 0 && 0 < eps < 1 && kW > 0 && kPhi0 > 0 && omega0 > 0;

sReq = FullSimplify[piTr/cMix, Assumptions -> $Assumptions];
(* Derive zetaReq by solving the lowest-lane support equation directly. *)
zetaSym = zetaSym;  (* fresh symbol *)
zetaSolList = Solve[(sReq - 1) - zetaSym (1 + eps (sReq - 2)) == 0, zetaSym];
If[Length[zetaSolList] =!= 1, fail["unique zeta solution"]];
zetaReq = FullSimplify[zetaSym /. First[zetaSolList], Assumptions -> $Assumptions];

Print["S_req = ", fmt[sReq]];
Print["zeta_req = ", fmt[zetaReq]];
expectZero["zeta_req at Pi_tr = C_mix", zetaReq /. piTr -> cMix];
expectZero["zeta_req at Pi_tr = 2 C_mix minus 1", (zetaReq /. piTr -> 2 cMix) - 1];

dZdPi = FullSimplify[D[zetaReq, piTr], Assumptions -> $Assumptions];
(* Alternative path: logarithmic differentiation of the rational form. *)
zetaTogether = Together[zetaReq];
numZ = Numerator[zetaTogether];
denZ = Denominator[zetaTogether];
dZdPiAlt = FullSimplify[
  zetaReq (D[numZ, piTr]/numZ - D[denZ, piTr]/denZ),
  Assumptions -> $Assumptions];
expectZero["dZdPi vs dZdPiAlt (independent path)", dZdPi - dZdPiAlt];
dZExpected = FullSimplify[cMix (1 - eps)/(cMix - eps (2 cMix - piTr))^2, Assumptions -> $Assumptions];
Print["d zeta_req / d Pi_tr = ", fmt[dZdPi]];
expectZero["dzeta_req/dPi - expected", dZdPi - dZExpected];

(* Build deltaZeta via Together-based algebra, independent of the SymPy form. *)
deltaZetaDerived = FullSimplify[Together[zetaReq - 1], Assumptions -> $Assumptions];
deltaExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(cMix - eps (2 cMix - piTr)), Assumptions -> $Assumptions];
Print["Delta_zeta = ", fmt[deltaZetaDerived]];
expectZero["Delta_zeta - expected", deltaZetaDerived - deltaExpected];

zetaPhys = FullSimplify[kW omega0^2/kPhi0, Assumptions -> $Assumptions];
omega0ReqSq = FullSimplify[zetaReq kPhi0/kW, Assumptions -> $Assumptions];
kPhi0Req = FullSimplify[kW omega0^2/zetaReq, Assumptions -> $Assumptions];

Print["zeta_0^(phys) = ", fmt[zetaPhys]];
Print["Omega_(0,req)^2 = ", fmt[omega0ReqSq]];
Print["K_(phi,0)^(req) = ", fmt[kPhi0Req]];
omegaSqSol = First[omega0Sq /. Solve[(kW omega0Sq/kPhi0) - zetaReq == 0, omega0Sq]];
expectZero["solve(zeta_phys = zeta_req) for Omega0^2 - expected",
           FullSimplify[omegaSqSol - omega0ReqSq, Assumptions -> $Assumptions]];

kPhi0Sol = First[kPhi0 /. Solve[(kW omega0^2/kPhi0) - zetaReq == 0, kPhi0]];
expectZero["solve(zeta_phys = zeta_req) for Kphi0 - expected",
           FullSimplify[kPhi0Sol - kPhi0Req, Assumptions -> $Assumptions]];

zetaTwin = FullSimplify[zetaPhys /. {omega0 -> 1, kPhi0 -> kW}, Assumptions -> $Assumptions];
omegaReqEqualStiff = FullSimplify[Sqrt[zetaReq], Assumptions -> $Assumptions];
kPhiReqEqualOv = FullSimplify[kW/zetaReq, Assumptions -> $Assumptions];
softFrac = FullSimplify[1 - kPhiReqEqualOv/kW, Assumptions -> $Assumptions];
(* Independent derivation of the softening fraction via Together[1 - 1/zetaReq]. *)
softFracDerived = FullSimplify[Together[1 - 1/zetaReq], Assumptions -> $Assumptions];
softExpected = FullSimplify[(1 - eps) (piTr - 2 cMix)/(piTr - cMix), Assumptions -> $Assumptions];
expectZero["softFrac vs Together[1 - 1/zetaReq] (independent path)", softFrac - softFracDerived];

expectZero["zeta_0^(twin) - 1", zetaTwin - 1];
Print["Required overlap boost at equal stiffness = ", fmt[omegaReqEqualStiff]];
Print["Required softened stiffness at equal overlap = ", fmt[kPhiReqEqualOv]];
Print["Exact softening fraction = ", fmt[softFrac]];
expectZero["softening fraction - expected", softFrac - softExpected];

Print[""];
Print["Stage 052 Mathematica audit passed."];

Exit[0];
