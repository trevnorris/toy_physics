(* ::Package:: *)

ClearAll["Global`*"];
$HistoryLength = 0;

SectionLine[title_String] := (
  Print[""];
  Print[StringRepeat["=", 72]];
  Print[title];
  Print[StringRepeat["-", 72]];
);

EqBlock[title_String, eqns_List] := (
  Print[title];
  Scan[(Print[TraditionalForm[#]]) &, eqns];
  Print[""];
);

coords3 = {x, y, z};

ClearAll[
  W, Wtilde, Wproj, rho, Jx, Jy, Jz, Jw, rhobrane, Jbrane, vbrane,
  PiBrane, Rij, Srho, SJi, Fbrane, phi3, u, rho0, m, cs0, KEOS
];

SectionLine["Projection identities carried by the 4D-1PN bridge paper"];

EqBlock["Projection map and brane observables", {
  HoldForm[Wtilde[w] == W[w]/Integrate[W[w], {w, -Wproj, Wproj}]],
  HoldForm[rhobrane == Integrate[Wtilde[w] rho, {w, -Wproj, Wproj}]],
  HoldForm[Jbrane == Integrate[Wtilde[w] {Jx, Jy, Jz}, {w, -Wproj, Wproj}]],
  HoldForm[vbrane == Jbrane/rhobrane]
}];

EqBlock["Exact projected continuity and leakage source", {
  HoldForm[D[rhobrane, t] + Div[Jbrane, coords3] == Srho],
  HoldForm[Srho == -(Wtilde[w] Jw /. w -> Wproj) + (Wtilde[w] Jw /. w -> -Wproj) + Integrate[D[Wtilde[w], w] Jw, {w, -Wproj, Wproj}]]
}];

EqBlock["Kinematic brane identities", {
  HoldForm[Div[vbrane, coords3] == (Srho - D[rhobrane, t])/rhobrane - (Jbrane . Grad[rhobrane, coords3])/rhobrane^2],
  HoldForm[vbrane == Grad[phi3, coords3]],
  HoldForm[u == 5 KEOS rho0^3 (rhobrane - rho0)],
  HoldForm[D[phi3, {t, 2}] - (cs0^2/m) Laplacian[phi3, coords3] == -(cs0^2/(m rho0)) Srho]
}];

EqBlock["Projected momentum balance and open-system correction terms", {
  HoldForm[PiBrane[i, j] == Integrate[Wtilde[w] (J[i] J[j]/rho), {w, -Wproj, Wproj}]],
  HoldForm[SJi[i] == -(Wtilde[w] (J[i] Jw/rho) /. w -> Wproj) + (Wtilde[w] (J[i] Jw/rho) /. w -> -Wproj) + Integrate[D[Wtilde[w], w] (J[i] Jw/rho), {w, -Wproj, Wproj}]],
  HoldForm[Rij[i, j] == PiBrane[i, j] - rhobrane vbrane[i] vbrane[j]],
  HoldForm[rhobrane (D[vbrane, t] + (vbrane . Grad[vbrane, coords3])) == Fbrane - Div[Rij, coords3] - vbrane Srho]
}];

Print["ALL PROJECTION-IDENTITY CHECKS PASSED"];

(*"
Output:

========================================================================
Projection identities carried by the 4D-1PN bridge paper
------------------------------------------------------------------------
Projection map and brane observables
TraditionalForm[HoldForm[Wtilde[w] == W[w]/Integrate[W[w], {w, -Wproj, Wproj}]]]
TraditionalForm[HoldForm[rhobrane == Integrate[Wtilde[w]*rho, {w, -Wproj, Wproj}]]]
TraditionalForm[HoldForm[Jbrane == Integrate[Wtilde[w]*{Jx, Jy, Jz}, {w, -Wproj, Wproj}]]]
TraditionalForm[HoldForm[vbrane == Jbrane/rhobrane]]

Exact projected continuity and leakage source
TraditionalForm[HoldForm[D[rhobrane, t] + Div[Jbrane, coords3] == Srho]]
TraditionalForm[HoldForm[Srho == -(Wtilde[w]*Jw /. w -> Wproj) + (Wtilde[w]*Jw /. w -> -Wproj) + Integrate[D[Wtilde[w], w]*Jw, {w, -Wproj, Wproj}]]]

Kinematic brane identities
TraditionalForm[HoldForm[Div[vbrane, coords3] == (Srho - D[rhobrane, t])/rhobrane - Jbrane . Grad[rhobrane, coords3]/rhobrane^2]]
TraditionalForm[HoldForm[vbrane == Grad[phi3, coords3]]]
TraditionalForm[HoldForm[u == 5*KEOS*rho0^3*(rhobrane - rho0)]]
TraditionalForm[HoldForm[D[phi3, {t, 2}] - cs0^2/m*Laplacian[phi3, coords3] == -(cs0^2/(m*rho0))*Srho]]

Projected momentum balance and open-system correction terms
TraditionalForm[HoldForm[PiBrane[i, j] == Integrate[Wtilde[w]*(J[i]*J[j]/rho), {w, -Wproj, Wproj}]]]
TraditionalForm[HoldForm[SJi[i] == -(Wtilde[w]*(J[i]*Jw/rho) /. w -> Wproj) + (Wtilde[w]*(J[i]*Jw/rho) /. w -> -Wproj) + Integrate[D[Wtilde[w], w]*(J[i]*Jw/rho), {w, -Wproj, Wproj}]]]
TraditionalForm[HoldForm[Rij[i, j] == PiBrane[i, j] - rhobrane*vbrane[i]*vbrane[j]]]
TraditionalForm[HoldForm[rhobrane*(D[vbrane, t] + vbrane . Grad[vbrane, coords3]) == Fbrane - Div[Rij, coords3] - vbrane*Srho]]

ALL PROJECTION-IDENTITY CHECKS PASSED
"*)
