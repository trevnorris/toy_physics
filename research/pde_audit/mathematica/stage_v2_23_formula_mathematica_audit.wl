Get[FileNameJoin[{DirectoryName[$InputFileName], "pde_v2_mathematica_common.wl"}]];

banner["Stage V2-23 formula and solver-packet Mathematica audit"];

Clear[w, D0, D2, D4, N0, N2, N4, B4, Z4, M, B2, Z2, G, cs, a, c, mhat, Sport];
$Assumptions = D0 != 0 && Element[{D0, D2, D4, N0, N2, N4, B4, Z4, M, B2, Z2}, Reals] &&
  G > 0 && cs > 0 && a > 0 && c > 0 && mhat > 0 && Sport > 0;

Dop = D0 + D2 w^2 + D4 w^4;
Num = N0 + N2 w^2 + N4 w^4;
Y = Expand[Normal[Series[D0/Dop, {w, 0, 5}]]];
Pref = Expand[Normal[Series[D0 Num/Dop^2, {w, 0, 5}]]];
u2 = responseU2[D0, D2];
u4 = responseU4[D0, D2, D4];
P0 = prefactorP0[D0, N0];
P2 = prefactorP2[D0, D2, N0, N2];
P4 = prefactorP4[D0, D2, D4, N0, N2, N4];

subbanner["Formula audit"];
checkEqual["response u2", Coefficient[Y, w, 2], u2];
checkEqual["response u4", Coefficient[Y, w, 4], u4];
checkEqual["prefactor P0", Coefficient[Pref, w, 0], P0];
checkEqual["prefactor P2", Coefficient[Pref, w, 2], P2];
checkEqual["prefactor P4", Coefficient[Pref, w, 4], P4];
A = M + B2 + Z2; Ciso = B4 + Z4;
Rpole = D0 Ciso - 3 A^2;
checkZero["one-pole equivalence", ((u4 - 4 u2^2) D0^2 /. {D2 -> -A, D4 -> -Ciso}) - Rpole];
N2const = 2 D2 N0/D0;
N4const = N0 (D2^2 + 2 D0 D4)/D0^2;
checkZero["constant P2", P2 /. N2 -> N2const];
checkZero["constant P4", P4 /. {N2 -> N2const, N4 -> N4const}];
Ptarget = 54 G cs^5/(5 a^5 c^5);
gammaEff = mhat^2 Sport Ptarget a^5/(27 cs^5);
gammaGR = 2 G/(5 c^5);
checkZero["normalization equivalence", (gammaEff /. {mhat -> 1, Sport -> 1}) - gammaGR];

subbanner["Solver residual packet"];
tolPath = artifactPath["stage_v2_23_tolerance_report.json"];
obsPath = artifactPath["stage_v2_23_observable_packet.json"];
tol = jsonImport[tolPath];
obs = jsonImport[obsPath];
checkTrue["tolerance report schema", tol["schema"] == "stage_v2_23_tolerance_report/v1"];
checkTrue["observable packet schema", obs["schema"] == "stage_v2_23_observable_packet/v1"];
checkTrue["open gate passes", TrueQ[tol["gates"]["open_gate_pass"]]];
checkTrue["stability gate passes", TrueQ[tol["gates"]["stability_gate_pass"]]];
checkTrue["target packet honestly fails", !TrueQ[tol["gates"]["target_packet_pass"]]];
checkTrue["dominant residual recorded", StringLength[tol["gates"]["failure_diagnosis"]["dominant_residual"]] > 0];

Print["All Stage V2-23 Mathematica checks passed."];
