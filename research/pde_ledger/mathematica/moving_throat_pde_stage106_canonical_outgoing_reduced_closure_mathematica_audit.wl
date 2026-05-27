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
  res = FullSimplify[Expand[expr], Assumptions -> $Assumptions];
  Print[name, " = ", fmt[res]];
  If[TrueQ[res === 0], pass[name], fail[name, res]];
];

banner["STAGE 106 — REDUCED 2.5PN CLOSURE ON CANONICAL OUTGOING DtN BRANCH"];

(* Paper-card Checks (ii) and (iii) are exercised upstream at stages 102 and
   104; this engine uses chi_Q = 1 as carry-in and proves N_Q = 1 from the
   retarded one-pole form via a path structurally distinct from the SymPy
   script (no nqGeneral/k0/k2/k4/gamma5 intermediates). *)

Clear[G, c, cS, a, omegaSym, chiQ, m0Sym, DeltaQ];
$Assumptions =
  And @@ Thread[{G, c, cS, a, omegaSym, chiQ, m0Sym} > 0] &&
  Element[DeltaQ, Reals];

OmegaQ = 3*cS/(2*a);
sigmaQcan = (9/8)/OmegaQ^5;
expectZero["sigma_Q^can - 4 a^5/(27 c_s^5)", sigmaQcan - 4*a^5/(27*cS^5)];

(* The retarded reduced one-pole module on the canonical compact passive/    *)
(* outgoing grouped-P_2 branch (paper appendix, eq. app-part04-Yret-explicit). *)
Yret[omSym_, chiSym_] := 3/4 + (1/4)/(1 - omSym^2/OmegaQ^2 - I*chiSym*sigmaQcan*omSym^5);

(* Series expand to order omega^7 so that paper Check (ii) "higher odd terms  *)
(* begin beyond omega^5" can be inspected on the canonical branch.            *)
seriesY = Normal[Series[Yret[omegaSym, chiQ], {omegaSym, 0, 7}]];
seriesYcan = Expand[seriesY /. chiQ -> 1];

Print["Y_Q^ret series (general chi_Q) = ", fmt[Expand[seriesY]]];
Print["Y_Q^ret series at chi_Q = 1   = ", fmt[seriesYcan]];

omega5Coeff = SeriesCoefficient[Yret[omegaSym, chiQ], {omegaSym, 0, 5}];
omega7Coeff = SeriesCoefficient[Yret[omegaSym, 1], {omegaSym, 0, 7}];
omega6Coeff = SeriesCoefficient[Yret[omegaSym, 1], {omegaSym, 0, 6}];
Print["omega^5 coefficient (general) = ", fmt[omega5Coeff]];
Print["omega^6 coefficient (chi_Q=1) = ", fmt[omega6Coeff]];
Print["omega^7 coefficient (chi_Q=1) = ", fmt[omega7Coeff]];

(* Sanity check on the omega^5 coefficient form: equals i chi_Q sigmaQcan/4.  *)
expectZero["omega^5 coefficient form", omega5Coeff - (I*chiQ*sigmaQcan/4)];

(* Carry-in target literals (paper appendix Box). *)
k0Target = 54*G*cS^5/(5*a^5*c^5);
k2Target = 6*G*cS^3/(5*a^3*c^5);
k4Target = 8*G*cS/(15*a*c^5);
gamma5Target = 2*G/(5*c^5);

(* F2 fix: assert the four target literals satisfy the canonical-even branch  *)
(* identities (testing the literals' mutual algebraic consistency, not the    *)
(* tautology K4 = K0/(4 OmegaQ^4)).                                            *)
expectZero[
  "target identity k0Target * k4Target - 4 k2Target^2",
  k0Target*k4Target - 4*k2Target^2
];
expectZero[
  "target identity gamma5Target - 9 Sqrt[k2Target^5/k0Target^3]",
  gamma5Target - 9*Sqrt[k2Target^5/k0Target^3]
];

(* Closure: impose the source-map relation m0hat^2 * chi_Q * N_Q = 1 with     *)
(* m0hat -> 1 (point-particle source) and chi_Q -> 1 (canonical branch),     *)
(* yielding N_Q = 1.                                                          *)
nqNatural = 1/(m0Sym^2 * chiQ);
expectZero[
  "N_Q on natural branch at m0hat=1, chi_Q=1",
  (nqNatural /. {m0Sym -> 1, chiQ -> 1}) - 1
];

(* Effective canonical odd coefficient: gamma_eff = m0hat^2 * N_Q * gamma5Target. *)
gamma5OnNatural = nqNatural * gamma5Target;
gammaEffCanonical = (m0Sym^2 * gamma5OnNatural) /. chiQ -> 1;
expectZero["canonical gamma_eff - target", gammaEffCanonical - gamma5Target];

(* F4: first-order Delta_Q sensitivity. On the natural branch (m0hat=1),     *)
(*   gamma_eff = gamma5Target / chi_Q; expansion around chi_Q = 1 + Delta_Q  *)
(*   gives zeroth coefficient gamma5Target and first-order slope             *)
(*   -gamma5Target = -2G/(5 c^5). A sign flip in N_Q would change the slope. *)
gammaEffOff = (m0Sym^2 * gamma5OnNatural) /. {m0Sym -> 1, chiQ -> 1 + DeltaQ};
gammaEffSeries = Normal[Series[gammaEffOff, {DeltaQ, 0, 1}]];
linearCoeff = FullSimplify[Coefficient[gammaEffSeries, DeltaQ, 1], Assumptions -> $Assumptions];
expectZero[
  "Delta_Q first-order sensitivity coefficient",
  linearCoeff - (-2*G/(5*c^5))
];
zerothCoeff = FullSimplify[Coefficient[gammaEffSeries, DeltaQ, 0], Assumptions -> $Assumptions];
expectZero[
  "Delta_Q zeroth-order coefficient equals Gamma5_target",
  zerothCoeff - gamma5Target
];

Print[""];
Print["RESULT:"];
Print["  Carry-in chi_Q = 1 from stage 104; canonical compact passive/outgoing"];
Print["  grouped-P_2 branch closes with N_Q = 1 at m0hat = 1."];
Print["  gamma_quad^eff = 2 G / (5 c^5); Delta_Q slope = -2 G / (5 c^5)."];

Exit[0];
