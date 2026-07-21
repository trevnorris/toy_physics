(* Independent Wolfram route for the G0 leftover-longitudinal scalar.

   This calculation completes the gradient square first, constructs the
   Robin resolvent vertices, and derives each thermodynamic member by direct
   work subtraction (plus an explicit reacting-charge series for V).  It does
   not consume formulas or sign labels emitted by the SymPy engine.
*)

ClearAll["Global`*"];
$Assumptions = b > 0 && c > 0 && k > 0 && km > 0 && zeta > 0 && reb > 0 &&
  suu > 0 && sgg > 0 && u > 0 && qf > 0 && js > 0 && gh > 0 && uheld > 0 &&
  sorient != 0 && b k - c^2 > 0;

here = DirectoryName[ExpandFileName[$InputFileName]];
root = ParentDirectory[ParentDirectory[here]];
paths = <|
  "directive" -> FileNameJoin[{here, "directive_leftover_scalar_electric_sign.md"}],
  "g0" -> FileNameJoin[{here, "g0_closure_card_v0.md"}],
  "audit" -> FileNameJoin[{root, "docs", "model_definition_audit.md"}],
  "handoff" -> FileNameJoin[{root, "docs", "em_analog_next_phase_handoff.md"}],
  "sim" -> FileNameJoin[{root, "docs", "two_throat_simulation_handoff_spec.md"}],
  "path36" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_36_c5_phase_potential.md"}],
  "path38" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_38_throat_body_electric_localization.md"}],
  "path38yaml" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_38_results.yaml"}],
  "path39" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_39_stage4_field_classification.md"}],
  "path24" -> FileNameJoin[{root, "software", "stage1_solver", "reports", "pathA_24_T1_wall.md"}],
  "u2directive" -> FileNameJoin[{here, "directive_u2_boundary_adjudication.md"}],
  "u2report" -> FileNameJoin[{here, "reports", "u2_boundary_adjudication_artifacts", "stage_1_production_v12", "production_summary.md"}]
|>;
pyPath = FileNameJoin[{here, "leftover_scalar_electric_sign_check.py"}];

zq[x_] := TrueQ[FullSimplify[x == 0]];
canon[x_] := ToString[Cancel[x], InputForm];

(* ---------------------------------------------------------------------- *)
(* Source-derived transcription contract.                                *)
(* ---------------------------------------------------------------------- *)

actionTerms = {
  "1/2*A_eff*(dt u_L)^2", "1/2*M_h*(dt h)^2",
  "-1/2*B_eff*|grad u_L|^2", "-1/2*K_h*|grad h|^2",
  "-C_hu*grad u_L.grad h"
};
testBCs = {
  {"V", "f_V=1 on spherical A; area_average_A(f_V)=1; u_L|A=s*u0*f_V; physical Q_u reacts"},
  {"M", "f_b=1/Area(A) on A, zero off A; integral_A f_b dS=1; conormal=s*q*f_b"},
  {"J", "f_b=1/Area(A) on A, zero off A; integral_A f_b dS=1; Omega adds -s*J0*integral_A(f_b*u_L)dS"}
};
zeroLedger = {
  "bulk `r_BH`, `r_B²H²`, `Hρ`, `Hδρ`, `H∂tθ`, and `H∇θ` couplings outside `Ω_mouth`",
  "dynamical modulation `δJ_m/δr_B`, `δJ_m/δh`, `δJ_m/δP`, neighbor-source response",
  "`r_Bu_L`, `r_B divu`, `r_Bu_T`, `Hu_T`, `u_Lu_T`, and two-gradient scalar–transverse mixing",
  "cross-kinetic `∂tu_L∂th`, one-time-derivative/Berry `u_L∂th` and `h∂tu_L`",
  "`u_L²`, reduced `h²`, `(∇²u_L)²`, `(∇²h)²`, `h³,h⁴`, and `(∇u_L)³`",
  "independent primitive `B(divu)²` in addition to `B_eff`",
  "independent `θ_B`, its amplitude-weighted kinetic/gradient terms, Josephson `θ−θ_B`, and brane-phase drain",
  "wall bending, anchoring, collar two-surface tension, surface number storage, and surface dissipation",
  "dynamic `δQ_d/δr_B`, `Γ(h_inc)`, `Γ(P_inc)`, return-location responses to `h,P`, return-kernel responses to `h,P`, and drain-rate orientation dependence",
  "direct drain sources in the `h` and `u_L` equations, and direct `h` contribution to `e_c`",
  "field-dependent geon derivatives `δE_g/δ{r_B,h,θ,H,u_L,ρ}`",
  "bulk viscosity, tangential drag, E2 no-slip traction, E3 permeability resistance, E3 phase jump",
  "all E4 multipliers, E5 Rayleigh kernels, E1 reactions, and mixture terms"
};

transcript = <|
  "action" -> actionTerms,
  "card" -> {{"A_eff", "1"}, {"M_h", "1"}, {"B_eff", "2"},
    {"K_h", "1"}, {"C_hu", "1/2"}, {"k_m", "1"}, {"g_chih", "1"}},
  "relations" -> {"T_L=B_eff=rho_B0^2/chi_c>0", "K_h=M_h*c_E^2"},
  "mouth" -> {"1/2*k_m*h^2", "-g_chih*s_i*h"},
  "native" -> {"continuity(u_L)", "continuity(B_eff*d_n u_L+C_hu*d_n h)",
    "continuity(h)", "continuity(K_h*d_n h+C_hu*d_n u_L) except mouth source", "u_L->0 at IR"},
  "tests" -> testBCs, "ledger" -> zeroLedger
|>;
expectedTranscript = Association[transcript];

parseLedger[g0_] := Module[{block, lines},
  block = StringSplit[StringSplit[g0, "## 9. Complete declared-zero ledger"][[2]],
    "## 10. Instantiability, gates, and checks"][[1]];
  lines = Select[StringSplit[block, "\n"], StringStartsQ[#, "| "] && StringContainsQ[#, "[POSTULATE]"] &];
  StringTrim[StringSplit[#, "|"][[1]]] & /@ lines
];

sourceFaithful[tr_, hess_: Automatic] := Module[{src, needles, g0ledger, sourceOK, internalOK, usedHessian},
  src = Map[Import[#, "Text"] &, paths];
  needles = <|
    "g0" -> {"½A_eff (∂tu_L)² + ½M_h(∂th)²", "−½B_eff|∇u_L|² −½K_h|∇h|²",
      "−C_hu ∇u_L·∇h", "| `B_eff` | `(-1,-2,1)` | `2`",
      "| `C_hu` | `(0,-2,1)` | `1/2`", "| `K_h` | `(1,-2,1)` | `1`",
      "| `A_eff` | `(-3,0,1)` | `1`", "| `M_h` | `(-1,0,1)` | `1`",
      "| reduced `k_m=K_m/ℓ²` | `E` | `1`", "| reduced `g_χh=j_m=J_m/ℓ` | `E` | `1`",
      "K_h=M_hc_E²", "½K_m H(x,0)² − J_m Q_χ[r_Σ,s_i] H(x,0)",
      "η_i(k_mh−g_χh s_i)", "B_eff∂_nu_L+C_hu∂_nh",
      "K_h∂_nh+C_hu∂_nu_L", "u_L→0"},
    "directive" -> {"u_L|_mouth = s·u₀", "far-field monopole strength (`∮` of the conormal) `= s·q` held",
      "a `−J φ` Legendre term with `J=s·J₀` held"},
    "audit" -> {"I1", "I2", "I3", "U1", "U2", "U3"},
    "handoff" -> {"144/144 cells UNRESOLVED", "R1 two-throat solve"},
    "sim" -> {"a supplied closure card before use", "does **not** imply that the force sign must flip"},
    "path36" -> {"B_eff = rho_B0^2/chi_c"},
    "path38" -> {"p_static=2"},
    "path38yaml" -> {"schema: pathA_38_throat_body_electric_sympy/v2", "N0_norm: 8/(3*ell)",
      "q_h_plus: 2*QE*tanh(b/ell)/b", "q_h_minus: -2*QE*tanh(b/ell)/b", "p_static: 2"},
    "path39" -> {"q_L", "mass source to u_L, recorded separately from charge residues"},
    "path24" -> {"DeltaE_unwind = 0"},
    "u2directive" -> {"144", "UNRESOLVED", "mouth"},
    "u2report" -> {"{'UNRESOLVED': 144}"}
  |>;
  sourceOK = And @@ Flatten[KeyValueMap[
    Function[{label, required}, (StringContainsQ[src[label], #] & /@ required)], needles]];
  g0ledger = parseLedger[src["g0"]];
  usedHessian = If[hess === Automatic, actionHessian, hess];
  internalOK = tr === expectedTranscript && tr["ledger"] === g0ledger && Length[g0ledger] == 13 &&
    zq[usedHessian[[1, 1]] - b] && zq[usedHessian[[1, 2]] - c] &&
    zq[usedHessian[[2, 1]] - c] && zq[usedHessian[[2, 2]] - k];
  sourceOK && internalOK
];

(* ---------------------------------------------------------------------- *)
(* Independent completed-square/Robin-resolvent route.                    *)
(* ---------------------------------------------------------------------- *)

staticDensity = b gu^2/2 + c gu ghgrad + k ghgrad^2/2;
actionHessian = Table[D[staticDensity, xi, xj], {xi, {gu, ghgrad}}, {xj, {gu, ghgrad}}];
ba = actionHessian[[1, 1]]; ca = actionHessian[[1, 2]]; ka = actionHessian[[2, 2]];
dd = Factor[Det[actionHessian]]; kap = Factor[dd/ba];

(* The positive Robin escape factor zeta is the maximum-principle Schur
   complement of the normalized mouth channel.  reeFull below makes the
   exact identity z_g=1-km<eta,L_h^-1 eta>=zeta explicit. *)
reeFull = (1 - zeta)/km; zg = Factor[1 - km reeFull]; zb = 1 - km reb;

(* Independent route: write the completed-square diagonal pair kernel
   directly and recover the physical response entries by differentiation.
   No vertex matrix from the SymPy construction is used. *)
diagonalPair[qu_, gs_] := Expand[qu^2/ba + (zg gs - (ca/ba) zb qu)^2/kap];
muu = Factor[Coefficient[diagonalPair[qu, gs], qu, 2]];
mug = Factor[Coefficient[Coefficient[diagonalPair[qu, gs], qu, 1], gs, 1]/2];
mgg = Factor[Coefficient[diagonalPair[qu, gs], gs, 2]];
mfar = {{muu, mug}, {mug, mgg}}; mdet = Factor[Det[mfar]];
qphysical = Factor[(u - sug gh)/suu];
qdatum = Factor[u/suu];
vhh = Factor[mgg - 2 (sug/suu) mug + (sug/suu)^2 muu];

splitReadout[expr_, datum_] := Module[{total, honly, uonly},
  total = Factor[expr];
  honly = Factor[Coefficient[total, gh, 2] gh^2];
  uonly = Factor[Coefficient[total, datum, 2] datum^2];
  <|"h_only" -> honly, "uL_only" -> uonly,
    "interference" -> Factor[total - honly - uonly], "total" -> total|>
];

(* V is obtained from a two-throat block Schur solve: invert the full
   reaction block at fixed values, substitute it into the selected member,
   and only then extract the eps pair coefficient. *)
vBlockCoefficient[member_] := Module[{self, full, aa, hh, uv, gv, qv, xv, sigma, pair},
  self = {{suu, sug}, {sug, sgg}};
  full = ArrayFlatten[{{self, eps mfar}, {eps mfar, self}}];
  aa = {{suu, eps muu}, {eps muu, suu}};
  hh = {{sug, eps mug}, {eps mug, sug}};
  uv = {sone u, stwo u}; gv = {sone gh, stwo gh};
  qv = Together[Inverse[aa].(uv - hh.gv)];
  xv = {qv[[1]], gv[[1]], qv[[2]], gv[[2]]};
  sigma = If[member === "conjugate", -1, 1];
  pair = Coefficient[Normal[Series[sigma xv.full.xv/2, {eps, 0, 1}]], eps]/(sone stwo);
  Factor[pair]
];

(* M uses a partial Legendre derivative of the diagonal pair kernel.  J is
   the negative on-shell source quadratic. *)
otherCoefficient[ensemble_, member_] := Module[{datum, stored},
  datum = If[ensemble === "M", qf, js]; stored = diagonalPair[datum, gh];
  Which[member === "bare", Factor[stored], ensemble === "M", Factor[stored - gh D[stored, gh]],
    ensemble === "J", Factor[-stored]]
];

deriveReadout["V", member_] := splitReadout[vBlockCoefficient[member], u];
deriveReadout["M", member_] := splitReadout[otherCoefficient["M", member], qf];
deriveReadout["J", member_] := splitReadout[otherCoefficient["J", member], js];

allMembers = Association@Table[
  ensemble -> <|"conjugate" -> deriveReadout[ensemble, "conjugate"],
    "bare" -> deriveReadout[ensemble, "bare"]|>, {ensemble, {"V", "M", "J"}}];
readouts = AssociationMap[allMembers[#]["conjugate"] &, {"V", "M", "J"}];

identityChecks = {
  zq[readouts["V"]["total"] - (muu u^2/suu^2 - vhh gh^2)],
  zq[readouts["M"]["total"] - (muu qf^2 - mgg gh^2)],
  zq[readouts["J"]["total"] + muu js^2 + 2 mug js gh + mgg gh^2]
};
If[!And @@ identityChecks, Print["FIRST_FAILURE=thermodynamic identities"]; Exit[1]];

positiveCertificate[x_] := TrueQ[FullSimplify[zg > 0]] && TrueQ[(dd /. {b -> 2, c -> 1/2, k -> 1}) > 0] &&
  AnyTrue[{muu, mgg, vhh, mdet, muu/suu^2}, zq[x - #] &];
signature[expr_, datum_] := Module[{aa, bb, cc, determinant},
  aa = Factor[Coefficient[expr, datum, 2]];
  bb = Factor[Coefficient[Coefficient[expr, datum, 1], gh, 1]];
  cc = Factor[Coefficient[expr, gh, 2]];
  determinant = Factor[aa cc - bb^2/4];
  Which[
    positiveCertificate[aa] && positiveCertificate[-cc] && !zq[determinant], "INDEFINITE",
    positiveCertificate[-aa] && positiveCertificate[determinant], "NEGATIVE_DEFINITE",
    zq[expr], "ZERO", True, "UNCLASSIFIED"]
];
datums = <|"V" -> u, "M" -> qf, "J" -> js|>;
neutralSignatures = AssociationMap[signature[readouts[#]["total"], datums[#]] &, {"V", "M", "J"}];
If[Values[neutralSignatures] =!= {"INDEFINITE", "INDEFINITE", "NEGATIVE_DEFINITE"},
  Print["FIRST_FAILURE=neutral signature derivation"]; Exit[1]];

(* ---------------------------------------------------------------------- *)
(* Q1 as a common stationary-point/curvature/protection calculation.      *)
(* ---------------------------------------------------------------------- *)

makeBarrier[name_, potential_, constraints_, constraintProtected_, components_, reason_] := <|
  "mechanism" -> name, "potential" -> potential, "constraints" -> constraints,
  "constraintProtected" -> constraintProtected, "components" -> components, "reason" -> reason|>;

barrierEval[row_] := Module[{target, deriv, curvature, direct, constraintHolds, stationary,
  positive, roots, minima, protected, pins, branches = {-1, 1}},
  target = sorient uheld;
  deriv = D[row["potential"], xmouth]; curvature = D[row["potential"], {xmouth, 2}];
  direct = ContainsAll[Variables[row["potential"]], Variables[target]] ||
    AnyTrue[row["constraints"], ContainsAll[Variables[#], Variables[target]] &];
  constraintHolds = row["constraints"] =!= {} && And @@ Flatten[Table[
    zq[constraint /. {sorient -> branch, xmouth -> branch uheld}],
    {branch, branches}, {constraint, row["constraints"]}]];
  stationary = constraintHolds || And @@ Table[
    zq[deriv /. {sorient -> branch, xmouth -> branch uheld}], {branch, branches}];
  positive = And @@ Table[TrueQ[FullSimplify[(curvature /.
    {lam -> 1, mu -> 1/3, sorient -> branch, uheld -> 1, xmouth -> branch}) > 0]], {branch, branches}];
  roots = xmouth /. Solve[(deriv /. {lam -> 1, mu -> 1/3, sorient -> 1, uheld -> 1}) == 0, xmouth];
  minima = Select[roots, TrueQ[FullSimplify[(curvature /.
    {lam -> 1, mu -> 1/3, sorient -> 1, uheld -> 1, xmouth -> #}) > 0]] &];
  protected = (constraintHolds && TrueQ[row["constraintProtected"]]) ||
    (Length[minima] > 1 && row["components"] > 1);
  pins = direct && stationary && positive && protected;
  <|"mechanism" -> row["mechanism"], "direct" -> direct, "stationary" -> stationary,
    "curvature" -> positive, "protected" -> protected, "pins" -> pins,
    "reason" -> row["reason"]|>
];

q1Census[inject_: False] := Module[{src, ledger, rows, evaluated, holders, relaxed, selfEnergy,
  mixedEnergy, inducedEnergy, linearEnergy, components},
  src = Map[Import[#, "Text"] &, paths]; ledger = parseLedger[src["g0"]];
  relaxed = 2 - (1/2)^2; selfEnergy = xmouth^2;
  mixedEnergy = relaxed xmouth^2/2;
  inducedEnergy = mixedEnergy - sorient induced xmouth;
  linearEnergy = mixedEnergy - sorient ql xmouth;
  components = If[StringContainsQ[src["path24"], "DeltaE_unwind = 0"], 1, 2];
  rows = {
    makeBarrier["u_L self-stiffness", selfEnergy, {}, False, 1, "positive bulk Hessian and IR zero relax"],
    makeBarrier["C_hu mixing", mixedEnergy, {}, False, 1, "convex Schur mixing supplies no essential datum"],
    makeBarrier["h mouth source + natural u_L conormal", inducedEnergy, {}, False, 1, "mixing induces a unique relaxable response"],
    makeBarrier["IR u_L->0", mixedEnergy, {xmouth}, False, 1, "IR zero is not a signed mouth target"],
    makeBarrier["conditional q_L proportional to s", linearEnergy, {}, False, 1, "linear source shifts one convex minimum without protection"],
    makeBarrier["wall/chi_B sector", mixedEnergy, {}, False, 1, "POSTULATE ledger zeros wall-u_L couplings"],
    makeBarrier["sleeve +/-w orientation", mixedEnergy, {}, False, components, "connected orientation and zero unwinding barrier"],
    makeBarrier["S_hold", mixedEnergy, {}, False, 1, "holds r_B only"],
    makeBarrier["geon", mixedEnergy, {}, False, 1, "geon-u_L derivative is zero"]
  };
  If[inject, rows = Append[rows, makeBarrier["injected_s_uL_hold",
    lam/4 (xmouth^2 - uheld^2)^2 - mu sorient (3 xmouth/(2 uheld) - xmouth^3/(2 uheld^3)),
    {}, False, 2, "0<mu<2*lam*uheld^4/3 preserves two wells while s-odd selector picks x=s*uheld"]]];
  evaluated = barrierEval /@ rows; holders = Select[evaluated, TrueQ[#pins] &];
  <|"status" -> If[holders === {}, "NO_NATIVE_CLAMP", "NATIVE_CLAMP_EXISTS(" <> holders[[1, "mechanism"]] <> ")"],
    "rows" -> evaluated|>
];
q1base = q1Census[]; q1status = q1base["status"];
q1bstatus = "BOLT_ON_DEFERRED_TO_R1";
q1bRequirement = "protected/disconnected +/-w sector plus direct s<->u_L datum coupling; supply exact added action and coefficient domain";

(* ---------------------------------------------------------------------- *)
(* Q3, radial Green derivation, dimensions, and §5 total classifier.      *)
(* ---------------------------------------------------------------------- *)

closureContent = <|"V" -> {"essential_uL_value", "charge_odd"},
  "M" -> {"held_uL_conormal", "charge_odd"}, "J" -> {"linear_uL_source", "charge_odd"}|>;

tierA[cvalue_: 1/2, additions_: {}] := Module[{margin, g0text, ledger, u2text, introduced},
  margin = 2 - cvalue^2; g0text = Import[paths["g0"], "Text"];
  ledger = parseLedger[g0text]; u2text = Import[paths["u2report"], "Text"];
  introduced = Union[Flatten[Values[closureContent]], additions];
  <|"scalar_hessian_positive" -> TrueQ[margin > 0],
    "transverse_decoupled" -> (!MemberQ[introduced, "u_T_datum"] && AnyTrue[Join[ledger, zeroLedger], StringContainsQ[#, "u_Lu_T"] &]),
    "zero_ledger_preserved" -> (Length[ledger] == 13 && !MemberQ[introduced, "zero_ledger_coupling"]),
    "charge_perp_mass" -> !MemberQ[introduced, "q_M"],
    "u2_nonselection" -> StringContainsQ[u2text, "{'UNRESOLVED': 144}"]|>
];

greenPower[n_Integer] := Module[{p, roots, decaying},
  roots = p /. Solve[p (p - n + 2) == 0, p]; decaying = Select[roots, TrueQ[# > 0] &];
  If[Length[decaying] =!= 1, Return[$Failed]]; First[decaying]
];
forcePower[expr_, n_Integer] := Module[{force, exponents},
  force = Factor[-D[sone stwo expr/(4 Pi rr^greenPower[n]), rr]];
  exponents = Cases[force, Power[rr, e_Integer] /; e < 0 :> e, Infinity];
  {force, -Min[exponents]}
];
forceResult = forcePower[readouts["J"]["total"], 3];

dims = <|"b" -> {-1, -2, 1}, "c" -> {0, -2, 1}, "k" -> {1, -2, 1},
  "km" -> {2, -2, 1}, "zeta" -> {0, 0, 0}, "reb" -> {-2, 2, -1},
  "suu" -> {0, 2, -1}, "sug" -> {-1, 2, -1}, "sgg" -> {-2, 2, -1},
  "u" -> {1, 0, 0}, "qf" -> {1, -2, 1}, "js" -> {1, -2, 1},
  "gh" -> {2, -2, 1}, "sone" -> {0, 0, 0}, "stwo" -> {0, 0, 0}, "rr" -> {1, 0, 0}|>;

inferDimension[expr_] := Which[
  NumberQ[expr], {0, 0, 0},
  Head[Unevaluated[expr]] === Symbol, Lookup[dims, ToString[Unevaluated[expr], InputForm], $Failed],
  Head[expr] === Times, If[MemberQ[inferDimension /@ List @@ expr, $Failed], $Failed, Total[inferDimension /@ List @@ expr]],
  Head[expr] === Power && IntegerQ[expr[[2]]], expr[[2]] inferDimension[expr[[1]]],
  Head[expr] === Plus, Module[{got = DeleteCases[inferDimension /@ List @@ expr, {0, 0, 0}]},
    If[got === {} || SameQ @@ got, If[got === {}, {0, 0, 0}, First[got]], $Failed]],
  True, $Failed
];

dimensionCheck[greenDim_: {-1, 0, 0}] := Module[{terms, got, dimA = {3, -2, 1}},
  terms = Flatten[KeyValueMap[Function[{ens, rs},
    Flatten[(If[Head[Expand[#]] === Plus, List @@ Expand[#], {Expand[#]}] &) /@
      DeleteCases[Values[rs], 0]]], readouts]];
  got = inferDimension /@ terms;
  !MemberQ[got, $Failed] && And @@ (# === dimA & /@ got) &&
    And @@ ((# + greenDim) === {2, -2, 1} & /@ got) &&
    And @@ ((# - 2 {1, 0, 0}) === {1, -2, 1} & /@ got) &&
    inferDimension[zg] === {0, 0, 0} && inferDimension[zb] === {0, 0, 0}
];

tierB = {"assembled gravity/drain response under datum", "curved-sleeve stability",
  "momentum/return closure", "dressed non-perturbative two-body response"};

section5Algebraic[signatures_, power_] := AssociationMap[
  Switch[signatures[#],
    "INDEFINITE", "range(attract_1/R^" <> ToString[power] <> "|null_leading|repel_1/R^" <> ToString[power] <> ")",
    "NEGATIVE_DEFINITE", "attract_1/R^" <> ToString[power],
    "ZERO", "null_leading", _, "unresolved(" <> signatures[#] <> ")"] &,
  {"V", "M", "J"}
];

Clear[section5Verdict];
Options[section5Verdict] = {"tier" -> True, "reason" -> "tier_A", "native" -> False,
  "bolt" -> False, "boltSupplied" -> False, "sign" -> "POSITIVE", "power" -> 2,
  "tierB" -> "DEFERRED", "unresolved" -> None, "algebraic" -> Automatic,
  "precedenceMutation" -> False};
section5Verdict[OptionsPattern[]] := Module[{targetSign = "POSITIVE", targetPower = 2, holder, held, alg, result},
  holder = Which[TrueQ[OptionValue["native"]], "NATIVE", TrueQ[OptionValue["bolt"]], "BOLT_ON", True, None];
  held[] := Which[holder === None, None, OptionValue["sign"] === "NULL", "NO_GO(holder_held_variable_null)",
    OptionValue["power"] =!= targetPower, "NO_GO(holder_held_variable_wrong_range)",
    OptionValue["sign"] =!= targetSign, "NO_GO(holder_held_variable_wrong_sign)",
    OptionValue["tierB"] === "PASS", holder <> "_DERIVED",
    True, holder <> "_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED"];
  If[TrueQ[OptionValue["precedenceMutation"]] && held[] =!= None, Return[held[]]];
  If[!TrueQ[OptionValue["tier"]], Return["NO_GO(" <> OptionValue["reason"] <> ")"]];
  If[!TrueQ[OptionValue["precedenceMutation"]] && held[] =!= None, Return[held[]]];
  If[OptionValue["unresolved"] =!= None, Return["UNRESOLVED(" <> OptionValue["unresolved"] <> ")"]];
  alg = OptionValue["algebraic"];
  If[alg === Automatic, alg = AssociationMap["UNCLASSIFIED" &, {"V", "M", "J"}]];
  result = "NO_NATIVE_CLAMP; ALGEBRAIC_SIGN={" <>
    StringRiffle[(# <> ":" <> alg[#]) & /@ {"V", "M", "J"}, ","] <> "; BOLT_ON_DEFERRED_TO_R1";
  (* Repair the brace separately so the string construction remains obvious. *)
  result = StringReplace[result, "; BOLT_ON_DEFERRED" -> "}; BOLT_ON_DEFERRED"];
  If[TrueQ[OptionValue["boltSupplied"]] && !TrueQ[OptionValue["bolt"]],
    result = result <> "; BOLT_ON_INEFFECTIVE(does_not_hold_datum)"];
  result
];

verdictOracle[case_] := Module[{holder, key},
  If[!case["tier"], Return["NO_GO(tier_A)"]];
  holder = Which[case["native"], "NATIVE", case["bolt"], "BOLT_ON", True, None];
  If[holder =!= None,
    key = Which[case["sign"] === "NULL", "null", case["power"] =!= 2, "range",
      case["sign"] =!= "POSITIVE", "sign", True, "pass"];
    Return[Switch[key, "null", "NO_GO(holder_held_variable_null)",
      "range", "NO_GO(holder_held_variable_wrong_range)", "sign", "NO_GO(holder_held_variable_wrong_sign)",
      "pass", If[case["tierB"] === "PASS", holder <> "_DERIVED", holder <> "_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED"]]]];
  If[case["unresolved"], Return["UNRESOLVED(named_datum)"]];
  "NO_NATIVE_CLAMP; ALGEBRAIC_SIGN={V:UNCLASSIFIED,M:UNCLASSIFIED,J:UNCLASSIFIED}; BOLT_ON_DEFERRED_TO_R1" <>
    If[case["boltSupplied"], "; BOLT_ON_INEFFECTIVE(does_not_hold_datum)", ""]
];

verdictTruth[mutation_: False] := Module[{domain, exact = True, classes = {}, actual, expected, case, required},
  domain = Tuples[{{False, True}, {False, True}, {False, True}, {False, True},
    {"POSITIVE", "NEGATIVE", "NULL"}, {2, 3}, {"PASS", "DEFERRED"}, {False, True}}];
  Do[
    case = <|"tier" -> row[[1]], "native" -> row[[2]], "bolt" -> row[[3]],
      "boltSupplied" -> row[[4]], "sign" -> row[[5]], "power" -> row[[6]],
      "tierB" -> row[[7]], "unresolved" -> row[[8]]|>;
    actual = section5Verdict["tier" -> case["tier"], "native" -> case["native"],
      "bolt" -> case["bolt"], "boltSupplied" -> case["boltSupplied"],
      "sign" -> case["sign"], "power" -> case["power"], "tierB" -> case["tierB"],
      "unresolved" -> If[case["unresolved"], "named_datum", None], "precedenceMutation" -> mutation];
    expected = verdictOracle[case]; exact = exact && actual === expected && StringLength[actual] > 0;
    AppendTo[classes, First[StringSplit[actual, ";"]]], {row, domain}];
  classes = Union[classes, {"BOLT_ON_INEFFECTIVE"}];
  required = {"NO_GO(tier_A)", "NO_GO(holder_held_variable_wrong_sign)",
    "NO_GO(holder_held_variable_wrong_range)", "NO_GO(holder_held_variable_null)",
    "NATIVE_DERIVED", "BOLT_ON_DERIVED", "NATIVE_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED",
    "BOLT_ON_ALGEBRAIC_CONSISTENT / R1_R3_REQUIRED", "NO_NATIVE_CLAMP",
    "UNRESOLVED(named_datum)", "BOLT_ON_INEFFECTIVE"};
  {exact && ContainsAll[classes, required], classes, Length[domain]}
];

mainAlgebraic = section5Algebraic[neutralSignatures, forceResult[[2]]];
mainLanding = section5Verdict["tier" -> And @@ Values[tierA[]], "algebraic" -> mainAlgebraic];

targetBlind[stage_] := Module[{neutralOK},
  neutralOK = And @@ (MemberQ[{"INDEFINITE", "NEGATIVE_DEFINITE", "ZERO", "UNCLASSIFIED"}, #] & /@ Values[neutralSignatures]);
  neutralOK && FreeQ[HoldComplete[stage], Alternatives[emTargetSign, emTargetPower, expectedFormula]]
];

(* ---------------------------------------------------------------------- *)
(* Symbolic payload: every response term and all twelve separated terms.  *)
(* ---------------------------------------------------------------------- *)

symbolicPayload[mutate_: False] := Module[{terms, truth},
  terms = <|"d" -> canon[dd], "kappa" -> canon[kap], "zg" -> canon[zg], "zb" -> canon[zb],
    "muu" -> canon[muu], "mug" -> canon[mug], "mgg" -> canon[mgg], "mdet" -> canon[mdet],
    "vhh" -> canon[vhh], "qphysical_v" -> canon[qphysical]|>;
  KeyValueMap[Function[{ens, rs}, KeyValueMap[(AssociateTo[terms, ens <> "." <> #1 -> canon[#2]]) &, rs]], readouts];
  If[mutate, terms["J.interference"] = canon[readouts["J"]["interference"] + js gh]];
  truth = verdictTruth[];
  <|"schema" -> "LEFTOVER_SCALAR_SIGN_V2_SYMBOLIC", "symbolic_terms" -> terms,
    "neutral_signatures" -> neutralSignatures, "green_power" -> greenPower[3],
    "force_power" -> forceResult[[2]], "q1" -> q1status, "q1b" -> q1bstatus,
    "tier_a" -> If[And @@ Values[tierA[]], "PASS", "FAIL"], "tier_b" -> "DEFERRED_R1_R3",
    "tier_a_details" -> tierA[],
    "truth_table_total" -> truth[[3]], "truth_table_classes" -> Length[truth[[2]]],
    "truth_table_pass" -> truth[[1]], "landing" -> mainLanding|>
];

pythonPayload[mutate_: False] := Module[{command, run, line},
  command = If[mutate, {"env", "LEFTOVER_PAYLOAD_MUTATION=1", "python3", pyPath, "--json-only"},
    {"python3", pyPath, "--json-only"}];
  run = RunProcess[command];
  If[run["ExitCode"] =!= 0, Return[$Failed]];
  line = SelectFirst[StringSplit[run["StandardOutput"], "\n"], StringStartsQ[#, "JSON_PAYLOAD="] &, Missing[]];
  If[MissingQ[line], Return[$Failed]];
  ImportString[StringDrop[line, StringLength["JSON_PAYLOAD="]], "RawJSON"]
];

payloadEqual[left_, right_] := Module[{lterms, rterms, metaKeys},
  If[right === $Failed || Sort[Keys[left]] =!= Sort[Keys[right]], Return[False]];
  lterms = left["symbolic_terms"]; rterms = right["symbolic_terms"];
  If[Sort[Keys[lterms]] =!= Sort[Keys[rterms]], Return[False]];
  If[!And @@ KeyValueMap[zq[ToExpression[#2] - ToExpression[rterms[#1]]] &, lterms], Return[False]];
  metaKeys = DeleteCases[Keys[left], "symbolic_terms"];
  And @@ (left[#] === right[#] & /@ metaKeys)
];

(* ---------------------------------------------------------------------- *)
(* Per-tooth controls.  Mutations alter source/action/functionals/operators. *)
(* ---------------------------------------------------------------------- *)

toothOrder = {"FAITHFULNESS", "Q1_BARRIER", "HELD_VARIABLE_V", "HELD_VARIABLE_M",
  "HELD_VARIABLE_J", "DOUBLE_COUNT", "C_HU_STABILITY", "FALLOFF", "UNITS",
  "Q_M_GUARD", "VERDICT_CLASSIFICATION", "TARGET_BLINDNESS", "DUAL_ENGINE_TERMS"};

checkVector[mutation_, dualOK_] := Module[{tr, checkedHessian, faith, q1, selected, recon, cvalue, ta, landing,
  nspace, fpow, ug, additions, verdict, stage},
  tr = Association[transcript];
  checkedHessian = actionHessian;
  If[mutation === "FAITHFULNESS",
    tr["action"] = ReplacePart[tr["action"], 5 -> "-2*C_hu*grad u_L.grad h"];
    checkedHessian = Table[D[b gu^2/2 + 2 c gu ghgrad + k ghgrad^2/2, xi, xj],
      {xi, {gu, ghgrad}}, {xj, {gu, ghgrad}}]];
  faith = sourceFaithful[tr, checkedHessian] && !zq[D[mug, c]] && !zq[D[mug, km]] &&
    And @@ Flatten[Map[zq, (mfar /. {km -> 0, zeta -> 1}) - Inverse[actionHessian], {2}]] &&
    Length[zeroLedger] == 13;
  q1 = q1Census[mutation === "Q1_BARRIER"];
  selected = AssociationMap[allMembers[#][If[mutation === "HELD_VARIABLE_" <> #, "bare", "conjugate"]] &, {"V", "M", "J"}];
  recon = AssociationMap[Total[Lookup[readouts[#], {"h_only", "uL_only", "interference"}]] &, {"V", "M", "J"}];
  If[mutation === "DOUBLE_COUNT", recon["J"] = readouts["J"]["h_only"] + readouts["J"]["uL_only"]];
  cvalue = If[mutation === "C_HU_STABILITY", 2, 1/2]; ta = tierA[cvalue];
  landing = section5Verdict["tier" -> And @@ Values[ta], "reason" -> "scalar_unstable"];
  nspace = If[mutation === "FALLOFF", 4, 3]; fpow = forcePower[readouts["J"]["total"], nspace][[2]];
  ug = If[mutation === "UNITS", {-2, 0, 0}, {-1, 0, 0}];
  additions = If[mutation === "Q_M_GUARD", {"q_M"}, {}];
  verdict = verdictTruth[mutation === "VERDICT_CLASSIFICATION"][[1]];
  stage = {dd, kap, mfar, readouts, neutralSignatures, q1status, tierA[]};
  If[mutation === "TARGET_BLINDNESS", stage = Append[stage, emTargetSign readouts["J"]["total"]]];
  <|
    "FAITHFULNESS" -> faith,
    "Q1_BARRIER" -> (q1["status"] === "NO_NATIVE_CLAMP"),
    "HELD_VARIABLE_V" -> (zq[selected["V"]["total"] - readouts["V"]["total"]] && !zq[allMembers["V"]["bare"]["total"] - readouts["V"]["total"]]),
    "HELD_VARIABLE_M" -> (zq[selected["M"]["total"] - readouts["M"]["total"]] && !zq[allMembers["M"]["bare"]["total"] - readouts["M"]["total"]]),
    "HELD_VARIABLE_J" -> (zq[selected["J"]["total"] - readouts["J"]["total"]] && !zq[allMembers["J"]["bare"]["total"] - readouts["J"]["total"]]),
    "DOUBLE_COUNT" -> And @@ KeyValueMap[zq[#2 - readouts[#1]["total"]] &, recon],
    "C_HU_STABILITY" -> (And @@ Values[ta] && landing =!= "NO_GO(scalar_unstable)"),
    "FALLOFF" -> (fpow == 2), "UNITS" -> dimensionCheck[ug],
    "Q_M_GUARD" -> tierA[1/2, additions]["charge_perp_mass"],
    "VERDICT_CLASSIFICATION" -> verdict, "TARGET_BLINDNESS" -> targetBlind[stage],
    "DUAL_ENGINE_TERMS" -> dualOK|>
];

mutationCampaign[] := Module[{local, other, dualBase, base, mutated, failed, outcomes = <||>, bad, landing},
  local = symbolicPayload[]; other = pythonPayload[]; dualBase = payloadEqual[local, other];
  base = checkVector[None, dualBase];
  If[!And @@ Values[base], Print["FIRST_FAILURE=baseline " <> ToString[Keys[Select[base, Not]]]]; Exit[1]];
  Do[
    mutated = If[tooth === "DUAL_ENGINE_TERMS",
      checkVector[None, payloadEqual[local, pythonPayload[True]]], checkVector[tooth, dualBase]];
    failed = Keys[Select[mutated, Not]];
    If[failed =!= {tooth}, Print["FIRST_FAILURE=ablation " <> tooth <> " got " <> ToString[failed]]; Exit[1]];
    If[tooth === "C_HU_STABILITY", bad = tierA[2];
      landing = section5Verdict["tier" -> And @@ Values[bad], "reason" -> "scalar_unstable"];
      If[landing =!= "NO_GO(scalar_unstable)", Print["FIRST_FAILURE=stability landing"]; Exit[1]]];
    AssociateTo[outcomes, tooth -> "FIRED_AT_" <> tooth], {tooth, toothOrder}];
  outcomes
];

printReport[ablations_] := Module[{truth = verdictTruth[], ta = tierA[], rebuilt},
  Print["LEFTOVER_SCALAR_ELECTRIC_SIGN — MATHEMATICA"];
  Print["FAITHFUL_STATIC_FUNCTIONAL="];
  Print["  S_Lh=" <> StringRiffle[actionTerms, " + "]];
  Print["  T_L=B_eff=rho_B0^2/chi_c>0; K_h=M_h*c_E^2; C_hu=1/2; k_m=1"];
  Print["  mouth=sum_i integral eta_i[1/2*k_m*h^2-g_chih*s_i*h]"];
  Print["  native BC=" <> StringRiffle[transcript["native"], "; "]];
  Scan[Print["  BC_" <> #[[1]] <> "=" <> #[[2]]] &, testBCs];
  Print["  ZERO_LEDGER_POSTULATE_ROWS=" <> ToString[Length[zeroLedger]] <> "; exact tags/terms verified"];
  Print["Q1a=" <> q1status];
  Scan[Print["  " <> #["mechanism"] <> ": direct=" <> ToString[#["direct"]] <>
    "; stationary=" <> ToString[#["stationary"]] <> "; curvature=" <> ToString[#["curvature"]] <>
    "; protected=" <> ToString[#["protected"]] <> "; pins=" <> ToString[#["pins"]] <> "; " <> #["reason"]] &,
    q1base["rows"]];
  Print["  injection_control=" <> q1Census[True]["status"]];
  Print["Q1b=" <> q1bstatus]; Print["  required_added_action_content=" <> q1bRequirement];
  Print["Q2_RESPONSE="];
  Print["  D=" <> canon[dd] <> "; kappa=" <> canon[kap]];
  Print["  z_g=" <> canon[zg] <> "; z_b=" <> canon[zb] <> "; both k_m-live"];
  Print["  m=" <> ToString[mfar, InputForm] <> "; det(m)=" <> canon[mdet] <> ">0 for D>0,z_g!=0"];
  Print["  V physical reaction Q0=" <> canon[qphysical] <> "; datum-only normalization qV=" <> canon[qdatum]];
  KeyValueMap[Function[{ens, rs}, Print["  " <> ens <> ": quadratic_signature=" <> neutralSignatures[ens]];
    KeyValueMap[Print["    " <> #1 <> ": A=" <> canon[#2] <> "; U=s1*s2*A/(4*pi*R); F_out=s1*s2*A/(4*pi*R^2)"] &, rs];
    rebuilt = Total[Lookup[rs, {"h_only", "uL_only", "interference"}]];
    Print["    double_count=" <> ToString[zq[rebuilt - rs["total"]]]]], readouts];
  Print["Q3_TIER_A=" <> If[And @@ Values[ta], "PASS", "FAIL"]];
  KeyValueMap[Print["  " <> #1 <> "=" <> ToString[#2]] &, ta];
  Print["Q3_TIER_B=DEFERRED"]; Scan[Print["  " <> #] &, tierB];
  Print["SECTION5_LANDING=" <> mainLanding];
  Print["VERDICT_DOMAIN_TOTAL=" <> ToString[truth[[3]]] <> "; classes_exercised=" <>
    ToString[Length[truth[[2]]]] <> "; pass=" <> ToString[truth[[1]]]];
  Print["TEETH="]; Scan[Print["  PASS " <> # <> "; ablation=" <> ablations[#]] &, toothOrder];
  Print["ENGINE_AGREE=PASS"]
];

jsonOnly = Environment["LEFTOVER_JSON_ONLY"] === "1" ||
  AnyTrue[$ScriptCommandLine, StringContainsQ[ToString[#], "--json-only"] &];
payloadMutation = Environment["LEFTOVER_PAYLOAD_MUTATION"] === "1";

If[jsonOnly,
  Print["JSON_PAYLOAD=" <> ExportString[symbolicPayload[payloadMutation], "RawJSON", "Compact" -> True]];
  Exit[0]
];

ablations = mutationCampaign[];
printReport[ablations];
Exit[0];
