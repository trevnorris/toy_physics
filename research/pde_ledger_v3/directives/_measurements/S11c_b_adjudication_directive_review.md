# S11c-b adjudication build directive — decision-list review, round 1 (two legs, NOT SOUND → folded)

**Artifact reviewed:** `directives/S11c_b_adjudication_build_directive.md` (orchestrator-written, first draft).
**Legs (orchestrator-written ⇒ Codex + Grok):** prompt `_legs/S11c_b_adjudication_directive_review.md`.
Raw: Codex `~/.s11_build/S11c_b_adjud_directive_codex.txt`; Grok `~/.s11_build/S11c_b_adjud_directive_grok.txt`.
**Both legs independently: NOT SOUND — do not launch the builder.** Bridge A verified sound by both.

## Verdict convergence
Both legs converged on the same core defect (the computed gate is ill-posed because the collapse lives in the
committed comparator, not the reconcile layer) and on the rule-5 leaks. Codex added the false-FLAG predicate
defect and the container-schema-mismatch defect. This is the rule-7 gate catching an architectural error in the
directive before any build round was spent.

## Orchestrator verification (rule 13 — I confirmed the load-bearing claim myself)
Grok/Codex claim: the comparator peels every WL `Derivative[...]` into a distinct jet token BEFORE the
reconcile layer sees an operand, so a post-canon "scan for `Derivative(f,x_i)`" gate finds an empty set for
every head and is extensionally a blanket applied→bare collapse.

```
$ grep -c 'Derivative(' ~/.s11_build/S11c_b_reconcile_run.out
0
$ grep -oE '(theta_d[123]|W_bg_d[123]|w1_profile_d[123]|mu_R_bg_d[123])' ~/.s11_build/S11c_b_reconcile_run.out | sort | uniq -c
    564 mu_R_bg_d1 ... 2582 w1_profile_d1 ... 972 W_bg_d1 ...   (jets already peeled to distinct symbols)
$ grep -nE 'def canon_wl_basic|Derivative\[|XJETX|FIELD *=' scripts/S11c_a_cross_engine_comparator.py
692: if text.startswith("Derivative[", index):  ... 714/721/730/736: token += "XJETX" + suffix
826: def canon_wl_basic(...)   534: FIELD = {   558: PROFILE = {"w1Jet": "w1_profile", ...}
```
CONFIRMED: 0 Derivative nodes in the post-reconcile residuals; the load-bearing Derivative→jet decoding is in
`C.s11ca` (L692, L826) + `C.EXTRA_HEAD` (S11c-b L315), all committed `5f01f9fa` and reviewed by the comparator
build legs. Codex's raw-transcript counts (1,272 `thetaWave`, 2,404 `eWave`, 11,019 `longitudinalTestPotential`
derivatives, all erased post-canon) are consistent with this.

## Blocking findings (all folded into directive v2)
1. **Bridge B computed gate is ill-posed** — scans already-canonicalized operands; real collapse is upstream in
   the legged comparator (`FIELD`/`EXTRA_HEAD`/`canon_wl_basic`), not `H.BARE_APPLIED`. **Fix:** DROP the
   computed-collapse gate. The never-blanket-collapse guard is already structural in the comparator's lossless
   jet decoding; the layer must introduce NO new applied→bare or jet fold.
2. **Wrong predicate → false FLAGs** (Codex) — "collapse iff never differentiated" would KEEP `thetaWave`/`eWave`
   (nonempty raw-jet set) as applied while PY uses bare base + separate jets — a manufactured representational
   FLAG. **Fix:** the correct criterion is lossless-jet decoding (every derivative → a distinct jet token, none
   discarded/mapped-to-base), which the comparator already implements; the layer VERIFIES jet-token
   conservation, it does not re-collapse.
3. **Container-schema families** (Codex) — ENERGY_BASIS / DIMENSIONS / HOMOGENEITY / ADMISSIBILITY_RESIDUAL /
   TERM_ORIGINS FLAGs are container/arity/key SHAPE mismatches, not algebraic; Bridges A/B cannot touch them.
   PY-selected reps `*_07`/`*_10` have no WL coefficient (one-engine physics) and were unprotected. **Fix:** a
   safe label-preserving leaf adapter only where it never identifies a quotient representative (DIMENSIONS
   values); classify the rest STRUCTURE_INCOMPLETE and OWE them; protect `*_07`/`*_10` explicitly.
4. **Ablation hooks contradictory/non-decisive** — `--force-collapse` of an applied head is a no-op on jets
   (symbols post-parse); `--drop-rename` must disable the comparator prepass (v1 does at `H:440`). **Fix:**
   redefine the jet ablation in RAW-TOKEN terms (`--collapse-jet <token>=<base>`, exact target named);
   `--drop-bridge-a` emits touched cases + pre/post operands; `--drop-rename` disables the prepass.
5. **Rule-5 leaks** — "every CORE family still FLAGs" (also FACTUALLY FALSE: ADMISSIBILITY_SUPPORT_OPERAND +
   ENERGY_BASIS_COUNT MATCH), "genuine non-self-adjoint Λ_X", "manufactures a false MATCH", "removal must move
   the residuals". **Fix:** strip every expected outcome from the directive; value-free structural DoD only;
   the load-bearing expectation lives ONLY in the build-leg prompt, not the build directive.

## Cleared by both legs
Bridge A (`B_rho_3 = bRho·W_0`) — coefficient alignment verified against WL:472/1625 + PY:1135 + spec:102.
CORE-family scope matches v1; control families stay NAMESPACE_INCOMPLETE; Bridge C integral-linearity reuse
sound.

## Disposition
Folded into `directives/S11c_b_adjudication_build_directive.md` v2 (one pass, per rule-7 "fold and go"), then
re-legged. If round 2 breeds further defects in the same material, delegate the directive rewrite to Codex
(rule 15).

---

# Round 2 — decision-list review of directive v2 (two legs, NOT SOUND → folded into v3)

Raw: Codex `~/.s11_build/S11c_b_adjud_directive_r2_codex.txt`; Grok `~/.s11_build/S11c_b_adjud_directive_r2_grok.txt`.
**Both legs: architecture now SOUND; remaining defects are PRECISION, not structural.** Both independently
re-verified Bridge A's factor and the one-engine protection (selected `{1,4,5,6,7,9,10,13}`; WL's 8 explicit
invariants → PY `{1,4,5,6,8,9,11,13}`, so WL-only DivGrad = 08/11, PY-only selected = 07/10; all four protected).

## Blocking findings (folded into v3)
- **Bridge A form** (both): stated as a PRODUCT rewrite `bRho·W_0 ↦ B_rho_3`, which misses `bRho·W_bg` and bare
  `bRho`. Fix: SYMBOL substitution `bRho ↦ B_rho_3/W_0`; non-ablated algebraic operands must carry no residual
  `bRho`.
- **Container adapter can manufacture agreement by subset-matching** (Codex): v2 forbade differently-labelled
  pairing but not partial/non-bijective alignment — a builder could match two agreeing leaves and drop the rest
  to MATCH. Fix: a container MATCH requires a source-cited TOTAL BIJECTION over ALL semantic leaves; any
  unmatched/duplicate/unlabeled leaf → STRUCTURE_INCOMPLETE; subset comparisons are diagnostics only. DIMENSIONS
  is NOT blanket leaf-safe (PY full dimension tree vs WL coarse schema).
- **A source-backed algebraic family would be wrongly deferred** (Codex): `SLAB_OPERATOR_TERM_ORIGINS/KINETIC`
  is tuple-vs-Association but source fixes the correspondence (PY tuple pos0=`U_MOMENTUM_ROWS`,
  pos1=`THICKNESS_ROW`, PY:1573; WL labels WL:851). Fix: the builder ENUMERATES each source-cited total-bijection
  label-restoration adapter (KINETIC is the pattern) and then compares algebraically; build legs verify each.
- **STRUCTURAL_ATOM_DISAGREE ≠ algebraic** (Codex): some admissibility operands are Boolean `Equivalent`/`Not`
  (PY:868-890, run:160-175), not arithmetic. Fix: algebraic route ONLY when both leaves are arithmetic `sp.Expr`
  (or same-shape matrices/tuples thereof); Boolean/Text/sort-or-shape-mismatch stays STRUCTURE_INCOMPLETE.
- **Jet conservation undefined across a legitimate rename** (both): v1 renames `theta_d1→grad_theta_1`; a raw
  `*_d1` counter false-reports JET_LOST. Fix: semantic jet ID `(canonical_base_after_rename, multiindex)`; the
  layer's SUBSTITUTION MAP must be jet-order-preserving (no substitution lowers a derivative order); before-set
  taken from comparator-captured operands pre-reconcile.
- **Ablation DoD admits no-ops** (Codex): "work and print" passes a hook that ignores an absent token. Fix: each
  hook selects its argument from the captured pre-transform inventory, reports a nonempty touched-case set,
  shows a syntactically changed operand, and prints recomputed before/after residuals; unknown/untouched arg =
  operational error.
- **Rule-5 leaks remain** (both — MY recurring defect, third occurrence this session): "mostly FLAG, two MATCH"
  (measured), "genuine agreement / genuine one-engine finding" (interpreted outcome), "the collapse ablation
  will force JET_LOST" (predicted), "absent an ablation every line is JET_CONSERVED" (expected classification).
  Fix: strip all; describe only computations; labels/classifications come from the run.

## Non-blocking (folded)
- Route `TextAtom('UNDEFINED_UNJOINED')` join atoms → COVERAGE/STRUCTURE_INCOMPLETE, never algebraic (Grok).
- Candidate-set citation is computed at ~PY:1273, not a literal at PY:1073 (both).

## Disposition — restructure, do not re-enumerate
Rounds 1-2 drew NEW schema cases each time (the pre-enumeration trap,
[[feedback_delegate_build_when_prose_directive_repeats]]). v3 restructures around INVARIANTS + ACCOUNTING +
VERIFIED REFS rather than a per-family schema list: lock Bridge A + the protected never-fold set (both
leg-verified), state routing as sort/bijection invariants that cover all cases uniformly, hand the builder the
verified refs (bRho identity, 07/10 & 08/11 sets, the KINETIC correspondence, H:440 prepass-disable, the
Boolean-atom locations), and let the BUILD legs gate the per-family adapter mechanics. If round-3 legs breed
NEW blocking defects (not nits), delegate the directive rewrite to Codex (rule 15).

---

# Round 3 — decision-list review of directive v3 (two legs, NOT SOUND → architecture confirmed; one comparator defect surfaced)

Raw: Codex `~/.s11_build/S11c_b_adjud_directive_r3_codex.txt`; Grok `~/.s11_build/S11c_b_adjud_directive_r3_grok.txt`.
**Both legs: architecture SOUND (sort invariants + total-bijection containers + Bridge A + protected one-engine
set), re-derived with evidence.** Remaining directive defects are precision; Codex additionally surfaced a
latent defect in the COMMITTED comparator.

## Directive precision defects (to fold)
- **Rule-5 leaks, FOURTH occurrence** (both legs; Codex found four): builder-report "every non-ablated case is
  JET_CONSERVED" (L127), "a case carrying 07/10 FLAGs" (L52), a build leg "breaks the reconcile" (L94),
  "--drop-bridge-a restores" (L109). Strip all; report computed labels only.
- **Incomplete taxonomy** (Codex): no `CONTAINER FLAG` for a total-bijection with nonzero leaf residual (KINETIC
  can be a faithful bijection AND nonzero); no route for joined scalar-sort mismatches (the 4
  COUPLING_KERNEL/ADJOINTNESS_RELATION cases: PY `NO_INDEPENDENT_SECOND_ROUTE` scalar PY:1980 vs WL relation
  WL:1130). Fix: ALGEBRAIC(MATCH|FLAG); CONTAINER(MATCH|FLAG|STRUCTURE_INCOMPLETE); STRUCTURE_INCOMPLETE for
  joined Boolean/TextAtom/relational/scalar-sort; COVERAGE.
- **Accounting denominator wrong** (Codex): must be `join + py_only + wl_only` (147 join + 84 one-engine = 231),
  not `join`; assert exact case-ID MULTISET equality (captured vs classified), operational failure on
  missing/dup/extra — a counter sum cannot detect one-drop-plus-one-dup.
- **Positive arithmetic-sort test** (Codex): in this SymPy, `Symbol` is both `sp.Expr` AND `Boolean`; define
  arithmetic leaf POSITIVELY as `isinstance(v, sp.Expr)` (Equivalent/Not/relational are not Expr; ordinary
  Symbols are). Fixtures for all five sorts.
- **Jet-ID alias** (both): `grad_theta_N` is not `_dN` grammar and `canon_jet_name` no-ops it; define
  `ID(theta_dN)=ID(grad_theta_N)=(theta,(N,))`, take before/after multisets AFTER applying the rename.

## ⛔ COMPARATOR DEFECT surfaced (Codex finding 3) — VERIFIED by orchestrator (rule 13)
`canon_jet_name` (`scripts/S11c_a_cross_engine_comparator.py:800-811`) records time differentiation as a
BOOLEAN `has_time`, emitting a single `_t` regardless of order:
```
has_time = False
for token in derivatives:
    if token == "t": has_time = True        # <-- boolean, not a count
suffix = "_t" if has_time else ""           # <-- one _t for any number of time derivatives
```
⇒ WL `D[...,{time,2}]` → `u_t_t` → `u_t` (order lost), while PY `u_tt` (single token `tt`, unrecognized as a
jet) passes through as `u_tt`. ASYMMETRIC time-order loss. Both engines emit second time derivatives only in
the KINETIC/inertia terms: PY `u_tt`/`e_W_tt` (`scripts/S11c_b_brane_operator_sympy_audit.py:248-252,1474,1490,
1578-1579`), WL `{time,2}` (`mathematica/S11c_b_brane_operator_mathematica_audit.wl:837-838,1658-1660`).

**Blast radius (computed):**
```
$ grep -oE '[a-zA-Z_]+_tt\b' ~/.s11_build/S11c_b_reconcile_run.out | sort | uniq -c   →  4 e_W_tt (only)
$ grep -n 'e_W_tt' ~/.s11_build/S11c_b_reconcile_run.out   →  all 4 in SLAB_OPERATOR/THICKNESS_ROW
$ grep '^A_minus_B COUPLING_KERNEL ' ... | grep -c '(u_[0-9]_t|e_W_t|theta_t)'   →  0
```
- **COUPLING_KERNEL (the first physics number) is tt-free** — built from the static energy catalogue; UNAFFECTED.
- The defect affects only the diagonal INERTIA of the SLAB brane operator (THICKNESS_ROW / U_MOMENTUM_ROWS) and
  the KINETIC term-origins bookkeeping — standard/inherited inertia, not the variable-coefficient coupling.
- The erased order CANNOT be recovered from comparator-captured operands (the info is lost pre-capture), so a
  correct fix requires repairing `canon_jet_name` and REGENERATING the comparator + v1 artifacts + re-legging.

## Disposition
Precision defects fold into the directive. The comparator time-order defect is a SCOPE/CORRECTNESS FORK (fix +
regenerate + re-leg the committed comparator, vs adjudicate the tt-free coupling-kernel number now and OWE the
canon_jet_name fix) — escalated to the user (rule 11: scaling is the user's call; reopening a committed
instrument is significant). Three rounds each drew new issues ⇒ per rule 15 /
[[feedback_delegate_build_when_prose_directive_repeats]], the directive rewrite is a candidate to delegate to
Codex once the scope is set.
