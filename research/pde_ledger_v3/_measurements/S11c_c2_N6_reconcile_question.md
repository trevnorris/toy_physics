# S11c-c2 N6 cross-engine reconcile — THE QUESTION (orchestrator-framed; ⛔ disposition NOT adjudicated)

⚠ This document frames the reconcile disposition of the SURFACED cross-engine residuals from the N6 T7
comparator run. It is **orchestrator-written** (the QUESTION is mine to own; the collapse instrument, if any,
is CODEX-written + G1-reviewed — CAS-authorship rule `6f8dbd34`). It **decides nothing**: per §N8/T7 and
[[feedback_reconcile_representational_bridge]], a NONZERO residual is ⛔ NOT yet a disagreement. Its role is to
be the orchestrator's frozen framing that the CODEX-written **instrument build directive** carries into its own
2 decision legs (G2).

⛔ Nothing here pre-decides that a surfaced residual reconciles or that it genuinely differs.

⭐ **Folded once (2026-09-09) from a convergent 2-leg question-vet** (Codex-sol + Grok, both QUESTION
NEEDS-WORK → A; record `_measurements/S11c_c2_N6_reconcile_question_vet.md`; leg reports in `_legs/`). The core
question (the OPERAND gap is the real open question, matched leaves are tractable) survived; the folds tightened
§1/§2 (sign + PIT qualification + relabel the trivial zeros), §3(b) (softened inference; BASELINE is a
nominal-control duplicate), §4 (whole-object "identities" → a predeclared PRIMITIVE bridge dictionary), §5
(remove σ_W binding), §3(d) (census is load-bearing control/premise, not dismissible metadata).

---

## 1. What the comparator surfaced (grounded — committed tally `run_tally.txt`, run record `run_data.md`)

The comparator (`scripts/S11c_c2_N6_cross_engine_comparator.py`, run `a094b284`) joins the blind WL `.out`
(`ae73b884`) against the 3 SymPy `.out` streams (`7c0790ab`) by object name, computes a symbolic residual per
matched key, and PRINTS it (three-valued; decides nothing). ⚠ **Sign:** the residual is
`A_minus_B = SymPy − WL` (operand A = the SymPy `.out`s, B = WL; `residual` returns `py_value − wl_value` —
`S11b_cross_engine_comparator.py:826`, `S11c_c2_N6_cross_engine_comparator.py:794`). Sign does not affect
zero/nonzero classification.

**(i) MATCHED residual-of-residuals is 0 (⛔ NOT itself the disposition — see §2 for why this is trivial):**
these are cross-engine differences of objects that are *already 0 within each engine*, so 0 cross-engine says
nothing about operand agreement. Per `run_tally.txt` SYMBOLIC channel (no-nonzero-found, PIT-qualified):
- `N6COV_R_COV`, `R_COV_BASELINE`, `R_COV_CONTROL_DELTA` — 160 ZERO / 0 NONZERO each.
- `N6COV_SOURCE_CONTROL_DELTA` — 160 / 0 (⚠ trivial by construction — see §3b: ACTUAL≡BASELINE at shipped
  settings).
- `N6RC_CARRIER_BRIDGE_RESIDUAL` — 320 / 0.
- `N6RC_ADVECTION_ABSENCE` 6/0; `N6RC_FROZEN_RELATIONS` 8/0.
- ⛔ The per-engine SymPy facts `R_cov=0` and carrier-bridge=0 (Reading B, `_measurements/S11c_c2_N6_RESOLVED.md`)
  are PER-ENGINE, PIT-qualified; the cross-engine 0 here is `(0)−(0)`, ⛔ not independent operand agreement.

**(ii) SURFACED computed residuals — matched keys, some NONZERO (the disposition's subject):**
- `N6RC_CARRIER_EULERIAN` / `CARRIER_MATERIAL` — 280 ZERO / **40 NONZERO** each.
- `N6COV_SOURCE_ACTUAL` / `SOURCE_BASELINE` / `SOURCE_PREDICTED` — 16 ZERO / **76 NONZERO** each.
- `N6COV_FROZEN_PHI` — 42 ZERO / **18 NONZERO**.
- control/premise (⛔ NOT dismissible bookkeeping — §3d): `N6COV_ACTUAL_CONTROL_PARAMETERS` 0/**4**;
  `N6COV_PHI_DOMAIN_CENSUS` 0/**4** (+44 BOOL rejected).

**(iii) ENTIRELY UNMATCHED cross-engine (different WL block/kernel/component vocabularies; also a `(face,grade)`
vs `{face,wave,grade}` schema non-join for the RC sources — `reconcile:105-111` vs `.wl:708-711` — this is a
schema non-join, ⛔ not a measured residual):** `N6RC_R_N6`, `SPLIT_CHECK`, `SPLIT_SUM`, `*_OPERAND`,
`*_CHANNEL`, `DIMENSIONS`, all guards, and the reconcile-engine `SOURCE_{EULERIAN,MATERIAL,BRIDGE_RESIDUAL}`.
All UNDECIDED (pairing table). ⇒ there is **NO direct cross-engine comparison of R_N6 itself** — these stay
UNDECIDED under both scope paths (schema, not heaviness).

---

## 2. The representational split (grounded — both engine sources)

Each surfaced residual is `SymPy_imported/native − WL_blind`. The two engines build the SAME physical object by
**different routes** by design (`S9_export_chain_rebuild_directive.md:16-18` — WL imports nothing and
re-derives; SymPy imports the frozen slab + c1 exports). The routes:

| object | SymPy — `scripts/S11c_c2_N6_{covariance,reconcile}_sympy.py` | WL (blind) — `mathematica/S11c_c2_N6_mathematica_audit.wl` |
|---|---|---|
| `C_E` (carrier, Eulerian) | `pressure_coefficients(expanded_rows(inputs.slab[α,ρ]))` (`reconcile:196-197`) — from the **IMPORTED** slab operator | `carrier[eulerianSlabFace[lawE,·]]` (`:327,336,866,868`); carrier = `D[rows,p]/.slotZero`; rows blind from `graphGeometry` |
| `C_M` (carrier, Material) | `pressure_coefficients(face_factory(…MATERIAL…))` (`reconcile:88-92,205`) — native material face | `carrier[rowsM]` (`:868`) — blind material rows |
| `SOURCE_{ACTUAL,PREDICTED,BASELINE}` | `source_terms(inputs,α,ρ,s, μ, V)` (`cov:191-193`); μ_pred = `mu_e.subs(Φ)`, μ_actual from `constitutive` | `sourceBind[sourceSolve, μ, V, density3]` (`:871-872`), μ from `materialAmplitude`/`muE`, V from `velocitiesE/M` |
| `Φ` (`FROZEN_PHI`) | `prolonged_phi` (`cov:63-129`) — N4 map + `DERIVATIVE_MAP` + `total_derivative`, jets to rank ≤2 | `phiMap[muE, ak·a, hk·h]` (`:236,849`) — blind graph advection a_ρ + thickness h |

⭐ **Structural fact (grounded, PIT-qualified).** Within EACH engine the carrier bridge `C_E − C_M` gives NO
nonzero found (SymPy `reconcile:210,254`; WL `:874,876`) and `R_cov = ACTUAL − PREDICTED` gives no nonzero
found (SymPy `cov:194`; WL `:878`); the *vanishing* is the committed PIT/tally, ⛔ not the subtraction line. Each
such object being 0 in each engine makes its cross-engine difference `(0)−(0)` — **0 trivially**. That says
NOTHING about whether the OPERANDS (`C_E`/`C_M`; ACTUAL/PREDICTED/BASELINE; Φ) agree cross-engine. That gap is
exactly why the operands are surfaced separately (§1.ii), and is the whole subject of this disposition. ⚠ Read
every "0" here as **"no nonzero found at the retained rectangle, under the adopted PIT-qualified disposition"**
(`cov:8`, `_measurements/S11c_c2_N6_RESOLVED.md:23`), ⛔ not an exact symbolic zero.

---

## 3. THE QUESTION — per surfaced channel (⛔ each is a question, not a claim)

**(a) CARRIER (40 nonzero, `N6RC_CARRIER_EULERIAN`/`MATERIAL`).** Under the frozen primitive bridge dictionary
(§4), do the matched carrier OPERAND residuals vanish coefficientwise per retained grade — SymPy's imported-slab
`C_E` vs WL's blind graph-geometry slab-face `C_E`, and SymPy's native `C_M` vs WL's blind `C_M`? The
within-engine and cross-engine bridge `C_E=C_M` already gives no nonzero (§1.i); the open question is the
carrier OPERANDS. ⛔ The bridge vanishing does NOT dispose of these 40. (c1 analogue: the DtN KERNEL closed
cross-engine while the raw `dtn_operator` whole-form stayed UNDECIDED — `[[project_s11c_c_state]]`.)

**(b) SOURCE + Φ together (76 + 18 nonzero) — THE substantive question.** Applying **one** frozen bridge
dictionary UNCHANGED, do the matched operand residuals vanish coefficientwise per grade **separately** for
`SOURCE_ACTUAL`, `SOURCE_PREDICTED`, and `FROZEN_PHI`? ⚠ `SOURCE_BASELINE` is a **nominal-control DUPLICATE** of
`SOURCE_ACTUAL` at shipped settings (`kappa_a=1, kappa_j=0` ⇒ μ_actual = μ_baseline; `cov:34-35,140-145`; WL
`:18,845`; confirmed cross-engine by `SOURCE_CONTROL_DELTA` 160/0), ⛔ NOT a third independent discriminator —
report it as the control duplicate. Source and Φ must be tested under the SAME chain (splitting them lets a map
error and its image cancel — `PREDICTED = μ_E.subs(Φ)`, `cov:132-134,185-191`; WL `:849,853,871`).
- If the operands collapse under the one unchanged dictionary ⇒ the engines AGREE **on those matched keys, under
  that frozen bridge** (⛔ NOT closure-wide; c1 over-claimed closure-wide AGREE and walked it back). This is
  strictly stronger than `R_cov=0` (a within-engine zero) and is the cross-engine corroboration that the
  source/map are right.
- If they do NOT collapse ⇒ **UNDECIDED** (the declared bridge did not reconcile them), ⛔ NOT a "genuine
  disagreement" — non-collapse establishes a genuine difference ONLY if the bridge's COMPLETENESS is
  independently established (§6). The `R_cov=0`-with-different-operands case is the two-errors-cancel pattern
  [[feedback_matching_number_is_not_evidence]] that this test EXISTS to expose, but exposing a candidate is not
  adjudicating it.

**(c) [merged into (b)]** — Φ is tested in the SAME chain as SOURCE; see (b). ⭐ This is the cross-engine
corroboration of per-engine premise caveat 1 (is Φ physically correct — R_cov cannot exclude an error SHARED by
the declared Φ and both of its own routes; `_measurements/S11c_c2_N6_RESOLVED.md`), because WL derived its map
BLIND and SymPy prolonged the supplied N4 map.

**(d) control/premise (`ACTUAL_CONTROL_PARAMETERS` 4, `PHI_DOMAIN_CENSUS` 4+44BOOL) — ⛔ LOAD-BEARING, not
dismissible metadata.** `ACTUAL_CONTROL_PARAMETERS` records the coefficients that SELECT the material amplitude
+ inserted junk (`cov:137,146-153`; WL `materialAmplitude`/`MATERIAL_NORMAL` `:845,977`); `PHI_DOMAIN_CENSUS`
records jet coverage/uncovered/max-rank and a failure ABORTS construction (`cov:105-115,125`; WL `:850`). Audit
the 8 nonzero leaves as **production-control equivalence** + **domain-coverage equivalence** (⛔ not "do they
carry physics"): e.g. is the `junk_symbol` diff pure spelling (`n6cov_J_mu` vs `junkMu`), and does WL emit a
THICKNESS actual-control that PY does not (PY `actual_amplitudes` retags only theta-advection + junk,
`cov:137-154`)? If a census leaf reflects differing μ-jet coverage, it routes to (b), not to bookkeeping.

---

## 4. The PREDECLARED PRIMITIVE BRIDGE DICTIONARY (⛔ frozen + reviewed BEFORE testing; ⛔ NO whole-object equality)

⛔⛔ The collapse test may apply ONLY **primitive, independently-derivable** conventions — each grounded in what
BOTH engines actually emit. ⛔ It may NOT apply any whole-carrier / whole-μ / whole-source equality: those ARE
the questions (§3a/b), and substituting them is a blanket collapse that assumes the answer
([[feedback_handcode_comparison_never_blanket_collapse]] / rule 5; the comparator already refused to
pre-register block/kernel/component equality — `_measurements/S11c_c2_N6_comparator_directive_review.md`).

**Admitted primitive maps (each cited in both engines):**
- **jet-vocabulary bridge** — `a.grad_theta[i] ↔ b.grad_theta[i]` via the common `wave_jet`
  (`reconcile:59-69`; WL `JET_VOCABULARY_BRIDGE` `:946-948`).
- **quadratic density Jacobian** `1 + tr(∇u)` and **wave-projection degree 2** (`reconcile:65-66`; WL
  `materialAmplitude` `:246`).
- **live-density rebind** — SymPy rebinds the imported c1 density before source extraction (`diagnostic:378`);
  WL `density4`, `density3 = density4·WBg` (`:839`). ⭐ c1 made the live-density issue MANDATORY for c2
  (`directives/_measurements/S11c_c1_comparator_reconcile.md:133`).
- **background/profile jets** — `WBg → W0(1+η·w1)`, `muRBg → muR(1+η·m1)` and their σ-scaled derivatives
  (WL `:127`; SymPy `S11c_c2_selfenergy_fold_sympy_audit.py:204,281`).
- **energy-basis map** — map coefficients BY CONTRACTION IDENTITY (WL quotient basis `:211`; SymPy termwise
  energy variation `diagnostic:318`). ⛔ Mapping coefficients by contraction identity is legitimate; equating
  the resulting μ objects is NOT.
- **face-velocity + source-solve factors** — both enter the emitted source operands directly (`cov:189`; WL
  `:359,364`). (Missing from the pre-fold draft; added.)
- **Φ spelling** as both engines construct it (`prolonged_phi` / `phiMap`) — ⛔ TEST equality of the maps, do
  NOT assume it.
- **c1 name/assumption maps** (`ε` placement, `omega` realness, Fourier-of-derivative) — ⛔ ONLY where those
  atoms ACTUALLY appear in these residuals; ⛔ do not import the whole c1 table unread.

**⛔ Excluded (each would make a collapse vacuous or leave the retained rectangle):**
- whole-carrier / whole-μ / whole-source equality (the §3 questions themselves);
- the **LAB_HELD ↔ MATERIAL_ADVECTED** cross-case chart identification — anchoring is a retained comparator key
  (`comparator:86`); a cross-case map is not a WL−SymPy bridge;
- **σ_W binding** `σ_W = η·W0/L_W` — it IDENTIFIES the two independent grade axes ⇒ leaves the retained
  rectangle (a cousin of the L-CAS `σ_W→0` proxy); see §5;
- **default on-shell / Fourier-kernel** identities — these are c1 KERNEL identities; SymPy removes the
  DtN/resolvent before source extraction (`diagnostic:379`) and WL's kernel construction begins only after
  `sourceBind` (`:380`), so admit an on-shell factor ONLY if an inspected residual actually carries a
  `q² − (ω²/c² − |k|²)` atom after the name maps.

⭐ The dictionary is applied UNCHANGED to ACTUAL and PREDICTED (and, for the carrier, to `C_E` and `C_M`). The
instrument's directive (Codex-written, 2 decision legs) fixes and freezes it before any test runs.

---

## 5. The retained order — coefficientwise, ⛔ NOT a proxy, ⛔ NOT a grade-combining binding

Both engines retain the **independent** rectangle of 4 coefficients `(η^i σ_W^j), i,j∈{0,1}`: SymPy
`GRADES = product((0,1),repeat=2)` (`diagnostic:43`), extracted independently (`diagnostic:245`); WL
`gradeIndices = Tuples[{Range[0,1],Range[0,1]}]`, `Series[…,{etaBg,0,1},{sigmaW,0,1}]` (`:123-124`). The
comparator joins `ETA, SIGMA ∈ {0,1}` as independent axes (`comparator:96-97`).
- ⛔ The collapse test must run **coefficientwise over all four grades**; ⛔ no grade-combining, no cross-grade
  cancellation, no turning `η·σ_W` into a discarded `η²`.
- ⛔ `σ_W→0` (the measured L-CAS over-clear) and any `σ_W`↔`η` binding are FORBIDDEN — retained order is WHERE
  the test runs, ⛔ not a map to apply. A physical σ/η relation may be REPORTED after the coefficientwise
  comparison, ⛔ never applied during it.
- The per-engine N6 residual lives in retained `σ_W^1` (`_measurements/S11c_c2_N6_reconcile_adjudication.md`);
  the cross-engine test must match that order.

---

## 6. What CANNOT be pre-adjudicated (surface, don't pre-decide)

- ⛔ Do NOT declare (a)/(b) reconcile before the collapse test runs and is reviewed.
- ⛔ Non-collapse is **UNDECIDED**, ⛔ NOT a disagreement — a genuine cross-engine difference is established ONLY
  if the frozen bridge's COMPLETENESS is independently proven. Collapse earns AGREE **only on the matched
  emitted operands, under that bridge** — ⛔ never closure-wide.
- ⛔ Do NOT re-litigate the per-engine N6 covariance resolution (`R_cov=0`, Reading B — stands, user-adopted
  `d21c8ff5`) or the cleared comparator/knife design.
- ⛔ Whatever does not collapse, and every UNMATCHED family (`R_N6`, channels, RC sources — §1.iii), stays
  SURFACED as UNDECIDED and is carried to the c2 step record — as c1 surfaced its UNDECIDED list. ⛔ NOT thereby
  a disagreement, ⛔ NOT thereby closed.
- ⛔ The 3 per-engine premise caveats (Φ physical-correctness; V face-velocity transform `V_E≡V_M`;
  extracted-block leakage) are carried OPEN regardless of the disposition outcome.

---

## 7. Disposition architecture (governing — ⛔ not negotiable)

1. **I own the QUESTION** (this document) and the final **adjudication** (G4).
2. ✅ **E1 question-vet DONE** — 2 convergent legs (Codex-sol + Grok), both QUESTION NEEDS-WORK → A; folded once
   (this document). Record `_measurements/S11c_c2_N6_reconcile_question_vet.md`.
3. **The collapse instrument is CODEX-authored (astra) against an orchestrator-written build directive** whose
    own **2 decision legs (G2)** freeze/verify the primitive bridge dictionary (§4) + the coefficientwise test +
   the extraction (a scoped `--family` comparator re-run regenerates the surfaced residuals — reproducible; the
   `.out` inputs are present on disk). The instrument PRINTS the collapsed residuals; it decides nothing. Its
   two re-review legs are a fresh Claude agent + Grok (Codex-written → not Codex-reviewed).
4. **Mechanical fact-lookups are the orchestrator's** — regenerating the reproducible comparator output and
   retrieving the surfaced residual EXPRESSIONS verbatim is a fact-lookup; the zero-test-after-substitution is
   the instrument.
5. **The disposition gets LEGS** — the c1 reconcile's correction-verify scoped my verdict TWICE
   ([[feedback_reconcile_representational_bridge]]); the c2 disposition adjudication is reviewed the same way.
6. **Then the c2 step record** — surface (⛔ not pre-adjudicate) the disposition outcome + the carried questions
   (2 S11c-b sign conventions, 6 §3d questions, c1 ENERGY); fix the misleading `I_{M→E}` "mapped-operand"
   terminology; preserve BOTH `R_N6`(=18/288, per-engine SymPy) and `R_cov`(=0, per-engine SymPy); carry the 3
   premise caveats; NO per-substep card.

## 8. Scope — path **A** (both legs, evidence-backed)

Build the collapse instrument to TEST §3/§4 collapse at retained order for the tractable channels, family-by-
family and grade-by-grade under bounded workers: the matched NONZERO leaves already ran (comparator materialized
`operand_A`/`operand_B` per matched key, `comparator:816-825`) at `peak_rss_kib 330412` / `runtime_s 10294` /
`deferred_oversize 0`, largest operand ~80 MB (`run_data.md:14-17`) — the same resource class as the finished
comparator, ⛔ not the c1 ≥64 GB regime. Only an actually-oversized individual channel falls back to surfacing
(B). The UNMATCHED `R_N6`/channels/RC-sources stay UNDECIDED under either path (schema, not heaviness). The 8
control/premise leaves (§3d) are a fact-lookup-sized semantic walk, part of A.

**Codex's exact operative reframing (adopted verbatim as the instrument's question):**
> For each fixed `(anchoring,density)` case and each retained coefficient `(η^i σ_W^j), i,j∈{0,1}`, does one
> predeclared, source-derived primitive bridge dictionary — containing field/profile-jet, live-density,
> energy-basis, face-velocity, source-solve and ε-placement conventions, but NO whole-carrier, whole-μ or
> whole-source equality — make the matched cross-engine residuals vanish separately for `C_E`, `C_M`,
> `SOURCE_ACTUAL`, `SOURCE_PREDICTED` and Φ? Apply the same dictionary unchanged to ACTUAL and PREDICTED; report
> BASELINE as the nominal-control duplicate. Collapse earns agreement only for those matched emitted operands.
> Failure to collapse remains UNDECIDED unless completeness of the bridge has independently been established.
