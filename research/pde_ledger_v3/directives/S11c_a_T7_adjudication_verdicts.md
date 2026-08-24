# S11c-a T7 — step-1 ADJUDICATION VERDICTS (which engine is correct, per deep family)

User chose FULL RECONCILE (2026-08-24). Step 1 of the plan: adjudicate each deep divergence from the
step-0 matrix (`S11c_a_T7_adjudication_matrix.md`, `3c7f9137`) against the CLOSED spec
`S11c_a_SHARED_PHYSICS.md` (`2926c71c`), deciding **which engine is correct** for each. Every verdict
rests on BOTH a decidable computation AND a spec citation (`[[feedback_one_engine_fix_is_a_spec_question]]`).
Commands + stdout: twin `_measurements/…_verdicts.md` and `_measurements/s11ca_t7_adjudication/`. This is
the DECISION LIST for the two engine patches ⇒ **two legs done (Codex + Grok), folded** — both reproduced
every engine ASSIGNMENT (no reversals) and corrected severity/scope/citations; each fold re-verified by
the orchestrator (rule 13). See "Round-1 leg fold" at end.

⛔ Fixes land on BOTH engines. ⚠ Downstream PATCH directives specify WHAT to compute, never an expected
value (rule 5); this document names the correct engine and the structural fix, not any target output.

## Inputs (hash-locked)
PY `~/.s11_build/S11c_a_sympy_engine.out` sha `6386471…`; WL committed `.out` sha `82062bd…`; spec
`2926c71c`; engines PY `9b6438fa`, WL `ddecdbc2`.

---

## VERDICT B — DENSITY axis: PY CORRECT on every divergent tag; WL wrong both ways
**Decidable test** (`_measurements/s11ca_t7_adjudication/density_adjudication.{py,stdout}`): within the
engine that emits both reps, pair cases by (branch,face,dof) and diff the two rep expressions.
- WL `KINEMATIC_BALANCE`, `RELATIVE_FLUX`: **IDENTICAL** across reps (0 differ / 8 same) ⇒ rep-INDEPENDENT.
- PY `PROJECTION_SHAPE_DERIV`, `_DYNAMIC_OPERAND`, `_RESIDUAL`, **`_TERM_ORIGINS`**: **DIFFER 4/4** ⇒
  rep-DEPENDENT. `PROJECTION_STATIC_OPERAND`: IDENTICAL (the dependence is dynamic).

**Spec.** §3b:351 `J_s^α ≡ ρ_m(v_bulk−v_face)·n̂` — flux/kinematic use the bound constant `ρ_m`
(`LEDGER['rho_m']`), NOT the density REPRESENTATIVES (ρ_4D vs ρ_br). Projection uses `ρ_4D` (§1b:76),
which is representative- and anchoring-dependent (§3a:315); §4:398 requires both reps "wherever density
enters." (⚠ fold: §2b:229 names only **T-g/T-h/T-i** — it does NOT cover T-f/projection; projection's
rep-dependence is supported by §1b/§3a/§4, not §2b:229.)

**Verdict.** PY's per-quantity density policy is correct. **WL is wrong both ways** — but the two WL errors
are of DIFFERENT depth:
- **WL kinematic/flux redundant axis → SHALLOW (fold, Codex).** WL's two rep expressions are identical and
  the physics is correct; the density axis is merely a redundant label. DROPPING it is schema
  normalization, ⛔ not a computation-layer fix. Handle at the emit/reconciliation stage (light review).
- **WL projection missing axis → COMPUTATION fix (full review).** WL emits no density axis on the
  rep-DEPENDENT projection objects (`PROJECTION_SHAPE_DERIV`, `_DYNAMIC_OPERAND`, `_RESIDUAL`,
  `_TERM_ORIGINS`) ⇒ WL is MISSING rep-dependent cases; WL must add the axis and emit both reps.
  `PROJECTION_STATIC_OPERAND` is rep-independent — adding a rep axis there is redundant (match PY's cluster
  or leave; not load-bearing).

## VERDICT C — VIRTUAL-WORK matrix: WL CORRECT (full grid); PY missing required pairing CASES
**Decidable test** (`_measurements/s11ca_t7_adjudication/virtualwork_offdiagonal.{py,stdout}`): WL's 8
off-diagonal `(physical DOF ≠ virtual DOF)` cases have NONZERO `SHAPE_DERIVATIVE.EXPRESSION`; PY (`:919-924`)
emits only the 8 diagonal cases. ⚠ **Fold (both legs, verified):** for fixed (branch,rep,virtual DOF) the
two WL physical-DOF rows are BYTE-IDENTICAL — the entire emitted record, not just `SHAPE_DERIVATIVE`. WL's
full grid is thus **physical-DOF-redundant** (2 distinct expressions per branch/rep, indexed by the virtual
DOF).

**Spec.** §4 T-d:419 — "…Which virtual-displacement pairings occur is **part of the computation**." The
engine must COMPUTE all physical×virtual pairings and let the result show which occur; PY pre-decides the
diagonal.

**Verdict.** **WL is correct** (computes the full grid). **PY is nonconforming** — it hardcodes the
diagonal, missing the REQUIRED PAIRING CASES (⛔ not "missing real physics" — the off-diagonals are
physical-DOF-redundant copies; PY's diagonal already captures both distinct virtual-indexed expressions).
This is a case-coverage / emission-conformance defect (§4:419), not deep physics.
**Fix (PY, emission/computation, review at full depth):** compute and emit `δ_v𝒲_bulk^α` for every
(physical DOF, virtual DOF) pairing. ⭐ **This also CROSS-CHECKS WL's physical-DOF-redundancy** — currently
a SINGLE-ENGINE result — turning it into a dual-engine measurement (rule 1/rule 6: keep the full grid in
both so they can disagree; PY's independent off-diagonals confirm or refute WL's redundancy).

## VERDICT H.1 — CONTROL coverage: PY CORRECT (covers all); WL missing 5 quantities
**Fact** (matrix fold + twin): WL form-control quantity set = 13 ⊂ PY = 18; WL-only = ∅; PY-only =
{`EVOLUTION_TERM_ORIGINS`, `PROJECTION_STATIC_OPERAND`, `PROJECTION_DYNAMIC_OPERAND`, `PROJECTION_RESIDUAL`,
`PROJECTION_TERM_ORIGINS`}. Same 5 absent from WL's uniform-limit.

**Spec.** §4:427-429/:438 enumerate those five as separate emitted objects; §5b:490 "For every T-object …",
§5c:497 "For each S11c-a object …". ⚠ **Fold (both legs): the spec is NOT under-specified** — the five are
separate T-objects and every one must be form-ablated + uniform-tested; no §5b/§5c addendum is needed
(my earlier "possible clarification" caveat is withdrawn).

**Verdict.** **PY correct; WL under-covers.** **Fix (WL, control-layer, full review):** extend WL's §5b
form ablation and §5c uniform limit to the 5 missing objects.

## VERDICT BG — BACKGROUND_STATE: WL CORRECT; PY missing BOUNDARY_LOADS only
**Fact** (fold, verified): WL `BACKGROUND_STATE` carries `BOUNDARY_LOADS`. PY `BACKGROUND_STATE` **already
emits** `θ⁰=0`, `V_0_*=0`, `J_0_*=0`, `A_0_*=0` (the §2d zero-conditions) but carries NO boundary loads
(holds appear only in `ADMISSIBILITY_PREMISE`). ⚠ **Fold: my earlier "PY carries none of the zero-condition
fields" was FALSE** — PY has the zeros; the gap is BOUNDARY_LOADS only.

**Spec.** §2d:251 `𝔅⁰ ≡ {…, boundary loads}`, :267-271 "Emit the state … `S11CA_BACKGROUND_STATE`."

**Verdict.** **WL correct; PY under-emits the boundary loads.** Supplied-field gap (bookkeeping), not
computed physics. **Fix (PY, emit-layer, light review under a coverage check):** add the supplied
boundary-load fields to `S11CA_BACKGROUND_STATE`; ⛔ do not duplicate the zeros it already emits.

---

## Shallow reconciliations (NOT engine-physics; resolved at the emit-patch / comparator stage)
No engine has wrong physics here; a comparator needs declared bridges. Handled after the engine patches,
under the provenance-manifest gate (plan §2-3):
- **Leaf representation (Family F):** PY graded `(bg, ε·deriv, total)` / baked-in `epsilon_shape` vs WL
  coefficient. Declared algebraic reconstruction identity (coefficient+order → graded value); ⛔ not
  relabeling. Only FACE_NORMAL/CONORMAL/FACE_MEASURE carry the explicit triple.
- **WL-only provenance leaves** (`EXACT_SOURCE`, `GRAPH_EVALUATED_SOURCE`, `ORIENTATION_OBJECT`,
  `SLAB_VELOCITY`, closure operands): WL-only, ⛔ not reconstructible ⇒ excluded from the residual, logged.
- **FIELD-explosion (FACE_SHIFT):** WL 10 field components ↔ PY 4 groups — lossless derivative-level
  decomposition; declared reconstruction identity.
- **ORIGIN repartition (EVOLUTION_TERM_ORIGINS 4↔3; PROJECTION_TERM_ORIGINS):** declared mapping of the two
  origin partitions; ⛔ not itemisation of one set.
- **WL kinematic/flux redundant density axis (Verdict B):** drop at the emit stage (schema normalization).
- **BACKGROUND_DENSITY_MAP branch axis:** ⚠ **adjudicated (legs split; §2b governs).** §2b:228-230
  constructs `Σ_E⁰(y)` on the pre-anchoring `y` and says "emit the **two** computed maps" (per
  REPRESENTATIVE); the map is branch-independent (PY LAB≡MATERIAL payloads, verified). So **2-per-rep is
  spec-faithful (WL correct); PY's branch axis is REDUNDANT** — a shallow PY over-emission, ⛔ NOT a WL gap.
  (§2c:243 "both branches emitted separately" governs the anchoring-DEPENDENT profiles, not this
  pre-anchoring map.) Reconcile by dropping PY's redundant branch axis or bridging.
- **Encoding:** FACE_MAP `±1`↔`PLUS/MINUS`; DIRECTION `1/2/3`↔`DIRECTION_n`; `VALUE↔EXPRESSION`,
  `MULTIGRADE↔MULTIGRADE_EPSILON_ETA_SIGMAW`; positional↔pipe keys.

## Patch plan (after this fold is committed)
1. **WL patch** (full review): Verdict B projection (add density axis to `PROJECTION_SHAPE_DERIV`,
   `_DYNAMIC_OPERAND`, `_RESIDUAL`, `_TERM_ORIGINS`) + Verdict H.1 (extend §5b/§5c to the 5 objects). The
   WL kinematic/flux redundant-axis DROP is schema normalization — do it in the same patch but it is not a
   physics change. Directive specifies WHAT to compute, ⛔ no expected values.
2. **PY patch** (full review): Verdict C (compute+emit the full virtual-work grid) + Verdict BG (add
   supplied boundary-load fields to BACKGROUND_STATE) + drop PY's redundant BACKGROUND_DENSITY_MAP branch
   axis (schema).
3. Re-run both engines (WL ~23 min, watch orphan kernel; PY ~3 min; ⚠ PY emit feeds `export_candidates` —
   decouple/preserve/review the regenerated export).
4. Emit-patch / comparator bridges for the shallow reconciliations under the provenance-manifest gate.
5. TRIVIAL join+residual comparator (frozen T7 contract) + 2 legs → freeze + run → cross-engine result.
6. Step record (carry all verdicts + bridges) → family card. ⭐ pin the schema FORWARD for b–e.

## Round-1 leg fold (provenance — both legs, orchestrator-verified)
Codex `~/.s11_build/s11ca_t7_independent_check.{py,stdout}` (+ doc audit); Grok
`~/.s11_build/grok_leg_verdict_*.{py,stdout}`. Both REPRODUCED every engine assignment (B PY-correct,
C WL-correct, H.1 PY-correct, BG WL-correct) — no reversal. Folded corrections (each re-verified against the
streams/spec): B adds `PROJECTION_TERM_ORIGINS` (4/4 differ) + reclassifies WL kinematic/flux drop as
shallow + fixes the §2b:229 miscite; C severity → "missing required pairing cases" (off-diagonals
physical-DOF-redundant, full record byte-identical) + rule-6 cross-check framing; H.1 confirmed, no
addendum; BG narrowed to BOUNDARY_LOADS (PY already emits the zeros); BACKGROUND_DENSITY_MAP adjudicated
(legs split → §2b governs, PY branch redundant). rule-2: renamed the B script to its cited path, added
executable H.1/BG commands to the twin. rule-5: removed the literal output-expression fragments both legs
flagged.

## Open items carried to the patches / comparator (⛔ do not assume — measure)
- Verdict C: PY's independently-computed off-diagonals will confirm or refute WL's physical-DOF-redundancy.
- Verdict B: `PROJECTION_STATIC_OPERAND` rep-axis is redundant (keep or drop — a schema choice).
