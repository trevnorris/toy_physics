# S11c-b slab_operator carrier — committed ablation-harness build directive (orchestrator-written knife-list)

## Role and authority
Build **two committed ablation harnesses** — one per S11c-b engine — that certify the **pressure-slot carrier**
is *genuinely computed* (nothing hand-typed / frozen / tautological). This directive is the **orchestrator-owned
knife-list** (per `docs/development_pipeline.md` §4). ⛔ The builder implements exactly these knives against the
live engines; ⛔ **the builder does not choose, add, or drop targets/knives.** A short manifest in each harness
*records* this list; it does not originate it.

**Governing standard:** `docs/development_pipeline.md` §4 "Committed ablation harnesses". Read it. The
non-negotiables it fixes: wrap the live engine (never re-derive); orchestrator-owned knives; construction site ≠
observation site; print the `(baseline, corrupted, diff)` triple then guard (⛔ no `PASS` tag); tagged-payload
drift guard (not full stdout); commit the transcript; ablate the harness itself.

## The target object (why this, why now)
c2's N6 binds the S11c-b slab operator's **pressure-slot carrier** `∂(slab rows)/∂(δp±, ∂_w δp±)|_{P=0}` (plus
the #90 `closure_shape_deriv` fold) as a **supplied/unfalsifiable** premise — it *becomes* c2's `C_E/C_M`. Two
review legs (codex-sol + Grok, 2026-09-07) flagged it as the highest-value early retrofit: c2 rests on it, yet
its genuineness has only ephemeral per-engine ablation and its cross-engine residual is deferred (≥64 GB). This
harness makes that premise a **proven-computed** one.

**⛔ Scope — per-engine GENUINENESS ONLY.** These harnesses prove each engine *computes* the carrier (controls
bite). They do ⛔ **NOT** resolve the deferred cross-engine questions (the kinetic ∓K and face-force ± sign
conventions, #90's closure-fold sign, the uniform-limit Λ survivor) — those remain ≥64 GB-deferred and are the
comparator's, not this harness's. A knife that "bites" says the object moved under a structural corruption; it
says nothing about which engine's sign is right.

## Engines (both — one harness each)
- **SymPy** `scripts/S11c_b_brane_operator_sympy_audit.py` (`PY_S11CB_*`).
- **WL** `mathematica/S11c_b_brane_operator_mathematica_audit.wl` (`WL_S11CB_*`, blind).

Note: the SymPy engine already ships `corrupt_virtual_constraint_source` (:2268) — an existing in-engine ablation
hook. The harness formalizes and commits systematic knives of that kind; reuse the engine's own hooks where they
match a knife below, ⛔ but never invent a hook the engine's construction doesn't support.

## Harness architecture (identical shape for both engines)
1. **Baseline = the canonical engine, unchanged.** Run the committed engine as-is; capture its emitted
   **tagged payloads** (the `PY_S11CB_*` / `WL_S11CB_*` objects), ⛔ not raw stdout.
2. **Each knife = a copy of the engine with exactly ONE construction site patched** (the FORM perturbation
   below), run, tagged payloads captured. ⛔ The patch is mechanically applied to the production source
   (copy-and-patch, or the engine's own hook); ⛔ never a reimplementation of the physics.
3. **Emit, per knife, the triple** `(baseline_payload, corrupted_payload, diff)` for every object in that knife's
   **impact cone** — both the must-MOVE objects and the must-NOT-move objects. **Print the triple, then guard**
   (E1). ⛔ No `PASS`/`FAIL` payload; ⛔ no assertion that a diff is zero/nonzero. A zero diff on a must-MOVE
   object, or a nonzero diff on a must-NOT-move object, is the *finding*, surfaced by the printed numbers.
4. **Drift guard:** the baseline tagged payloads must match the committed engine's tags (copy-identity, ⛔ not
   full stdout — exclude progress/RSS/timing/banners). Emit the guard as a printed comparison, ⛔ not an assert.
5. **Single-case** (one anchoring × density, e.g. `LAB_HELD` / `RHO4_CONSTANT`) is sufficient — a bite is a
   bite. ⛔ Do not run all 4 cases (the heavy ≥64 GB territory is the cross-engine residual, not this).
6. **Construction site ≠ observation site:** every patch is on the object's *construction*; ⛔ never on the emit
   or diff layer (mutating the emitter makes every payload "bite").

## The knives — carrier-only (K1–K4). Each: SITE → FORM perturbation → impact cone.
Line numbers are current-tree anchors; the builder confirms the exact site by reading, ⛔ never trusts the number blindly.

**K1 — Freeze knife (rule-17 / the 26=26 test).** SITE: the background-gradient spurion — SymPy
`DERIVATIVE_MAP[i][W_bg]=grad_W[i]` (:754), `spurion=("BACKGROUND_FIRST_JET",…)` (:1519), `basis_second` (:1453);
WL `widthSpurion=gradient[anchoredWidth]/WZero` (:721), `modulusSpurion` (:722). FORM: freeze the live background
first-jet (`grad_W → 0` / spurion → 0) — leaves the variable-coefficient family for the uniform one.
- MUST-MOVE: the 15 ∂W_bg + 15 ∂μ_R,bg spurion basis rows, and every carrier object built on them (§3a basis
  40 → 10).  MUST-NOT-move: the 10 uniform-sector terms.

**K2 — Constraint-fold knife (pin B; the reaction is computed, not typed).** SITE: SymPy
`corrupt_virtual_constraint_source` (:2268) / `constraint_fold_from_source` (:2444); WL
`virtualConstraintSource[route,branch,density,…]` (:947). FORM: structurally alter the virtual-constraint source
(a form change, ⛔ not a rescale).
- MUST-MOVE: the constraint-reduced U/e_W slab rows.  MUST-NOT-move: the bulk θ mass-balance base.

**K3 — Closure / face-response knife (the direct δp carrier, #90).** SITE: SymPy `closure_shape_deriv`
(:389; folded into the θ-row at :2214 / :2820), `face_generalized_force_rows` (:2135), `Lambda_A_0`/`Lambda_V_0`
(:408-409); WL closure/face route (`dimensionLambdaA/V` :115-116 and its face-response construction). FORM: alter
the closure shape-derivative *structure* feeding the slot carrier.
- MUST-MOVE: the θ-row closure term + the Λ_A/Λ_V face-response channels.  MUST-NOT-move: the bulk mass-balance
  base, and **Λ_X** (traction — a DISJOINT source; that disjoint channels stay put is the #90 test).

**K4 — Slot-structure knife (the carrier's own δp± index wiring).** SITE: the δp±/∂_wδp± slot dependence —
SymPy `delta_p_{face}` (:486), `d_w_delta_p_{face}` (:489), and where the rows contract the closure over these
slots. FORM: break the slot structure (swap the +/− face pairing, or the ∂_w-slot ↔ value-slot mapping) — a
structural mis-wiring, ⛔ not a coefficient.
- MUST-MOVE: the extracted carrier `∂(rows)/∂(δp±, ∂_wδp±)`.  MUST-NOT-move: the P-independent operator base.

**Coefficient exceptions:** none. All four are genuine FORM. If, during the build, a channel proves FORM-blind
(a FORM knife cannot reach it), ⛔ do not invent a coefficient knife — **stop and report it to the orchestrator**;
a named coefficient exception is added to THIS list by the orchestrator, never by the builder.

## Deliverables
- `scripts/S11c_b_carrier_ablation_harness_sympy.py` — the SymPy harness.
- `mathematica/S11c_b_carrier_ablation_harness.wl` (with a small Python/driver if needed) — the WL harness.
- Each opens with a **manifest** header mirroring K1–K4 (site, object, why load-bearing, the FORM change, the
  must-MOVE / must-NOT-move cone). ⛔ The manifest records this list; it does not originate it.
- Each run's exact invocation + source/output digests + the literal `(baseline, corrupted, diff)` transcript go
  in `_measurements/S11c_b_carrier_ablation_harness_{sympy,wl}.md` (E1 — the harness is evidence only once run).

## Three script clauses (verbatim, non-negotiable)
1. The harness may **PRINT** computed objects (payloads, diffs). It may ⛔ NOT state conclusions — no `PASS`,
   no verdict, no "bites"/"fails" tag. Interpretation is the orchestrator's, from the printed triples.
2. **Print operand and residual, then guard.** Emit `baseline`, `corrupted`, and their `diff`; a residual
   asserted zero/nonzero carries no information.
3. Interpretation belongs to the review/step record, ⛔ not the script.

## Build & review discipline
- **Author:** a fresh `gpt-6-astra` (code; WL + Python). The knives are orchestrator-owned and G2-cleared, so
  astra implements without choosing targets. Each engine's harness may be a separate invocation; ⛔ neither
  receives the sibling engine's output or the comparator.
- **Reviewed** (astra-written script → G1): a **fresh Claude agent + Grok**, review-until-clear, serialized if
  both ablate Mathematica. **Verify the harness by ablating the harness itself** — a coefficient-rescale or a
  dead-path mutation of a knife must **fail to report a bite** (else the harness is a self-report).
- Commit the reviewed harness baseline before any repair overwrites it; commit the run transcript with it.
- ⛔ Never run two memory-heavy CAS jobs concurrently (30 GB box; serialize).
