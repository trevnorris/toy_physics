# S11c-b slab_operator carrier — committed ablation-harness build directive (orchestrator-owned knife-list)

## Role and authority
Build **two committed ablation harnesses** — one per S11c-b engine — that certify the **pressure-slot carrier**
of the slab operator is *genuinely computed* (nothing hand-typed / frozen / tautological). This directive is the
**orchestrator-owned knife-list** (per `docs/development_pipeline.md` §4). ⛔ The builder implements exactly these
knives against the live engines; ⛔ **the builder does not choose, add, drop, or re-target knives.** ⛔ The builder
receives ONLY this directive — ⛔ no comparator, no sibling engine's output, no other context.

**Governing standard:** `docs/development_pipeline.md` §4 "Committed ablation harnesses". Read it. The
non-negotiables it fixes: wrap the live engine (never re-derive); construction site ≠ observation site; print the
`(baseline, corrupted, diff)` triple then guard (⛔ no `PASS`); tagged-payload drift guard; commit the transcript;
ablate the harness itself. ⛔ **The harness makes no selection and no interpretation** — it prints the COMPLETE
carrier (every row × every slot) for the baseline and for every knife, and their diffs; deciding what "should"
move is not the harness's job and no expected outcome appears in it or in this directive.

## The target object
The **pressure-slot carrier** of the slab operator:

> **carrier = ∂( emitted slab-operator rows ) / ∂( native pressure atoms ) , then evaluate at (all pressure atoms → 0).**

Pressure enters the emitted rows through the permeable-face laws. A downstream step binds this carrier as a
supplied premise; these harnesses exercise each engine's construction of it.

**⛔ Scope — per-engine GENUINENESS ONLY.** Exercise each engine's construction of the carrier under the listed
knives. Do ⛔ **NOT** resolve any cross-engine sign / convention question — those are the comparator's and are
deferred.

## Engines, pinned case, and the operator-only entrypoint
- **SymPy** `research/pde_ledger_v3/scripts/S11c_b_brane_operator_sympy_audit.py` (`PY_S11CB_*`).
- **WL** `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (`WL_S11CB_*`, blind).

Line numbers are current-tree anchors; the builder **confirms each site by reading**, ⛔ never trusts the number
blindly. If a named function/symbol is not where stated, stop and report — ⛔ do not retarget.

**Pinned case (both engines): route `EULERIAN`, branch `MATERIAL_ADVECTED`, density `RHO4_CONSTANT`.** One case
suffices. ⛔ Do not run all 4 cases; ⛔ do not touch the `COUPLING_KERNEL` / tower / heavy-controls path (that is
the ≥64 GB OOM territory).

- **SymPy — import and call the production function.** The 4-case emit lives in `run()` under
  `if __name__ == "__main__":` (`:5585`), so **`import`ing the module is side-effect-free.** Call
  `build_operator("MATERIAL_ADVECTED", "RHO4_CONSTANT", "EULERIAN")` (signature `:2723`). ⚠ It returns
  `(casify(operator), origins, mu_theta_value)` — `casify` (`:603-613`) turns the operator dict into a **nested
  `sp.Tuple`**, so `operator["U_BODY_BALANCE"]` does **not** work; access rows with `named_tuple_row(operator,
  "…")` (`:2584`) and the sub-payload with `named_tuple_row(row, "EXPANDED")`. Mirror the emit transform
  `retained_grade(operator)` (`:4121`) so the carrier is on the object that is actually emitted. ⛔ Do not call
  `build_kernel` or any `COUPLING_KERNEL` path.
- **WL — load DEFINITIONS ONLY, then call `evaluatedModel` directly.** ⚠ The emit driver `Do[…]` (`:2190`) calls
  `extractCouplingData` (`:2199`) and `kernelOriginsFromOrigins` (`:2203`) **before** it tests
  `S11CB_PRIMARIES_ONLY` (`:2206`) — that is the ≥64 GB OOM path, and it runs even under `PRIMARIES_ONLY` and even
  for one case. So ⛔ **do not run the emit `Do` at all.** Instead load only the **definitions** (everything up to
  the `(* Main variable-coefficient objects. *)` marker at `:1878`; `evaluatedModel` is defined at `:1317`, before
  it) — e.g. copy the engine to a temp file truncated just above `:1878`, then `Get` it — and call
  `evaluatedModel["EULERIAN", "MATERIAL_ADVECTED", "RHO4_CONSTANT"]["OPERATOR"]` (3-arg form ⇒ `corrupted=False`).
  ⛔ Never call `extractCouplingData`, `kernelOriginsFromOrigins`, or `frozenEvaluatedModel` (`:1230`, a different,
  frozen object).

## The observation object (pinned — do NOT reimplement N6)
Take the operator rows **from the live engine** (above) and form the carrier on them. ⛔ Do not reimplement N6's
row re-derivation (`face_factory`); ⛔ do not patch the emit layer. The generic per-slot differentiation is trivial
and may be written directly (it mirrors N6's `pressure_coefficients`,
`research/pde_ledger_v3/scripts/S11c_c2_N6_diagnostic_sympy.py:410-412` — read **only** as the reference for the
atoms and the reduction; ⛔ not imported).

- **Rows** (observation targets):
  - SymPy: `named_tuple_row(operator, "U_BODY_BALANCE")`, `…"E_W_BALANCE"`, `…"THETA_BALANCE"`, each reduced to its
    `"EXPANDED"` sub-payload. ⚠ These are the **post-fold, post-face** rows written at `:2998` / `:3028` / `:3042`.
    ⛔ Not the energy templates at `:2367` / `:2372` / `:2377` (`THETA_BALANCE` there is `epsilon·mu_theta`, later
    **overwritten** at `:3042`).
  - WL: from `evaluatedModel[…]["OPERATOR"]`, the `U_MOMENTUM_ROWS`, `MASS_EVOLUTION_ROW`, `THICKNESS_ROW` (the
    `evaluatedModel` operator, `:1344-1351`). ⛔ Not `frozenEvaluatedModel` (`:1256-1263`) nor the `rawModel`
    face-bundle (`:1189-1190`).
- **Native pressure atoms** (the differentiation slots; pin this exact order):
  `(delta_p_plus, d_w_delta_p_plus, delta_p_minus, d_w_delta_p_minus)`.
  - SymPy: the symbols of those four names (registered `:484-493`), retrieved by name from the live engine's atom
    table.
  - WL: the **applied** pressure functions `pressureUpper[xOne,xTwo,xThree,time]` and `pressureLower[…]`
    (`pressureField[1]`/`[-1]`, `:1014-1015`). ⚠ Differentiate the **applied** function (`D[row, pressureUpper]` on
    the bare head is 0); there is no in-plane `∂_w` pressure slot in WL, so its `d_w_*` carrier columns are
    absent — print them as not-applicable, ⛔ do not fabricate them.
- **Carrier** = `∂(row)/∂(atom)` **first**, **then** substitute all pressure atoms → 0. ⚠ Order matters (a
  self-test below checks it). ⛔ WL `truncateBackground` (`:165-170`) is an `(η,σ_W)` projector, **not** a P→0 —
  do the `→0` explicitly. In WL the emitted row is a `modelRecord` wrapper (`:1572-1581`); differentiate the
  `"EXPRESSION"` slot (or the raw `["OPERATOR"]` row), ⛔ not the wrapper Association.
- **Print the complete carrier** — every row × every slot — for the baseline and for each knife, plus the diff.
  ⛔ No selection; ⛔ no "expected"/"should move" annotation.

## Harness architecture (identical shape for both engines)
1. **Baseline = the canonical engine, unchanged.** Obtain the single-case operator via the entrypoint above; form
   the baseline carrier over all rows × all slots.
2. **Each knife = a copy of the engine with exactly ONE construction site patched** (the FORM below). Re-obtain the
   operator, re-form the carrier. ⛔ Copy-and-patch the production source; ⛔ never re-implement the physics; ⛔
   never patch the carrier extractor, the emit, or the diff layer.
3. **Emit, per knife, the triple** `(baseline_carrier, corrupted_carrier, diff)` over the complete carrier. **Print,
   then guard** (E1). ⛔ No `PASS`/`FAIL`, ⛔ no assertion that a diff is zero/nonzero, ⛔ no verdict.
4. **Drift guard.** The baseline carrier must match a fresh canonical run: compute the retained single-case object
   (or its digest) from the imported `build_operator` / definitions-loaded `evaluatedModel`, and from an
   **unablated temporary copy** of the engine, and print their comparison. ⚠ There is no committed single-case
   `SLAB_OPERATOR` tag to compare against (emission happens only inside the guarded `run()` / emit `Do`), so
   compare the **object**, not a tag run. ⛔ Not full stdout (exclude progress/RSS/timing/banners).
5. **Construction site ≠ observation site:** every patch is on the object's *construction*; ⛔ never on the carrier
   extractor / emit / diff.

## The knives — one SITE + one FORM each, both engines
⛔ Exactly ONE function and ONE FORM per knife — ⛔ no "or", no "i.e."-glued alternative, no builder discretion. ⛔
State no expected motion anywhere.

### K_A — closure face-response (Λ_A) knife
- **SymPy.** The response fold builds `closure_residuals` from `closure_shape_deriv` (`:2815-2828`), sums them to
  `closure_residual_sum` (`:2829`), and folds them into the source at `:2834`. ⚠ Each face's `closure_shape_deriv`
  is an **unexpanded `Mul`** in which `Lambda_A_0` and `Lambda_V_0` share the top-level term, so it must be
  **expanded before the term filter** or the filter deletes the whole expression. **FORM = for each face, expand
  and delete the `Lambda_A_0`-bearing additive terms, keeping every other term:**
  `patched = sp.Add(*(t for t in sp.Add.make_args(sp.expand(closure_shape_deriv)) if not t.has(Lambda_A_0)))` — ⛔
  NOT a coefficient substitution `Lambda_A_0 → 0`, ⛔ NOT the `:408` registration bind. Name `build_operator` as
  the single function.
- **WL.** `flux = lambdaAResponse affinity + lambdaVResponse normalVelocity` (`:1080`). **FORM = drop the
  `lambdaAResponse affinity` term** (leaving `flux = lambdaVResponse normalVelocity`). ⛔ Do not patch
  `lambdaVResponse normalVelocity`, ⛔ not the shared `affinity` (`:1079`, also feeds traction), ⛔ not the
  tower-only flux at `:2862`. Name `faceSources` as the single function.

### K_T — traction / virtual-work face knife (whole channel)
- **SymPy.** `face_generalized_force_rows` (`:2135-2178`) extracts `u_face`/`e_face`, added inside `build_operator`
  at `:2998-3040` (`+ face_u[a]` `:3005`/`:3021`, `+ face_e` `:3031`/`:3039`). **FORM = remove all four
  `face_u`/`face_e` additions** at `:2998-3040`. Name `build_operator` as the single function.
- **WL.** `tractionPressure = pressureField[sign] + lambdaXResponse affinity` (`:1081`) →
  `virtualWork = -graphMeasure … tractionPressure … virtualNormalDisplacement` (`:1082-1083`). **FORM = set the
  whole `virtualWork` product to `0`** in `faceSources` (`:1082-1083`) — ⛔ not only the bare `pressureField[sign]`
  term. Name `faceSources` as the single function.

### K_W — native-pressure face **collapse** (the +/− face identification)
- **SymPy.** **FORM = extend `substrate_substitutions` (`:1995`) to add a minus→plus face collapse to its returned
  substitution dict**: `delta_p_minus → delta_p_plus` and `d_w_delta_p_minus → d_w_delta_p_plus` (same-dimension).
  That dict is consumed only by `filtered_substrate`, which applies it (`:2045-2048`,
  `expression.subs(..., simultaneous=True)`) to every `SHAPE_SUBSTRATE_KEYS` substrate — including the delta_p-bearing
  `closure_shape_deriv` and `virtual_work_shape_deriv` — before they become rows. ⛔ Not the registration `:484-493`
  (declares symbols only). ⛔ Not a value-slot ↔ ∂_w-slot map (dimensionally invalid: `delta_p` is DIM_PRESSURE
  `:486`, `d_w_delta_p` is DIM_PRESSURE/DIM_L `:492`).
- **WL.** **FORM = redefine the lower face from the upper**: `pressureField[-1] := pressureUpper[…]` (replace
  `:1015`), a construction-site face collapse propagating through `affinity`/`tractionPressure`.

**No builder-chosen knives.** ⛔ Implement exactly the three knives above. If a listed site cannot be patched as
specified — the named function/line does not carry what is described — **stop and report to the orchestrator**; ⛔
do not substitute a different site, a different FORM, or a coefficient knife. Any change to this list is the
orchestrator's, never the builder's.

## Self-tests the harness must RUN and PRINT (⛔ do not state an expected result — print and stop)
For each engine, in addition to the knives:
1. **Extractor-order self-test.** Re-run the whole carrier extraction with the `→0` (P=0) applied **before** the
   `∂/∂p` instead of after; print the resulting carrier triples.
2. **Dead-path self-test (one per engine, at a pressure-free site).** Structurally delete a pressure-free kinetic
   addend — SymPy: the e_W kinetic addend at `:2970`; WL: the analogous pressure-free kinetic addend at `:1347` —
   and print the carrier diff.
3. **Live-rescale contrast (×2 at each knife's own site).** *Rescale* (multiply by 2 — ⛔ do not remove): K_A —
   ×2 the same expanded `Lambda_A_0`-bearing addends the K_A FORM deletes (keeping the others); K_T — the
   `face_u`/`face_e` additions (SymPy) / `virtualWork` (WL); K_W — the minus pressure slot
   (`delta_p_minus → 2·delta_p_minus`, ⛔ without collapse) / `pressureLower → 2·pressureLower` (WL). Print the
   carrier diff for each.
Print all transcripts. ⛔ No assertion, no verdict — interpretation is the orchestrator's.

## Deliverables
- `research/pde_ledger_v3/scripts/S11c_b_carrier_ablation_harness_sympy.py`.
- `research/pde_ledger_v3/mathematica/S11c_b_carrier_ablation_harness.wl` (with a small Python/driver if needed).
- Each opens with a **manifest** naming, per knife, the construction site it patches and the FORM applied — ⛔ **no
  expected motion / must-MOVE / must-NOT / target-row wording** (that is the orchestrator's, kept out of this
  packet). The manifest records the list; it does not originate it.
- Each run's exact invocation + source/output digests + the literal transcript go in
  `research/pde_ledger_v3/_measurements/S11c_b_carrier_ablation_harness_{sympy,wl}.md`.

## Three script clauses (verbatim, non-negotiable)
1. The harness may **PRINT** computed objects (carriers, diffs). It may ⛔ NOT state conclusions — no `PASS`, no
   verdict, no "bites"/"fails". Interpretation is the orchestrator's, from the printed triples.
2. **Print operand and residual, then guard.** Emit `baseline`, `corrupted`, and their `diff`; a residual asserted
   zero/nonzero carries no information.
3. Interpretation belongs to the review / step record, ⛔ not the script.

## Build discipline — YOUR TASK ENDS AT: write the two harnesses, run each once, report
- Write the two harness scripts (implementing the three knives + the three self-tests above) and run **each once**
  single-case to emit its transcript, then **report** what you produced. ⛔⛔ **That is the whole task.** ⛔ Do NOT
  review the harnesses, ⛔ do NOT launch or spawn any review agent / second engine / "independent" leg, ⛔ do NOT
  iterate "review-until-clear," ⛔ do NOT ablate-the-harness as a review pass. The harness's own three self-tests
  ARE the deliverable — implement and run them once; their **interpretation** and any independent review happen
  **downstream, not here.** ⛔ Do NOT `git commit`.
- Each engine's harness may be a separate run; ⛔ neither receives the sibling engine's output.
- ⛔ Never run two memory-heavy CAS jobs concurrently (30 GB box; serialize). ⛔ Wrap every WL kernel run in
  `timeout 600`. ⛔ Never run the emit `Do` / `COUPLING_KERNEL` / 4-case path (≥64 GB).
