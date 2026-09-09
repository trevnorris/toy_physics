I'll review the collapse-instrument directive as an independent decision leg: first the engines and framing, then each admitted map, retained-order design, and controls.I'll read the framing, engines, and cited line ranges next so each map can be derived from the actual code rather than the directive's prose.The c1 path missed; I'll locate it and start reading the engine sources at the cited lines.I'll pull the cited engine sections next—grades, jets, energy basis, source maps, and comparator extraction.Map 4 is the load-bearing physics map; I'll derive the energy-basis change from both engines next.I'll enumerate first-jet contractions from the engine code and check whether the table's index pairing is grounded.I'll fact-lookup which first-jet coefficients the engines actually emit, and whether WL's retained `energyCoefficient` order is on disk.# Independent decision-leg review

**Artifact:** `research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_collapse_build_directive.md` (v3, Codex-authored)  
**Role:** decision list + physics-bearing directive for a collapse instrument that does not yet exist. Reviewed by reading the engines first; no script ablation.

---

## 1. Frozen dictionary — primitive vs whole-object, grounding

Maps 1, 2, 3, 5, 7, 8 are primitive and grounded. Map 6 is correctly absent from production. Map 4 is a legitimate coefficient-basis change **as an ungraded constructor identity**, and the 19-row pairing is grounded in both engines — but several rows are the **wrong object to apply to already-graded numeric leaves**. That is the only MUST.

### Map 1 — grad-θ jet duals — primitive, grounded

Atom map, not a field/carrier equality. SymPy emits `(a.grad_theta[i], b.grad_theta[i])` at `S11c_c2_N6_reconcile_sympy.py:59-69`; operator atoms at `S11c_b_brane_operator_sympy_audit.py:227-233`; `wave_jet` rewrites `grad_theta_<i>` → `theta_d<i>` at `S11c_c2_selfenergy_fold_sympy_audit.py:139-169`. WL constructs `jet["theta",{i}]` at `.wl:84-89` and the same physical jet in `JET_VOCABULARY_BRIDGE` at `.wl:945-948`.

### Map 2 — source-wave point-preserving map — primitive, grounded, scoped

Bare `jet[f,I]` ↔ applied `∂_I s11cc2Field_f(Y,t)` at the same point. WL `sourceMap` keeps registered bare wave atoms (`.wl:707-711`); SymPy `source_value` applies `wave_jet(..., Y)` (`S11c_c2_N6_diagnostic_sympy.py:392-396`); comparator preserves applied heads (`S11c_c2_N6_cross_engine_comparator.py:287-289,990-993`). Scope `S_A,S_P,S_B` only, Φ excluded: Φ is an abstract map with no `Y` evaluation (`S11c_c2_N6_covariance_sympy.py:63-129`; `.wl:236-243`). Not source equality.

### Map 3 — profile jets / constructor scales — primitive, grounded

Spelling `w1ProfileJet<I> → w1_profile_d<I>` (and m1) plus the r≥1 scales. WL `profileRules` at `.wl:127-132`; SymPy `background_jet_expression` at `S11c_b_brane_operator_sympy_audit.py:767-778` and `profile_definitions` at `:900-908`. Numeric path is leftover spelling only; zero-jet `WBg→W0(1+η·w1)` is Φ-only, before coefficient extraction. That split matches both engines’ `finish`/ `grades` (`.wl:136-141`; `diagnostic_sympy.py:245-248`).

### Map 4 — energy-basis change — primitive as a constructor table; MUST on stage of `R_W`

This is **not** whole-μ / whole-carrier equality. It is a 19-row coefficient identification by contraction, with μ explicitly not a pairing oracle (`diagnostic_sympy.py:318-349` is correctly barred).

**Pairing is grounded**, independently:

- WL retained list is the named object `THETA_DEPENDENT_CONTRACTIONS` in `WL_S11CC2_N6RC_PROVENANCE` (retrieved from `mathematica/out/S11c_c2_N6_mathematica_audit.out`). Order is exactly the table’s E01–E19, including the gW shear swap E18/E19 vs E10/E11. `energyCoefficient<i>` is `retained[[i]]` for i≥3 (`.wl:222-226`).
- SymPy first-jet indices 04, 06, 07, 08, 12, 13, 14 match an independent enumeration of `enumerate_new_candidates` (`brane_operator_sympy_audit.py:1510-1527,1334-1353,1313-1323`). The same seven indices per source are the `gamma_s11cb_{w_bg,mu_r_bg}_{04,06,07,08,12,13,14}` atoms in `scripts/S11c_c1_exports.py`.
- Uniform scales E01/E02/E04/E07/E09 match `uniform_coefficient` (`:1584-1594`). Local thickness `E=W_0 e_W/W_bg` and `D_i E = R_W(r_i − e_W gW_i/W_bg)` are `:1631-1641`. Two-background-jet remainders of that chain rule sit outside `σ_W^1`; the one retained mixing is E14’s `q·D(E)` contribution to `eW·(gW·q)`.

**MUST — `R_W` images are ungraded constructor identities, but production applies them to already-graded numeric leaves.**

The table lives in `Q(W_0,W_bg)` and the production direction is “WL coefficient atom → last column” on “coefficient atoms in materialized numeric operands” (scope table, map 4). Numeric leaves are already grade-keyed (`§ Retained order`; WL `gradePart` after `finish` at `.wl:136-141,713-714`; SymPy `Compiler.grades` after `xreplace(profiles)` at `diagnostic_sympy.py:245-248`). The directive also forbids expanding `R_W` in `η` inside a graded leaf.

Those two requirements do not commute for every `R_W` row (E02, E03, E04, E05, E13, E14).

Worked example, E02. Ungraded constructors:

- WL: `cCoupling · WBg · θ · eW` (`.wl:222`)
- SymPy: `C · W_bg · θ · E` with `E = W_0 e_W/W_bg` (`:1589-1590,1631`), which **cancels** to `C · W_0 · θ · eW` **before** grading.

Ungraded identification `cCoupling → R_W·C = (W_0/W_bg)·C` is correct. After each engine grades:

- WL `η^0`: `cCoupling · W_0 · θ · eW` (`WBg` already expanded by `profileRules`)
- SymPy `η^0`: `C · W_0 · θ · eW` (no leftover `W_bg`)

Substituting the frozen image into the graded WL leaf yields `(W_0/W_bg)·C·W_0·θ·eW`, which equals the SymPy leaf only if `W_bg` is frozen to `W_0`. Leaving `W_bg` live reintroduces a zero-jet field into a grade. Expanding `R_W` in `η` mixes grades, which the tripwire forbids.

The same stage error hits every `R_W` row. Rows without `R_W` (E01, E06–E12, E15–E19) **do** commute with grading (E01’s `bRho → B_rho_3/W_0` is a constant `W_0`, and both engines still carry the `WBg` factor in that contraction). Rank/invertibility of the 19×19 matrix over `Q(W_0,W_bg)` does not catch this: that check is on the ungraded table, not on the object applied to graded leaves.

This changes every residual that uses E02/E03/E04/E05/E13/E14. It is not a collapse prediction; it is the dictionary applying the wrong stage of a correct pairing.

### Map 5 — Φ spelling — primitive, grounded

Domain/multi-index/derivative syntax only. SymPy `prolonged_phi` (`covariance_sympy.py:63-129`); WL `phiMap` (`.wl:84-94,236-243`). Images stay operands. No whole-Φ equality, no `Y` evaluation.

### Map 6 — correctly removed

Neither engine emits a leftover density **atom**. SymPy substitutes `inputs.density[(rho,)][1]` inside `source_terms` (`diagnostic_sympy.py:378-383`); WL builds `density3=density4·WBg` and `sourceBind` does `rhoFace→density` (`.wl:380-381,839-840,869-872`). A production rewrite would be either inert or a constant↔live-field fold (the c1 rule-17 hazard). Census-only `LIVE_DENSITY_PREMISE` is the right remainder.

### Map 7 — leftover Jacobian spelling — primitive, grounded, occurrence-gated

`1+tr(∇u)` ↔ `1+Sum_i jet["u"<>i,{i}]` only if that factor still sits in a scoped operand. Emitted at `reconcile_sympy.py:59-68`; used inside `materialAmplitude` at `.wl:244-247`. Degree-2 projection is a construction stage (`brane_operator_sympy_audit.py:1970-1981`; `.wl:244-247`) and is excluded. Scope `C_M,S_A,S_B` (not `C_E,S_P,PHI`) matches where the factor is built.

### Map 8 — ε / ω occurrence-conditioned — primitive, not a default on-shell map

ε path cited at `diagnostic_sympy.py:378-389` and `.wl:364-385`; `omega` declared at `.wl:102-113`. Activate only if the atom survives maps 1–7. No Fourier-of-derivative identity.

### Excluded list

Whole-carrier/μ/source/Φ equality, construction-stage projectors, `LAB_HELD↔MATERIAL_ADVECTED`, `σ_W` binding / `σ_W→0`, default on-shell/Fourier: all correct. Face velocity and source-solve factors stay in the residual (`.wl:274,304,366-381,862-872`; `reconcile_sympy.py:82-92`; `diagnostic_sympy.py:378-389`). Those are independently computed subobjects, not primitives.

---

## 2. Completeness vs over-reach

Needed primitives that the surfaced operands actually use are present, except that map 4’s `R_W` rows are the wrong **stage** of a needed scale (finding above). No extra whole-object map. Comparator mechanical folds (pressure-slot spelling, jet camelCase) stay available and do not enlarge physics scope (`comparator.py:287-326`). Inactive maps are reachability records, not a licence to broaden.

---

## 3. Retained order / no proxy

Forced coefficientwise on `(η^i σ_W^j), i,j∈{0,1}`. SymPy `GRADES = product((0,1), repeat=2)` at `diagnostic_sympy.py:43`, extracted at `:245-248`. WL `gradeIndices = Tuples[{Range[0,1],Range[0,1]}]` and nested `Series` at `.wl:123-141`. Comparator keys `ETA,SIGMA ∈ {0,1}` independently (`comparator.py:96-97`). Grade combining / `σ_W→0` / `σ_W↔η` / `η` re-introduction into a graded leaf are forbidden and tripwired. Φ is the one family that is ungraded in the comparator (`extract_meta` at `:688-765`; `FROZEN_PHI` at `covariance_sympy.py:116-124`; `.wl:975-976`) and is explicitly re-graded by profile expansion then four-coefficient extraction. That is the retained rectangle, not a proxy.

---

## 4. Collapse witness

Three-valued `{collapsed, residual_remains, undecided}`. Residual printed first. Fresh PIT, primes/seeds disjoint from both engines (`diagnostic_sympy.py:669` → `(1000000009, 998244353, 1000000033)`; `.wl:733` → `{1000000009, 998244353, 1004535809}`). Denominator rejection, derived bound, attempt/time/RSS caps, exhaustion → `undecided`, provenance emitted. No verdict path; native Booleans not coerced to a CAS zero. Controls print deltas and do not assert a production residual.

---

## 5. Mandatory controls

By design they can fail:

| Control | Can bite? |
|---|---|
| One-sided corruption + locality | Yes: sentinel image, `control_delta` nonzero on the predeclared use-set and zero outside. Catches an over-broad rule. |
| Dropped primitive | Yes: delete one active primitive, require a residual on the **defining-relation** pair, not a production leaf. |
| Reachability | Honest inert census, not a bite test (correct). |
| Blanket-collapse | Yes: whole-operand identification on the same pipeline must yield `collapsed` on every ablation leaf; unreachable from production. |
| Grade-combining tripwire | Yes: structural reject of zero-jet re-expansion / `η` insertion / `σ_W` binding / cross-grade access. |

They are not tautological **except** that they will not catch the map-4 stage error: the defining-relation for map 4 **is** the ungraded table row, so dropped-primitive and rank checks can pass while production still applies `R_W` to graded leaves.

---

## 6. BASELINE + census

`SOURCE_BASELINE` is the nominal-control duplicate at `kappa_a=1, kappa_j=0` / `actualAdvectionCoefficient=1, actualJunkCoefficient=0` (`covariance_sympy.py:34-35,137-154`; `.wl:18-20,845-849`). Not a third discriminator.

Census is load-bearing, not metadata. Closed six-row crosswalk (`kappa_a↔ADVECTION`, `kappa_j↔JUNK`, `imported_theta_e_jets↔DOMAIN`, `coverage↔COVERAGE`, `uncovered↔UNCOVERED`, `max_present_jet_rank↔MAX_RANK`) matches the engines (`covariance_sympy.py:87-115,137-153`; `.wl:236-243,850-852,977-980`). WL-only `THICKNESS` / `MATERIAL_NORMAL` stay one-sided; no synthesized SymPy zero. `JUNK_SYMBOL` may still join by exact `field_name` lowercasing — that is allowed, not pre-adjudicated. No expected census residual.

---

## 7. Value-free / leak

No expected collapse outcome, witness count, or production pass condition. Prior-run 40/76/18 counts are not reproduced. Map-4 `W_0/W_bg` and `1/2` are constructor facts. “The expected table is exactly the one printed above” specifies the dictionary, not a collapse target.

---

## 8. Fence + DoD + extraction reuse

Astra fence is explicit: build → verify → run → report → STOP; no self-review; no comparator edit; no edits outside deliverable/test/`.out`. Reusing `load_py_jsonl`, `load_wl`, `verified_spelling_maps`/`ACTIVE_NAMES`/`POINT_NAMES`, leaf decoders, `extract_*`, `extract_meta`, `materialize` is correct: `compare_family` materializes then **releases** (`comparator.py:794-851`) and `object_work` returns accounting only (`:862-898`). DoD can detect a missing deliverable, a reachable whole-object map, missing map-6 absence, typed-scope violations, and a verdict/`assert` path. It cannot detect the map-4 graded/`R_W` stage error, because it requires “map 4 is the exact 19-row table above.”

---

## Finding (MUST)

**Map 4 `R_W` rows (E02, E03, E04, E05, E13, E14) are ungraded constructor identifications specified as production rewrites on already-graded numeric leaves, with expansion of `R_W` in `η` forbidden.** That does not commute with both engines’ grade extraction (`diagnostic_sympy.py:245-248`; `.wl:127-141,713-714`). Production residuals on those contractions would carry a `W_bg` vs already-cancelled/`W_0` mismatch that is an artifact of the dictionary stage, not the §3 question. Repair is to freeze either (i) apply-map-4-then-grade on ungraded coefficient atoms, or (ii) an explicit per-grade image of those six rows with no live `W_bg` reintroduced into a grade. Rank/invertibility of the ungraded 19×19 matrix is not that freeze.

All other admitted maps are primitive, cited, and not §3 in disguise. Witness, retained rectangle, exclusions, BASELINE, census, leak discipline, fence, and extraction reuse are sound.

---

**DIRECTIVE NOT-SOUND** — MUST: map 4 `R_W` images applied to already-graded numeric leaves.
