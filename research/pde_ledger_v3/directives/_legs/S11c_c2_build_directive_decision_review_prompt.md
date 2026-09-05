# Decision-list review — the S11c-c2 SymPy build directive (pre-builder, two legs)

## Artifact
`directives/S11c_c2_sympy_build_directive.md` — an orchestrator-written **pre-builder decision list** for the
SymPy-engine build of step S11c-c2 (the self-energy fold). Its physics authority is the **already-cleared** spec
`directives/S11c_c2_SHARED_PHYSICS.md` (v2, committed `16849fc6`; 2 decision legs + fold — ⛔ do NOT re-review the
spec's physics decisions, they are settled). This directive is **thin**: it POINTS at the spec for physics and fixes
only the **build-mechanical** layer. Your job is to check that layer + that the directive carries the spec faithfully,
in **ONE pass** (this is a decision list, not a spec — fold-and-go, ⛔ not iterate-to-green).

## What to check — the build-mechanical DECISIONS and their evidence

### 1. Import wiring (§2) — VERIFY BY RUNNING `ledger_fold`, ⛔ not by reading
Run `load_model("scripts/S11c_b_exports.py","scripts/S11c_c1_exports.py")` (from
`/var/projects/toy_physics/research/pde_ledger_v3`; `scripts/ledger_fold.py:102`). Confirm or refute, each with your
command + literal stdout saved to a named `/tmp` path:
- The fold: base 2441 + c1 delta 44 → 2485, empty exact-key intersection, `overwrites==[]`.
- The **18 closure-covering object roots** the directive lists (13 from S11c-b + 5 from c1): does
  `check_consumer(fold, <the 18>)` resolve, and does its closure **cover the entire §7 consume-set** (the 18 + the 8
  `s11cc1_response_resolvent_*` + the 6 `s11cc1_k_*` + `s11cc1_q_out_{input,output}` + `s11cc1_w1_profile_hat_transfer`
  + 3 `s11cc1_w1_profile_jet_hat_*` = 38)? Is any of the 18 actually reachable from the others (i.e. is 18 truly
  minimal), or is a needed root missing?
- The **`IMPORT_KEYS` rule** the directive states: is it correct that `assert_lookups_equal_manifest` requires
  `IMPORT_KEYS` = **exactly** the build's direct `fold[key]` lookups (undeclared **and** declared-but-unused both
  fail)? Is declaring `v_bulk_normal_0` genuinely wrong here (in the fold but not referenced by c2's consume-set)?
- The **bare-vs-prefixed `face_response`** hazard: do BOTH `face_response` (step S11b) and `s11c_c1_face_response`
  (step S11c-c1) exist with **different** values, and do `check_consumer`/`assert_lookups_equal_manifest` pass on
  **either** (so the guard will NOT catch a wrong-provenance bind)? Is the directive right to make binding the
  prefixed one a hard requirement?
- Is `coupling_kernel` correctly restricted to the §5a ordering-ablation **regression** operand (⛔ never a §3a/§3c
  construction operand)?

### 2. The fold symbol map (§3) — VERIFY against the real rows
- Are the four δp_s slots **`delta_p_plus, delta_p_minus, d_w_delta_p_plus, d_w_delta_p_minus`** bare Symbols present
  in **`slab_operator`, `closure_shape_deriv`, AND `traction`**, and is `Lambda_X_0` present in `traction` but ABSENT
  from `closure_shape_deriv` (confirming Λ_X-in-traction-only)? Show it.
- Is the two-stage `subs` correct: (Stage 1) replace the 4 slots with c1's closed `δp_s`(+w-jets) per case
  `(anchoring, ±1, density)` from `s11c_c1_face_response{,_coeffs}`; (Stage 2) identify the c1 symbols with slab
  fields? Are the DISJOINT names (slots absent from `s11c_c1_face_response`) correctly handled as an explicit map
  rather than a silent identification?
- Are the **two asymmetries** correctly identified and resolved: (a) **μ_θ** — slab `mu_theta_L/M` is per-anchoring
  only, c1 `s11cc1_mu_theta_{α}_{s}` is per anchoring×face, so both faces' c1 μ_θ map to the single slab
  `mu_theta_{L/M}` (μ_θ face-independent, spec §1d); (b) **V_s** — the slab has NO bare `s11cc1_V_*` symbol (V_s is
  kinematic), so the map is a substitution to the slab's kinematic V_s (from `kinematic_balance`/`face_velocity`, spec
  §1d `n̂_s·v_bulk,s=V_s+J_s/ρ_m`), ⛔ not a rename? Is anything in the map wrong, ambiguous, or missing (e.g. a fourth
  needed identification)?

### 3. Faithfulness to the cleared spec (§1 pointers) — ⛔ do NOT re-derive the physics; check the POINTERS don't drift
Confirm the directive's HELD-PHYSICS carries match the spec `16849fc6` (⛔ report a finding only if the directive
**contradicts, weakens, or omits** a spec requirement, ⛔ not if you'd word the spec differently): close-then-extract
ordering (`extract(close(SLAB))`); substitute the closed **δp_s (+w-jets), NOT a closed J_s** (no J_s slot);
**Λ_X in traction t_s only**; the increment **SURFACES** (⛔ not cancels) the two cross-engine-UNVALIDATED sign
conventions and the **§3d.4 mechanical-power pairing adjudicates the face-force sign**; the three rule-17/UNDECIDED
carry-ins (§3d.1 density field-vs-field with the pre-fold `background_density_map` bind and `∇ρ→0` barred; §3d.2 `t_s`
native covector; §3d.3 DtN kernel-vs-whole-form). Does the directive correctly forbid treating any c1 §1b UNDECIDED
item as cross-engine-closed?

### 4. Build-skill compliance + leak discipline
Are the three script clauses present verbatim + the structural rule + the "no tautological residual / increment is an
export representation" corollary? Does the directive leak any c2 **output** value, sign, order, parity, or grade, or
any acceptance criterion referencing an expected value? (The supplied slab-row / operator-inverse **input** structures
are legitimately supplied — flag only an expected **output**.)

## ⛔ Method boundaries (decision-list review)
- ⛔ **Do NOT ablate a fictional script.** The SymPy deliverable does **not exist yet** (rule 14: a control written
  against a nonexistent script is reviewed by reading — the weakest instrument — and belongs in the **build** review,
  ⛔ not here). Defer every executable script-control test (FORM ablation, one-sided corruption of the built script,
  emit-before-guard) to the build's own two legs. Here, check only the **decisions + their evidence**.
- ✅ **DO run `ledger_fold`** to verify §1–2 wiring/map claims — those are decisions about the REAL files, and a prose
  claim about the fold is worth nothing (show the command + literal stdout). Save scripts + stdout to named `/tmp`
  paths; ⛔ do not modify the working tree.
- This is **ONE pass**: report the findings that change what the builder will compute or may claim; ⛔ do not iterate.

## Output
Ranked findings (MUST-FIX / SHOULD-FIX / NIT) — each with the quoted directive text, the quoted spec/real-file
evidence, and your command+stdout paths for §1–2. End with an explicit verdict: **READY TO BUILD** (the wiring, map,
carries, and clauses are correct and evidence-backed) or the exact list to fix before the astra build.
