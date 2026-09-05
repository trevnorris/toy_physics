# S11c-c2 SymPy build directive — decision-list gate record (rule 7 / G2: one two-leg pass, fold once, go)

The c2 SymPy build directive (`directives/S11c_c2_sympy_build_directive.md`) is a **pre-builder decision list** —
orchestrator-written, thin (physics = pointers to the cleared spec `16849fc6`), fixing only the build-mechanical
layer. Per G2 it gets **one** two-leg pass (Codex sol xhigh + Grok), verify, **fold once**, go to the build — ⛔ not
iterated to green. Findings that change what is computed route to the **build** gate (the astra build's own two legs),
⛔ not to another decision-review round.

## Commands (identical prompt to both legs; prompt `_legs/S11c_c2_build_directive_decision_review_prompt.md`)
```
codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
  "$(</abs/_legs/S11c_c2_build_directive_decision_review_prompt.md)" > …_decision_review_codex_sol.txt 2>&1
grok --prompt-file /abs/…_decision_review_prompt.md --cwd .../pde_ledger_v3 --model grok-4.6 --effort high \
  --permission-mode bypassPermissions --output-format plain > …_decision_review_grok.txt 2>&1
```
Reports: `S11c_c2_build_directive_decision_review_{codex_sol,grok}.txt` (Codex clean-extract; raw 2.1 MB transcript
outside the repo at `/var/projects/toy_physics_ext_logs/…_codex_sol_RAW.txt`, tree hygiene).

## What both legs CONFIRMED sound (computation-backed, ran `ledger_fold`)
Import wiring §2: fold 2441+44→2485 (empty intersection, no overwrites); the **18 closure-covering object roots** are
minimal and their `check_consumer` closure (355 keys, no ambiguity) covers the whole 38-key consume-set; the
`IMPORT_KEYS` = exact-`fold[key]`-lookup rule is right (undeclared **and** declared-but-unused both fail);
`v_bulk_normal_0` must not be declared (in fold, not in closure/consume-set); the bare-vs-prefixed `face_response`
hazard is real (both exist, values differ, guard passes on either — so the guard will NOT catch a swap). Faithfulness
to the cleared spec: no contradiction/weakening/omission (both legs). Leak: clean (only meta-prohibitions). The four
δp_s slots are bare COORDINATE Symbols in `slab_operator`/`closure_shape_deriv`/`traction`; `Λ_X` in `traction` only.

## Findings (verdict NOT READY) — verified by me (rule 13) against the real fold, then FOLDED once
1. **Closed δp_s source = `DELTA_P`, ⛔ NOT `PRESSURE`** (both legs). The operator-valued closed pressure is
   `s11c_c1_face_response → CASES → (α,±1,ρ) → VALUE → DELTA_P`; the `PRESSURE` object lives on
   `s11c_c1_face_response_coeffs` and is the **scalar flat/coefficient** regression object (the "scalar division" the
   directive forbids). Verified: face_response display has DELTA_P not PRESSURE; coeffs has PRESSURE not DELTA_P.
   → §3.1.
2. **w-jets are COMPUTED, not imported** (both legs). c1 exports no `d_w_delta_p_*`; compute them as the shape-deriv
   of `DELTA_P` (S11c-a interface geometry). Verified: `C1_RESPONSE_EXPLICIT_D_W_DELTA_P_REFS=[]`. → §3.2.
3. **ε-normalization** (Codex). c1's `DELTA_P` is ε-scaled (inputs = ε·symbol, c1 `.py:834`) and the slab slot carrier
   is multiplied by an outer ε ⇒ direct subs double-counts to **O(ε²)** vs the correct O(εη). Verified: face_response
   carries `epsilon_shape`; the coeff·DELTA_P product factors `epsilon_shape**2`. → §3.3 (reconcile to one ε; EMIT
   the (ε,η,σ_W) order as the check).
4. **`V_s` → `face_velocity`, ⛔ NOT `kinematic_balance`** (both legs); pin the `DELTA_W`/`ZETA_C` representation.
   Verified: `face_velocity` is the 8-case V_s Tuple; `kinematic_balance` is a `delta_v_bulk_*` residual identity. → §3.4.
5. **Kernel bridge** (Codex). A 4th c1 identifier class `s11cc1_dtn_operator_*` (whole-form `Z`) must be bridged to
   the AGREE'd two-momentum `dtn_kernel` per case, else the whole-form silently becomes the construction operand
   (spec §3d.3). Verified: face_response references `s11cc1_dtn_operator_*` + resolvent. → §3.4 table row + §3d.3.
6. **`coupling_kernel` = mandatory 19th regression import** (Codex). Not in the 18 nor closure-pulled; the design
   PROMOTES a declared regression operand into `IMPORT_KEYS` (`export_ledger_bind_closure_design.md:148`). Verified:
   `coupling_kernel` in fold, not in 18-closure. → §2.
7. **Three script clauses byte-verbatim** (both legs). → §0 restored to the build-skill text verbatim.

## Disposition
All seven folded into the directive in ONE pass (§0, §2, §3.1–3.4). ⛔ No decision-review re-leg (G2 one-pass). The
directive is the reviewed baseline committed before the astra build runs against it; the build's own two re-review
legs (Codex-written → fresh Claude agent + Grok) re-verify the built script + the corrected map + the ε-order check +
the kernel bridge against the real files.
