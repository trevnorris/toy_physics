# Compact-prep verify — are the c2 N6 state docs accurate + in the clear?

You are Codex, verifying the session's state records are accurate before a compaction. ⛔ Document-only; ⛔ do NOT
modify the tree; run no CAS. Read the sources + `git log`, confirm each claim, flag any over-claim/error with
`file:line`. Working dir `/var/projects/toy_physics`; paths under `research/pde_ledger_v3/` unless absolute.
Relevant commits (`git log --oneline -14`): `a7770aff` STATUS, `4eb1a9ac` N6 BUILD CLEAR, `a24acfb6` N6 build,
`595a5024` CLEAR TO BUILD, `c24a9ade` astra spec folded, `be542cb2` r3 change-of-author, `ead5588a` r2, `d0666953` r1,
`77d812a0` compact-prep (prior session).

## The claimed arc (verify each against the commits + records)
1. **N6 route-2 took 4 directive rounds + a rule-15 CHANGE OF AUTHOR.** My hand-written route-2 was wrong 3× (r1
   ordering `d0666953`; r2 `material_pullback`-on-rows `ead5588a`; r3 legs split + the `T`-on-increment framing is the
   wrong TYPE `be542cb2`) ⇒ **astra authored the route-2 spec** `_measurements/S11c_c2_N6_route2_spec_astra.md`,
   reviewed **fresh Claude + Grok BOTH SPEC CLEAR** → folded into §5c + directive `c24a9ade` → final Codex+Grok pass
   CLEAR TO BUILD `595a5024`. Verify the adjudication records `_measurements/S11c_c2_N6_directive_review_r{1,2,3}_adjudication.md`.
2. **The cleared construction:** `I_{M→E}=extract(close(SLAB_M)−SLAB_M)`, **NO `T` on the increment** (native S11c-a
   MATERIAL face sources into the same δp symbols, differenced directly — the parent `task_rep_invariance` pattern);
   **material μ_θ** at BOTH face + c1-source (route 1 Eulerian μ_θ); **reverse-u GRADE-SUPPRESSED** (the mandatory-
   survival requirement REMOVED). Verify §5c (`directives/S11c_c2_SHARED_PHYSICS.md` ≈303-393) + the directive
   `directives/S11c_c2_N6_sympy_build_directive.md` faithfully carry this, and both point at the cleared astra spec.
   ⛔ Flag any residual stale `T_{M→E}`-on-increment / "reverse channel must survive" language in either.
3. **The BUILD:** `scripts/S11c_c2_N6_diagnostic_sympy.py` (918 lines) + the audit `ANCHORING_L_MINUS_M` edit committed
   `a24acfb6`; **2 build legs (fresh Claude + Grok) BOTH BUILD CLEAR** `4eb1a9ac`. Verify
   `_measurements/S11c_c2_N6_build_clearance.md` (FORM ablation + one-sided route corruption + AST 0-asserts;
   guards=0; PIT sound). Confirm the audit edit in `scripts/S11c_c2_selfenergy_fold_sympy_audit.py` is raw
   `LAB−MATERIAL` with no `representation_pullback`, density loop kept
   (`git diff a24acfb6~1 a24acfb6 -- …selfenergy_fold_sympy_audit.py`).
4. **The DISPOSITION** (the one to scrutinize hardest): `R_N6` is **certified-nonzero in ~18 forward-block columns** and
   recorded as a **computed FINDING, NOT pre-adjudicated** — ⛔ NOT stated as "representation invariance FAILS" and
   ⛔ NOT as "N6 passes". Its disposition (real failure vs reconcilable remainder) is the reconcile/step-record stage.
   Verify the build-clearance record + STATUS + memory all frame it this way (raw nonzero ≠ disagreement).
5. **STATUS.md** new top clause `a7770aff` + the memory
   `/home/trevnorris/.claude/projects/-var-projects-toy-physics/memory/project_s11c_c_state.md` (body "NEXT" section,
   ≈ lines 215-235) + `MEMORY.md` pointer line for `project_s11c_c_state.md` — accurate, consistent with 1–4, prior
   STATUS clause marked [SUPERSEDED], no stale "NEXT = N6 BUILD correction".

## Flag specifically
- Any **over-claim**: N6 stated as "DONE"/"passed" (it is per-engine SymPy built+cleared; the **blind Wolfram N6 +
  cross-engine comparator + reconcile + step record are still owed**); the `R_N6` nonzero stated as an invariance
  failure or a pass (it is an un-adjudicated finding); the reverse-u removal or the μ_θ-material decision misstated;
  the audit edit misdescribed.
- Any **factual error** vs the commits/records (SHAs, what cleared, what is owed).
- Any **owed item omitted** anywhere: the reconcile adjudication of R_N6; the blind Wolfram N6 engine; numeric F/G
  re-grounding; `MEMORY.md` compaction (~21 KB → <17 KB); the ephemeral 499 MB `.out`.

## Output
End with **CLEAR TO COMPACT** (docs accurate + open items correctly recorded) or the exact fix list. Keep it brief.
