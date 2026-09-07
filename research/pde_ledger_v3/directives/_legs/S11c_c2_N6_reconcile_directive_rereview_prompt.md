# Re-review (round 2) — the S11c-c2 N6 RECONCILE build directive, after folding round-1 findings

## Artifact + what changed
`/var/projects/toy_physics/research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_directive.md` (v2) — an
orchestrator-written directive for an astra-authored reconcile diagnostic. Round 1 (2 legs) reported FOLD REQUIRED;
the five findings were folded. Working dir `/var/projects/toy_physics`; paths under `research/pde_ledger_v3/`.

This is a **physics-bearing pre-builder directive re-review** (review-until-clear). ⛔ No CAS, ⛔ no build, ⛔ no
fictional-script ablation (that is the build legs' job). Your task: **independently re-verify each fold against the
sources** (⛔ not merely that the fold text is present), and **scan for any NEW defect or regression the folds
introduced**.

## The five folds to re-verify (the adjudication record is `_measurements/S11c_c2_N6_reconcile_directive_review_adjudication.md`)
Re-derive each from the sources; confirm it landed correctly AND is faithful:
1. **A1 — separated controls.** `a_ρ` lives only in the constitutive material μ (`constitutive` :328-342); open
   pressure-slot coefficients carry neither θ nor μ (route-2 §4 :130-150). The directive's "Controls" must now bite
   the **carrier** via a covector/normal-map corruption (`ms` fixed) and the **source** via `a_ρ→0` (RHOBR-only,
   carrier + Eulerian byte-identical), OBSERVE (not require) `R_N6`, and emit computed absence for RHO4. Confirm no
   residual "must move the carrier under `a_ρ`" language and no A−A remains.
2. **A2 — three-way split.** `R_N6 = I(ΔC,S_M) + B(C_M,ΔS) + B(ΔC,ΔS)` (carrier/source/cross). Verify the algebra
   (`I(C,S)=−C·p+B(C,S)`, `build_increment` affine in S :491-506), that CARRIER uses the **material** source `ms`
   (not `S_E`), and that SOURCE/CROSS are signatures-{6,9,12}-only contractions with no re-added bare term.
3. **A3 — source pinned to `es`/`ms` + in-process single `pit()`.** Confirm `SOURCE_{EULERIAN,MATERIAL}` are the
   diagnostic `es`/`ms` circuits (`source_terms` :378-389, :841-843), the `b_{r,s}` formula is a slot identification
   only (⛔ not a re-coded expression; `V̄=V/ε` already inside), and that `E,M,R_N6` are recomputed in-process under
   one `pit()` (⛔ no cross-run join of stored PIT tables).
4. **A4 — SPLIT_CHECK contract.** Confirm the guard is **samplewise-zero numerators on shared PIT samples** (⛔ not a
   structural zero node — `plus`/`minus` :137,:163 don't reduce `x−x`; ⛔ not satisfiable by emitting `number(0)`).
5. **A5 — jet-vocabulary bridge.** Confirm `a.grad_theta[i] (theta_d1) ↔ b.grad_theta[i] (grad_theta_1)` (route-2
   :128; diagnostic `N6_JET_BRIDGE` :351-356) is now in the frozen + emitted relation set, with broader renaming
   forbidden.

## Also re-confirm the round-1 "held sound" survived the folds
The corrected zero-modulo-definitions question; no nonzero offset `J`; fixed-anchoring restriction; coefficient-level
carrier test; combined-source formula kept outside opaque `Z[b]`; one-sided PIT interpretation (nonzero = certificate,
all-zero = "no nonzero found" at conditional δ); no expected value / no residual-zero exit; the three script clauses +
four corollaries.

## Sources (read first, form your own view)
`scripts/S11c_c2_N6_diagnostic_sympy.py` (:318-360, :378-389, :476-513, :841-853, :668-773);
`_measurements/S11c_c2_N6_route2_spec_astra.md` (:39-83, :99-150, :152-178, :128, :209); §5c of
`directives/S11c_c2_SHARED_PHYSICS.md`.

## Physics filter
Report a finding only if it catches a way the reconcile could be wrong, vacuous, leak-prone, mislabel a real failure
as invariance (or vice versa), or intractable — or a NEW defect/regression the folds introduced.

## Output
Findings each with `file:line` + minimal fix. End with **DIRECTIVE SOUND — CLEAR TO BUILD** or the exact remaining
fold list. Brief, evidence-first.
