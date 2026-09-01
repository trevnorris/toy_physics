# S11c-b #89a WL §3a basis-repair engine — build-review clearance (2 legs, both CLEAR on physics)

⚠ ANSWER-BEARING (states the completed count 40 / 15). For the orchestrator + review record only.

Artifact reviewed: `research/pde_ledger_v3/mathematica/S11c_b_brane_operator_mathematica_audit.wl` (Codex build,
+605/−73) + its regenerated output `…/mathematica/out/S11c_b_brane_operator_mathematica_audit.out` (156 MB).
Directive: `directives/S11c_b_89_wl_3a_repair_directive.md`. Build log `~/.s11_build/S11c_b_89_wl_build_codex.log`.

## Legs (rule 7: engine is Codex-written ⇒ legs = fresh Claude agent + Grok; Mathematica SERIALIZED)
- **Leg A — fresh Claude general-purpose agent.** VERDICT 1 finding (0 blockers). Evidence
  `~/.s11_build/S11c_b_89a_wl_buildleg_claude/` (derive_family, ablate_family, operator_freeze, span_check2.wl,
  ablation_harness.wl + .stdout).
- **Leg B — Grok (grok-4.6, high).** VERDICT 1 finding (0 blockers). Evidence
  `~/.s11_build/S11c_b_89a_wl_buildleg_grok/` (derive_complete_family.py, span_engine_vs_derivation.py,
  harness_enumeration_rank.wl + .stdout). ⚠ Grok's first launch failed on an orchestrator path typo
  (`toy_projects/..`, exit 1, no work) — relaunched with the correct path; this run is the leg of record.
- Prompt (both legs, byte-identical): `directives/_legs/S11c_b_89_wl_basis_build_review.md`.

## Physics claims and the computations behind them (rule 2 — command + literal output)
Both legs derived the complete family INDEPENDENTLY from first principles (own SymPy enumerators), then matched
the engine; two independent derivations + the engine all agree — ⛔ not a matching-integer coincidence.

- **Basis is 40 = 10 uniform + 15(∂W_bg) + 15(∂μ_R,bg), restricts to 26.** Independent enumeration: one
  spurion (vector) × a quadratic in `{u,∇u,θ,∇θ,e_W,∇e_W}`, all δ-only (O(3), Kronecker) contractions, ε
  (Levi-Civita) EXCLUDED by parity, exact thickness map imposed ⇒ **15/source, rank 15**. Committed-8 ⊊ 15
  (strict). Engine emits `ENERGY_BASIS_COUNT.COUNT_OPERAND = 40` both anchorings (`RANK_SELECTED_INDICES
  {1..40}`), `RESTRICTED_FORM_OPERAND = 26` (residual 0). (Leg A derive_family.stdout; Leg B
  derive_complete_family.stdout.)
- **⭐ SPAN identity (the decisive check — verify SPANS not COUNTS).** Both legs rebuilt the 15 forms with
  their own contraction code: `rank(engine ∪ independent) = 15 = each` per source; Grok found the engine's
  pairings instantiate to BYTE-IDENTICAL scalars; label sets equal. The engine's 15 span the SAME 15-dim space
  as the independent derivation, not merely the same integer. (Leg A span_check; Leg B span_engine_vs_derivation.)
- **Count is COMPUTED not typed.** No `40`/`15` literal in source; `COUNT_OPERAND = Length[MatrixRank-selected]`;
  the sole `26` literal is the restrict-switch `REFERENCE_OPERAND`. (Both legs grep.)
- **FORM ABLATION (mandatory) bites.** On a `/tmp` copy / minimal harness: restrict-to-8 → 40→26 (Δ14, family
  back to 8/source); +1 ε pseudoscalar → per-source 15→16 (parity exclusion load-bearing); drop one shear
  Kronecker pairing → 40→38 / 15→14 (completeness load-bearing). Byte-identical under (a)/(c) would mean the
  repair did nothing; it moved. (Both legs, ablation harness.)
- **Thickness map imposed.** New invariants use the mapped `localEw`; `THICKNESS_MAP_RESIDUAL` nonzero at grade
  {η¹,σ_W¹} for the e_W-bearing records, zero otherwise. (Both legs.)
- **Operator-freeze diagnostic HONEST (scope split verified).** `PRODUCTION_FROZEN_EL_RANK=26`,
  `LIVE_BACKGROUND_EL_RANK=40`, difference 14, NOT asserted equal; `evaluatedModel` still `applyProfile`s
  before `constrainedRows` — the operator freeze is DOCUMENTED, not repaired (deferred to #89b). Both legs'
  abstract EL tests independently confirm the Hessian is non-absorbable operator content, ⚠ by DIFFERENT
  measures (not the same total-rank number): the engine's total frozen-EL 26 vs live-EL 40; Grok's abstract
  total frozen 7 vs live 12; Leg A's abstract Hessian-SECTOR rank 2→0 (its total abstract rank was 4=4 — the
  drop is confined to the Hessian sector). All three agree the freeze removes the Hessian.
- Existing controls (`REP_INVARIANCE`, `INDEPENDENCE`, `FORM`, `UNIFORM`, `HOMOGENEITY`) still emit; per-
  anchoring ranks computed separately with span-equivalence residual 0; coefficient control count-invariant;
  Hessian-in-energy guard Δ 0.

## The one converged finding — CONTROL-QUALITY, deferred to #89b (not a physics defect)
Both legs (SHOULD-FIX): the §5.E new-invariant dimension residual `RESIDUAL_DERIVED_MINUS_STORED`
(`…_audit.wl:2269–2290`; STORED set at enumeration ~L482–484) is a **copy-check** — `deriveInvariant
DimensionFromFactorContent` and `STORED_INVARIANT_DIMENSION` are the same `dimensionGradient + Σ(factor dims)`
over the same `dofFactorSpecifications`, so the residual is structurally 0 for any factor-dimension assignment
(Grok proved: mutating both copies → residual 0; only desyncing the lookup after freeze → `{1,0,0}`). It cannot
catch a wrong invariant dimension — a check auditing its own input ([[feedback_a_check_cannot_audit_its_own_input]]).
⚠ Verified benign: the OLD hand-typed `newInvariantDimensions` list IS removed, and the derived dimensions ARE
correct (both legs checked against factor content: U_THETA→dimensionless, GRADU_GRADEW→L⁻²). ⇒ The dimensions
are validated OUT-OF-BAND (the legs), not by the in-band residual. **FOLD into #89b**: replace §5.E with a
genuinely independent dimension route (dimensional analysis of the emitted invariant EXPRESSION vs the factor-
metadata sum), or delete the tautological residual. ⛔ Do NOT claim the in-band residual validates the dims.

## Disposition
#89a WL §3a basis repair (completed O(3)-Kronecker enumeration + thickness map) is **CLEARED on physics** by
two independent legs — completed count 40, span identical to independent derivation, computed-not-typed, form
ablation load-bearing, operator-freeze honestly deferred (scope split verified). Both engines (PY #89, WL #89a)
now emit 40 on the CORRECT completed basis (formal cross-engine comparator = integration step). One converged
control-quality nit (§5.E tautology) folded to #89b. #89a is COMMITTED + PUSHED (basis `d4adbd99`, STATUS
`6f5c8c38`; origin + GIN, `.out` has a gin annex copy). NEXT = #89b WL operator unfreeze: ⛔ spec-confirm §3b
WL-side FIRST → orchestrator directive → 2 decision legs (Codex+Grok, rule 7) → Codex build → 2 build legs
(fresh Claude + Grok, Mathematica SERIALIZED; the live-operator rank WITHHELD, basis 40 already public,
frozen operator = public regression) → integration/re-adjudicate KINETIC+θ + fold 2 owed #88 hardenings →
#90 PY §3c content → close.
