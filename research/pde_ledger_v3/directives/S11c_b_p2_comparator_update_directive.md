# P2 — cross-engine comparator update (extract_slab fix + §3a bridge) — decision list

Orchestrator-written ⇒ 2 decision legs (Codex + Grok) before the builder (rule 7). States what must be TRUE, not
the recipe. Physics-bearing (a mis-map manufactures false cross-engine agreement OR disagreement) ⇒ ⛔ never
blanket-collapse; a name-grep across engines is the bug.

## Objective
Update the cross-engine comparator so it correctly JOINS the CURRENT engines' emitted structure and a genuine PY↔WL
`row_residual` can be read. Targets: `scripts/S11c_b_cross_engine_comparator.py` (`extract_slab`, `extract_coupling`,
name tables) + the §3a name-bridge tables in `S11c_b_handcoded_comparison.py` / `S11c_b_adjudicated_comparison.py`.
Inputs are the fresh 4-case primaries `.out` from BOTH engines' `S11CB_PRIMARIES_ONLY` (WL = P1-WL gate; PY =
existing) — matching case-sets, so `row_residual`'s union+align consumes both. The comparator must neither
MANUFACTURE agreement (a wrong bridge / blanket-collapse) nor MANUFACTURE disagreement (a stale bridge joining
unlike objects).

## Why (context — verified, do not re-litigate)
Since the comparator was written, PY folded (constraint-reduced U/E_W rows carrying `−∇μ_θ`; θ-row =
`evolution_mass_balance − Σ closure_shape_deriv`; μ_θ separate; face rows) and #90 (face+response into the operator),
and WL renamed its §3a coefficients to `gammaWJET*`/`gammaMUJET*` at `d4adbd99` — AFTER the bridge tables
(`70164909`). Verified (rule 13): `extract_slab` still maps `THETA_BALANCE → ("MU_THETA","THETA")`
(`comparator.py:766-770`) — post-fold that row is mass-evolution, not μ_θ; and the 12 `gammaWidth*` bridge entries
(`handcoded_comparison.py:107-118`) match NOTHING post-rename ⇒ currently 0 of ~30 §3a coefficients bridged.

## Governing invariants
- ⛔ Every cross-engine identification is justified by the OBJECT it names, not the name. Never blanket-collapse.
  [[feedback_handcode_comparison_never_blanket_collapse]] [[feedback_name_grep_across_engines_is_the_bug]]
- ⛔ Emit operands + residual, then guard; NO expected-value acceptance (rule 5). The residual is OUR finding.
- ⛔ Every bridge/identification must be TESTABLE by one-sided corruption: breaking one MOVES the residual.

## extract_slab structural correctness — what must be TRUE
D1. PY `THETA_BALANCE` is the θ-row = the sourced mass-evolution balance (`evolution_mass_balance − Σ
    closure_shape_deriv`), NOT μ_θ; its cross-engine partner is WL `MASS_EVOLUTION_ROW`. The comparator must NOT
    label PY `THETA_BALANCE.EXPANDED` as the canonical `MU_THETA`, and must handle its new sub-keys
    (`SOURCE_OPERAND`, `EVOLUTION_TERM_ORIGINS`, `EXPANDED`).
D2. The μ_θ↔μ_θ comparison flows ONLY through the `MU_THETA_OPERATOR` primary (`extract_mu`). After D1 no slab-row
    path may also emit a canonical `MU_THETA` object (else the join double-counts μ_θ).
D3. The PY mass-evolution content under test is the full sourced sum (from `THETA_BALANCE`), not the single
    `ADVECTIVE_MASS_OPERAND` summand (which stays an origin/diagnostic, not the row's partner).
D4. WL `CENTER_FACE_GENERALIZED_ROW` and PY's `FACE_GENERALIZED_FORCE_ROWS.CENTER_FACE_GENERALIZED_ROW` join to each
    other — or both are explicitly excluded with a recorded reason. No engine-only orphan.
D5. U/E_W stay joined `U_BODY_BALANCE↔U_MOMENTUM_ROWS`, `E_W_BALANCE↔THICKNESS_ROW`; both `EXPANDED` values are the
    constraint-reduced + face-folded rows (shape aligns; CONTENT is what the residual tests).
D6. Every PY slab sub-object and every WL slab sub-object is JOINED to its partner or explicitly recorded as
    engine-only with a stated reason — NO silent reclassification (the failure mode: THETA_BALANCE's new sub-keys
    silently becoming unjoined PY-only).
D6b. `extract_coupling` is RE-VERIFIED against the post-#90 kernel structure (face-force / closure changes in
    `build_kernel`); the slots it reads still name the same objects on both engines, or it is updated. (Open — must
    be checked, not assumed clean.)

## §3a coefficient bridge — what must be TRUE
D7. Each WL `gammaWJET<SUFFIX>`/`gammaMUJET<SUFFIX>` bridges to the PY coefficient `gamma_s11cb_{w_bg|mu_r_bg}_NN`
    that multiplies the SAME O(3)-Kronecker invariant (both engines enumerate the same 15/source family since
    #89/#89a).
D8. The pairing is established by matching the INVARIANT each coefficient multiplies (PY `ENERGY_BASIS_NEW_INVARIANTS`
    ↔ WL new-invariant records, under the already-complete profile-jet/DOF renames), COMPUTED — ⛔ NEVER by assuming
    PY's positional index `NN` corresponds to a suffix by position. A wrong index manufactures false agreement.
    [[feedback_name_binding_is_unpoliced]]
D9. The bridge covers every coefficient appearing in the compared rows/kernel. The stale `gammaWidth*`/`gammaModulus*`
    entries (matching nothing) are removed/re-keyed; the never-bridged suffixes are added.
D10. The completed bridge is TESTABLE: corrupting one coefficient pairing MOVES the cross-engine residual.

## scope + provenance — what must be TRUE
D11. The change is validated by RUNNING `row_residual` against the FRESH 4-case primaries `.out` from both current
     engines (P1-WL gate output + PY `S11CB_PRIMARIES_ONLY`), not the stale committed `.out`. Case-sets match (both
     emit all branches×densities), so no case-intersection is needed; if a scope mismatch surfaces, it is a finding
     to record, not to silence by collapsing.
D12. S11c-a is NOT regressed (the bridge additions are S11c-b-scoped; S11c-a joins/renames still hold).
D13. Output states, per joined object: operand_PY, operand_WL, residual — no verdict token, no assert on the
     residual.

## Legs
Orchestrator-written ⇒ 2 decision legs (Codex + Grok). Build (Codex) ⇒ 2 build legs (fresh Claude agent + Grok),
FORM ablations mandatory: corrupt a coefficient bridge → residual moves; corrupt the THETA_BALANCE→MASS_EVOLUTION
identification → residual moves; a stale/no-op bridge is byte-identical and must be caught. The build validates
against the fresh `.out` (available after P1-WL clears + the production emit).
