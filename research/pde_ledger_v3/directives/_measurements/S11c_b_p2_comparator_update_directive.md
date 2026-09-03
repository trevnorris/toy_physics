# P2 comparator update — decision-leg review record (2 legs, computation-backed, convergent → v1 REJECTED)

Directive `directives/S11c_b_p2_comparator_update_directive.md` (orchestrator-written ⇒ Codex + Grok decision legs,
rule 7). Prompt `directives/_legs/S11c_b_p2_comparator_decision_review.md`. Raw logs
`~/.s11_build/S11c_b_p2/decision_{codex,grok}.log`. **v1 REJECTED by BOTH legs (10/9 defects, file:line + literal
probe output). Both independently CLOSED D6b.** These reshape P2's scope; a v2 fold is owed before any builder.

## The scope-reshaping findings (both legs, orchestrator rule-13 verified)
1. **`row_residual.py` is IN P2's blast radius — the objective wrongly omitted it.**
   - (Codex) `row_residual.py:427` computes `residual = py_trunc − (wl_complete − face_attributed)`, removing face
     from WL ONLY (`:721`). Post-#90 PY ALSO carries face in its complete U/E rows (`sympy_audit.py:2980/2998/3028`;
     WL `…audit.wl:1345`) ⇒ this now manufactures a FALSE disagreement equal to the face contribution. The stale
     one-sided face subtraction — not `extract_coupling` — is the #90 structural break.
   - (Both) `row_residual` never loads `MU_THETA_OPERATOR` (its `_family_cases` calls are only SLAB_OPERATOR,
     SLAB_OPERATOR_TERM_ORIGINS, COUPLING_KERNEL, ADMISSIBILITY_OPERATOR_OPERAND — rule-13 verified `:566-602`). Once
     D1 removes slab's false canonical `MU_THETA`, μ_θ SILENTLY vanishes from the residual D11 uses to validate.
   - (Both) `SLAB_OBJECTS` (`row_residual.py:43`) still contains `MU_THETA`; a one-sided SLAB object RAISES
     (`:616-620`). The object remap (D1/D2/D4) can make even a matching 4-case `.out` fail to load.
2. **§3a bridge needs SCALE FACTORS, not string renames (the manufactured-agreement trap in its exact form).**
   (Codex measured) `EXACT_UNIQUE 0 / SCALED_UNIQUE 30` per branch: every W-pair satisfies `I_PY = W_0·I_WL`, every
   μ-pair `I_PY = μ_R·I_WL` (WL spurion `gradient[anchoredWidth]/WZero`, `…/muR` `…audit.wl:722`; PY raw jets
   `grad_W` `sympy_audit.py:182`; rule-13 verified WL uses `(muR/WZero)`). A direct `gammaWJET*→gamma_s11cb_w_bg_NN`
   rename leaves `1/W_0` ⇒ false disagreement. The bridge must be an EXPRESSION bridge
   (`gammaWJET*→W_0·gamma_s11cb_w_bg_NN`), but the committed bridge only accepts string→string renames
   (`handcoded_comparison.py:72,279`) ⇒ the bridge REPRESENTATION must change. ⛔ Folding `W_0` OUT of the invariant
   to force a match is a blanket collapse of background-jet normalization — the scale must be EMITTED.
   `COEFFICIENT_STANDARD_NAME` is `gamma_W_<suffix>` (`…audit.wl:629-635`), NOT `gamma_s11cb_*_NN` — a no-op bridge.
3. **A bijection CERTIFICATE is required — D10 ("corrupt one pairing → residual moves") passes a WRONG pairing.**
   (Codex installed a wrong first pairing) `DROP_WRONG_PAIRING_MOVES_RESIDUAL True` while `WRONG_PAIR_INVARIANT_MATCH
   False`. Need a per-coefficient printed certificate: both invariant operands, the adjudicated normalization factor,
   exact residual, one unique partner, 30/30 coverage, both branches — a swapped target must FAIL it.

## The other convergent defects (fold into v2)
4. **D1 θ sub-keys unclassified** (both): `SOURCE_OPERAND ≡ EXPANDED` (same `mass_balance`, `sympy_audit.py:3043/3045`);
   `EVOLUTION_TERM_ORIGINS` is provenance (also in the origin family). Join ONLY `EXPANDED`; `SOURCE_OPERAND` =
   duplicate diagnostic; origins → origin family. Leaving both risks the `_family_cases` DUPLICATE-AXIS raise
   (`adjudicated_comparison.py:1091-1096`, rule-13 verified) if `ADVECTIVE→MASS_EVOLUTION_ROW` (`comparator.py:777`)
   is also left.
5. **D4 CENTER_FACE must be JOINED (has a real partner), NOT optionally excluded**; but ⛔ do NOT re-join
   `FACE_GENERALIZED_FORCE_ROWS.{U,E_W}` (already inside the D5 rows = double-count) or `THETA_FACE_FLUX` (already in
   the mass row); watch the virtual-work SIGN (PY `+∂/∂δ_vζ_c` vs WL `−linearVirtualVariation` — join same virtual
   displacement or manufacture a sign residual). CENTER_FACE also not in `row_residual.SLAB_OBJECTS`.
6. **D6 needs an EXHAUSTIVE path-by-path disposition table** (both), not "engine-only with reason": many post-fold
   keys are unnamed orphans (`FACE_GENERALIZED_FORCE_ROWS.*`, `MU_THETA_FACE_BINDING`, `U/E_W_BALANCE.{LOCAL,
   DIVERGENCE_FLUX}`, WL `DIVERGENCE_FORM_SOURCE.*`/`BACKGROUND_FIRST_JETS`/`THETA_VARIATION_RULE`, origin-name
   mismatches FACE/FACE_VIRTUAL_WORK, FLUX/FACE_FLUX). `FACE_SHAPE_SUBSTRATE` whole-container join is unsound (PY 20
   categories vs WL 8 fields) — needs an adapter + cardinality guard.
7. **D9 must RETIRE `PROTECTED_ATOM_NAMES`** (both): still holds PY `07/10` + dead `gammaWidth*`; these divert
   nonzero residuals to `PROTECTED_UNREDUCED` (`adjudicated_comparison.py:953`) and are excluded from ablation-touch
   (`:1478`) ⇒ HIDES disagreement once bridged. D10 must exercise formerly-protected 07/10.
8. **D9/`_extra_basic` blanket applied→bare erases arguments** (Codex): `widthBackground[x…] −
   widthBackground[chi…] PARSED 0` (`comparator.py:379/396`; `handcoded:159/300`) — erases a real base-vs-shifted
   material/background-anchor distinction ⇒ manufactured agreement. Need an argument-sensitive fixture (only inert
   evaluation args removable).
9. **D11 validation too weak** (both): union-then-align does not require the expected Cartesian set (a case absent
   from BOTH vanishes silently); energy families are 2-branch not 4-case; row_residual hardcodes default `.out`
   paths (`:1084/1095`) — need expected-set guards per family/object + selectable paths / emitted provenance so a
   builder cannot validate against stale defaults.
D6b (both): `extract_coupling` needs NO structural update — both θ channels are already mass-evolution
(PY `weak_operator_blocks` `:3728`; WL forward density `…audit.wl:1690`). ⛔ Do not retarget the θ channel at μ_θ.

## Consequence for the plan
P2 (v1) treated the comparator update as `extract_slab` + §3a name-bridge. The legs show it is a LARGER,
physics-bearing job spanning `extract_slab`, the §3a bridge (with a NEW scale-factor representation + a bijection
certificate + an argument-sensitive applied→bare guard), AND `row_residual.py` (remove the stale one-sided face
subtraction, load `MU_THETA_OPERATOR`, expected-set guards, retire `PROTECTED_ATOM_NAMES`). ⇒ recommend splitting
into two focused builders (shorter decision lists = lower defect rate): **P2a** slab-row join + `row_residual`
structural fixes; **P2b** §3a coefficient bridge (scale factors + certificate + guard). Each own decision list + 2
legs. Surfaced to the user before proceeding (scope grew materially beyond the resume prompt).
