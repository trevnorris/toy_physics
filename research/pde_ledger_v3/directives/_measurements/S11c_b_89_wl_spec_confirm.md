# S11c-b #89 WL side — spec-confirm (§1d/§3a binds WL; no spec change) + defect localization

Prereq to drafting the WL repair directive (resume step 1: "spec-confirm §1d/§3a applies WL-side →
draft"). This records the claims the directive's §1–§2 rest on, each with the command/lines behind it
(rule 2). No engine changed; no spec changed.

## Claim 1 — §3a (244–256) is ENGINE-NEUTRAL and names the OBJECT; it binds WL
`research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md:244–256` (read verbatim) says construct
"**every** scalar bilinear in the DOF data allowed by the S11b symmetry group with those spurions carrying
their transformation" and "Independence is judged as field bilinears with B1's constraint not applied;
carry **every** independent invariant with its own free symbolic coefficient." §1d:163–168 adds that a
variable-coefficient representative difference (`c∇·F ≡ −(∇c)·F`) is "physics in the operator/kernel, not a
representational identity … a computed object adjudicated after the run."
⇒ The spec names an OBJECT (the complete independent field-bilinear family), NOT a mechanism — it prescribes
neither PY's divergence quotient nor WL's hand-coded list. Both engines owe the SAME object by their own
route. **No spec change is required; the WL fix is a spec-compliance repair of an incomplete enumeration**,
parallel to the PY #89 repair, DIFFERENT mechanism.

## Claim 2 — WL implements a hand-picked SUBSET (8/source), not "every" bilinear (code-verified)
The WL new-invariant family is sourced from FOUR parallel length-8 hand-coded lists, keyed by position:
- `newInvariantExpressions` (`mathematica/S11c_b_brane_operator_mathematica_audit.wl:417–436`) — 8 forms.
- `newInvariantLabels` (L410–415), `newCoefficientSymbols` (L390–398), `newCoefficientNames` (L400–408) — 8 each.
Both `constructEnergyData` (L501–520) and `constructFullFieldBackgroundEnergy` (L570–576) consume these,
so completing them propagates to BOTH the emitted basis AND the operator input (via
`basisRepresentativeIndices`, L622). The count is COMPUTED, not asserted: `independentRepresentativeIndices`
(L599–618) takes incremental `MatrixRank` (L612) over `rankVariables` (L590–597). WL uses genuine symbolic
`gradient[anchoredWidth]` (L568), not a frozen lookup-map — so unlike PY it has no separate §3b frozen-map
defect; its ONLY §3a defect is the short enumeration.

## Claim 3 — WL's committed public count = 26; new family = 8/source (measured)
`grep -a -o 'ENERGY_BASIS_COUNT[^}]*}' mathematica/out/S11c_b_brane_operator_mathematica_audit.out`
→ `"LAB_HELD" -> <|"COUNT_OPERAND" -> 26, "RANK_SELECTED_INDICES" -> {1,…,26}|>` (all 10 uniform + 8 WJET
+ 8 MUJET independent, full rank 26). `grep -a -o 'WJET_[A-Z_]*\|MUJET_[A-Z_]*' … | sort -u` → 8 WJET + 8
MUJET labels. ⇒ **26 is the WL public regression value** (the frozen coincidence both engines currently
share). The completed count is WITHHELD (rule 5).

## Claim 4 — the correct family is 15/source of FIRST-JET field-bilinears (raw rank 15) — #87
`_measurements/S11c_b_87_wl_subspace_result.md`: WL's 8 ⊊ span(correct 15), raw rank(correct 15)=15 AND
quotient rank(correct 15)=15 (nullity 0). ⇒ The correct basis lives in the RAW first-jet field-bilinear
space; with the first jet live the divergence quotient is trivial (nullity 0), so WL's plain `MatrixRank`
over its first-jet `rankVariables` is ADEQUATE to reach the correct rank once the family is complete. No
Hessian atom is a basis element. (This is the DIFFERENT-mechanism claim: WL needs no quotient and no
Hessian un-freeze — only a complete enumeration.)

## ⚠ Rule-17 residual — must become a COMPUTED control, not an assertion (Claim 4 is "only if")
`independentRepresentativeIndices` (L602) uses the standard `applyProfile` (L363), whose
`profileDerivativeRules` (L297–310) zeroes every 2nd/3rd jet (`higherRules`, L304–307). WL also has a
live-higher-jet path `applyBackgroundProfileWithGeneratedJets` (L360) via
`profileRulesRetainingGeneratedJets` (L339–358). Claim 4 predicts the completed basis carries no higher-jet
atom, so the two paths give the SAME rank — but that must be COMPUTED and emitted (frozen-path rank vs
live-path rank); a difference is a finding (a completed-family form carries a higher jet the freeze
truncates → rule 17 / [[feedback_basis_independence_must_not_freeze_spurion]]). This is the WL analog of
"verify SPANS not COUNTS."

## Disposition
Spec-confirm PASSES: §3a binds WL, no spec change, object named. Directive obligations: (1) systematic
COMPLETE enumeration of the symmetry-allowed first-jet field-bilinear family (⛔ not a hand-extended list,
⛔ not a transliteration of PY's `enumerate_new_candidates` — WL builds its own, rule 1); (2) count computed
by rank, WITHHELD (public regression = restrict-to-original-8 ⇒ 26); (3) the frozen-vs-live-path rank
control above. NEXT = draft `S11c_b_89_wl_3a_repair_directive.md` → 2 decision legs (Codex + Grok).
