# S11c-d SymPy build directive — decision review (G2 TRIGGER), round 1

**Artifact:** `directives/S11c_d_sympy_build_directive.md` (orchestrator-written thin SymPy build directive).
**Physics authority:** cleared spec `directives/S11c_d_SHARED_PHYSICS.md` v10 (`399a8516`).
**Legs (orchestrator-written → Codex + Grok):** Codex-sol `gpt-5.6-sol` xhigh; Grok `grok-4.6` high. Identical rendered
prompt `_legs/S11c_d_sympy_build_directive_decision_review_prompt.md`. Both verified the pinned facts computationally
(scripts + literal stdout under each leg's own `/tmp` copies; working tree untouched). Raw transcripts (outside repo):
`scratchpad/S11cd_dir_{codex,grok}.txt`.

## Verdicts
- **Grok: NOT CLEAR** — one spec-bar must-fix (F1: census incomplete). Import wiring, VALUE-level census, leak
  discipline, HELD-PHYSICS, script clauses, builder-lane, EMIT/EXPORT otherwise sound.
- **Codex: NOT SOUND AS-IS** — three must-fix (Q3 map-in-prohibition, Q4 named withheld order, Q7 Riesz export).
  Q1/Q2/Q5/Q6 independently verified sound.
- **Complementary + disagreeing:** Grok caught F1 (Codex missed it); Codex caught Q3/Q4/Q7 (Grok cleared them). Two
  legs each caught what the other missed — the two-engine design working as intended.

## Findings — orchestrator G4 verification + disposition (ALL FOUR ACCEPTED, must-fix)

**F1 (Grok) — the §2/§3 census was taken against the `VALUE` slot only; each `(α,ρ)` case of BOTH closed rows is a
FIVE-slot payload.** Slots: `VALUE`, `MULTIGRADE`, `DIMENSION_L_T_M`, `COMPUTED_BRANCH_BINDINGS`,
`FOURIER_PROFILE_BINDINGS`. The last two carry convention-bearing 3-D content the census never names:
`FOURIER_PROFILE_BINDINGS` = 4 Equalities per case defining the c2 hats as 3-D `d³Y` integrals with a numeric `(2π)`
measure **at the transfer argument only** (a ready-made 3-D convention map sitting on the consumed row);
`COMPUTED_BRANCH_BINDINGS` = 3 Equalities defining `s11cc2OutgoingNormalMomentum` as a 3-D `Piecewise` root at all
three momentum sites (`k_out`, `k_in`, `k_mid`). A builder that extracts only `VALUE` (what §2 says the payload is)
either leaves these 3-D definitions unreduced OR `subs`'s c2's 3-D convention as the reduction — both violate spec
§1c. **VERIFIED (mechanical):** `grep -oE "'(FOURIER_PROFILE_BINDINGS|COMPUTED_BRANCH_BINDINGS|VALUE|MULTIGRADE|
DIMENSION_L_T_M)'" scripts/S11c_c2_exports.py | sort | uniq -c` → each slot **8×** (= 4 cases × 2 rows). CONFIRMED.
→ **must-fix; changes what is computed.**

**Q3 (Codex) — the Fourier-reduction candidate map leaks through prohibitions.** §3 types `[L_W/(2π)]`, `(2π)²L_W`,
`δ²(Q_∥)` and the assembled relation `A_3D ≡ [L_W/(2π)]δ²(Q_∥)A_edge` (as prohibitions, copied from spec
§250–266). ⚠ Grok read these as "spec's own negative instruction, not a finding." **Adjudication:** the governing
build-skill rule (`.claude/skills/build/SKILL.md:186–193`, MEASURED) — *"a PROHIBITION leaks the answer as surely as
an assertion; write forbidden-pattern examples with PLACEHOLDER content, never the step's real content"* — governs the
**builder-facing** directive (the spec is the orchestrator-facing authority and may be more explicit). F1 shows the
3-D map factor literally sits on the row, so repeating the candidate in a prohibition doubly risks anchoring the
engine's normalization. Codex is correct. → **must-fix; changes what may be computed.** Fix: genericize (no specific
factors/relation).

**Q4 (Codex) — the withheld criterion is named as "the `O(1)` grating reductio"** (§6), disclosing the expected order
while claiming no expected order is handed to the builder. ⚠ Grok read leak discipline as clean. **Adjudication:**
`O(1)` is the strong-edge (S11c-e) order, out of S11c-d scope, so it is not strictly an S11c-d builder target — but
the fix (a value/order-free "the falsification acceptance criterion is withheld, orchestrator-side") is free and
strictly safer, and the build-skill spirit forbids stating the withheld thing's order. Accept the conservative fix.
→ **fix (leak discipline).**

**Q7 (Codex) — export membership contradicts the declared downstream boundary.** `S11c_decisions.md:52` hands
S11c-e "scattering amplitudes / **resonances / local spectrum**." The spec's canonical resonance/local-spectrum
outputs are `S11CD_BOUND_POLE_SET_AND_RIESZ_DATA` + `S11CD_BOUND_SPECTRAL_OVERLAP` (§3b). The directive marked "the
bound Riesz data" **EMIT-only** (§5), which excludes a declared handoff. ⚠ Grok read the export set as sound
(Riesz-as-emit-only consistent with "no capture protocol → spectral overlap only"). **Adjudication:** the decision
list literally lists "resonances / local spectrum" as what d hands e; the bound pole set + Riesz + overlap IS that
object. Codex's reading is grounded in the declared boundary; Grok's excludes it. With no S11c-e manifest yet, the
canonical pole/Riesz/overlap data (+ recursive closure) must remain EXPORTABLE under the declared scope, or the
exclusion must await an actual S11c-e manifest. Codex is correct. → **must-fix; can prevent S11c-e binding the
canonical spectrum/confinement operands.** Fix: move bound pole set + Riesz + overlap into the export list; remove
from emit-only.

## Independently verified SOUND (both legs, computational, agreeing)
- **Q1 import wiring:** `2441 + 44 + 70 → 2555`, additive, `overwrites==[]`, pairwise intersections empty; closed-root
  closure 213 / 273 symbol / 45 dimension edges; `assert_lookups_equal_manifest` fails on undeclared AND
  declared-unused; the open-vs-closed provenance hazard is real (both row-pairs pass the existence guard); casing
  `s11cc2Fieldtheta`; no `*Dimension` companion on the closed rows.
- **Q2 VALUE-level census:** hat/Integral/DiracDelta counts reproduce exactly (transfer 88/550, jet 100/136 each,
  Integral 218/526, atoms 26/90, 0 delta in closed rows, 3 delta on dtn_kernel); c1 snake_case hats absent; ONM at
  out/in/mid. (Completeness fails only via F1's sibling slots.)
- **Q5 HELD-PHYSICS:** reduced-representation rule, computed `K₀/A₀`, full vertex, surfaced reduced-block/kernel
  residual, no global `ω(k)`, two distinct photon-kill channels + no-weak-well, conditional N12, live c2 debt — all
  faithfully carried.
- **Q6 builder-lane + script clauses:** three clauses, structural rule, four corollaries, non-tautological §5a,
  FORM ablation, value-independent emission, build→verify→report→stop — all present.

## Disposition (round 1)
NOT CLEAR round 1. Fold all four verified must-fix findings (F1, Q3, Q4, Q7), then re-review (physics-bearing content
meets the spec bar → review-until-clear). Reviewed baseline preserved before the fold (`86b71819`).

## Round 2 — re-review of the folded directive → CLEAR (both legs SOUND)
Fresh legs, identical prompt `_legs/S11c_d_sympy_build_directive_decision_review_prompt_r2.md` (Codex-sol xhigh +
Grok-4.6 high), both verifying the folds AND re-scanning the whole directive computationally against the real rows.

- **Grok r2: CLEAR (sound)** — all four folds correct + complete; whole-directive re-scan found no new must-fix. One
  cosmetic **nit** (a leftover "§F1 slot below" label with no matching heading; does not change computation/claims) —
  fixed before commit.
- **Codex r2: CLEAR (sound)** — no must-fix, no nit. Independently confirmed: five-slot payload across all 8 cases;
  `FOURIER_PROFILE_BINDINGS` = transfer-only 3-D definitions, not valid reduced answers; ONM bindings cover
  `k_out`/`k_in`/`k_mid`; no convention-bearing element missed; no typed reduction relation or withheld acceptance
  order remains; bound pole/Riesz + spectral overlap correctly exported; import provenance, script clauses, leak
  discipline, HELD-PHYSICS, builder lane sound.

Both legs independently re-derived the census against the real rows (scripts + literal stdout under each leg's own
`/tmp`). Raw transcripts: `scratchpad/S11cd_dir_r2_{codex,grok}.txt`.

**CLEARED** (round-2 gate: both legs SOUND). Stopping rule met — nothing outstanding changes what is computed or may
be claimed (the fixed cosmetic nit did neither). The build directive is the governing build-mechanical authority for
the S11c-d SymPy engine; the cleared spec `399a8516` remains the physics authority. NEXT = the SymPy build
(`gpt-6-astra` high) → 2 build legs (fresh Claude + Grok) → blind WL → T7 → reconcile → step record.
