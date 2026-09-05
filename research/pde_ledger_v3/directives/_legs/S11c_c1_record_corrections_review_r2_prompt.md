# Independent physics review (ROUND 2) — the S11c-c1 RECORD corrections, after folding round-1 findings

## Context
Round 1 (two legs, Codex-sol + Grok) reviewed a set of record-only corrections to the S11c-c1 records (⛔ the
ENGINES/EXPORTS are NOT edited and are NOT under review — c1's computed objects STAND; both round-1 legs
independently re-confirmed NO REOPEN by re-deriving the density recovery, residual 0). Round 1 returned exactly two
defects, now FOLDED with the legs' own prescribed wording. Your job: verify the two folds are physically CORRECT and
introduce NO new defect, and that the rest of the corrections remain accurate + complete. This is a DOCUMENT review.

## The two round-1 findings that were folded (verify each fix)
1. **Grazing inverse mis-attribution.** The *Exact grazing* bullet in `steps/S11c_c1_curved_bulk_closure.md`
   previously named `[I+(Λ_A/ρ_m²)Z]⁻¹` as the `~1/η` nonanalytic object. It now names the **DtN inverse `N⁻¹`**
   (hence `Z=iρ_mω·N⁻¹`) as the singular object, and states the permeable face-closure resolvent
   `[I+(Λ_A/ρ_m²)Z]⁻¹` is **regular** at grazing for `Λ_A≠0` (`O(η)`; `[I+(Λ_A/ρ_m²)Z]⁻¹·Z → ρ_m²/Λ_A`). Confirm:
   is the corrected attribution physically right (`N⁻¹`/`Z` carry the `~1/η` pole; the permeable resolvent does
   not), and is the retained scope (`NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER`, `Z₁` non-grazing only, `‖N₀⁻¹N₁‖≪1`)
   still correct and sufficient?
2. **Energy-residual parenthetical sign.** The *Carry-forward › Energy-residual orientation* item previously glossed
   `P_face + P_∞ = 0` as "(traction work + **minus** outgoing far-field Poynting = 0)", which is self-contradictory.
   It now reads "(traction work + **positive outgoing** far-field Poynting = 0, with `P_∞` the outgoing flux)" and
   keeps the separate, correct instruction that to form a **vanishing `A−B` residual** the subtraction operand `B` is
   minus outgoing Poynting (`A−B = P_face−(−P_∞) = 0`), or emit `A+B`. Confirm the identity and the residual-operand
   instruction are both now correct.

## Also re-confirm (should be unchanged from round 1, where both legs cleared them)
- Density channel is O(εη), first-order-in-η; the earlier "harmless O(η²)" was wrong. (seal 5, both records)
- No-reopen / recoverability: PY carries `rho_br_bg_rho4_constant` opaquely (0 derivatives, no `w1_profile` baked
  in); re-binding to `background_density_map` recovers `−μ_θ·w₁/ρ_br` exactly (residual 0). c1 STANDS.
- Independence scoping is accurate (supplied composition/expected-values = fidelity; the two-momentum kernel is the
  independently-confirmed object).
- The 7 lower-severity carry-forward items are each accurately dispositioned and the set is complete vs the two
  source reports.

## What you are handed (read SOURCES first, form your own view, THEN read the corrected records)
- The two ROUND-1 leg reports (your predecessors' full computation-backed findings):
  `research/pde_ledger_v3/directives/_measurements/S11c_c1_record_corrections_review_{codex_sol,grok}.txt`.
- The two ORIGINAL source retro reports: `…/_measurements/S11c_c1_spec_retro_review_{grok,codex_sol}.txt`, and the
  adjudication `…/_measurements/S11c_c1_spec_retro_review_adjudication.md`.
- The corrected records: `steps/S11c_c1_curved_bulk_closure.md` and
  `directives/_measurements/S11c_c1_comparator_reconcile.md` (the full edit set is
  `git -C /var/projects/toy_physics diff` against HEAD `f8509b7a`).
- The real engine + fold if you want to re-check recoverability: `scripts/S11c_c1_exports.py`,
  `scripts/S11c_b_exports.py`, `scripts/ledger_fold.py` (module is under `scripts/`, `load_model` at line 102).

## Required method (DOCUMENT branch)
Read the sources first, form your own view, THEN read the corrected records; quote both sides for every finding.
⛔ A prose "I re-derived it and agree" is worth nothing — for finding 1 (which inverse is singular) save your SymPy
script AND its literal stdout to named `/tmp` paths and cite them. ⛔ Do not modify the working tree.

## Physics filter
Report a finding ONLY if it catches a way the corrected record (a) misstates the physics, (b) over- or under-scopes
what may be claimed, or (c) drops/misattributes a claim-level item a source report raised. ⛔ Do not re-litigate the
settled no-reopen unless you have a concrete computation showing it is wrong.

## Output
Ranked findings (MUST/SHOULD/NIT) with quoted record text + quoted source + (for finding 1) script+stdout paths.
End with an explicit verdict: **CLEAR** (the two folds are correct and the corrections are accurate + complete) or
the exact list to fix. If you conclude no-reopen is actually wrong, say **REOPEN c1** and show the irrecoverable channel.
