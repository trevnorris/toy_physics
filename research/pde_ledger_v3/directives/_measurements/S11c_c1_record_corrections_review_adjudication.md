# S11c-c1 RECORD corrections — review-until-clear gate record (rule 2 + rules 8–10, 13)

The retrospective spec review (`S11c_c1_spec_retro_review_adjudication.md`) earned **record-only** corrections to
two committed c1 RECORDS (⛔ NOT the engines/exports — c1's computed objects STAND, no reopen). Those corrections
are physics-bearing prose ⇒ **review-until-clear** (not the decision-list one-pass). This file is the gate record:
two rounds, both legs computation-backed, converged to CLEAR. Orchestrator-written ⇒ legs = **Codex sol xhigh +
Grok**.

## The corrected records (engines untouched)
- `steps/S11c_c1_curved_bulk_closure.md` — rule-17 density note (O(εη)), new *Exact grazing* bullet, new
  *Independence is SCOPED* bullet + a scope pointer on the blindness claim, new *## Carry-forward corrections*
  subsection (7 items).
- `directives/_measurements/S11c_c1_comparator_reconcile.md` §3.5 — "harmless O(η²)" → O(εη)-recoverable
  representational framing.

## Commands (identical prompt to both legs each round)
```
# Round 1 — prompt _legs/S11c_c1_record_corrections_review_prompt.md
codex exec -m gpt-5.6-sol -c model_reasoning_effort=xhigh --sandbox danger-full-access \
  "$(</abs/_legs/S11c_c1_record_corrections_review_prompt.md)" > _measurements/S11c_c1_record_corrections_review_codex_sol.txt 2>&1
grok --prompt-file /abs/_legs/S11c_c1_record_corrections_review_prompt.md --cwd .../pde_ledger_v3 \
  --model grok-4.6 --effort high --permission-mode bypassPermissions --output-format plain > …_grok.txt 2>&1
# Round 2 — prompt _legs/S11c_c1_record_corrections_review_r2_prompt.md → …_r2_{codex_sol,grok}.txt
```
⚠ Launch-discipline notes (all caught by deliverable-verification, ⛔ never by the masked background exit-0): a first
launch died because the outside-repo log dir did not exist (redirect failed, leg never ran); a later Codex launch
got an EMPTY prompt because a persisted `cd` made a RELATIVE `"$(<…)"` path resolve to nothing (~3.2k-token
"What would you like me to work on?" tell) — fixed by an ABSOLUTE prompt path. Both round-1 legs were once KILLED
simultaneously by a phone-interrupt (mobile back-button ≈ Stop); relaunched (stateless).
⚠ The Codex report files here are the **clean extracted report blocks**; the full raw transcripts (3.4 MB + 2.5 MB,
each carrying the whole folded export verbatim) live OUTSIDE the repo at
`/var/projects/toy_physics_ext_logs/S11c_c1_record_corrections_review{,_r2}_codex_sol_RAW.txt` (tree hygiene, ⛔ not
a blindness claim).

## Round 1 — two findings (both legs computation-backed; both independently NO-REOPEN)
Both legs re-ran the density recovery against the REAL export (residual 0) and independently confirmed **no reopen**
— my override now verified a 3rd and 4th time. They converged on exactly two RECORD defects:
1. **Grazing inverse mis-attribution** (Codex MUST / Grok SHOULD). The *Exact grazing* bullet named
   `[I+(Λ_A/ρ_m²)Z]⁻¹` as the `~1/η` object. That is wrong: the singular object is the **DtN inverse `N⁻¹`** (hence
   `Z=iρ_mω·N⁻¹`); the permeable face-closure resolvent `[I+(Λ_A/ρ_m²)Z]⁻¹` is **regular** at grazing for `Λ_A≠0`
   (`O(η)`, `[…]⁻¹·Z → ρ_m²/Λ_A`). Both printed `Z=C·N⁻¹ 1/η-pole = True`, `D⁻¹ 1/η-pole = False`.
2. **Energy-residual parenthetical sign** (Codex SHOULD). "`P_face+P_∞=0` (traction work + **minus** outgoing …)"
   was self-contradictory; with `P_∞` = positive outgoing flux the identity is `P_face+P_∞=0`, and the "minus"
   belongs only to the subtraction operand `B=−P_∞` in an `A−B` residual.

**Folds** (the legs' own prescribed wording): (1) name `N⁻¹`/`Z` as the singular object + state the permeable
resolvent is regular; keep scope `NOT_ESTABLISHED_AT_FIRST_SHAPE_ORDER`, `Z₁` non-grazing only, `‖N₀⁻¹N₁‖≪1`.
(2) "(traction work + **positive outgoing** far-field Poynting = 0, with `P_∞` the outgoing flux)"; keep the
`A−B` operand-`B=−P_∞` / emit-`A+B` instruction. My verification (rule 13): `N=ηB ⇒ N⁻¹~1/η`; `D=I+aZ ⇒ D⁻¹≈(aZ)⁻¹
=O(η)`, `D⁻¹Z→1/a=ρ_m²/Λ_A` — both folds physically correct.

## Round 2 — both legs CLEAR (zero findings)
Re-legged on the folded records. **Codex sol: CLEAR** (MUST none / SHOULD none / NIT none — `1/η`-pole verdicts as
above, `RESOLVENT*Z → (ρ_m²/Λ_A)I` residual 0; energy `P_face+P_out=0`, `A−B=0` for `B=−P_∞`; 24 opaque density
leaves residual 0; independence + 7 carry-forwards accurate). **Grok: CLEAR** (same computations independently; plus
an output-name audit confirming FLUX/TRACTION density leaves are opaque too — 24 leaves = 4 keys × 3 outputs).
Reports: `…_record_corrections_review_r2_{codex_sol,grok}.txt`.

## Verdict — STEP 1 review-until-clear COMPLETE
Both legs CLEAR, nothing outstanding that changes what is computed or may be claimed (G4 stopping rule met). c1's
**engines and exports STAND — NO REOPEN**; the density O(εη) channel is recoverable in c2 by re-binding
`rho_br_bg_rho4_constant` to `background_density_map` (⚠ c2 not yet built). The seven carry-forward items are
committed as owed (energy orientation, `h_s`/`N`–`Z` terminology, density-multiplication-operator naming **for the
c2 directive**, `K_a` Hermitian, evanescent grades η²/η·σ_W/σ_W², drain-tilt vs convection, flat-`Z₀` leakage).
⚠ For c2 wiring: the fold module is `scripts/ledger_fold.py` (`load_model` at line 102), ⛔ not `reduction/`.
