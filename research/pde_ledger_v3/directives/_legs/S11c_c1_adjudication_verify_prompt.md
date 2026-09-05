# Verify — the S11c-c1 retro-review adjudication + the compact-prep state docs

The orchestrator ran a retrospective review of the S11c-c1 SHARED_PHYSICS spec (it had been folded once, never
re-legged). The two legs **disagreed**: Grok CLEAR, Codex-sol **BLOCK — reopen c1**. The orchestrator **overrode
the reopen verdict** and ruled **c1's engines/exports STAND** (no reopen), with only 3 claim-level corrections owed
to the c1 *records*. Your job is to verify that override and the state docs are correct before a context compaction.
⛔ Do not rubber-stamp; if the override is wrong, say so plainly and say c1 must reopen.

## The load-bearing question (settle it independently, with computation)
Codex-sol's reopen rested on: the RHO4_CONSTANT face response dropped a genuine **O(εη)** background-density channel
(`d(μ_s)/dη|₀ = −μ_θ·w₁/ρ_br ≠ 0`), so "δp_s/J_s agree" is invalid → regenerate the response. The orchestrator's
override rests on: PY carries the density **opaquely** (`1/rho_br_bg_rho4_constant`, 0 derivatives, never combined
with `w1_profile`), so **c2 recovers the exact channel by re-binding** `rho_br_bg_rho4_constant` to
`background_density_map` — the channel is NOT irrecoverably dropped, so c1 need not reopen.

**Verify with your own inspection + SymPy:**
1. In `scripts/S11c_c1_exports.py` (fold it via `ledger_fold.load_model("scripts/S11c_b_exports.py",
   "scripts/S11c_c1_exports.py")`), does `s11c_c1_face_response_coeffs` for RHO4_CONSTANT carry
   `rho_br_bg_rho4_constant` as an **opaque symbol** (a bare `1/ρ`), with **no** `w1_profile` combined with it and
   **no** derivative taken on it? Quote the actual coefficient.
2. Does `background_density_map` in the fold carry the live relation `rho_br_bg_rho4_constant = W_bg·ρ_br/W_0 =
   ρ_br(1+η·w₁)`?
3. Does substituting that relation into the opaque `μ_θ/rho_br_bg_rho4_constant` and expanding to first order in η
   **recover** `−μ_θ·w₁/ρ_br` (the dropped channel)? Run it.
4. ⇒ **Is the orchestrator's override correct** (opaque ⇒ recoverable ⇒ no reopen, c2 recovers per its §3d.1), or is
   Codex-sol's original reopen right (the channel is baked in / lost / the "agree" is a blanket-collapse of PY's
   constant onto WL's field)? The distinguishing test: is anything about PY's exported RHO4_CONSTANT response
   **irrecoverable** by a downstream re-bind of the opaque density symbol? If yes → reopen; if no → no reopen.

## Also verify (accuracy)
- The adjudication record `directives/_measurements/S11c_c1_spec_retro_review_adjudication.md` — are its claims and
  its transcript of both legs accurate against the two reports (`_measurements/S11c_c1_spec_retro_review_{grok,codex_sol}.txt`)?
- The **3 owed corrections** (seal-5 O(εη)-recoverable-representational; grazing NOT_ESTABLISHED; independence
  scoping) — are they correctly stated, and are they the COMPLETE set of claim-level defects the retro-review found,
  or does either report carry a MUST-FIX the adjudication dropped? In particular: does Codex-sol's grazing MUST-FIX
  (the inverse is nonanalytic at exact grazing) and its rule-5/3-leakage MUST-FIX map to corrections 2 and 3, and
  is downgrading Codex-sol's reopen (correction "none — c1 stands") justified by the recoverability finding?
- The STATUS.md top clause (2026-09-05) and its consistency with the adjudication.

## Method
This is a document + computation verification. ⛔ A prose claim about the export is worth nothing — inspect the real
`scripts/S11c_c1_exports.py`/`ledger_fold.py` and run your SymPy; save script+stdout to named paths. Do NOT modify
the working tree.

## Output
1. Your independent verdict on the override: **NO REOPEN (override correct)** or **REOPEN (override wrong)**, with the
   quoted coefficient + the recovery computation. 2. Any inaccuracy in the adjudication record or STATUS clause.
3. Any MUST-FIX the adjudication dropped. End: **CLEAR TO COMPACT** or the exact list to fix first.
