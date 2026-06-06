---
unit_id: 087
batch: III.5
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T00:07:48Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 087

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what is specified. Do NOT touch the load-bearing reduction asserts (`unblocked zeta_req`, `d zeta_req exact formula`). Do NOT change the three window literal VALUES — they are correct against the notes. Do NOT touch paper.tex or notes/.

After editing, RUN both scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:69-75`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:46-63`

**Issue:**
The `rho_suff/rho_fail/rho_max` "cross-check against upstream stage-086" assertions compare each window literal against the SAME literal re-typed in the same file (sympy py:58 vs py:73; mathematica wl:57 vs wl:61). Output records confirm self-comparison: sympy diff ~1e-16 (15-vs-30-digit parse rounding only), mathematica `diff = 0.`. The comments (py:69-72, wl:46-55) claim these guard against renumber/transcription drift, but a hollow self-comparison cannot detect a mistyped literal because both sides change together. The Mathematica `zeta_*` numeric checks (wl:73-75) are similarly self-referential (target = `rho_X - 1` re-typed). The three literal values are CORRECT versus the notes — this is a verification-mechanism defect, not a value mismatch.

**Required change:**
Replace the per-literal self-comparisons with structural relations among the three literals that genuinely fail if any one is mistyped, and stop the comment overclaiming an upstream cross-check. Apply symmetrically to BOTH engines.

SymPy (`...stage087..._sympy_audit.py`):
1. Reword the comment block at py:69-72 to drop the "cross-check ... against the upstream stage 086 quoted values" framing; state plainly that these are the Family-1 window literals carried from the notes and that the ordering/gap relations below are the sanity check.
2. Replace the three `expect_close(... vs stage-086 ...)` self-comparisons (py:73-75) with ordering/gap asserts that can fail under bad transcription. Add a helper or inline checks asserting:
   - `rho_suff < rho_fail` (strict ordering)
   - `rho_fail < rho_max` (strict ordering)
   - `0 < rho_max - rho_fail < sp.Float("1e-6")` (tight constructive-ceiling gap)
   Print each comparison result; raise `AssertionError` if any relation is violated. Keep printing the three literal values so they remain in the transcript.

Mathematica (`...stage087..._mathematica_audit.wl`):
1. Reword the comment block at wl:46-55 the same way (drop the "anchors them against the upstream stage-086 paper values" overclaim).
2. Replace the three `expectApprox[... vs stage-086 ...]` self-comparisons (wl:61-63) with the same ordering/gap checks via an `If[!(rhoSuff < rhoFail < rhoMax) || !(0 < rhoMax - rhoFail < 10^-6), fail[...], pass[...]]` pattern (or three separate `If`/`fail` checks). Keep printing the three literal values.
3. The `zeta_*` numeric checks at wl:73-75 are ALSO self-referential (target = `rho_X - 1` re-typed as a 25-digit decimal). Re-anchor each to the carried `rho_X` literal so they genuinely test the `epsBlk->0` substitution of `zetaReq`: change the targets from the hardcoded decimals to `rhoSuff - 1`, `rhoFail - 1`, `rhoMax - 1` respectively, e.g. `expectApprox["zeta_suff = rho_suff - 1", zetaSuff, rhoSuff - 1, 10^-20];` (and likewise for fail/max). A bug in `zetaReq` that did not vanish at `epsBlk=0` would now fail. Keep the load-bearing `unblocked zeta_req` (wl:44) and `d zeta_req exact formula` (wl:43) asserts untouched.

Do NOT modify the `unblocked zeta_req` asserts (py:55, wl:44) or the `d zeta_req exact formula` assert (wl:43). Do NOT change the literal values `3.46622291347846 / 3.46752913273870 / 3.46752922945601`.

**Self-test (already performed by auditor):**
- `rho_suff = 3.46622291347846 < rho_fail = 3.46752913273870 < rho_max = 3.46752922945601` → ordering asserts pass with the correct literals.
- `rho_max - rho_fail = 3.46752922945601 - 3.46752913273870 = 9.671731e-8`, which is `> 0` and `< 1e-6` → gap assert passes with the correct literals; it would FAIL if `rho_max` or `rho_fail` were mistyped to differ by more, so the check is genuinely load-bearing.
- These relations introduce no new constant and keep the three literals unchanged → no new paper_misalignment.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 087` and `redteam exec-mathematica 087` and confirms: (a) the new ordering/gap checks appear and pass, (b) the self-comparison `vs stage-086` lines are gone or reworded, (c) the `unblocked zeta_req` (and `d zeta_req exact formula`) asserts still print and pass, (d) both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl`
- summary: Replaced the self-referential ratio literal comparisons with structural ordering/gap checks and re-anchored Mathematica zeta checks to the carried rho literals.
- deviation: none

## F1-followup — finish the over-claim removal (top docstring)

The inline comment at py:69-72 was correctly reworded, but the SAME over-claim survives in the
SymPy TOP docstring at `scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:17-20`:
> "This script restates the unblocked one-ratio criterion `zeta_req(rho_alpha; 0) = rho_alpha - 1`
>  as a downstream-consistency probe and **cross-checks the Family-1 window literals against the
>  upstream stage-086 quoted values to catch renumber or transcription drift.**"

That now mis-describes the script (the checks are ordering/gap structural relations on the carried
literals, not an upstream literal cross-check). Reword lines 17-20 to match the applied fix, e.g.:
> "This script restates the unblocked one-ratio criterion `zeta_req(rho_alpha; 0) = rho_alpha - 1`
>  as a downstream-consistency probe and sanity-checks the Family-1 window literals (carried from
>  the Stage-086 notes) via their ordering and constructive-ceiling gap."

Docstring-only, non-functional — the committed transcript is unchanged. Do NOT touch any code,
assertion, or literal value. After editing, RUN `python3 <path>` to confirm it still exits 0 with
the identical transcript. Append an `## Applied: F1-followup` block (`files_changed`, `summary`,
`deviation`).

## Applied: F1-followup

- files_changed:
  - `scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py`
- summary: Reworded the SymPy top docstring to describe the carried Family-1 literals as ordering/gap sanity checks rather than an upstream literal cross-check.
- deviation: none
