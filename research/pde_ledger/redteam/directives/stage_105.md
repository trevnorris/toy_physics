---
unit_id: 105
batch: IV.2
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 105

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:26-66`

**Issue:**

The `.wl` script is a line-by-line port of the `.py` script: the same intermediate variables (`omegaQ`, `sigmaCan`, `yRet`, `chiQ`, `lamDef`, `yDef`) are introduced in the same order, the same `Series`/`Coefficient` extraction is used at the same order, and the same closed-form right-hand sides are compared. The deformed-branch normalization `Y_def = -3/Lam_def` used by both engines is not stated in the notes — it is a script-side choice replicated on both sides. The second-engine policy requires Mathematica to derive `chi_Q = 1` from an independent algebraic path so that a transcription error in either the module form or the closed-form RHS cannot pass silently on both engines.

**Required change:**

Rewrite the `.wl` audit (current lines 26-66) so that the verification of `chi_Q = 1` follows a different algebraic path than the `.py`. One acceptable approach (the residue path):

1. Keep the banner line (after F2 lands) and the `$Assumptions` block (`.wl:28-31`).
2. Define `omegaPole = 3 cSound / (2 aThroat)` as the (squared-)pole location of the retarded module.
3. Build the retarded module in factored form `yRetFactored = (3 (omegaPole^2 - omega^2) + 4 I chiQ sigmaQ omegaPole^2 omega^5 + (omegaPole^2)) / (4 (omegaPole^2 - omega^2 - I chiQ sigmaQ omegaPole^2 omega^5))` (or any algebraically distinct presentation that does not reuse the `3/4 + (1/4)/(...)` decomposition from the `.py`). State `sigmaQ = (9/8)/omegaPole^5` and assert `FullSimplify[sigmaQ - 4 aThroat^5/(27 cSound^5)] === 0` to re-confirm D1.
4. Replace the series-coefficient choreography with a residue/limit-based extraction. Concretely: compute `imagPart5 = Coefficient[ComplexExpand[Im[yRetFactored /. omega -> omega]], omega, 5]` *only after* an algebraically distinct simplification path (e.g. `Series[yRetFactored, {omega, 0, 5}] // Normal // ComplexExpand // Im // Expand`), then assert `FullSimplify[imagPart5 - chiQ aThroat^5/(27 cSound^5)] === 0`. The point is to avoid the `Coefficient[ySeries, omega, k]` extraction with the exact same intermediate `ySeries` variable.
5. For \(\chi_Q = 1\): instead of `Solve` over a linear equation, evaluate the imaginary-part residual at the outgoing-branch target `aThroat^5/(27 cSound^5)` and `Reduce` (or `FindInstance`) for `chiQ` over the reals, then assert the unique solution equals 1.
6. For the deformed branch, derive `yDefDirect` from the outgoing operator without going through `-3/lamDef`. Acceptable: invert `Λ_2^def · Ŷ_2^def = -3` symbolically using `SolveAlways` or `Reduce` over the series coefficients up to order 5, recovering `Ŷ_2^def = 1 + z^2/9 + 4 z^4/81 + I xiQ z^5/27`, then assert that match.

The bottom-line assertion `(chiQ /. chiSol) - 1 == 0` (or its `Reduce`-based equivalent) must remain. Variable names must differ from the `.py` (no `yRet`, `lamDef`, `yDef`, `omegaQ`, `sigmaCan`).

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 105` and confirm:
- the script exits 0;
- the substrings `(1/4)/(1 - omega^2/omegaQ^2`, `lamDef`, and `yDef` no longer appear in the `.wl`;
- the final transcript still asserts `chiQ - 1 == 0` (or equivalent) and prints `Stage 105 Mathematica audit passed.` (the latter conditional on F2 being applied first).

## F2 — paper_misalignment

**Subtype:** target_mismatch

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_105.tex:2` quote: `\label{stage:105}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md:1` quote: `# Moving-Throat PDE — Stage 105: Exact Fixing of chi_Q`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:3` quote: `Stage 88 SymPy audit.`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:28` quote: `banner("STAGE 88 — EXACT FIXING OF chi_Q")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:26` quote: `banner["STAGE 088 — EXACT FIXING OF chi_Q"];`

## Resolve before fix_loop

The script's docstring and banner labels declare "Stage 88" / "Stage 088" while the file path, notes, paper card label (`\label{stage:105}`), and the `.wl`'s own closing print line (`"Stage 105 Mathematica audit passed."`) all identify this as Stage 105. The paper card's display section also names it "Stage~122" (the public display number from the global stage→display map). Which label should the script's banner/docstring carry?

Possible directions (the user picks one):
- (a) Script labels should read "Stage 105" (matching the file path, notes, and label-id) → Codex updates `.py:3` to `"""Stage 105 SymPy audit."""`, `.py:28` to `banner("STAGE 105 — EXACT FIXING OF chi_Q")`, `.wl:26` to `banner["STAGE 105 — EXACT FIXING OF chi_Q"];` — no paper change.
- (b) Script labels should read "Stage 122" (matching the displayed section title) → Codex updates the same three lines using `Stage 122` instead — no paper change.
- (c) The stale "Stage 88" label is intentional historical provenance and should remain → Codex leaves the docstring/banner untouched; no script change.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
