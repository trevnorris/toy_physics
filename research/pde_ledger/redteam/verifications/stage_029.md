---
unit_id: 029
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:35:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 029 (batch II.1 v2)

This verification supersedes the prior batch II.1 v1 verification (2026-05-22) on the same unit. The findings inventory in this round (F1-F4) is a new set raised by the re-audit; it is not a continuation of the prior round's F1-F3.

## Per-finding outcomes

### F1 — paper_misalignment (legacy stage-number drift in docstrings/banner)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:3` — docstring filename header now reads `moving_throat_pde_stage029_dynamic_loading_sympy_audit.py` (was `..._stage12_...`).
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:5` — `SymPy audit for Stage 029 of the moving-throat PDE program.` (was `Stage 12`).
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:33` — `banner["STAGE 029 — DYNAMIC LOADING"];` (was `STAGE 012`).
- Mathematica exec_log line 8 confirms `STAGE 029 — DYNAMIC LOADING` at runtime.

**Assessment:**
The relabel matches the directive exactly. The retained `Stage-11 loading parameter` mention at sympy docstring line 21 is a physical-content reference (not a file-identification label), correctly preserved per directive. The grep gate `Stage 12|STAGE 012|stage12_|stage_12` returns no matches in either script. Resolved during the batch II.1 v2 paper_alignment Codex apply session (Q1).

### F2 — insufficient_verification (selected odd coefficient not asserted as combined identity)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:291-303` — added a new `expect_zero("delta D_-^odd (script) - delta D_-^odd (paper formula)", delta_D_script - delta_D_paper)` block. `delta_D_paper` is built from the paper's explicit `(Omega_U**2 * lambda_W + lambda_R * lambda_U)**2 / (Omega_U**2 * Omega_W**2 - lambda_R**2 * sigma)**2` denominator with `Gamma_port * kappa_sel_sq * omega**5`; `delta_D_script` is `-sp.I * beta5 * kappa_sel_sq * omega**5` with the script-derived `beta5`.
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:222-232` — added the mirror `expectZero`, with `deltaDPaper` (paper closed form using `delta0^2`) versus `deltaDScript` (`-I*beta5*kappaSelSq*omega^5`).

**Assessment:**
The assertion binds the paper's eq `selected-odd` to the script-computed combination. Non-tautological: in sympy, `beta5` is derived as `beta_clean.subs(omega, 0) * Gamma_port` from the upstream Schur-complement chain (line 216), while `delta_D_paper`'s denominator is independently constructed from `Omega_U**2 * Omega_W**2 - lambda_R**2 * sigma`; equality collapses to zero only because `beta_clean(omega=0)` actually equals the paper closed form — a genuine algebraic check that exercises the M7 identity inside sympy too. The Mathematica twin similarly compares `beta5` (built from `betaClean /. omega -> 0`) against `delta0^2`-denominator form.

Sympy exec_log line 916: `delta D_-^odd (script) - delta D_-^odd (paper formula) = 0`.
Mathematica exec_log lines 54-55: `delta D_-^odd (script) - delta D_-^odd (paper formula) = 0` followed by `PASS:`.

Codex's reported deviation ("Parenthesized the Mathematica multi-line deltaDPaper RHS to preserve the requested formula safely in .wl script form") is purely syntactic — wrapping a multi-line expression in `(...)` so the line-continuation parses in Wolfram Language. Mathematical content unchanged.

### F3 — insufficient_verification (sympy lacks direct eigenvector projection cross-check)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:269-285` — inserted a SymPy direct nullspace eigenvector projection block before the existing weak/strong limit checks. Constructs `K_eff_al` symbolically, takes `(K_eff_al - lam_minus_template * I).nullspace()`, projects `v` onto the resulting eigenvector, and asserts `expect_zero("kappa_sel^2 closed-form vs eigenvector projection", sp.simplify(kappa_sel_template - kappa_sel_sq_direct_template))`.

**Assessment:**
This mirrors Mathematica's existing M8 cross-check in sympy. Non-tautological: `kappa_sel_template` (HF derivative `-d lambda_- / d alpha`) and `kappa_sel_sq_direct_template` (`(v.vec_lo)^2 / |vec_lo|^2`) are computed via genuinely independent symbolic paths, and the equality holds by Hellmann-Feynman for the specific `K_eff_al = M_0 - al v v^T` structure. The directive's `nullspace()` route is preserved verbatim, sidestepping `eigenvects()`'s occasional `CRootOf` issues. Sympy exec_log line 913: `kappa_sel^2 closed-form vs eigenvector projection = 0`. The existing weak/strong limit checks remain in place at lines 288-289 (log lines 914-915), so the new assertion adds coverage without removing prior anchors. Codex declared `deviation: none`; confirmed.

### F4 — paper_misalignment (alpha_crit verified in scripts, not in paper card)

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:17` — removed the `• the refined softening threshold alpha_crit,` bullet from the docstring Scope section.
- `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:189-194` (original line range) — removed the `alpha_crit` definition + `expect_zero("alpha_crit - expected", ...)` block from `conservative_profile_selection`.
- `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:159-166` (original line range) — removed the `al = Symbol["alphaLoad"]; detTemplate = ...; alphaCrit = ...; expectZero["det(alpha_crit)", ...]` block.

**Assessment:**
User approved direction (b) — trim `alpha_crit` from stage 029 since stage 031 owns the refined threshold. Independently verified: `alpha_crit`/`alphaCrit` lives natively at `scripts/moving_throat_pde_stage031_selected_branch_reachability_sympy_audit.py:74,77,81,89,98-99` and `mathematica/moving_throat_pde_stage031_selected_branch_reachability_mathematica_audit.wl:97,102,106-108,110-111,120,128,142-143`, including the `expectZero` PASS bindings. Grep for `alpha_crit|alphaCrit|alphaLoad` in stage 029 returns only one hit at sympy:260, which is the local Hellmann-Feynman parameter `al = sp.symbols('alpha_load', real=True)` used by the F3 nullspace block — a different, correctly retained variable from the removed solver. Stage 029 exec logs contain no `alpha_crit` line and both scripts exit 0. Resolved during the batch II.1 v2 paper_alignment Codex apply session (Q2).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Line 916: `delta D_-^odd (script) - delta D_-^odd (paper formula) = 0` (new F2 anchor)
- Line 913: `kappa_sel^2 closed-form vs eigenvector projection = 0` (new F3 cross-check)
- Lines 914-915: pre-existing weak/strong-limit `kappa_sel` checks still PASS
- Line 825: `extracted beta_5 - expected beta_5 = 0` (upstream chain intact)
- Line 547-548: `beta - clean transfer factor = 0` and the `alpha - (alpha_cons + beta*Pi) at O(Pi) = 0` (D3 anchor still in place)
- Line 364: `DeltaK_tilde - DeltaK_bare (Xi_0 cancellation) = 0` (typo guard from prior batch unchanged)

**Mathematica:** exit=0. Notable lines:
- Line 8 (banner): `STAGE 029 — DYNAMIC LOADING` (F1 fix observed at runtime)
- Lines 14-15: `Sigma_seq - (Xi I + alpha vv^T) = {{0, 0}, {0, 0}}` / `PASS`
- Lines 47-48: `kappa_sel^2 closed-form vs eigenvector projection = 0` / `PASS` (pre-existing M8 still holds)
- Lines 54-55: `delta D_-^odd (script) - delta D_-^odd (paper formula) = 0` / `PASS` (new F2 anchor)
- Line 57: `Stage 029 Mathematica audit passed.`

**Output freshness:** The saved `.txt` outputs (`scripts/output/moving_throat_pde_stage029_dynamic_loading_sympy_audit.txt` mtime 2026-05-25 23:24:03 and the mathematica twin mtime 23:24:57) are OLDER than the script edits (sympy 23:38:01 and wl 23:38:22) and do NOT contain the new F2/F3 assertion lines (a grep confirms 0 occurrences). The fresh runs are captured in `redteam/exec_logs/stage_029_{sympy,mathematica}.log` (timestamps 2026-05-26 00:16/00:17), both showing the new lines and `exit_code: 0`. Per this verifier prompt's scope (exec logs are the canonical source for assertion verification), all assertions pass; the `.txt` mirrors should be regenerated before any downstream consumer treats them as canonical. Flagged as a side observation, not a blocker.

## Material-change assessment

`material_change`: false.

Reasoning:
- F1 is label-only; no mathematics affected.
- F2 adds a new asserted equality; the equality holds because the underlying objects (`beta5`, `kappa_sel_sq`, `Delta_0`) were already correct upstream. No derived value changes.
- F3 adds a new asserted equality for `kappa_sel^2` in sympy; no value or formula is altered, only an independent verification path is now in place in sympy too (matching the existing Mathematica M8).
- F4 removes an extra assertion that is now owned by stage 031. The remaining stage 029 surface area is exactly the paper-card-stated D1-D4 deliverables; downstream consumers of stage 029 do not import `alpha_crit` from this script.

No downstream unit's input shape changes. The orchestrator's automatic `upstream_stale: true` marker for units > 029 is conservative and can remain; no specific re-audit narrowing is needed beyond confirming stage 031 still owns `alpha_crit` (verified above).

## Side observations (non-blocking)

1. Saved `.txt` outputs under `scripts/output/` and `mathematica/output/` predate the script edits and do not yet contain the F2/F3 PASS lines (they also still contain the F4-removed `alpha_crit` lines). The exec_logs do reflect the fresh state. Re-running `$RT exec-sympy 029` and `$RT exec-mathematica 029` after this verification (subject to the user's Mathematica single-seat rule) would re-canonicalize the saved transcripts.
2. The retained `al = sp.symbols('alpha_load', real=True)` at sympy:260 is the Hellmann-Feynman loading parameter local to `selected_mode_projection` (different role from the removed `alpha_crit` solver). It is correctly kept and used by the new F3 nullspace block.
3. Codex's F2 Mathematica deviation note (parenthesizing the multi-line RHS) is purely syntactic; mathematical content unchanged.
4. The prior batch II.1 v1 verification at this path (verifier 2026-05-22) addressed a different findings set (tautological_check + mathematica_transliteration). This new batch's findings target subsequent issues raised by the re-audit and are independent of those resolutions.

## Verdict justification

All four findings are legitimately resolved. F1 (docstring/banner relabel Stage 12 → 029) and F4 (alpha_crit trim with stage 031 destination) were applied during the batch II.1 v2 paper_alignment Codex apply session per the user's direction; both edits are conservative, surgical, and externally verified against the destination (stage 031) and the runtime banner (mathematica log line 8). F2 (sympy + Mathematica selected-odd combined-identity assertion) and F3 (sympy direct nullspace eigenvector projection) were applied via fix_loop iter1; diffs match the directive verbatim except for a benign Mathematica F2 parenthesization. Both engines run clean to exit 0, with all new assertions emitting `= 0` / `PASS`. No regressions; no tautological assertions; no collateral edits beyond the directive scope. The only operational note is that the saved `.txt` mirrors are stale relative to the script edits — flagged as a side observation, not a verification blocker since exec_logs are the canonical evidence source.
