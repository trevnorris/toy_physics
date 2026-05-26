---
batch_id: II.1-v2
label: Part II.1 — Overlap isotropy through continuum kernel (paper-grounded re-audit)
started: 2026-05-25
completed: 2026-05-26
stages_total: 13
stages_verified: 13
stages_blocked: 0
material_change_any: false
---

# Batch II.1 v2 — Final summary (paper-grounded re-audit)

Third batch processed under the v2 paper-grounded auditor (after I.1 and I.2). All 13 stages reached `verified`. **No `material_change` flagged for any stage** — all script-side edits were additions (new assertions, lane-collapse checks, paper-anchor blocks) or removals of tautological self-checks; no downstream-visible closed forms changed. The v2-introduced finding set was dominated by `insufficient_verification` (script verifies the right object but the assertion is thin) and `tautological_check` (assertion algebraically guaranteed by upstream definitions) — distinct from I.2's duplication pattern.

## Per-stage outcomes

| Stage | Audit findings | Resolution path | Verifier verdict | material_change |
|---|---|---|---|---|
| 024 | F1 mathematica_transliteration (Sections III/V port), F2 insufficient_verification (anchor Z/N rationals to 2x2 Maxwell inverse), F3 tautological_check (C_alpha self-equality + equal-lane substitutions), F4 insufficient_verification (explicit O(3)-collapse + lane-breaking witness) | fix_loop iter1 F3/F4 + F1/F2 blocked iter1 (sign convention) + Codex remediation pass for F1/F2 (`-R` off-diagonal derived from `paper/parts/part01_parent_geometry.tex:956` `+R_l A_l W_l` Lagrangian) | verified | false |
| 025 | F1 tautological_check (N0 self-substitution), F2 insufficient_verification (P=0 corollary not exercised) | fix_loop iter1 | verified | false |
| 026 | F1 insufficient_verification (K_req structural form anchor) | fix_loop iter1 | verified | false |
| 027 | 0 (clean) | n/a | verified | false |
| 028 | 0 (clean) | n/a | verified | false |
| 029 | F1 paper_misalignment (label drift Stage 12 → 029), F2 insufficient_verification (selected odd coefficient combined identity), F3 insufficient_verification (SymPy nullspace eigenvector cross-check), F4 paper_misalignment (alpha_crit verified in scripts but not in paper card) | Q1=a apply (relabel) + Q2=b apply (trim, destination stage 031) + fix_loop iter1 (F2, F3) | verified | false |
| 030 | F1 insufficient_verification (HF identity (v.e_-)^2 = -dλ_-/dα not directly exercised) | fix_loop iter1 | verified | false |
| 031 | F1 mathematica_transliteration (full PART I-II port), F2 tautological_check (generic quotient/HF scaffold) | fix_loop iter1 | verified | false |
| 032 | F1 insufficient_verification (eigensystem-based (v.e_-)^2 anchor), F2 mathematica_transliteration (Stage 15.4-15.5 port), F3 insufficient_verification (universal-identity scaffolding block) | fix_loop iter1 | verified | false |
| 033 | F1 tautological_check + redundant gate block | fix_loop iter1 | verified | false |
| 034 | 0 (clean) | n/a | verified | false |
| 035 | F1 paper_misalignment (target_mismatch — paper polynomial coefficients 206/138 vs scripts' correct 189/121) | Q3=a apply (paper + notes fix; scripts unchanged) | verified | false |
| 036 | F1 tautological_check (M_mix witness construction was arithmetic identity), F2 tautological_check (definitional self-consistency labelling) | fix_loop iter1 | verified | false |

## Findings breakdown

- `paper_misalignment`: 3 substantive (Q1=a, Q2=b, Q3=a across stages 029, 029, 035) — 2 of 3 trim or relabel; 1 paper-side coefficient fix
- `insufficient_verification`: 8 (024 F2/F4, 025 F2, 026 F1, 029 F2/F3, 030 F1, 032 F1/F3) — the dominant v2-introduced category
- `tautological_check`: 6 (024 F3, 025 F1, 031 F2, 033 F1, 036 F1/F2)
- `mathematica_transliteration`: 4 (024 F1, 031 F1, 032 F2, plus implicit at 024 Section IV which was fixed alongside)
- `stop_cold`: 0
- iter-2 / orchestrator hot-fix: 1 (stage 024 F1/F2 remediation pass — see below)

**V2-introduced findings missed by v1: ~18 substantive across 10 stages.** This is the heaviest v2 yield to date (I.1: ~10 substantive; I.2: ~10 substantive). The bulk is `insufficient_verification`: scripts that verify the right object but with thin or missing direct anchors. Pattern differs from I.2: I.2 was duplication-dominated; II.1 was script-thinness-dominated. Likely the first sweep wrote terse anchors and v2's paper-grounded reading exposed gaps that v1 (script-only) couldn't see.

## Resolution methodology

Same v2 pipeline as I.1/I.2:

1. **Audit wave**: 13 paper-grounded auditor agents (10 in parallel, then 3) read paper + notes + appendix BEFORE scripts. 3 paper_misalignment flagged (029 F1, 029 F4, 035 F1).

2. **Questions session**: `redteam/resolutions/batch_II1_paper_alignment.md` + `codex_prompt_batch_II1.md`. Codex answered all 3 — Q1=(a) cosmetic relabel, Q2=(b) trim with destination stage 031, Q3=(a) paper-side polynomial fix with full quotient-rule expansion.

3. **User review**: User approved all three directions in a single decision (no redirection this batch — Codex's first-pass recommendations all held up; cross-verification of Q2's destination claim by orchestrator independently confirmed stage 031 owns `alpha_crit`).

4. **Apply session**: `codex_apply_batch_II1.md` with per-question scope and destination-verification guardrail. Applied 3/3 cleanly. Notable: Q3 was the first batch-level paper-side edit since the v2 sweep began (I.1 had paper expansions per (a) direction; I.2 was all script-side trims).

5. **Output refresh**: stage 029 sympy + mathematica serial.

6. **fix_loop**: 9 stages with script-side findings (024, 025, 026, 029, 030, 031, 032, 033, 036). Processed sequentially per `fix_parallelism: 1` in two waves of 5+3 background bash. 8/9 iter1 clean; stage 024 F1/F2 blocked on sign convention (Codex correctly refused directive-prescribed `+rPair` matrix that contradicted paper's `+2 G_U G_W R` mixed term).

7. **Codex remediation pass for stage 024**: combined prompt with explicit math-authority delegation. Codex read `paper/parts/part01_parent_geometry.tex:956` (`+R_l A_l W_l` Lagrangian), derived `-R` off-diagonal for the conservative spring matrix, rewrote Section III/V independently, and instrumented Section IV with symbol-context resets and memoized sphere integrals (i4/i6) — see "Notable per-stage detail" below.

8. **Verifier wave**: 9 agents in parallel. All 9 returned `verified`, `material_change: false`. 0 needs_rework, 0 blocked.

## Notable per-stage detail

### Stage 024 — sign-convention remediation + Section IV performance fix

**Iter1 block**: directive prescribed `mPair = {{omegaU^2 - omega^2, +rPair}, {+rPair, omegaW^2 - omega^2}}` for the F1 matrix-inverse rewrite. Codex correctly detected that inverting this matrix gives `-2 g_U g_W R` in `g^T M^{-1} g`, which contradicts the paper's `Q_r = G_U² Ω_W² + 2 G_U G_W R + G_W² Ω_U²` (`paper/stages/stage_024.tex:108`). Blocked both F1 and F2 pending sign resolution.

**Remediation pass**: Codex read `paper/parts/part01_parent_geometry.tex:956`, found the upstream `+R_l A_l W_l` Lagrangian coupling, and derived: positive `R U W` Lagrangian → `+R` Hamiltonian off-diagonal → `−R` conservative spring matrix off-diagonal (since spring matrix sits on the LHS in frequency space). Independently verified: with `M = [[Ω_U² − ω², −R], [−R, Ω_W² − ω²]]`, `g^T M^{-1} g` at ω=0 gives `(g_U² Ω_W² + 2 g_U g_W R + g_W² Ω_U²)/Δ_r` ✓ matches paper's Q. Same chain confirms `(M^{-1} g)_W = (R g_U + Ω_U² g_W)/det` ✓ matches paper's P. Applied F1 and F2 with `-R` off-diagonal, both anchor assertions pass.

This is the second clean instance of "Codex as math authority" working as designed (first was stage 023 F2 in batch I.2). The pattern is now established: when the directive prescribes a sign or convention that contradicts the paper, Codex defers to the paper and rederives.

**Section IV performance hang** (orchestrator-detected, separate from the F1/F2 blocks): during the post-fix-loop exec wave, the `Table[tripleOverlap[basis[[i]], qMat, basis[[j]]], {i, 1, 5}, {j, 1, 5}]` at `mathematica/moving_throat_pde_stage024_..._mathematica_audit.wl:202` hung at >18 min CPU with no output past the Section IV banner. The original v1 script ran in <30s; symbol context contamination from the F3/F4 additions in earlier sections was the suspect. Remediation prompt instructed Codex to (i) add a `ClearAll[gU, gW, rPair, omegaU, omegaW, mPair, zFromMatrix, nFromMatrix, qRef, hRef, pRef, deltaRef, sRef, zRefRational, nRefRational, ...]` reset at the top of Section IV, and (ii) memoize i4/i6 sphere integrals via `i6[i_, j_, k_, l_, m_, nn_] := i6[i,j,k,l,m,nn] = Integrate[...]`. After applying both, full Mathematica runtime dropped from >1080s (killed) to 25.05s. This is the first confirmed instance of **symbol-leakage induced hang in a long Mathematica script**; documented for the pitfalls list (potential pitfall #6).

### Stage 029 — combined paper_misalignment + script-side finding

Stage 029 was the only stage in this batch with both paper_misalignment items (F1, F4) and script-side items (F2, F3). The paper_misalignment items were applied via the Q1/Q2 Codex apply session (Q1=a relabel; Q2=b trim with destination stage 031 verified). The script-side items were then applied via fix_loop iter1. The directive `redteam/directives/stage_029.md` carries an annotated `## Resolve before fix_loop` block (now marked RESOLVED) plus four `## Applied: F<n>` blocks. Verifier confirmed all four findings legitimately resolved, alpha_crit absent from stage 029 audits, present in stage 031 audits. Stage 031 owns the refined threshold from this batch onward.

### Stage 035 — first batch-level paper-only edit since v2 began

Q3 was a target_mismatch: paper boxed `dF/dξ`'s bracket polynomial coefficients as `206 δ² ξ + 138 ξ³`; both engines (and Codex's independent quotient-rule expansion) gave `189 δ² ξ + 121 ξ³`. Numerical sanity at δ=0 confirmed scripts are correct (`121/[81(1-ξ)²]` matches scripts; paper form gives wrong `138/[81(1-ξ)²]`). Applied as paper-side edit to `paper/stages/stage_035.tex:71` and `notes/stages/moving_throat_pde_stage035_...md:86`; scripts unchanged. Stage 035 went straight to `verified` without needing a Codex fix_loop pass.

## Cascade observations

**No material_change downstream cascades to flag.** All 13 stages have `material_change: false` (no derived closed forms changed; only assertion additions and tautological removals). This is in contrast to I.1 v2 (1 material_change: stage 004) and I.2 v2 (4 material_change: stages 013, 014, 015, 018). With II.1 v2 closed at depth, the entire range 004–023 + 024–036 is paper-aligned at v2 depth.

## Process observations / pitfalls

1. **Codex-as-math-authority pattern is stable.** Two instances now (023 F2 in I.2; 024 F1/F2 in II.1) where the directive's prescribed sign convention contradicted the paper. Codex correctly deferred and read the upstream Lagrangian. Pattern should be retained for III.1+ batches.

2. **Symbol-leakage induced Mathematica hang (NEW pitfall #6 candidate).** Stage 024 Section IV would have been unaudible without the `ClearAll[...]` reset added by Codex. The trigger: large additions to earlier sections that introduce many global symbols, then a later section using FullSimplify with `$Assumptions = True`. Specifically: a `Table[...sphere integral...]` over 5×5 entries with a 6-fold sum inside FullSimplify will explore astronomical simplification paths if symbol context is contaminated. Pre-emptive defense: any time a directive adds many new global symbols to a Mathematica section before a heavy FullSimplify/Series block, also prescribe a `ClearAll[<symbol-list>]` reset.

3. **Cross-verification of Codex destination claims worked first-try.** For Q2 (029's alpha_crit destination), orchestrator independently grepped stage 031's paper card and scripts before user gate. Confirmed `paper/stages/stage_031.tex:43,65` and `scripts/moving_throat_pde_stage031_..._sympy_audit.py:87,94,116` + `mathematica/moving_throat_pde_stage031_..._mathematica_audit.wl:59,61,71`. Pattern from I.2 transferred cleanly.

4. **User-redirection rate dropped to 0 this batch.** I.2 had load-bearing user redirection (Q1/Q2 flipped from (c) to (b)). II.1 had no redirection — Codex's first-pass recommendations all held up. Difference: II.1's paper_misalignment items were less ambiguous (label drift, missing-from-card with clear downstream destination, polynomial coefficient with independent derivation). Going forward batches with similar character (script-thinness rather than duplication) should also have low redirection rates.

5. **All 13 verifiers flagged stale canonical `scripts/output/*.txt` / `mathematica/output/*.txt` as non-blocking observation.** The fix-loop run via codex updates scripts but not the canonical .txt mirrors; `$RT exec-mathematica` writes to `redteam/exec_logs/` not to canonical paths. Orchestrator added a final canonical-refresh wave at end of batch. Could be incorporated into the SKILL.md pipeline as an automatic step after fix_loop, before verifier wave.

## Cosmetic follow-ups (non-blocking)

- **Stage 024**: Section II.1 subbanner is now empty after F3 deletions. Header with no body. Could remove the empty subbanner in a cosmetic cleanup pass.
- **Stage 026**: Mathematica final banner still reads `Stage 9 Mathematica audit passed.` (legacy renumbering carryover). Same template-drift class as stage 022/028's legacy labels.
- **Stage 032**: `expectZero` label at `mathematica/.../stage032_..._mathematica_audit.wl:172` is malformed (SymPy + Mathematica labels concatenated with `;`). Residual passes but label is ugly.

None of these block downstream consumption. Could be batched into a cosmetic-cleanup pass at the end of the full v2 sweep.

## Stage 023 follow-up (from stage 024 verifier note)

Stage 024 verifier raised a side-observation: "Stage 023 Mblock sign question that left stage 023 `blocked_unfixable` is the same matrix structure resolved here for stage 024." This is **incorrect** — stage 023 was verified in I.2 v2 close (commit f16a6f4) via the Codex remediation pass that derived the sign from `paper/stages/stage_023.tex` `eq:app-stage023-full-lagrangian` (`+R U W` → `−R` off-diagonal). Stage 023 is currently `verified`, not blocked. The verifier likely confused the older v1 stage 023 directive with current state. No action needed; included here for record.

## Coverage status post-II.1 v2

| Range | Status |
|---|---|
| 001–036 | verified at v2 depth (I.1 + I.2 + II.1 v2 sweeps complete) |
| 037–084 | verified at v1 depth (III.1, III.2, III.3, III.4) |
| 085–253 | pending |

Cumulative: 97 of 253 stages verified (84 stages at the v1 boundary; 84+13 = 97 with II.1 v2 batch added at v2 depth — but the count is the same 84 stages, just at deeper v2 verification for batches I.1, I.2, II.1).

Wait — count check: I.1 (12) + I.2 (11) + II.1 (13) + III.1 (12) + III.2 (12) + III.3 (12) + III.4 (12) = 84. Of these, batches I.1, I.2, II.1 (36 stages) are at v2 depth; III.1-III.4 (48 stages) at v1 only. Cumulative stages verified: 84 (unchanged from previous batches — same stages, deeper verification for the 36 v2-completed).

## Next batch

III.1 (stages 037–048, 12 stages, "Continuum kernel, generalized branch, rank-2") — first non-EM-projected, non-Maxwell-bridge batch under v2. Expected character: less duplication (post-I.2's trim spree), possibly more script-thinness (post-II.1's pattern). User-approval gate required between batches per `[[feedback-sequential-audit-chunks]]`.
