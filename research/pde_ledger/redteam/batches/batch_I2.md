---
batch_id: I.2
label: Part I.2 — Maxwell bridge, parent throat action, reduced one-port
started: 2026-05-21
completed: 2026-05-21
stages_total: 11
stages_verified: 11
stages_blocked: 0
material_change_any: false
---

# Batch I.2 — Final summary

All 11 stages reached `verified` status. Both engines (SymPy + Mathematica) exit 0 on every unit; all findings closed substantively. No `material_change` flags, so no downstream cascade — batch II.1 may proceed without re-auditing I.2.

## Per-stage outcomes

| Stage | Audit findings | Codex iterations | Verifier verdict | Notes |
|---|---|---|---|---|
| 013 | 3 (missing_math, taut×1, hardcoded) | 1 | verified | M3 primitive form deviated from directive; codex's choice was correct and produced matching M5 deltaP4 coefficient |
| 014 | 4 (taut×3, missing_math) | 1 | verified | replaced `A+(-A)=0` tautology; new `Z2_slot` linkage asserts |
| 015 | 3 (missing_math, taut, hardcoded) | 1 | verified | removed `if m == 0: return sp.Integer(1)` Gaunt short-circuit; closed Gaussian overlap evaluations |
| 016 | 1 (missing_math, 11 claims) | 1 | verified | dependent-symbol `R[t,w,u,v]`, symbolic IBP product-rule check |
| 017 | 3 (taut×2, missing_math) | 1 | verified | replaced wall-only specialization tautology with 6 lane cross-checks against generic K1/H_even |
| 018 | 4 (taut×2, missing_math, insuff_verif) | 1 | verified | new `Xi1_from_expected` block via closed-form `expected_dK/expected_dM` |
| 019 | 1 (missing_math, 12 claims) | 1 | verified | codex deviated from directive's ansatz; deviation reproduced SymPy's P2/P4 so M4-M8 still pass |
| 020 | 2 (missing_math, taut) | 1 | verified | own `GauntIntegral` from `ThreeJSymbol`; m=0 lane now exercises Wigner machinery |
| 021 | 2 (insuff_verif, math_transliteration) | 2 (iter1 EL bug) | verified | iter1 manual EL pattern `D[expr, qFun[t]]` returned 0 for `qFun[t]`-containing terms; iter2 switched to `EulerEquations` from `VariationalMethods` AND incidentally caught the same `lRed` multi-line continuation defect seen at I.1 stage 003 |
| 022 | 2 (math_transliteration, taut) | 1 | verified | Sections I/II/IV switched from `Series` extraction to `Solve[coeffEqs, ...]`; `SphericalHankelH1[2, z]` replaced `j2+i*y2`; round-trip tautology removed |
| 023 | 4 (taut×3, math_transliteration) | 1 | verified | additivity check on `grouped_parts`; closed-form `N2/N4_target_closed` substitutions; numerical-substitution and direct-Bessel-expansion cross-checks |

## Process bugs discovered and fixed

1. **`$RT exec-sympy` / `$RT exec-mathematica` race condition on MANIFEST.yaml:** When the orchestrator ran 11 parallel `$RT exec-sympy` calls to refresh saved outputs, the atomic-write-temp-then-rename pattern collided. Symptoms: `FileNotFoundError: MANIFEST.yaml.tmp -> MANIFEST.yaml` and `yaml.scanner.ScannerError` on a corrupted final key (stray `-05-20T01:00:00-06:00'` line). Fixed manifest by hand. Memory saved at `feedback_no_parallel_exec_sympy.md`. Workaround: refresh outputs via direct `python3 script.py > output.txt` (no manifest write) in parallel; reserve `$RT exec-*` for stages where exec metadata matters and run sequentially.
2. **Recurring `lRed = ...` multi-line continuation defect:** Now caught twice — I.1 stage 003 and I.2 stage 021. The Mathematica parser truncates a multi-line `lRed = a + b\nc + d` to just the first line; codex's fix is to parenthesise the RHS. Task `P4-32` opened in `PAPER_CLEANUP_TRACKER.md` to sweep all stage `.wl` files for the same class of defect.
3. **Codex EL pattern `D[expr, qFun[t]]` returns 0:** When variables are bound to head[arg] forms (e.g. `q = qFun[t]`), Mathematica's `D` treats them as not-a-variable and returns 0 instead of the expected partial derivative. Codex's manual EL operator pattern `D[D[lRed, D[q,t]], t] - D[lRed, q]` therefore silently produces a wrong result. Fix: use `EulerEquations[L, {f1[t], f2[t], ...}, t]` from the `VariationalMethods` package. Recorded in stage 021's directive iter 2 and in the `MATHEMATICA_MIRROR_POLICY` Stage 021 entry.

## Prompt hardening landed this batch

None this batch — the auditor self-test step added at the end of I.1 appears to have caught most directive-level bugs preemptively. Codex caught two engine-level bugs (stage 021 EL pattern + lRed continuation) that no auditor would have spotted without running scripts.

## Auditor error rate observed

0 of 11 directives had math errors that codex caught via the Block-and-delta loop. (Compare I.1: 3 of 12.) The auditor self-test prompt-hardening from the end of I.1 is doing real work.

1 of 11 codex iter-1 outputs regressed due to a Mathematica-engine quirk (stage 021's `D[..., qFun[t]]` bug) — not an auditor error, not a codex misread of the directive, but a Mathematica behaviour codex didn't anticipate. The verifier wrote a clean delta directive recommending `EulerEquations` and codex iter 2 landed it (plus the bonus `lRed` parenthesization catch).

## Tracker updates landed in this commit

- `notes/MATHEMATICA_MIRROR_POLICY.md`: 11 new Independent-Mirror Set entries (013-023). The existing stage 022 entry was extended with the batch-I.2 transliteration-rewrite notes.
- `notes/CHECKPOINT_TRUST_AUDIT.md`: 022 and 023 rows extended; stage 003 row cross-references the recurrence of the `lRed` defect class at stage 021.
- `notes/CHECKPOINT_CONSTANT_PROVENANCE.md`: no edits (no new constants surfaced; 021 is not a checkpoint).
- `notes/PAPER_CLEANUP_TRACKER.md`: new P4-31 (I.2 batch), P4-32 (`lRed` defect class sweep), P3-04 (paper-side propagation); change-log entry dated 2026-05-21.
- `notes/EM_PROJECTED_INTEGRATION_TRACKER.md`: verification-policy bullet updated to reflect Stages 013-020 now have mirrors and Stage 021 was substantially rewritten; new completed-checks bullet for batch I.2.
