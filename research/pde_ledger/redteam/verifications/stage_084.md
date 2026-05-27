---
unit_id: 084
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 084 (v2)

## Per-finding outcomes

The auditor's v2 report on stage 084 returned `verdict: clean` with `findings_count: 0`. There were no formally-filed findings to verify. The auditor noted in the "Cosmetic note" paragraph (report line 90) and again in the "What the script claims to verify" paragraph (report line 44) that the script banner string at line 32 still read `"STAGE 067 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON"` despite the unit being renumbered to 084 during the great-reorder commit `0d09ef6`. The auditor explicitly declined to file this as a math finding ("cosmetic only, not a math finding"; "deferred to a renumbering sweep per the doc-alignment exclusion"). The orchestrator nonetheless promoted this observation to a `paper_misalignment` / `notes_contradicts_script` style relabel and applied it directly.

Treating the orchestrator-applied relabel as F1 (banner string `STAGE 067` -> `STAGE 084`), classification follows.

### F1 — paper_misalignment (banner relabel, orchestrator-applied)

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:32` — the literal `banner["STAGE 067 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON"];` was replaced by `banner["STAGE 084 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON"];`. This is the sole hunk in `redteam/exec_logs/stage_084_diff.patch`: 1 `-` line, 1 `+` line, 3 lines of unchanged context above and below. No other file was touched.

**Assessment:**
The relabel matches the orchestrator's stated change verbatim. No collateral edit: the diff is precisely 1 line modified, with the rest of the banner string (the en-dash, the descriptive title) preserved character-for-character. The output file `mathematica/output/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.txt:3` now also reads `STAGE 084 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON`, confirming a re-run after the edit. No assertion was touched, so there is no tautology risk to evaluate. Script algebra (the five substantive checks A1-A7 from v2: inverse demand map, cross-route `Xi_F1`, Family-1 `Pe -> oo` limit, four ordering inequalities) is byte-identical to the v2 audit's reviewed state.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script exists for unit 084 (status-only carve-out: `is_status_only_candidate: true`, `is_checkpoint: false`). The orchestrator deposited no `stage_084_sympy.log` — correct, since there is no SymPy artifact to run.

**Mathematica:** exit=0 (inferred from the closing `Stage 084 Mathematica audit passed.` line and the complete output through the implicit `Exit[0]`). The orchestrator did not deposit a `stage_084_mathematica.log` under `redteam/exec_logs/`, but the canonical script-output file is fresh and complete. Notable lines from `mathematica/output/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.txt`:

- Line 3: `STAGE 084 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON` (relabel propagated into output, confirming re-run).
- Lines 10, 16, 21, 28, 30, 32, 34: seven `PASS:` lines — same set and count as the v2 audit reviewed (`inverse demand map`; `Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta)`; `zeta_phys at Family-1 (Pe->oo limit) matches upstream zeta_max^(F1)` with `diff = 2.23e-15`; chi-window ordering; hard-ceiling gap; J-window ordering; fail-side J vs chi).
- Line 18: `Limit::alimv: Warning: Assumptions that involve the limit variable are ignored.` — same benign warning the v2 auditor investigated and confirmed harmless; the limit value on line 19 still prints to 15-digit agreement with upstream `zeta_max^(F1)`.
- Line 36: `Stage 084 Mathematica audit passed.` followed by the implicit `Exit[0]`.

**Output freshness:** confirmed. Script mtime is `2026-05-27 02:05`; output mtime is `2026-05-27 02:19` — output is 14 minutes newer than the script, demonstrating the orchestrator re-ran Mathematica after applying the banner relabel.

## Material-change assessment

`material_change`: false.

The applied diff is a string literal inside a `banner[]` print call (a cosmetic header) on line 32 of the `.wl`. It does not change any symbolic definition, any assertion, any numeric value, or any control-flow path. The output PASS set is unchanged (same 7 assertions, all PASS, same residuals to 15 digits). No downstream unit consumes the banner string. Any orchestrator-policy `upstream_stale: true` flag for units > 084 is not warranted by anything in this diff and may be cleared without re-audit.

## Side observations (non-blocking)

- The orchestrator did not write a dedicated `redteam/exec_logs/stage_084_mathematica.log`; only the canonical Mathematica `output/*.txt` file is present. The verifier prompt allows for this. For consistency with the per-stage log convention used elsewhere, the orchestrator may wish to mirror the canonical output into `exec_logs/` going forward. Non-blocking.
- This verification file overwrote a v1 verification (`verify_date: 2026-05-25`, `findings_resolved: 3`) covering the original three-finding directive. The v1 verification's content is now superseded; if the orchestrator wants to retain v1 history, consider archiving prior verification snapshots before overwrite. Non-blocking.
- The banner string is a print-only side effect; a future renumbering sweep could derive the stage number from the filename or a top-of-script constant to prevent recurrence. Out of scope for this verification.

## Verdict justification

`verified`. The orchestrator applied the banner relabel (`STAGE 067` -> `STAGE 084`) exactly as described, the diff contains nothing beyond that single string change, the post-fix Mathematica re-run produced an output file with the new banner and the same seven PASS lines from the v2 audit, and the script exits cleanly. No assertions, definitions, or numeric values were touched; `material_change` is false. The auditor's underlying v2 `verdict: clean` carries forward unchanged.
