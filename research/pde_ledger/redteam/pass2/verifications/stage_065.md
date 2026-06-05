---
unit_id: 065
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 065

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
No script-logic edit (none required). The committed `.txt` transcripts were re-generated post-fix. The refreshed SymPy transcript (`exec_logs/stage_065_sympy.log`, dated 2026-06-05T13:52:31) now emits `STAGE 65 — THIN-WALL CONFINEMENT BRANCH` (line 8) and `STAGE 65 AUDIT PASSED` (line 53). The refreshed Mathematica transcript (`exec_logs/stage_065_mathematica.log`, dated 2026-06-05T13:54:24) now emits `STAGE 065 — THIN-WALL CONFINEMENT BRANCH` (line 21) and `Stage 065 Mathematica audit passed.` (line 54).

**Assessment:**
Correct. The stale `STAGE 48 / STAGE 048` self-labels in the prior committed transcripts are gone, replaced by the canonical `STAGE 65` (SymPy) / `STAGE 065` (Mathematica) banners that the banner code already emitted. All math content is unchanged: `g_phi = V0/ell`, the `(1,2,1)` `I1` structure, `J1_num = sqrt(2)*sqrt(pi)/2 = Sqrt[Pi/2]`, `J3_num = 3*sqrt(2)*sqrt(pi)/8 = (3*Sqrt[Pi/2])/4`, and every residual `= 0` / `PASS`. No collateral content change. The directive scoped F1 to the orchestrator-side refresh (no Codex script edit), which is what happened.

### F2 — numbering self-label (stale docstring/comment labels)

**Classification:** resolved

**What changed:**
The diff (`exec_logs/stage_065_diff.patch`) shows exactly one changed line in `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py`: the docstring title line, `-Moving-Throat PDE — Stage 48 SymPy audit` → `+Moving-Throat PDE — Stage 65 SymPy audit`. Confirmed against the live file: line 3 now reads `Moving-Throat PDE — Stage 65 SymPy audit`.

**Assessment:**
Correct and strip-the-number identical: only `48`→`65`, no format change, no padding, no other line touched. The directive deliberately narrowed F2 to the single SELF-label (line 3) and explicitly GUARDED line 22 (`5. Inserting the Stage-44 thresholds ...`, a CROSS-reference to upstream stage 061) and the `STAGE 65` banner as DEFERRED. I confirmed line 22 of the live file still reads `5. Inserting the Stage-44 thresholds and kappa = K_X L^2/T_X gives` — correctly LEFT UNTOUCHED. The original report's F2 also listed the py:22 cross-ref under "required change," but the orchestrator/directive intentionally deferred it to the dedicated numbering plan; that scoping is consistent with the in-loop Reading-2 policy (self-labels in-loop; cross-refs deferred). No `.wl` label edit needed (the `.wl` header self-label was already canonical). The docstring is not printed, so the transcript is unaffected by this edit — confirmed (the refreshed SymPy transcript carries no `Stage 48`/`Stage 65` docstring echo).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- L8: `STAGE 65 — THIN-WALL CONFINEMENT BRANCH` (canonical banner)
- L24: `thin-wall remainder after dropping O(ell/a) correction = 0`
- L40/41: `K_X cancellation in V0_fail^2 = 0` / `... suff ... = 0`
- L53: `STAGE 65 AUDIT PASSED`; L60: `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- L21: `STAGE 065 — THIN-WALL CONFINEMENT BRANCH` (canonical banner)
- L16/18: `PASS: independent: I1 polynomial expansion matches direct integral` / `PASS: independent: J1 = I_f / H_w under constant compressibility`
- L54: `Stage 065 Mathematica audit passed.`; L55: `# exit_code: 0`

All `expect_zero`/`expectZero` residuals are `0` / `PASS` in both engines; the two agree on every symbolic form and on the Gaussian-moment values (`J1_num = sqrt(2)*sqrt(pi)/2 = Sqrt[Pi/2]`, `J3_num = 3*sqrt(2)*sqrt(pi)/8 = (3*Sqrt[Pi/2])/4`).

**Output freshness:** The exec logs are dated 2026-06-05T13:52:31 (SymPy) and 2026-06-05T13:54:24 (Mathematica), both post-fix (after the docstring edit landed at the directive's 2026-06-05T19:46:42Z apply timestamp's local-equivalent run window) and both emit the canonical banners, confirming the committed `.txt` transcripts were regenerated post-fix. Both runs exit 0.

## Material-change assessment

`material_change`: false.

Both edits are label-only: F1 refreshes transcript stage self-labels (`STAGE 48/048` → `STAGE 65/065`) with no math-content change, and F2 changes a single docstring NUMBER (`Stage 48` → `Stage 65`) that is not printed and never enters any assertion. No derived symbolic result changed — `g_phi`, `I1`, `G_eq`, the thresholds, and all residuals are byte-for-byte identical in content. No downstream unit can depend on a docstring title or a transcript banner. No downstream staleness.

## Side observations (non-blocking)

- The original report's F2 "required change" listed the py:22 `Stage-44`→`Stage-061` cross-ref edit, but the directive correctly re-scoped it as a CROSS-reference and deferred it to the dedicated SCRIPT/OUTPUT-band numbering plan. The cross-ref remains in the file by design; it is tracked, not lost. Not a blocker.
- The directive cites the self-label as py:3 while the original report cited it as py:2; this is an off-by-one in line counting of the same docstring-title line. The diff and live file confirm the actual edited line is line 3. No substantive discrepancy.

## Verdict justification

Both findings are resolved. F2 is a strip-the-number-identical single-line docstring edit (`48`→`65`) confirmed by the one-line diff and the live file, with the py:22 cross-reference and the `STAGE 65` banner correctly left untouched per the directive's scope guard. F1 is satisfied by the regenerated transcripts, which now carry the canonical `STAGE 65 / STAGE 065` banners with all math content unchanged. Both engines exit 0 with every residual `= 0` / `PASS` and agree on all symbolic forms. No regressions in the diff, no collateral edits, no math change. `material_change: false`.
