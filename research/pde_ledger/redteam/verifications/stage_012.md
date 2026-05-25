---
unit_id: 012
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 012 (v2 re-audit follow-up)

This supersedes the v1 verification. The v2 paper-grounded re-audit (report `redteam/reports/stage_012.md` dated 2026-05-25) raised a single low-severity `tautological_check` finding against the Mathematica `M1` block; Codex applied the directive's prescribed fix. This verification confirms the fix.

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
At `mathematica/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.wl:50-56`, Codex deleted the six tautological `expectZero` calls (former lines 50-66, e.g. `expectZero["M1 Z0 primitive one-port", Z0form - Q/Delta]` where line 41 had just bound `Z0form = Q/Delta;`) and replaced them with the exact Print group prescribed by the directive:

```mathematica
Print["M1 primitive one-port forms (carried from Stage 4 / Stage 5):"];
Print["  Z0 = ", fmt[Z0form]];
Print["  Z2 = ", fmt[Z2form]];
Print["  Z4 = ", fmt[Z4form]];
Print["  N0 = ", fmt[N0form]];
Print["  N2 = ", fmt[N2form]];
Print["  N4 = ", fmt[N4form]];
```

The form bindings on lines 41-48 (`Z0form = Q/Delta;` through `N4form = ...;`) are preserved — they remain needed by the substantive M2..M9 blocks downstream.

**Assessment:**
- The diff at `redteam/exec_logs/stage_012_diff.patch` is exactly the prescribed surgical edit (17 lines removed, 7 lines added); no collateral changes elsewhere in the file, no other files touched.
- The new block is documentation-only (no assertion), so it cannot reintroduce a tautological PASS — it simply emits the carried-forward forms for transcript readability.
- All M2..M9 substantive checks are intact: `grep -n "^expectZero|^expectNonZero" ...wl` shows 24 `expectZero` calls at lines 156-265 and 4 `expectNonZero` calls at lines 269-282, matching the original M2..M9 inventory from the report (12 M2 checks, 1 M3, 2 M4, 1 M5, 2 M6, 2 M7, 4 M8 zero + 2 M8 nonzero, 2 M9 nonzero).
- The directive's `## Applied: F1` block reports `deviation: none`, consistent with the diff.

## Exec log assessment

**SymPy:** exit=0 (inferred from `STATUS: PASS` transcript end-state). No `redteam/exec_logs/stage_012_sympy.log` was captured this iteration — consistent with the directive's instruction "Do NOT run python or mathematica. Only edit files." The SymPy script itself was not touched (mtime 2026-05-21 11:37 unchanged), but the saved output `scripts/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_sympy_audit.txt` has mtime 2026-05-25 17:15 — the orchestrator regenerated it post-edit. Final line: `STATUS: PASS`. Since `assert_zero`/`assert_nonzero` raise on failure, reaching that line confirms all 15 sympy assertions still pass.

**Mathematica:** exit=0 (inferred from `STATUS: PASS` transcript end-state). No `stage_012_mathematica.log` captured this iteration (same reason). The saved transcript `mathematica/output/moving_throat_pde_stage012_projected_maxwell_primitive_bridge_mathematica_audit.txt` has mtime 2026-05-25 17:15 (script mtime 17:14) — regenerated post-fix. Notable transcript content:

- Lines 2-8: the new header and six `Z0 = ... N4 = ...` documentation lines, exactly as the directive prescribed.
- Lines 9-32: all twelve M2 PASS lines (six `series closed form` + six `partial route`) preserved, each with `residual = 0`.
- Lines 33-64: M3..M9 PASS lines unchanged (static Xi1, two K round-trips, fixed-target compat shift, transported normalization K + round-trip, transported compat surface + shift, four z0-channel cancellation PASSes, two z0-channel retention PASSes with explicit non-trivial residuals `ell/Delta` and `-(ell*(2*P^2 + Delta*pGoal*Q))/(Delta^3*pGoal)`, two z4-flip mutation PASSes with explicit non-trivial residuals).
- Line 65: `STATUS: PASS`.
- The six pre-fix `PASS: M1 Z0/Z2/Z4/N0/N2/N4 primitive one-port` lines are absent, as intended.

**Output freshness:** confirmed. Mathematica transcript mtime (17:15) > script mtime (17:14). SymPy transcript mtime (17:15) > script mtime (May 21 11:37). Both outputs are post-edit.

## Material-change assessment

`material_change`: false.

The removed checks were identically-zero residuals by construction (`Z0form - Q/Delta` with `Z0form` bound to `Q/Delta` one line earlier; analogous for Z2..N4). They verified nothing beyond `FullSimplify[0] === 0`. Their removal:
- Shrinks the transcript PASS count by 6 (from 34 to 28), but
- Does not alter any derived closed form, solved K surface, compatibility-shift expression, partial-derivative cancellation outcome, or mutation residual.
- Does not change any quantity that downstream stages might consume (the substantive provenance of the M1 primitive forms lives in Stages 4 / 5, per the SymPy docstring at line 55).

No downstream re-audit warranted by this edit.

## Side observations (non-blocking)

- The orchestrator did not produce `stage_012_sympy.log` or `stage_012_mathematica.log` this iteration. This is consistent with the directive's "do not run scripts" clause, and output transcripts (regenerated post-edit) provide sufficient pass evidence. If the orchestrator's contract normally requires exec logs for every iteration, that's a process gap worth noting but does not block verification.
- The pre-existing v1 verification file (dated 2026-05-21, covering four findings F1-F4 from the v1 audit) was overwritten by this v2 verification. The v1 fixes (second-engine `.wl` creation, K solve round-trips, derived-z0 q1/d1 cancellation, z0..n4 closed-form anchors) remain in the scripts and continue to PASS — they are the substantive M2..M9 / A1..A15 content this v2 verification confirms is intact.

## Verdict justification

The single v2 finding is fully resolved with a clean surgical edit that exactly matches the directive. The six tautological M1 `expectZero` calls became Print documentation; the upstream form bindings and all downstream substantive M2..M9 checks are untouched; the regenerated Mathematica transcript ends in `STATUS: PASS` with the prescribed new header and the prescribed six PASS removals; the SymPy script and its transcript are unchanged (still `STATUS: PASS`). No regressions, no new tautologies, no collateral edits. Removing inert tautologies does not alter any derived content, so `material_change` is false.

stage 012: verified, material_change: false
