---
unit_id: 093
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 093

## Per-finding outcomes

No findings were raised by the auditor (verdict: clean, findings_count: 0). No directive file exists for this unit; there is nothing to verify against findings.

The auditor report's "Findings" section explicitly states "(no findings)" and documents three concerns that were considered and rejected (potential weak `c_pole=1/4` check, missing scriptable coverage for hygiene Checks items 2/3 and the omega-dependent module, and a banner mislabel). All three were correctly classified as either status-only carve-out cases or cosmetic-only.

## Banner-relabel cross-check (orchestrator carve-out)

The auditor flagged one cosmetic side observation in §Findings #3: line 26 originally printed `banner["STAGE 076 — GROUPED-P2 STATUS UPDATE"]`. Per orchestrator note, the IV.1 batch banner sweep was applied to relabel this to STAGE 093.

Verification:

- `mathematica/moving_throat_pde_stage093_grouped_p2_status_update_mathematica_audit.wl:26` now reads `banner["STAGE 093 — GROUPED-P2 STATUS UPDATE"];` — relabel landed cleanly.
- The exec log header (`output/...stage093...txt:3`) likewise shows `STAGE 093 — GROUPED-P2 STATUS UPDATE`, consistent with a post-fix rerun.
- The pre-existing terminal banner at script line 46 (`"Stage 093 Mathematica audit passed."`) was already correct and is unchanged.
- No other lines in the script were touched by the relabel; the obstruction formula at wl:30, the four `expectZero` assertions at wl:40–43, and the surrounding scaffolding are byte-identical to what the auditor reviewed.

The two other numbering legacies the auditor noted (paper section header "Stage~110", notes references to "Stages 74–75") are in prose, out of scripts-only verifier scope, and explicitly out of the IV.1 banner sweep's mandate.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script exists for this unit. Per the status-only carve-out (`is_status_only_candidate: true` in MANIFEST), missing_sympy is not a finding.

**Mathematica:** exit=0 (inferred — terminal line `"Stage 093 Mathematica audit passed."` is reached, and `Exit[0]` follows at wl:48; no `Exit[1]` fail-path triggered). Notable lines from the log:

- L9–L10: `c_pole - 1/4 = 0` / `PASS: c_pole - 1/4`
- L11–L12: `c_geom - 3/4 = 0` / `PASS: c_geom - 3/4`
- L13–L14: `rho_alpha - 4/3 = 0` / `PASS: rho_alpha - 4/3`
- L15–L16: `zeta_req - 1/3 = 0` / `PASS: zeta_req - 1/3`
- L18: `Stage 093 Mathematica audit passed.`

All four `expectZero` checks resolve to literal `0` after `FullSimplify`; no `FAIL` lines appear.

**Output freshness:** confirmed. The Mathematica `.txt` output mtime is 2026-05-27 14:29; the `.wl` script mtime is 2026-05-27 11:12. The log post-dates the script (and post-dates the banner relabel that landed at 11:12), so the saved output reflects the current, relabeled script.

## Material-change assessment

`material_change`: false.

The only edit in this unit's scope was a banner string relabel (display text in a `banner[...]` call). No formula, assertion, variable, or numeric constant was touched. Downstream units cannot depend on a print-side label, so no upstream-stale propagation is warranted for unit 093 in isolation. (The orchestrator's batch-wide IV.1 banner sweep may still mark units >093 stale for its own bookkeeping, but stage 093's edit specifically introduces no derived-result change.)

## Side observations (non-blocking)

1. Prose-side numbering legacies remain: the paper card's section header reads "Stage~110" and the notes file references "Stages 74–75" for the upstream obstruction-formula provenance. These are out of scripts-only verifier scope and are documented in the original auditor report (§Findings #3). No action requested here; flagging for the prose-cleanup tracker, not for this verification.
2. The four `expectZero` assertions A2–A4 (`c_geom - 3/4`, `rho_alpha - 4/3`, `zeta_req - 1/3`) are algebraic consequences of A1 (`c_pole - 1/4`) given the deterministic `cGeom = 1 - cPole`, `rhoAlpha = 1/cGeom`, `zetaReq = cPole/cGeom` definitions at wl:31–33. The auditor already noted this redundancy and accepted it as legitimate carry-forward depth for a status-only unit. Not a finding; recorded only because verifier independently confirmed the dependency structure on the script.

## Verdict justification

Audit unit 093 was rated clean by the auditor with zero findings. The single side-observation banner relabel (STAGE 076 → STAGE 093) landed correctly at `.wl:26`, the exec log was regenerated after the edit (log mtime 14:29 > script mtime 11:12), the rerun shows all four `expectZero` checks passing with residual `0`, and no derived numerical content was modified. `material_change` is false. Verdict: verified.
