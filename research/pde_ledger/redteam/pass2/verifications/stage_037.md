---
unit_id: 037
batch: III.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 037

The auditor recorded two items (F1 `stale_output`, F2 deferred SCRIPT/OUTPUT-band stale self-labels). The directive carries a single in-loop finding (`findings_count: 1`): the UNAMBIGUOUS stale SymPy SELF-labels (docstring filename + header `Stage 20`, where `20` is THIS stage's pre-renumber number → canonical `037`). The `.wl` was correctly excluded (already canonical) and the F2 deferred cross-refs (`Stage-17/19`) plus banner zero-padding were explicitly routed away from this fix_loop. The F1 output-freshness component is the orchestrator's standard fresh-run step. I verify the directive's single applied finding.

## Per-finding outcomes

### F1 — stale_output (self-label) [label-only]

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py` docstring:
- L3 `moving_throat_pde_stage20_continuum_kernel_sympy_audit.py` → `moving_throat_pde_stage037_continuum_kernel_sympy_audit.py`
- L5 `Stage 20 SymPy audit:` → `Stage 37 SymPy audit:`

The captured diff (`stage_037_diff.patch`) contains exactly and only these two lines; no other hunk is present. Confirmed against the live file: L3/L5 now read the canonical stage-number tokens.

**Assessment:**
The edit is strictly LABEL-ONLY. Only the stage-number tokens (`20`→`037` in the filename, `20`→`37` in the header) changed; the surrounding text (`moving_throat_pde_..._continuum_kernel_sympy_audit.py`, `... SymPy audit:`) is byte-identical otherwise. No equation, value, variable name, or assertion byte was touched.

The directive's "DO NOT TOUCH" items were correctly left in place, confirmed by reading the live source:
- L6 `the Stage-17/19 reduced branch data` — cross-ref to upstream sources, UNTOUCHED (correct).
- L224 `Stage-17/19 branch data exactly.` — cross-ref, UNTOUCHED (correct).
- L222 banner `STAGE 37 AUDIT COMPLETE` — already-canonical banner, UNTOUCHED (correct; the `37`→`037` zero-pad normalization is deferred to the dedicated numbering pass, not this loop).

The `.wl` source was not touched (the diff patch references only the `.py`; the `.wl` mtime is unchanged at 2026-06-03 15:59). Codex's `## Applied: F1` block records `files_changed` = the `.py` only, `deviation: none` — consistent with the captured diff. No collateral edit.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 37 — CONTINUUM-KERNEL EXTRACTION AUDIT` (now canonical, was `STAGE 20`)
- `STAGE 37 AUDIT COMPLETE` (now canonical, was `STAGE 20 AUDIT COMPLETE`)
- All residuals zero and all closed forms emitted: `A continuum formula = 0`, `M_mix continuum formula = 0`, `delta continuum formula = 0`, `Sigma_wall - [Xi I + alpha v v^T] = [[0,0],[0,0]]`, `kappa0 - 2 sqrt(2)/pi = 0`, etc.

**Mathematica:** exit=0. Notable lines:
- `STAGE 037 — CONTINUUM-KERNEL EXTRACTION` (canonical; was `STAGE 020` in the stale committed output)
- `Stage 037 Mathematica audit passed.`
- Every `PASS:` line present, including the load-bearing `PASS: Sigma_wall - [Xi I + alpha v v^T]`, `PASS: Sigma_wall (2,2) consistency with ansatz`, `PASS: A numerator matches Schur form`, `PASS: delta numerator matches closed form`, and all 7 closed-form / 2 gate-numerator checks.

The `.wl` was not edited, so its banner was already canonical (`STAGE 037`); the prior stale `STAGE 020` in the committed output was the F1 staleness, now eliminated by the fresh run.

**Output freshness:** confirmed. Source mtimes: `.py` 2026-06-05 07:40:20, `.wl` 2026-06-03 15:59:11. Refreshed output mtimes: sympy `.txt` 2026-06-05 08:09:33, mathematica `.txt` 2026-06-05 08:09:33 — both newer than their respective sources. The prior `STAGE 20/020` banners are gone; both transcripts now show the canonical banners and every prior PASS/`= 0` residual remains.

## Material-change assessment

`material_change`: false.

The only edit is two stale stage-number tokens in a SymPy docstring (non-executing comment text). No derived result, assertion, closed-form, or printed value changed — confirmed by the result lines in both refreshed transcripts being identical (modulo the banner number) to the audit's described pre-fix content. No downstream unit can depend on a docstring label.

## Side observations (non-blocking)

- The F2 deferred-pass items remain open by design: the SymPy banners L51/L222 still read `STAGE 37` (un-padded) vs. the `.wl`'s `STAGE 037`, and the `Stage-17/19` cross-refs (L6, L224) are intact. These are correctly routed to the dedicated content-keyed numbering pass (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`) per the directive; not a blocker here.
- The committed mathematica `.txt`'s prior internal inconsistency (stale `STAGE 020` header but current `Stage 037` footer) is resolved by the fresh run — header now `STAGE 037`.

## Verdict justification

The single in-loop finding (F1, stale SymPy self-labels) is resolved: the captured diff is strictly label-only (exactly L3/L5 stage-number tokens), the deferred cross-refs and `.wl` were correctly left untouched, both refreshed `.txt` outputs now show canonical banners with newer mtimes than their sources, every prior PASS/`= 0` residual is intact, and both engines exit 0. No regression appears in the diff or logs. No math/value/assertion changed, so `material_change: false`. Verdict: verified.
