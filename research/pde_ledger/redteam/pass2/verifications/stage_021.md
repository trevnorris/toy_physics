---
unit_id: 021
batch: I.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-04T21:30:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 021

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
A single label-only edit in `mathematica/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.wl:195`:
`Print["Stage 004 Mathematica audit passed."];` → `Print["Stage 021 Mathematica audit passed."];`
plus regeneration of the committed transcript
`mathematica/output/moving_throat_pde_stage021_reduced_one_port_normal_form_mathematica_audit.txt`.
The captured diff (`exec_logs/stage_021_diff.patch`) shows exactly this one line changed and nothing else — the surrounding `expectZero[...]` lines (192–193) appear only as unchanged context, and the `Exit[0]` is untouched.

**Assessment:**
The edit is exactly what the directive's required change asked for (step 1, the optional cosmetic closing-label fix). The earlier line-35 banner fix (`STAGE 004`→`STAGE 021`, commit e2a4780) was already in the live script; this resolves the matching leftover closing-print literal that the audit flagged. The changed string is a non-load-bearing `Print` argument — not an equation, variable, assertion, or `expectZero` target — so it carries zero math impact. No collateral edit is present in the diff. Both halves of the required change (label fix + output regeneration) are in place. Resolved.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script was in scope for this finding (F1 is Mathematica-output staleness only); no `stage_021_sympy.log` action required.

**Mathematica:** exit=0. The orchestrator's independent re-run (`exec_logs/stage_021_mathematica.log`) confirms:
- Header (log line 8): `STAGE 021 — MAXWELL + MIXED-SECTOR REDUCTION` — banner now current.
- Closing (log line 87): `Stage 021 Mathematica audit passed.` — leftover label now fixed.
- `# exit_code: 0` (log line 88).
- Every `PASS:` line and `= 0` residual is present and unchanged from the prior transcript: gauge variations (E_w, C_a); EL residuals (Q/A/W); `sigmaCons from LinearSolve matches closed form`; `A/W exact solution residual`; `z0/z2/z4 formula` and `Sigma z0/z2/z4`; `N(omega) compact formula`; `N(0) positive-square form`; `delta D_2^(odd) composed …`; `Y2_hat minimal branch`; `Gamma5_port - a^5/(27 c_s^5)`; `N_scalar leading term`; `scalar odd order`. All 18 PASS flags identical.

**Output freshness:** The committed `…_mathematica_audit.txt` matches the orchestrator exec log line-for-line in content (banner `STAGE 021`, closing `Stage 021 Mathematica audit passed.`, all closed forms — Σ, D_cons, z0/z2/z4, N(ω), N(0), Λ2, Ŷ₂, Γ₅, N_scalar — and all PASS lines identical). The committed output now reflects the post-fix script and is no longer the ~9-day-stale `STAGE 004` transcript the auditor flagged. Output is fresh.

## Material-change assessment

`material_change`: false.

The only edit is a cosmetic closing-label print string (`Stage 004`→`Stage 021`). No equation, variable, assertion, residual, or derived result value changed — every `= 0` residual and every PASS flag is byte-identical in content to the prior transcript. No downstream unit can depend on this label string, so no unit > 021 is invalidated by this change.

## Side observations (non-blocking)

None. The diff is minimal and exactly scoped; the regenerated transcript is internally consistent and matches the orchestrator's independent re-run.

## Verdict justification

The single low-severity `stale_output` finding (F1) is fully resolved: the diff is label-only with zero math/equation/assertion change, the regenerated committed output reads `STAGE 021` in both the banner (line 3) and the closing line (`Stage 021 Mathematica audit passed.`), the orchestrator's independent Mathematica re-run exits 0 with every PASS/`= 0` residual unchanged, and the change is cosmetic with no material downstream effect. Verdict: verified, material_change false.
