---
unit_id: 044
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

# Verification — unit 044

## Per-finding outcomes

### F1 — stale_output (self-label) [label-only]

**Classification:** resolved

**What changed:**
SymPy source `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`,
6 lines (per `stage_044_diff.patch`):
- line 3 docstring: `Moving-throat PDE — Stage 27 SymPy audit.` → `… Stage 44 SymPy audit.`
- line 53 subbanner: `27.1 — Exact continuum-selected branch equation and quadratic theorem` → `44.1 — …`
- line 80 subbanner: `27.2 — Exact continuum-selected normalization function` → `44.2 — …`
- line 111 subbanner: `27.3 — Minimal-kernel source-tied surface` → `44.3 — …`
- line 134 subbanner: `27.4 — Interference-matched tracking surface` → `44.4 — …`
- line 148 subbanner: `27.5 — Exact mismatch penalty` → `44.5 — …`

In every case ONLY the leading numeric token changed; the trailing subbanner prose is verbatim.

**Assessment:**
Correct and complete against the directive. The change is strictly label-only:
- No equation, value, variable name, or assertion appears anywhere in the diff hunks — the
  only modified text is the stage-number token. The diff context lines (the `q = sqrt(lambda0)*RU`
  computation, the `D_cont` definition, the `expect_zero` calls, etc.) are unchanged context, not edits.
- CROSS-refs were correctly LEFT untouched: I read the live file head — line 8 `Stage-24` and line 11
  `Stage-25` (upstream stages 041/042) are intact, exactly as the directive's DO-NOT-TOUCH list required
  (lines 8, 11, 55, 82, 121). The already-canonical `STAGE 44` banner (line 51) and `STAGE 44 THEOREM
  LEDGER` (line 156) were not disturbed.
- The `.wl` was not touched (no `.wl` hunk in the diff), consistent with its labels already being canonical
  (`STAGE 044` / `Stage 044 …`).
- No assertion was added, removed, or weakened, so there is no tautology risk introduced by this fix.

This was the known user-deferred SCRIPT/OUTPUT-band drift; the directive applied only the UNAMBIGUOUS
self-labels (`27`→`044`) and deferred the cross-refs to the dedicated pass. Resolved as scoped.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 44 — CONTINUUM-SELECTED RANK-2 CLOSURE` (banner now canonical)
- subbanners `44.1`–`44.5` present and in order; `STAGE 44 THEOREM LEDGER`
- residuals preserved: `quadratic branch equation = 0`, `zero-load root = 0`,
  `third-slice F at Rphi=2 = 0`, `source-tied n = 0`, `source-tied F = 0`,
  `tracking collapse of n_req = 0`, `tracking F collapse = 0`,
  `mismatch penalty in C coefficient = 0`; `# exit_code: 0`

**Mathematica:** exit=0. Notable lines:
- `STAGE 044 — CONTINUUM-SELECTED RANK-2 CLOSURE` (untouched, already canonical)
- `PASS: quadratic branch equation`, `PASS: zero-load root`, `PASS: xiPhys solve match`,
  `PASS: third-slice F at rPhi=2`, `PASS: source-tied n`, `PASS: source-tied F`,
  `PASS: tracking collapse of n_req`, `PASS: tracking F collapse`,
  `PASS: mismatch penalty in C coefficient`
- `Stage 044 Mathematica audit passed.`; `# exit_code: 0`

Every PASS / `= 0` residual present in the auditor's pre-fix inventory is still present post-fix,
and the printed symbolic forms (n_req^(cont), xi_phys, F_cont, the Solve cross-route root) are unchanged
in content. No math moved.

**Output freshness:** Both refreshed `.txt` logs carry `# date: 2026-06-05T08:11:…-06:00` headers
(post-fix orchestrator re-run), newer than the applied edit; the SymPy transcript now shows the canonical
`STAGE 44` / `44.x` labels (previously `STAGE 27` / `27.x`), confirming regeneration. Fresh.

## Material-change assessment

`material_change`: false. The edit changed only display strings (docstring + subbanner numeric tokens).
No derived result, residual, symbolic form, or assertion changed; both engines still exit 0 with identical
PASS/`= 0` lines. No downstream unit can depend on a banner label, so no `upstream_stale` propagation is
warranted on math grounds.

## Side observations (non-blocking)

- The SymPy ledger banner `STAGE 44 THEOREM LEDGER` and the Mathematica `STAGE 044` banner use different
  zero-padding (`44` vs `044`); this is a pre-existing cosmetic inconsistency in the canonical labels,
  not introduced here, and belongs to the deferred dedicated script/output-band pass. Non-blocking.
- The remaining `Stage-24` / `Stage-25` cross-refs (upstream 041/042) are intentionally still pre-renumber
  per the directive's deferral; they are out of scope for this unit and were correctly left.

## Verdict justification

Codex applied exactly the six label-only token changes the directive specified, made no collateral edits,
left all cross-refs and the `.wl` untouched, and altered no equation, value, or assertion. The refreshed
SymPy output now carries the canonical `STAGE 44` / `44.1`–`44.5` labels while every prior PASS / `= 0`
residual is preserved, and both engines exit 0. The single low-severity stale-label finding is resolved
with no math impact, so the unit verifies with `material_change: false`.
