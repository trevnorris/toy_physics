---
unit_id: 059
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 059

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:73` — inserted `expect_positive("Xi_suff - Xi_fail (ordered)", Xi_suff_ordered - Xi_fail_ordered)` directly after the `zeta_req_branch = ...` line (originally line 72).
- `mathematica/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.wl:73` — inserted `expectPositive["Xi_suff - Xi_fail (ordered)", xiSuffOrdered - xiFailOrdered];` directly after the `xiSuffOrdered = ...` line (originally line 72).

**Assessment:**
The edits match the directive's required-change byte-for-byte. The diff capture (`redteam/exec_logs/stage_059_diff.patch`) confirms only these two single-line insertions; no collateral edits. Both insertions sit inside the existing positive-real assumption blocks (`Pe_req > 0`, `Delta0 > 0`, `delta_gap > 0` for SymPy; `peReq > 0`, `delta0 > 0`, `deltaGap > 0` for Mathematica), so the helpers operate under the right hypotheses. The new assertion is non-tautological: it tests positivity of `Pe_req*delta_gap/(Delta0*(Delta0+delta_gap))`, which is positive only because the structural relation `DeltaInf_ordered = Delta0 + delta_gap` is correct and `Pe_req` carries the prescribed sign — a sign or role flip upstream would break this check. The SymPy transcript line `Xi_suff - Xi_fail (ordered) = Pe_req*delta_gap/(Delta0*(Delta0 + delta_gap))` matches the auditor's prescribed factored form exactly, and Mathematica emits the expected `PASS: Xi_suff - Xi_fail (ordered)`. The unused `zeta_req_branch` line was explicitly permitted to remain, and it does.

## Exec log assessment

**SymPy:** exit=n/a. The canonical exec log `redteam/exec_logs/stage_059_sympy.log` was not produced by the orchestrator (only the `stage_059_diff.patch` is present in `exec_logs/`). Falling back to `scripts/output/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.txt` (mtime 2026-05-26 11:05, newer than the script's 11:04). Notable lines:
- `Xi_suff - Xi_fail (ordered) = Pe_req*delta_gap/(Delta0*(Delta0 + delta_gap))` (new, matches directive)
- `Omega^2 linear coefficient = 0` (pre-existing, passes)
- `Xi_fail*DeltaInf saturates at Pe_star diff = 0` (pre-existing, passes)
- `Stage 42 audit passed.`

**Mathematica:** exit=n/a. The canonical exec log `redteam/exec_logs/stage_059_mathematica.log` was not produced. Falling back to `mathematica/output/moving_throat_pde_stage059_operator_branch_residual_bounds_mathematica_audit.txt` (mtime 2026-05-26 11:05, newer than the script's 11:04). Notable lines:
- `Xi_suff - Xi_fail (ordered) = (deltaGap*peReq)/(delta0*(delta0 + deltaGap))`
- `PASS: Xi_suff - Xi_fail (ordered)` (new)
- `PASS: Xi_fail*DeltaInf saturates at Pe_star` / `PASS: Xi_suff*Delta0 saturates at Pe_star` (pre-existing)
- `Stage 059 Mathematica audit passed.`
- Benign `Limit::alimv` warning is unchanged from the audit transcript (acknowledged in the original report; the limit value `-1 + 4/Pi` is correct because `Omega_Pe^2` is analytic at zero).

**Output freshness:** Both `.txt` outputs (mtime 11:05) were regenerated after the script edits (mtime 11:04). Both transcripts terminate at their `passed` banners, indicating clean exits despite the missing dedicated exec-log captures.

## Material-change assessment

`material_change`: false.

The edit only adds a positivity assertion on previously-defined symbolic quantities. No derived result, threshold value, expansion coefficient, or downstream-visible constant changes. The exposed downstream invariants (`Xi_fail = Pe_req/DeltaInf`, `Xi_suff = Pe_req/Delta0`, the `(4 - pi)/pi` weak-coupling slope, and the threshold-saturation relation) are all unaffected. No downstream unit depends on the new assertion line.

## Side observations (non-blocking)

- The SymPy and Mathematica banners still say "Stage 42" / "STAGE 042" while the file/paper say Stage 059. The auditor flagged this as cosmetic-only and did not raise it as a finding; Codex correctly left it alone.
- The `zeta_req_branch = sp.simplify(A_K * Omega(Pe_req) ** 2)` line on SymPy:72 remains unreferenced. The directive explicitly noted this as collateral and instructed leaving it; this is by design.
- An older verification file from the 2026-05-22 first-pass iteration (covering F1–F4 from that earlier audit) was present at this path and has been overwritten by this 2026-05-26 second-pass verification, which corresponds to the single-finding audit in `redteam/reports/stage_059.md`.

## Verdict justification

The single finding was applied exactly as prescribed: one-line insertions in each engine, no deviations, no collateral edits per the diff. The new assertion exercises previously-dead `_ordered` scaffolding via a non-tautological positivity check — a quotient whose numerator carries `delta_gap` makes the check sensitive to the upstream `DeltaInf_ordered = Delta0 + delta_gap` construction and to the sign of `Pe_req`. Both transcript outputs include the new line in the form the auditor prescribed and terminate at their `passed` banners. Dedicated exec-log files were not captured by the orchestrator, but the regenerated `.txt` outputs (mtimes newer than the script mtimes) and the explicit `PASS:` lines provide equivalent evidence of clean runs. Verdict: verified.
