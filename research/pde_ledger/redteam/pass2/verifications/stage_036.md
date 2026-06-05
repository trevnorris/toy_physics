---
unit_id: 036
batch: II.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-04T23:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 036

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
Two label-only edits to the SymPy source, per the diff patch (`redteam/pass2/exec_logs/stage_036_diff.patch`):
- `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:3` — docstring `Moving-throat PDE Stage 19 SymPy audit.` → `... Stage 36 ...`
- `scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:163` — `print("All Stage 19 checks passed.")` → `print("All Stage 36 checks passed.")`

Both output transcripts were then regenerated (orchestrator's exec re-run, both engines exit 0). The `.wl` was not modified (mtime unchanged at 2026-06-03 15:59:11, consistent with the directive's "make NO change to .wl logic").

**Assessment:**
The edit matches the directive exactly. The diff is purely two string substitutions of `19`→`36` in a docstring comment and a final print statement — zero changes to any equation, value, symbol, or assertion. No collateral edits: the diff hunks touch only lines 3 and 163; no `expect_zero`/`expect_true` call, no closed-form, and no constant is altered. The fix is label-only and non-functional, so there is no tautology risk and no new assertion to evaluate.

The refreshed committed outputs now read canonical banners:
- SymPy `.txt`: `STAGE 36.1 — EXACT SUPPORT-FEASIBILITY FUNCTION` (L3), `STAGE 36.2`, `STAGE 36.3`, `STAGE 36.4`, closing `All Stage 36 checks passed.` (L36). The old `STAGE 19.x` / `All Stage 19 checks passed.` labels are gone.
- Mathematica `.txt`: `STAGE 036 — SUPPORT-FEASIBILITY FRONTIER` (L3), `STAGE 036.3`, `STAGE 036.4`, closing `Stage 036 Mathematica audit passed.` (L56). The old `STAGE 019` labels are gone.

Every residual / `= 0` line and every `PASS:` line is preserved unchanged from the auditor's recorded values: SymPy `g_B,req^2/varpi^2 - (pi^2 A/8)(G-Mmix) = 0`, `dG/dxi - manifestly positive form = 0`, `G(0,delta)=0`, `G_max - closed form = 0`, final-test identity `= 0`, kappa anchor `F - R_target_sym = 0`, near-onset series `= 0`; numeric witnesses `M_mix=0.0279506713496104 / 0.810569469138702`, `G=0.465517241379310`, `R_target=1414562/558009`. Mathematica mirrors all of these as `PASS:` lines including the independent discriminant check `-72 delta^2` and the per-coefficient near-onset checks. No FAIL anywhere; both engines exit 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 36.1 — EXACT SUPPORT-FEASIBILITY FUNCTION`
- `symbolic kappa derivation: F(xi,delta) - R_target_sym = 0` (the load-bearing non-tautological anchor, unchanged)
- `G near-onset series through O(xi^2) = 0`
- `All Stage 36 checks passed.`

**Mathematica:** exit=0. Notable lines:
- `STAGE 036 — SUPPORT-FEASIBILITY FRONTIER`
- `PASS: dG/dxi numerator discriminant equals -72 delta^2`
- `PASS: symbolic kappa derivation: F(xi,delta) - R_target_sym`
- `Stage 036 Mathematica audit passed.`

**Output freshness:** confirmed. Both saved `.txt` outputs have mtime 2026-06-04 23:13:28, newer than the SymPy script (23:05:48) and the Mathematica script (2026-06-03 15:59:11). The committed `.txt` files match the exec-log transcripts line-for-line.

## Material-change assessment

`material_change`: false.

No derived result changed. The only edits are two human-readable label strings (docstring + closing print) that do not enter any computation. Every numeric and symbolic value emitted is identical to what the auditor recorded. No downstream unit can depend on a banner/print label, so no unit > 036 is rendered stale by this fix.

## Side observations (non-blocking)

None.

## Verdict justification

The single finding (F1, stale_output) is fully resolved. The diff is purely label-only (`Stage 19`→`Stage 36` at `.py:3` and `.py:163`) with zero math/value/assertion change; the `.wl` was correctly left untouched. Both refreshed transcripts now carry the canonical `STAGE 36.x` / `STAGE 036` banners and the `All Stage 36 checks passed.` / `Stage 036 Mathematica audit passed.` closers, every `= 0` residual and `PASS:` line is preserved, there are no FAILs, both engines exit 0, and the saved outputs are freshly regenerated (newer mtimes). Verdict: verified, material_change false.
