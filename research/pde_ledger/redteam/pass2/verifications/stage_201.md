---
unit_id: 201
batch: VI.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 201

The single finding in the original report is non-script (F1 = stale committed SymPy `.txt` carrying a pre-renumber "STAGE 184" banner, refreshed by the orchestrator's independent re-run). No Codex directive was written (the source already prints `STAGE 201`), so `redteam/pass2/directives/stage_201.md` and `..._diff.patch` correctly do not exist. Verification confirms (A) the output refresh landed clean, (B) the audit's clean/independent disposition still holds on the refreshed artifacts, and (C) `material_change: false`.

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 184" banner)

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent re-run regenerated the committed SymPy output. `scripts/output/moving_throat_pde_stage201_explicit_realization_compiler_and_canonical_orbit_projection_sympy_audit.txt` now reads `STAGE 201 — EXPLICIT REALIZATION COMPILER AND CANONICAL ORBIT PROJECTION` (line 3) and `STAGE 201 LEDGER` (line 350); the stale "STAGE 184" banner (orig report cited `:11`/`:358`) is gone — `grep "STAGE 184"` returns nothing. The `.txt` mtime is Jun 9 16:51, newer than the `.py` (Jun 3 15:59). The Mathematica output was already current and was likewise refreshed (mtime Jun 9 16:51), banner `STAGE 201 -- EXPLICIT REALIZATION COMPILER AND CANONICAL ORBIT PROJECTION` (line 3) and `STAGE 201 MATHEMATICA LEDGER` (line 110).

**Assessment:**
Correct and complete. The orchestrator re-run is the prescribed remedy for a stale capture whose live source string is already correct. The refreshed SymPy transcript carries every result as an equality: section identity `M_* S - I_3` (line 57), chart consistency `exp(q_tr/nt/eta) - R = 0` (lines 101–103), mismatch chart `m_T/m_K/m_mu = 0` (lines 104–106), uniqueness triple `dT/dKeta/dmu` residuals `= 0` (lines 175–177), and the fixed-point closure `Delta_real^int = 0` (line 371). No FAIL / AssertionError / Traceback anywhere. The Mathematica output prints 22 `PASS:` lines covering M1–M8 (lines 17, 20, 26, 28, 30, 37, 40, 43, 50, 52, 61, 64, 73, 75, 77, 79, 81, 83, 86, 93, 104, 107) with no FAIL. Output refresh is clean.

## Disposition re-confirmation (post-refresh)

- **Genuinely independent `.wl`:** confirmed still holding on the refreshed artifacts. The load-bearing canonical right-section `S` (and `dx_rep = -S q`, the projection, fixed-point, and linearized compiler all derive from it) is extracted by two different methods — SymPy inverts the dependent 3×3 block of `M_*` and scatters it into the 8×3 matrix (`.py:66–74`); Mathematica stacks `M_*` on five unit-row free-coordinate constraints into an 8×8 system and `LinearSolve`s it column-by-column (`.wl:120–144`). The refreshed Mathematica output independently re-confirms `Det[dependent block] - (1+chi0_star) = 0` (line 73), a diagnostic the `.py` lacks. Both engines land on the identical closed forms for `S`, `dx_rep`, and the uniqueness triple. Not a transliteration.
- **Non-tautological assertions:** confirmed. Every `expect_zero`/`PASS` subtracts an independently-built object from a separately-defined closed form or contracts it with `M_*` (e.g. SymPy A4 hand-typed `expected_dx_rep`, A5 `M_* dx_rep + q`; Mathematica M6 `Solve`-derived triple vs closed form). The trivial-case probes are real can-fail statements: at unit ratios `dx_rep = 0` and `Pi^can(x) = x` (SymPy line 371 / Math line 50, 64), while the eta-only perturbation gives a nonzero `Keta` repair (Math line 52).
- **Reconciliation:** confirmed. The report's 11 symbolic deliverables (`M_*`, `S`, `dx_rep`, `m_T/m_K/m_mu`, intrinsic packet, `Pi^can`, uniqueness triple, `Det = 1+chi0_*`, linearized compiler) all reconcile MATCH against the notes/appendix; 0 misaligned. The lone `Det = 1+chi0_*` literal is correctly classed INTERNAL (invertibility diagnostic, not a boxed deliverable).

## Exec log assessment

**SymPy:** exit=0 (`stage_201_sympy.log:377` `# exit_code: 0`). Notable lines: the captured ledger tail reproduces the fixed-point criterion `x in Z_* <=> chi_Q(x)=1 and Pi_O^can(x)=x  equivalently  Delta_real^int(x|Z_*) = 0`. All §I–§VII residuals print `= 0`; no FAIL/AssertionError/Traceback.

**Mathematica:** exit=0 (`stage_201_mathematica.log:122` `# exit_code: 0`). Notable lines: `PASS: M1 residual: M_* . S - I_3` (txt:17); `PASS: M6 residual: Det[dependent block] - (1+chi0_star)` (txt:73); `PASS: M6 check: Solve returns a unique dependent triple` (txt:77); `PASS: M8 residual: M_* . (Delta x + Delta x_rep^lin)` (txt:104). 22 PASS, 0 FAIL.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime Jun 9 16:51, newer than the `.py` (Jun 3 15:59) and `.wl` (Jun 2 00:36). The SymPy banner is now canonical "STAGE 201" (the F1 stale-184 banner is cleared).

## Material-change assessment

`material_change`: false. No source code changed; the only edit was the regenerated committed SymPy `.txt` (banner relabel 184→201 + transcript refresh). No derived result changed, so no downstream unit is affected.

## Side observations (non-blocking)

None beyond the single reported finding. The `exp(q)-R` chart-consistency checks (A2) are definitional rather than deep, but the auditor flagged them as such and the load-bearing independence lives in the `S`/`dx_rep` route, which is genuinely cross-checked. I concur and add nothing.

## Verdict justification

`verified`. The sole finding is a non-script stale_output resolved by the orchestrator re-run: the refreshed SymPy `.txt` now carries the canonical "STAGE 201" banner (line 3 / ledger line 350), the stale "STAGE 184" is gone, its mtime postdates the `.py`, and every residual reads `= 0`; the Mathematica output prints 22 PASS with no FAIL; both exec logs report `# exit_code: 0`. The audit's clean/independent disposition holds on the refreshed artifacts — the `.wl`'s augmented-`LinearSolve` route is genuinely independent of the `.py`'s dependent-block-inversion route, assertions are non-tautological, and reconciliation is 11/11 MATCH with 0 misaligned. No regressions; `material_change: false`.
