---
unit_id: 217
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

# Verification — unit 217

The single finding in the original report is non-script (F1 = stale committed SymPy `.txt` carrying the pre-renumber `STAGE 200` banner, refreshed by the orchestrator's independent re-run). No Codex source edit was directed and none was needed — the directive has no `## Applied:`/`## Blocked:` block and no `stage_217_diff.patch` exists, consistent with a pure output-refresh. `git status --porcelain` shows the only modified tracked file is the SymPy `.txt`; the `.py` and `.wl` sources are untouched. Verification confirms (A) the output refresh landed clean (SymPy now `STAGE 217`), (B) the audit disposition holds on the refreshed artifacts, and (C) the published-value `162` resolution still HOLDS.

## Per-finding outcomes

### F1 — stale_output (committed SymPy `.txt` carried pre-renumber "STAGE 200" banner)

**Classification:** resolved

**What changed:**
No source change. The orchestrator's independent re-run (per context: sympy exit 0, mathematica exit 0, the mathematica re-run 148s under the 600s cap) regenerated the committed SymPy output. `scripts/output/moving_throat_pde_stage217_..._sympy_audit.txt` now reads `STAGE 217 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER` at line 3 (the stale `STAGE 200` banner is gone). Its mtime is Jun 9 16:51, newer than the `.py` (Jun 3 15:59). The Mathematica output was already fresh and was likewise re-captured (mtime Jun 9 16:51), banner `STAGE 217 -- FIVE-COORDINATE SIMPLEX OPTIMIZER MATHEMATICA AUDIT` (line 8).

**Assessment:**
Correct and complete. Re-run is the prescribed remedy; no `.py` logic edit was asked for or made. The refreshed transcript carries every load-bearing result: `Lifted Bezout bound = 162` (L34), `Projected one-chart Bezout bound = 750` (L43), degree pattern `(3,3,3,3,2)` (L34-38), `(5,5,5,6)` (L44-47); all four stationary-numerator identities `= 0` (L26-29); all 12 diag-reduction residuals `= 0` (L53-64) and 8 equal-mix residuals `= 0` (L69-76); `exit_code: 0` (L88). No content disagreement — only the banner/date label refreshed.

## 162 RE-CONFIRMATION (published value resolution still holds)

The pass-1 user-resolved `179`/`230` → `162` correction HOLDS. Independently re-confirmed (read of committed outputs + card/appendix/notes, no script execution):

- **Scripts DERIVE 162 (not hardcoded):** SymPy output L34-39 prints the per-factor degrees `(3,3,3,3,2)` then `Lifted Bezout bound = 162` (the value is the product of computed `total_degree`s, asserted `!= 162`, not a posited literal). Mathematica output prints `M3 lifted degrees {F_r,F_s,F_t,F_u,F_delta} = {3, 3, 3, 3, 2}` (L37), `M3 lifted Bezout product = 162` (L38), and the independent re-statement `M3 literal 3^4*2 evaluates to 162` (L39) with `PASS: M3 lifted product equals 3^4*2` (L43). Both engines DERIVE 162.
- **Arithmetic chain intact:** appendix carries `3^4\cdot2=162` (L1200) → `2\times162=324` (L1205) → `1140+324=1464` (L1261); fallback `750` per envelope → `1500` (L1207) → `1140+1500=2640` (L1266). The `1140` itself is `600+540=1140` (L1079). Card states the per-envelope `162` bound (L13), preferred budget `324` and fallback `1500` (L15). Notes carry boxed `3\cdot 3\cdot 3\cdot 3\cdot 2 = 162` (L409), `\boxed{162}` (L616), `\boxed{324}` (L620), and the candidate-set theorem `lifted Bézout bound 162` (L31, L634).
- **Zero surviving 179/230:** `grep -E "\b(179|230)\b"` across stage_217.tex, the notes `.md`, part06 appendix, both scripts, and both outputs returned ZERO matches. The wrong pass-1 values would break `2×=324` and `1140+=1464`; only `162` is consistent with the `(3,3,3,3,2)→324→1464` chain.
- **`.wl` independence (corroborating, unchanged):** the `.wl` reaches `L_x` via the symbolic `S D[delta,#] − D[S,#] delta` Jacobian and clears the M1 stationary numerators with `Cancel[Together[...]]` (the long `M1 cleared derivative numerators = {...}` block, output L14), whereas the `.py` substitutes the closed form `2*r*Delta`. Degree/Bézout checks are engine-agnostic canonical facts, legitimately cross-checked. Not a port.

## Exec log assessment

**SymPy:** exit=0. Notable lines: `r-/s-/t-/u-stationary numerator identity = 0` (committed-out L26-29); `Lifted Bezout bound = 162` (L34/log L39); `Projected one-chart Bezout bound = 750` (L48 in log); all diag and equal-mix residuals `= 0`. `# exit_code: 0`.

**Mathematica:** exit=0. Notable lines: `PASS: M1 stationary numerator identity r/s/t/u` (log L16-22); `PASS: M3 lifted product equals 3^4*2` (L43) with product `= 162` (L38); `PASS: M4 projected product equals 5*5*5*6` (L54) with product `= 750` (L49); every M5 (12) and M6 (8) check prints `PASS` with `= 0`; `All Stage 217 Mathematica audit checks passed.` (L104). `# exit_code: 0`. No FAIL anywhere.

**Output freshness:** confirmed. Both `.txt` outputs carry mtime Jun 9 16:51, newer than the `.py` (Jun 3 15:59) and `.wl` (Jun 2 12:16). SymPy banner is now canonical `STAGE 217` (the F1 stale-`STAGE 200` banner is cleared). `git status` shows the SymPy `.txt` as the only modified tracked source — no `.py`/`.wl` logic change.

## Material-change assessment

`material_change`: false. No source code changed; the only edit is the regenerated committed SymPy `.txt` (banner relabel `STAGE 200`→`STAGE 217` + transcript date refresh). No derived result changed — 162, 750, the (3,3,3,3,2)/(5,5,5,6) degree patterns, and all residual zeros are identical to the prior transcript. No downstream unit is affected.

## Side observations (non-blocking)

None beyond the single reported finding. The 324/1500/1464/2640 budget figures are prose-level arithmetic combinations of the script-derived 162/750 (and the imported 1140 support-≤4 ledger), correctly classed as out-of-stage / downstream-arithmetic in the report; the load-bearing script deliverables are the symbolic identities + the integer degree/Bézout facts, all of which hold. I concur and add nothing.

## Verdict justification

`verified`. The sole finding is a non-script `stale_output` resolved by the orchestrator's independent re-run: the refreshed SymPy `.txt` now carries the canonical `STAGE 217` banner with `Lifted Bezout bound = 162`, `Projected one-chart Bezout bound = 750`, every stationary/diag/equal-mix residual `= 0`, and `exit_code: 0`; the Mathematica output prints `PASS` on all 30 checks with no FAIL (`exit_code: 0`). No source edit was directed or made (no diff patch; git shows only the `.txt` modified), so there is no regression surface. The published-value resolution HOLDS: both engines DERIVE 162 as the product of computed degrees `(3,3,3,3,2)` (never hardcoded), the `162→324→1464` / `750→1500→2640` arithmetic chain is intact across card/appendix/notes, and a grep for `179`/`230` across all seven 217 artifacts returns zero matches. `material_change: false`.
