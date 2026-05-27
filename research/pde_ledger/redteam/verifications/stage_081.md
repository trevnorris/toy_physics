---
unit_id: 081
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 081

## Per-finding outcomes

### F1 — paper_misalignment (target_mismatch)

**Classification:** resolved

**What changed:**
Orchestrator-direct edit (no Codex invocation) applied direction (a) — script-side relabel only — to the SymPy script:

- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:3` docstring changed from `SymPy audit for Stage 64.` to `SymPy audit for Stage 081.`.
- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:28` banner changed from `banner("STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS")` to `banner("STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS")`.

The diff captured at `redteam/exec_logs/stage_081_diff.patch` shows exactly these two hunks and nothing else (one-line edits at lines 3 and 28 of the `.py`). The companion `.wl` file already carried the correct `STAGE 081` banner per the auditor report (auditor cited `mathematica/...wl:38` as already-correct) and the diff confirms the `.wl` was not touched.

**Assessment:**
The edit precisely matches the auditor's "Verification" prescription in `redteam/reports/stage_081.md:95`:
> After resolution, the SymPy script's line 3 docstring reads "SymPy audit for Stage 081." and line 28 banner reads `"STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS"`. A re-run of `redteam exec-sympy 081` produces a transcript whose third line reads `STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS`.

Reading the current files:

- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:3` now reads `SymPy audit for Stage 081.` (verified).
- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:28` now reads `banner("STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS")` (verified).
- `scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt:3` now reads `STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS` (verified).

No algebra, assertion, or symbol-definition lines were touched. The two `expect_zero` assertions (`Q(0)-1` at line 40 and `Q(1)-2` at line 41) are identical to v1, and the refreshed transcript at `scripts/output/...txt:7-8` confirms both still produce `Q(0)-1 = 0` and `Q(1)-2 = 0`. There is no tautological-assertion concern because no assertion changed.

No collateral edits beyond the two-line label change.

## Exec log assessment

**SymPy:** exit=0 (per orchestrator note; no `redteam/exec_logs/stage_081_sympy.log` file was captured this iteration, only the diff patch — see note below). The saved transcript `scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt` was refreshed post-fix (mtime 2026-05-27 02:15:47 > script mtime 2026-05-27 02:05:20). Notable lines from the refreshed transcript:

```
STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS
Q(0)-1 = 0
Q(1)-2 = 0
Blocking ceiling eps_blk < 0.40526368971137149977
```

The `expect_zero` helper raises `AssertionError` on any non-zero residue; both `= 0` prints together with the script reaching the `FINAL LEDGER` banner imply exit 0. The PASS set is identical to v1 (the two anchor `= 0` results and the same `Pi_xxx/C_mix` numerical values to all printed digits).

**Mathematica:** exit=n/a. The `.wl` was not touched this iteration (the auditor explicitly noted the `.wl` banner at line 38 was already correctly `STAGE 081`). No `stage_081_mathematica.log` is present, which is correct given no Mathematica re-run was required.

**Output freshness:** confirmed.

- `scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py` mtime: 2026-05-27 02:05:20.
- `scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt` mtime: 2026-05-27 02:15:47 (newer than script — refreshed post-fix).
- `mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl` mtime: 2026-05-23 10:36:58.
- `mathematica/output/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.txt` mtime: 2026-05-25 00:21:34 (newer than the `.wl`; untouched this iteration, as expected).

Note (non-blocking): no `redteam/exec_logs/stage_081_sympy.log` file exists; only `stage_081_diff.patch` is present in `redteam/exec_logs/`. The orchestrator's contextual statement plus the refreshed transcript timestamps and assertion prints are sufficient to confirm a passing run.

## Material-change assessment

`material_change`: false.

The only change is a two-line string-literal relabel in the SymPy script's docstring and banner. No symbol definitions, no assertion expressions, no numerical thresholds, no print statements with derived results, and no Mathematica content were modified. No downstream unit depends on the banner text of an upstream sympy script.

## Side observations (non-blocking)

- The SymPy script's `dQ/dzeta` print (line 42) and blocking-ceiling print (line 73) remain print-only, not assertions, consistent with the auditor's row "partial (print-only)". This is the v1 design choice the auditor explicitly tolerated (the load-bearing Mathematica side carries the corresponding assertions M9 and M4–M8). Not a finding for this iteration.
- A prior version of this verification file (dated 2026-05-25) on disk covered the earlier four-finding iteration (F1–F4 against the Mathematica mirror). That file is overwritten by this current verification for the new single-finding iteration.

## Verdict justification

The single F1 paper_misalignment was resolved by orchestrator-direct two-line relabel of the SymPy docstring (line 3) and banner (line 28) from "Stage 64" to "Stage 081", matching the auditor's prescribed verification text exactly. The diff at `redteam/exec_logs/stage_081_diff.patch` shows no collateral edits, no algebra changes, and no touch of the `.wl`. The refreshed SymPy transcript now reads `STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS` on its banner line and still emits both assertion prints (`Q(0)-1 = 0`, `Q(1)-2 = 0`) — meaning the script still exits 0 with the v1 PASS set, no tautology, no regression. `material_change: false`. Verdict: `verified`.
