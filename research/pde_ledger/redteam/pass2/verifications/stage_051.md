---
unit_id: 051
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T18:05:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 051

## Per-finding outcomes

### F1 — stale_output (unambiguous filename self-label; number-only)

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py:3` — the
docstring filename line was corrected from the pre-renumber stem
`moving_throat_pde_stage34_lowest_twin_criterion_sympy_audit.py` to the canonical
on-disk 3-digit form `moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py`.
The captured diff (`stage_051_diff.patch`) shows exactly one substantive hunk: line 3,
`stage34` → `stage051`. Both committed `.txt` transcripts were refreshed by the
orchestrator re-run.

**Assessment:**
The change matches the directive's required change item 1 verbatim (number-only,
3-digit because the on-disk filename is 3-digit). It is strip-the-number identical to
HEAD apart from the stage index — no assertion, symbol, or numeric expression touched.
I confirmed the diff contains nothing beyond the docstring line. The DO-NOT-TOUCH
deferred references were correctly left intact:
- py:20-21 still read "the Stage 050/034 product law" and "the Stage 047 coherent map"
  (compound dual-epoch + cross-ref) — unchanged.
- py:126 still reads "Stage 047/030 coherent forward map" — unchanged.
- wl:87 still reads "Stage 047/030 coherent forward map" — unchanged.
- The `STAGE 51` banners (py:63, py:149) retain their existing 2-digit format — unchanged.

As a checkpoint (`\StatusExactClosure{}`, anchor MTDC-T6), no certified constant or
cross-stage consistency claim was disturbed: `C_mix=8Λ(1-eps)/π²`, the boxed criterion
`Pi_tr ≤ 2 C_mix = 16Λ(1-eps)/π²`, `Λ_twin,req`, `M_mix^(twin,req)=G_tr/2`, the
`Z_W` forward-map consistency, and the `xi_(2x)` root all remain exactly as audited and
all still evaluate to `= 0` / `PASS` in the refreshed transcripts. Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `STAGE 51 — EXACT TRACKING-BRANCH PRODUCT` (corrected banner now in transcript).
- `Pi_tr - expected closed form = 0`
- `zeta_req at Pi = C_mix = 0` / `zeta_req at Pi = 2 C_mix minus 1 = 0`
- `Z_W^(twin,req) - forward-map(M_mix=G_tr/2) = 0`
- `G_tr(xi_(2x)) - 2 M_mix = 0`; `STAGE 51 THEOREM LEDGER`.

**Mathematica:** exit=0. Notable lines:
- `STAGE 051 — LOWEST TWIN CRITERION` (banner already canonical, unchanged).
- `PASS: Pi_tr - expected closed form`
- `PASS: zeta_req at Pi = C_mix`, `PASS: zeta_req at Pi = 2 C_mix minus 1`, `PASS: zeta_req - 1 at Pi = 2 C_mix`
- `PASS: Z_W^(twin,req) - forward-map(M_mix=G_tr/2)`
- `xi_(2x): Solve vs claim = 0` → `PASS`; `PASS: G_tr(xi_(2x)) - 2 M_mix`; `Stage 051 Mathematica audit passed.`
(The two `Limit::alimv` warnings on the endpoint limits are pre-existing and benign — the limits still evaluate correctly to 0 and Infinity.)

**Output freshness:** Confirmed. Both `.txt` outputs have mtime 2026-06-05 12:22, newer
than the SymPy script (11:49) and the Mathematica script (2026-06-03 15:59). Both
transcripts now carry the corrected stage banners (`STAGE 51` / `STAGE 051`), closing
the stale-banner half of F1.

## Material-change assessment

`material_change`: false. The only source edit is a docstring filename self-label
(number-only). No derived result, constant, threshold, or assertion changed; the
refreshed transcripts reproduce the identical math the auditor verified. No downstream
unit can depend on a docstring filename string. Nothing for downstream re-audit.

## Side observations (non-blocking)

- The compound/cross labels py:20-21, py:126, wl:87 carry stale pre-renumber indices
  (".../034", ".../030") but are correctly deferred to the dedicated SCRIPT/OUTPUT-band
  numbering plan per the directive; not blocking here.

## Verdict justification

The single low-severity `stale_output` finding is fully resolved: the unambiguous
docstring filename self-label was corrected `stage34`→`stage051` (3-digit canonical),
the diff shows no collateral edits, all deferred compound/cross references were left
untouched, and both committed transcripts were refreshed (mtimes newer than scripts,
both engines exit 0 with all checks `= 0`/`PASS` and corrected banners). As a checkpoint
stage, no certified constant or cross-stage consistency claim was perturbed.
`material_change` is false. Verdict: verified.
