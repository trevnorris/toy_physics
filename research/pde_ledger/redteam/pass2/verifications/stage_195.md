---
unit_id: 195
batch: V.3
verifier_model: Opus 4.8 (1M context)
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 195

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
No source edit (the directive required only a re-run, per report F1 "Required change"). The committed transcripts were regenerated:
- `scripts/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.txt` — mtime now 2026-06-09 14:06, ≥ `.py` (2026-06-03 15:59).
- `mathematica/output/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_mathematica_audit.txt` — mtime now 2026-06-09 14:06, ≥ `.wl` (2026-06-01 11:49).

**Assessment:**
The stale STAGE-178 banner that triggered F1 is gone. The refreshed SymPy `.txt` line 3 reads `STAGE 195 — EXACT SOURCE-MAP REDUCTION OF THE CANONICAL OUTGOING BRANCH` and line 122 reads `STAGE 195 LEDGER`. A `grep` for `178` across both outputs returns nothing. Every SymPy residual line ends `= 0` (Section I lines 30–31; II lines 56–58; III lines 73–75; IV lines 98–99, 110–111; V lines 116–119). The Mathematica `.txt` carries `STAGE 195 - SOURCE-MAP REDUCTION…` (line 3) and `STAGE 195 MATHEMATICA LEDGER` (line 73), with all 16 checks emitting `= 0` and a matching `PASS:` line; no `FAIL`/`error`/inequality line exists in either file. Refresh landed; the disposition (genuinely independent `.wl`, X−X de-taut held, 0 reconciliation misalignments) is undisturbed since no source changed.

## Exec log assessment

**SymPy:** exit=0 (inferred from refreshed committed `.txt`; orchestrator re-run). Notable lines:
- `STAGE 195 — EXACT SOURCE-MAP REDUCTION OF THE CANONICAL OUTGOING BRANCH` (line 3) — canonical banner, no longer STAGE 178.
- `observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1) = 0` (line 56) — the headline factorization derived from the observable residual, reduces to 0.
- `linearized N_Q - 1 = 0` (line 110) and `Delta_norm(canonical) = 0`-family closes V (lines 116–119).

**Mathematica:** exit=0 (inferred from refreshed committed `.txt`). Notable lines:
- `chi_Q from native DtN deformation = (3*(betaStretch^5*scaleS + 9*sigma5))/(3*scaleS - sigma0)` (line 46) — χ_Q DERIVED from the Hankel-operator route, matching the `.py`'s posited closed form.
- `PASS: observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1)` (line 25).
- 16 `PASS:` lines total; `N_Q deformation law` / `linearized N_Q - 1` / canonical-closure all PASS.

**Output freshness:** Confirmed. Both `.txt` mtimes (2026-06-09 14:06) exceed their script mtimes (`.py` 06-03 15:59; `.wl` 06-01 11:49). Stale STAGE-178 transcript fully replaced.

## Material-change assessment

`material_change`: false. No source code was edited — only the committed `.txt` transcripts were regenerated. No derived result changed; downstream units are unaffected.

## Side observations (non-blocking)

The paper-side card-text lag noted by the audit context is USER-DEFERRED to paper-cleanup and is out of scope for the scripts-only verifier; it does not block verification.

## Verdict justification

The single low-severity F1 (`stale_output`) is resolved: both committed transcripts were re-run post-fix, now carry the canonical `STAGE 195` banner (no surviving `178` reference), and every check emits `= 0` / `PASS:` with no failures. The disposition holds unchanged — the `.wl` Section IV remains a genuinely independent re-derivation of χ_Q from the `SphericalHankelH1[2,x]` operator (not a port), the F1 X−X de-taut held with no surviving definitional echo (the headline `m̂₀²χ_Q N_Q=1` is built from the observable odd residual in both engines), and reconciliation is 0-misaligned. Since no source changed, `material_change: false`. Verdict: verified.
