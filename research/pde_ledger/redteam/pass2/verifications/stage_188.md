---
unit_id: 188
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 188

## Per-finding outcomes

### F1 — stale_output

**Classification:** resolved

**What changed:**
No source edit was required (the finding explicitly states the `.py` already prints `STAGE 188`; the stale banner lived only in the committed `.txt` transcript). The orchestrator's standard independent re-run regenerated both committed outputs. The empty `stage_188_diff.patch` (0 bytes) confirms no script source was touched this pass, exactly as expected for a transcript-refresh-only finding.

**Assessment:**
The refresh landed correctly.
- SymPy output `…_sympy_audit.txt:3` now reads `STAGE 188 — BRANCH-OBSERVABLE COMPLETION AND THE EXACT FIRST-ORDER OBSERVABLE COMPILER`, and the ledger banner at `:272` reads `STAGE 188 LEDGER`. The pre-renumber `STAGE 171` banner is fully gone (a `grep -E 'STAGE [0-9]+'` over both files returns only `STAGE 188` hits — no surviving stale label anywhere).
- Mathematica output `:3` reads `STAGE 188 -- BRANCH-OBSERVABLE COMPILER (MATHEMATICA AUDIT)` (was already current; still current).
- Freshness: both `.txt` mtimes are `2026-06-09 14:03:46`, newer than the `.py` (`2026-06-03 15:59`) and the `.wl` (`2026-06-01 11:20`). Output now post-dates source. The exact staleness condition the finding flagged (`.txt` older than `.py`) is cleared.
- All check lines pass: every SymPy residual prints `= 0` / a zero matrix / zero vector (sections I–VII); the Mathematica transcript prints 26 `PASS:` lines and the closing `Stage 188 Mathematica audit passed.` No `FAIL`/`error`/`traceback` in either output or either exec log. Both exec logs end `# exit_code: 0`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `C_* - 1/C_tr,* = 0` and `A_tr,* - B_* C_tr,* = 0` (the two coherent-branch coefficient identities, including the hinge identity).
- `factorized compiler - expected compiler =` → zero matrix; `det(C_obs->def) = epsilonetaₛₜₐᵣ/(epsilonetaₛₜₐᵣ - 1)`.
- `observable roundtrip - Delta_obs^(1) =` → zero vector; bijection-on-generic-packet residuals (sec VII) print zero.

**Mathematica:** exit=0. Notable lines:
- `PASS: observable-to-quotient compiler agrees with SymPy`, `PASS: quotient-to-defect compiler agrees with SymPy`, `PASS: direct compiler agrees with SymPy expected compiler` (cross-engine consistency layer, applied after independent Solve/Series/D derivation).
- `PASS: C_def->obs * C_obs->def - I` and `PASS: observable roundtrip - Delta_obs^(1)` (inverse + round-trip).
- `PASS: 1/det(C_obs->def) well-defined (nonzero det)` and the def/quot bijection PASSes (sec VII iff is exercised on the generic packet, not the zero vector).

**Output freshness:** confirmed. Both `.txt` mtimes (Jun 9 14:03) post-date both the `.py` (Jun 3 15:59) and `.wl` (Jun 1 11:20). Regenerated post-refresh.

## Material-change assessment

`material_change`: false. No source-code change occurred (empty diff patch); the only change is a refreshed transcript banner. No derived result changed, so no downstream unit's inputs are affected by this stage's refresh.

## Side observations (non-blocking)

- There is no `directives/stage_188.md` with `## Applied:` / `## Blocked:` blocks for this unit, which is consistent with a refresh-only finding that needs no Codex source edit. Nothing to verify there.
- On a fresh read, the audit's independence assessment holds: the `.wl` derives each compiler from finite branch-ratio first-log-drifts via `Series`/`Log`/`Solve` + `Outer[D,…]` and only cross-checks the hand-typed SymPy reference matrices afterward (the `sympy*` arrays are a post-hoc consistency layer, not the derivation). The clearest tell — the third defect row derived from `SeriesCoefficient[Log[(1 - eta*Exp[small*etaObs])/(1 - eta)], …]` rather than the `.py`'s pre-collapsed `-epsetas/(1-epsetas)` coefficient — is borne out by the transcript (`delta ln(1 - epseta) = (eta*etaObs)/(-1 + eta)` at mathematica `:83`). Genuinely independent, not a transliteration.
- Value reconciliation in the report (16 deliverable values, 0 misaligned) is consistent with the refreshed transcripts; the printed final forms (`C_tr,*`, `B_*`, `A_tr,*`, both determinants, the triangular compilers) match across both engines and the report's table.

## Verdict justification

The single low-severity `stale_output` finding is resolved: the orchestrator's independent re-run refreshed both committed transcripts, the SymPy banner now reads `STAGE 188` (the pre-renumber `STAGE 171` label is gone everywhere), output mtimes post-date the scripts, every check line reads `= 0` / `PASS` / zero-matrix, and both exec logs exit 0 with no failures. The empty diff confirms no source change, matching the finding's "no source edit required" note. The audit's disposition holds on a fresh read — the `.wl` is genuinely independent (finite-ratio Series/Solve/autodiff derivation, SymPy cross-check applied afterward) and value reconciliation shows 0 misalignments. Verdict: verified.
