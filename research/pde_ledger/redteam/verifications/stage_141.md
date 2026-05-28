---
unit_id: 141
batch: IV.5
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 0
findings_total: 0
material_change: false
---

# Verification — unit 141

## Per-finding outcomes

No findings in the original report (`findings_count: 0`, `verdict: clean`). Nothing to resolve.

## Exec log assessment

**SymPy:** n/a — no `.py` audit script exists for stage 141 (confirmed by directory listing: `scripts/moving_throat_pde_stage141_*.py` returns no matches). No exec log expected.

**Mathematica:** n/a — no `.wl` audit script exists for stage 141 (confirmed by directory listing: `mathematica/moving_throat_pde_stage141_*.wl` returns no matches). No exec log expected.

**Output freshness:** n/a — no scripts means no outputs to refresh.

## Status-only confirmation

- The auditor's report (`stage_141.md`) declared the unit clean under the status-only carve-out (`is_status_only_candidate: True` per manifest), citing the paper card's verbatim `Verification` field "SymPy audit: none yet.  Mathematica audit: none yet."
- I re-confirm via filesystem checks that no SymPy or Mathematica script has been added since the audit ran: neither `scripts/moving_throat_pde_stage141_*.py` nor `mathematica/moving_throat_pde_stage141_*.wl` exists.
- No directive was written (consistent with zero findings); `redteam/directives/stage_141.md` is absent.
- The carry-forward attribution chain documented in the auditor report (mouth-gain formulas attributed to Stages 188–191, plus the self-matched susceptibility closure block) is the upstream anchor for this ledger entry and was not affected by any edits in this verification pass — no edits were attempted.

## Material-change assessment

`material_change`: false.

This is a status-only ledger stage with no script edits and no new identities introduced. Downstream units are not affected because nothing was changed.

## Side observations (non-blocking)

None.

## Verdict justification

Nothing has changed since the auditor's clean verdict: no scripts have been added, no directive exists, and the paper card's "none yet" declaration remains consistent with the absent script files. The status-only carve-out applies as documented in the original report, the attribution chain to upstream stages (188–191 plus self-matched susceptibility closure) is intact, and there is nothing to rework. Verdict: `verified`.
