---
unit_id: 135
batch: IV.4
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-06T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 135

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
Codex deleted the single tautological assertion line
`expectApprox["closure residual", residual, 0, 10^-14];` from
`mathematica/moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl`.
The captured diff (`stage_135_diff.patch`) is a one-line removal at the former
`.wl:78` and touches nothing else. No replacement assertion was added.

**Assessment:**
The edit exactly matches the directive's "required change" (delete the line, do
not replace it, leave everything else as-is). Verified against the current
file state:

1. **No "closure residual" assertion remains.** Grep over the `.wl` for
   `expect*` returns lines 56, 67, 74, 75, 76, 77, 78 — the names are
   `outlet-consistent reduction` (56), `0 < S_q(Pi_star) < 1` (67), the four
   numeric checks `Sigma_m^*`/`M_s^*`/`M_q^*`/`mixed-lane` (74-77), and
   `3 Sigma_m^* < Pi_* < 4 Sigma_m^*` (78, now an `expectTrue`, not the deleted
   `expectApprox`). None is "closure residual". The tautology is gone.
2. **Residual PRINT retained.** `.wl:72`
   `Print["Pi_star - Sigma_star*(4 - S_star) = ", residual];` is present and
   unchanged; `residual` is still defined at `.wl:64`
   (`residual = N[piStar - sigmaStar*(4 - sStar), 30];`). The value is still
   emitted as a sanity probe.
3. **Refreshed Mathematica output has no `PASS: closure residual` line.** The
   transcript contains exactly 7 PASS lines (output lines 12, 15, 22, 24, 26,
   28, 30): outlet-consistent reduction, `0<S_q<1`, the four numeric checks, and
   the `3σ*<Π_*<4σ*` bound. The residual still appears as a bare value at output
   line 20 (`Pi_star - Sigma_star*(4 - S_star) = 0...28.68818091139617`), i.e.
   printed, not asserted.
4. **Parity with SymPy is restored.** `.py:84-86` keeps the print with the
   comment "Numerical sanity probe (was the original tautological closure
   residual)" and no assert; the `.wl` now matches that exact posture.

The removed assertion was provably tautological: `sigmaStar` is solved at
`.wl:60` from `piStar == sigmaVar*(4 - sStar)`, so
`piStar - sigmaStar*(4 - sStar) ≡ 0` by construction. Removing it (rather than
re-anchoring) is correct — the kernel value is independently exercised by the
numeric anchors at lines 74-77 and the range checks at 67/78. No collateral
edits, no regressions.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Sigma_m^*, M_s^*, M_q^* anchored to notes values within tolerance.`
- `M_q^* * S_q(Pi_*) = -0.297111597463198745446779533971`
- `Pi_star - Sigma_star*(4 - S_star) = -7.22456073765234902012354489563e-17`
  (printed-only residual, no assert — as intended).

**Mathematica:** exit=0. Notable lines:
- `PASS: outlet-consistent reduction`
- `PASS: Sigma_m^* numeric check` / `PASS: M_s^* numeric check` /
  `PASS: M_q^* numeric check` / `PASS: mixed-lane correction numeric check`
- `PASS: 3 Sigma_m^* < Pi_* < 4 Sigma_m^*`
- `Pi_star - Sigma_star*(4 - S_star) = 0``28.68818091139617` (printed, no PASS).
- No `PASS: closure residual` line anywhere. `# exit_code: 0`.

All remaining `expect*` checks PASS; both engines exit 0.

**Output freshness:** confirmed. Both saved `.txt` outputs have mtime
2026-06-06 21:47:39, newer than the edited `.wl` (mtime 2026-06-06 17:14:34) and
the `.py` (2026-05-27). Outputs were regenerated post-fix.

## Material-change assessment

`material_change`: false.

The edit only removed an assertion that verified nothing (a `0 == 0` identity by
construction). No deliverable value changed: the residual is still computed and
printed, and every numeric anchor (`Σ_m^*`, `M_s^*`, `M_q^*`, mixed-lane
correction, kernel `S_q(Π_*)`, the imported `Π_*`) is byte-for-byte identical to
the pre-fix outputs. Nothing downstream of unit 135 can depend on a removed
no-op check, so no downstream re-audit is warranted on account of this change.

## Side observations (non-blocking)

None. The diff is a clean single-line deletion with no spillover.

## Verdict justification

The sole finding (F1, tautological `closure residual` assertion in the
Mathematica script) is fully resolved: Codex deleted exactly the named line,
added no replacement, retained the residual print, brought the `.wl` into the
parity the SymPy side already established, and the refreshed transcript no longer
emits a `PASS: closure residual` line while all remaining `expect*` checks PASS
and the script exits 0. No regressions in the diff, outputs are fresh, no
deliverable value changed. Verdict: `verified`.
