---
unit_id: 194
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-09T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 194

## Per-finding outcomes

### F1 — paper_misalignment (paper_missing_script_claim — card-text lag, paper-side)

**Classification:** resolved (USER-DEFERRED, non-blocking)

**What changed:**
No script change. The finding is a paper-side prose lag — the stage card's
`\stagefield{Verification}` still reads "Mathematica audit: none yet" while a
complete, passing `.wl` exists. The audit report correctly routed this to the
USER (Codex does not edit `paper/`), and the orchestrator handoff confirms the
USER has DEFERRED it to the paper-cleanup tracker. No edit is in scope for this
red-team unit.

**Assessment:**
Correct disposition. This is explicitly the paper-cleanup class flagged in the
V.3 context; it is non-blocking for script verification (scripts-only scope).
The underlying technical claim — that a genuine second engine exists and passes —
is independently confirmed below from the fresh `.wl` output (all `PASS`). Holds.

### F2 — stale_output

**Classification:** resolved

**What changed:**
No source edit. The orchestrator re-ran the SymPy script and recommitted the
refreshed `.txt`. The committed
`scripts/output/moving_throat_pde_stage194_..._sympy_audit.txt` now carries banner
`STAGE 194 —` (line 3) and `STAGE 194 LEDGER` (line 143); the stale `STAGE 177`
banner is gone. File mtime is Jun 9 14:06, newer than the `.py` (Jun 3 15:59).

**Assessment:**
Refresh landed correctly. Banner now matches the current script. Every SymPy
residual prints `= 0` (lines 27–29, 50–52, 65, 92, 103–104, 110–111, 138–140) —
no result drifted from the prior (stale-bannered) capture, consistent with the
auditor's note that only the banner, not the result values, was stale. Resolved.

## Exec log assessment

**SymPy:** exit=0 (re-run by orchestrator; refreshed `.txt` is the artifact).
Notable lines:
- `STAGE 194 — EXACT OUTGOING l=2 DTN FINGERPRINT, …` (line 3) — canonical banner present.
- `Lambda_out series - expected = 0`, `Yhat_out series - expected = 0`, `static DtN slot + 3 = 0` (27–29).
- `sigma_Q^can - 4 a^5/(27 c_s^5) = 0` (50), `chi_Q - 1 deformation law = 0` (111), `Gammabar_5 - 2 G/(5 c^5) = 0` (140).

**Mathematica:** exit=0 (fresh; banner `STAGE 194`, line 3). Notable lines:
- `h_2^(1)(x) = (I*E^(I*x)*(-3 + (3*I)*x + x^2))/x^3` (9) — native `SphericalHankelH1` representation, a distinct route from SymPy's trig-closed-form `j2+i*y2`.
- `positive pole from even matching = (3*soundSpeed)/(2*radius)`, `positive sigma from odd normalization = (4*radius^5)/(27*soundSpeed^5)` (22–23) — `Solve`-derived, vs SymPy pin-and-verify.
- Every check carries an explicit `PASS:` (lines 13–84); no `FAIL`.

**Output freshness:** Confirmed. Both `.txt` files are mtime Jun 9 14:06, newer
than the SymPy `.py` (Jun 3 15:59) and the `.wl` (Jun 1 11:42). Both banners read
`STAGE 194`. The stale `STAGE 177` SymPy capture is fully replaced.

## Material-change assessment

`material_change`: false. No source-code change occurred — F2 was an output
refresh (banner-only correction; all result values unchanged) and F1 is a
USER-deferred paper-side prose lag with no script edit. No derived result that
downstream units depend on was altered. No narrow re-audit of units > 194 is
warranted on account of this unit.

## Side observations (non-blocking)

- The audit's independence verdict is corroborated directly by the fresh `.wl`
  output: native special-function Hankel route (line 9) + `Solve`-with-positive-root
  pole/σ extraction (lines 22–23) vs SymPy's explicit-trig + pin-and-verify. Genuine
  independent re-derivation, not transliteration.
- Value-reconciliation: the audit checked 14 deliverable values with 0 misalignments;
  all fresh-output values (Λ₂^out, Ŷ₂^out, Ω_Q, σ_Q^can, χ_Q=1, Σ₂/Σ₄, χ_Q map,
  invariant tuple K̄₀/K̄₂/K̄₄/Γ̄₅) match the report's reconciliation table. Nothing to flag.

## Verdict justification

Both findings are resolved. F2's stale-output refresh landed: the SymPy `.txt`
now carries the canonical `STAGE 194` banner with every check `= 0` and an mtime
newer than its script; the Mathematica `.txt` is likewise fresh with every check
`PASS`. F1 is a paper-side card-text lag, correctly routed to the USER and
DEFERRED to the paper-cleanup tracker — non-blocking under the scripts-only
verifier scope, and the substantive claim it concerns (a real, passing second
engine) is independently confirmed from the `.wl` output. The disposition holds:
the `.wl` is a genuinely independent re-derivation (native `SphericalHankelH1` +
`Solve`-based pole/σ/K̄₀ extraction), engines agree, and reconciliation shows 0
misalignments. No source change, no material/downstream impact. Verdict: verified.
