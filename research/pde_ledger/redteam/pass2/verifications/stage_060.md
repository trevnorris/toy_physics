---
unit_id: 060
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T12:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 060

## Per-finding outcomes

### F1 — stale_output (unambiguous self-label; number-only)

**Classification:** resolved

**What changed:**
`scripts/...sympy_audit.py:3` docstring `Stage 43 SymPy audit ...` → `Stage 60 SymPy audit ...` (diff hunk lines 18-25). Both committed `.txt` outputs were regenerated.

**Assessment:**
Correct and number-only. The diff touches exactly one source line (py:3); the banner format is preserved (no "43" remains, no spurious padding). Confirmed by sed: py:3 reads "Stage 60"; py:22/py:156 banners still read "STAGE 60" (not padded, per directive). The DEFERRED cross-ref at py:159 ("Stage-39 exponential family", a cross-reference to stage 056) was correctly LEFT UNTOUCHED — sed confirms it still reads "Stage-39", and it is absent from the diff. Output `.txt` mtimes (1780683732) postdate both script mtimes (1780682731/1780682732), so the saved transcripts are fresh. The refreshed SymPy output banner reads "STAGE 60 — ENTROPIC SOURCE MICROCLOSURE" with no "43". No result value changed.

### F2 — tautological_check (Mathematica side; the genuine math fix)

**Classification:** resolved

**What changed:**
At `wl:140` the old `expectZero["Xi_micro - Lambda^2 L^2/(Theta T_X)", xiMicro - lambdaPhi^2*ell^2/(theta*tX)]` (expr − its own wl:132 definition ≡ 0) was replaced with `expectZero["Xi_micro chi-route equals D/M-route", xiMicroFromChi - xiMicroFromDM]` (diff hunk lines 5-12).

**Assessment:**
The new check is genuinely NON-tautological. Neither operand references the literal wl:132 definition `xiMicro = lambdaPhi^2*ell^2/(theta*tX)`. Instead:
- `xiMicroFromChi` (wl:133) is built from the susceptibility relation `chiSigma*lambdaPhi^2*ell^2/tX /. chiSigma -> 1/theta`.
- `xiMicroFromDM` (wl:136) is built from the Einstein relation `mSigma*lambdaPhi^2*ell^2/(dSigma*tX) /. dSigma -> mSigma*theta`, where `mSigma` must cancel (mSigma/(mSigma·theta)=1/theta).
These are two independently-constructed routes; their difference is 0 only because both physical relations (`chiSigma=1/theta` and the Einstein relation `dSigma=mSigma·theta` with mSigma cancellation) are correctly stated. Mis-stating either would make the residual nonzero, so the assertion can genuinely fail. The old `expr − expr` form is fully removed (grep for `xiMicro - lambdaPhi^2*ell^2/(theta*tX)` returns NONE in the .wl). The `Xi_micro` deliverable remains verified: by this new cross-route check plus the retained wl:134-138 chi/D-M consistency checks and wl:141-142 susceptibility/phenomenological-form checks. The headline RESULT value is UNCHANGED — math output line 37 `Xi_micro = (ell^2*lambdaPhi^2)/(theta*tX)`, identical to the audit-recorded value and to the SymPy value `L**2*Lambda_phi**2/(T_X*Theta)` — so material_change is false. The fix is exactly as the directive's F2 block described; no collateral edits (the susceptibility/phenomenological lines below were untouched).

## Exec log assessment

**SymPy:** exit=0. Banner "STAGE 60 — ENTROPIC SOURCE MICROCLOSURE" (no "43"); `Xi_micro = L**2*Lambda_phi**2/(T_X*Theta)`; all residuals 0 and `PASS: dissipation density nonnegative`. (SymPy was untouched by F2 except its banner already correct; its `Xi_micro` line is its own pre-existing derived check.)

**Mathematica:** exit=0. New line present and passing:
```
Xi_micro chi-route equals D/M-route = 0
PASS: Xi_micro chi-route equals D/M-route
```
plus retained `xiMicro consistency via chi substitution`, `... via D/M substitution`, `Xi_micro susceptibility form`, `Xi_micro phenomenological form` all = 0 / PASS. `Xi_micro = (ell^2*lambdaPhi^2)/(theta*tX)` unchanged. Final line "Stage 060 Mathematica audit passed."

**Output freshness:** Confirmed. `.txt` mtimes (1780683732) are newer than both script mtimes (py 1780682731, wl 1780682732). Outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false. The F2 edit only removes a tautology and replaces it with an independent-route consistency assertion over the SAME already-correct value `Xi_micro = Lambda_phi^2 L^2/(Theta T_X)`. No verified RESULT value moved; both engines still agree on `Xi_micro`. F1 is label-only. No downstream unit is affected by a value change.

## Side observations (non-blocking)

None blocking. The new check is a cross-route equality between the two substitution-derived expressions rather than tying back to the Pe×support-drop route literally suggested in the report's F2 prose; the directive explicitly authorized the implementer to design the route and accepted any non-tautological independent construction, so this satisfies the acceptance criteria. The deliverable is additionally derived independently on the SymPy side.

## Verdict justification

Both findings are resolved. F1 is a clean number-only self-label fix (py:3) with the deferred cross-ref (py:159) correctly preserved and both outputs refreshed. F2, the batch's only genuine math fix, removes a true `expr − expr` tautology and replaces it with a real chi-route-vs-D/M-route equality whose two sides are independently constructed from distinct physical relations and could fail if either were mis-stated; the `Xi_micro` deliverable stays verified by this plus the retained checks, the value is unchanged, and both engines exit 0 with all PASS lines and fresh outputs. No regressions in the diff. Verdict: verified, material_change false.
