---
unit_id: 056
batch: III.2
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T12:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 056

## Per-finding outcomes

### F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Classification:** resolved

**What changed:**
The captured diff (`exec_logs/stage_056_diff.patch`) shows exactly two line edits to `scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py`:
- `py:3` docstring: `Moving-Throat PDE — Stage 39 SymPy audit.` → `Moving-Throat PDE — Stage 56 SymPy audit.`
- `py:97`: `print("\nStage 39 audit passed.")` → `print("\nStage 56 audit passed.")`

Both verified against the live source (read post-fix): py:3 reads "Stage 56 SymPy audit.", py:97 reads "Stage 56 audit passed."

**Assessment:**
The change exactly matches the directive's required change (number-only, 2-digit format preserved). The diff contains zero math/assertion/symbol/numeric edits — only the two label strings and the refreshed `.txt`. Confirmed all the must-not-touch / leave-as-is points hold:
- `py:7` cross-reference "reproduces the Stage-36 boost formula" (cross-ref to stage 053, 36+17=53) is LEFT UNTOUCHED — verified in live source, still reads "Stage-36". It is correctly absent from the diff.
- `py:31` banner "STAGE 56 — ..." was already canonical and was NOT padded — verified unchanged.

Stripping the digits from the two edited labels makes the file byte-identical to HEAD apart from the number — pure label normalization, not a semantic change. Non-tautology of the underlying assertions is unaffected (no assertion touched); the auditor's A1–A7/A8–A14 inventory carries forward intact.

## Exec log assessment

**SymPy:** exit=0. Notable lines (`exec_logs/stage_056_sympy.log`):
- L8 `STAGE 56 — TRANSPORT ORIGIN OF THE SOURCE-SHAPE ASYMMETRY` (banner now 056-aligned)
- L10–L18 all residuals 0: `zero-flux transport residual = 0`, `normalization ... - L = 0`, `Omega_Pe - expected formula = 0`, `dOmega/dPe - Cov/I_W = 0`
- L24 `Stage 56 audit passed.` (closing label fixed)
- L25 `# exit_code: 0`

**Mathematica:** exit=0. Notable lines (`exec_logs/stage_056_mathematica.log`):
- L8 `STAGE 056 — ...` (already canonical; .wl untouched, as expected)
- L11/L16/L18/L26/L28/L30/L39/L41 all `PASS:` lines present; residuals 0
- L43 `Stage 056 Mathematica audit passed.`; L44 `# exit_code: 0`
- The `Limit::alimv` warnings (L20/22/32) are benign (assumption-on-limit-variable ignored), already noted by the auditor; limits still computed correctly.

Engines agree on every result (I_W, I_Pe, Omega_Pe, endpoints 1 and π/2, small/large series, covariance residual 0) — unchanged from the pre-fix transcripts the auditor inventoried, confirming label-only.

**Output freshness:** confirmed. Both `.txt` outputs have mtime 2026-06-05 12:22, newer than the SymPy script (11:56) and the Mathematica script (06-03 15:59). The committed sympy `.txt` now reads "STAGE 56" (line 3) and "Stage 56 audit passed." (line 19), replacing the prior "Stage 39" banners; the mathematica `.txt` reads "STAGE 056" (line 8). All numeric result lines identical to the auditor's transcript.

## Material-change assessment

`material_change`: false.

The edit touches only two docstring/print label strings plus the regenerated transcript. No derived quantity, assertion, symbol, or numeric expression changed; `Pe`, `Ω_Pe`, `I_W`, the asymptotics, and the covariance identity are all byte-identical apart from the label number. No downstream unit can depend on a print label. No `upstream_stale` propagation warranted on substance.

## Side observations (non-blocking)

None. The deferred SCRIPT/OUTPUT-band cross-ref handling (py:7 "Stage-36") is consistent with the dedicated numbering plan and was correctly excluded.

## Verdict justification

The sole finding (F1, low-severity stale_output / self-label numbering drift) is fully resolved: the two stale "Stage 39" self-labels are now "Stage 56" in the exact format requested, the deferred cross-reference at py:7 was correctly left untouched, the banner was not over-padded, and no math changed. Both engines re-ran to exit 0 with all assertions passing and all result values unchanged, and both `.txt` outputs were refreshed (mtimes newer than scripts) with the corrected banners. material_change is false. Verdict: verified.
