---
unit_id: 168
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-08T20:05:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 168

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex edited only the `.wl` (`mathematica/moving_throat_pde_stage168_off_bundle_slippage_mathematica_audit.wl`), per the diff at `exec_logs/stage_168_diff.patch`:

- wl:71-74 — the hand-typed literal `epsPerp = g*epsT + (g + bCoeff)*epsv + cCoeff*epsL;` is removed. In its place the three slippage weights are DERIVED:
  `wT = Coefficient[Expand[-deltaPerpSlip], epsT];` (and `wv`, `wL`), then `epsPerp = wT*epsT + wv*epsv + wL*epsL;`. `deltaPerpSlip` is itself built upstream (wl:64-70) by substituting the Stage-165 transport laws + slippages into `deltaPerp`.
- wl:75-76 — new `expectZero["epsPerp weights match boxed form (g, g+1/(2s), 2g+3/(4s))", (wT - g)^2 + (wv - (g + bCoeff))^2 + (wL - (2*g + 3/(4*s)))^2];`.
- wl:78 — pre-existing `expectZero["delta_perp + eps_perp", deltaPerpSlip + epsPerp]` left unchanged (now compares against the derived `epsPerp`).
- wl:105-106 — `rNum = SetPrecision[1.77799353547498, 30];` replaced by `rExact = Sqrt[4107 - 100*Pi^2]/(10*Pi);` then `rNum = N[rExact, 30];`. `gNum` left as decimal per directive (no closed form available).

**Assessment:**
The edit matches the directive exactly; no collateral changes (mouth-bias, outlet, carry-forward blocks untouched). The `.py` is NOT in the diff — untouched, as required.

(a) **Derivation genuinely independent (not a relabel):** `wT/wv/wL` are read off `-deltaPerpSlip`, which was assembled by substitution into `deltaPerp` — a path distinct from the `.py`'s hand-typed `eps_perp = g*epsT + (g+B)*epsv + C*epsL` (py:77). The `delta_perp + eps_perp` check is therefore no longer transcription-vs-transcription: it compares the substituted coordinate against a coefficient-read of that same coordinate, so a wrong upstream `deltaPerp` coefficient or transport law would surface (it would shift `wT/wv/wL` and break the check) rather than being silently re-typed identically into both engines.

(b) **New `expectZero` non-tautological:** the sum-of-squares confirms the DERIVED weights equal the boxed Family-1 forms `(g, g+1/(2s), 2g+3/(4s))`. The RHS targets (`g`, `g+bCoeff`, `2g+3/(4s)`) are typed independently of how `wT/wv/wL` are computed (coefficient extraction), so it is a real cross-check, not `X-X`. Log line 24-25 shows residual 0 / PASS.

(c) **Decimal de-transcription landed:** the bare `1.77799353547498` no longer appears as the radius source; it is now `Sqrt[4107 - 100*Pi^2]/(10*Pi)`. Numeric coefficient prints are unchanged (log lines 51-56 match the pre-fix report values to 20 digits).

All seven original `expectZero` checks still PASS, plus the new one (8 PASS total), exit 0.

## Exec log assessment

**SymPy:** exit=n/a. No sympy exec log exists for this unit (`ls` confirms only `stage_168_mathematica.log` and the diff). The F1 re-author touched only the `.wl`; the `.py` was not modified, so no sympy re-run was warranted.

**Mathematica:** exit=0. Notable lines:
- L19 `PASS: bundle tangency (delta_perp on exact lower branch)`
- L24-25 `epsPerp weights match boxed form (g, g+1/(2s), 2g+3/(4s)) = 0` / `PASS`
- L26 `delta_perp with slippages = -2*epsL*g - epsT*g - epsv*g - (3*epsL)/(4*s) - epsv/(2*s)` then L28 `PASS: delta_perp + eps_perp`
- L34, L40, L42, L44, L46 — `PASS` on deltaPi + all four outlet defects
- L65 `Stage 168 Mathematica audit passed.` / L66 `exit_code: 0`

**Output freshness:** `.wl` mtime 1780957970 < committed output `.txt` mtime 1780960353 — output is newer than the post-fix script, so it was regenerated after the edit. Fresh.

## Material-change assessment

`material_change`: false.

The edit is method-only. `epsPerp` is now derived but reduces to the identical boxed weights `(g, g+1/(2s), 2g+3/(4s))`; the radius is computed from a closed form that evaluates to the same `1.77799353547498...`. No numeric coefficient or symbolic identity changed (all residuals still 0; printed coefficients unchanged to 20 digits). No downstream unit can observe a different result.

## Side observations (non-blocking)

- The directive's `## Applied: F1` block is present and accurate (files_changed, one-sentence summary, deviation: none).
- The `gNum` decimal `0.758035078944663` was deliberately left as a literal per the directive (no closed form available); not a defect.

## Verdict justification

The single low-severity transliteration finding is fully resolved. The Mathematica path to `ε_⊥` now genuinely derives the slippage weights by coefficient extraction from the substituted normal coordinate, breaking the verbatim-port symmetry with the `.py`; the new boxed-weight `expectZero` is a real, non-tautological independent cross-check that passes; the decimal radius is de-transcribed to its Family-1 closed form; all seven original checks still pass; the `.py` is untouched; and the script exits 0 with a freshly regenerated output. No regressions in the diff, no material change. Verdict: verified.
