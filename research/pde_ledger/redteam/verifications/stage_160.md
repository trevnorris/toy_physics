---
unit_id: 160
batch: IV.6
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 160

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl:28-49` was rewritten. (a) Symbols renamed per directive: `rc→rStar`, `drc→deltaR`, `dk0→dKappa0`, `dg0→dGamma0`, `k0Star→kappa0Canon`, `g0Star→gamma0Canon`, `dkW→deltaKappaW`, `dgW→deltaGammaW`; `eps` removed entirely. (b) The `Series[..., {eps,0,1}]` + `Coefficient[..., eps, 1]` recipe was replaced by direct chain-rule total-differentials `deltaKappaW = Together[dKappa0/(1+rStar) - (kappa0Canon/(1+rStar)^2)*deltaR]` and similarly for `deltaGammaW` (lines 34-39). The identity, gate substitution, and both Print labels were retargeted to the new names. Tangential block (lines 51-64) untouched per directive item 4. SymPy script untouched (verified via mtime 1779945093 < .wl mtime 1779984284, and grep of new symbol names returns NONE in `.py`).

**Assessment:**
Directive followed exactly. The chain-rule derivation is structurally distinct from SymPy's `series(...).coeff(eps,1)` recipe — it computes the differential analytically at the canonical point without introducing a perturbation parameter. Forbidden substrings `Series[`, `Coefficient[`, `k0Star`, `g0Star`, `\bdk0\b`, `\bdg0\b`, `, eps,` all return NONE in the `.wl`. The two `expectZero` calls remain non-tautological: the residual `deltaGammaW - (1/3)*deltaKappaW - (dGamma0 - (1/3)*dKappa0)/(1+rStar)` is a nontrivial four-symbol cancellation (rStar, deltaR, dKappa0, dGamma0), and pure-scale residual is `(dKappa0/3 - dKappa0/3)/(1+rStar)` after substitution — both verified `= 0` in the exec log (lines 7-8, 10-11). No collateral edits.

## Exec log assessment

**SymPy:** exit=0 (script untouched; log re-run confirms no regression). Notable:
- `dκ_W = (-dr_c/3 + dκ0)/(r_c + 1)`
- `exact compensated-branch slippage identity = 0`
- `pure-scale harmlessness = 0`

**Mathematica:** exit=0. Notable:
- `dκ_W = -((deltaR - 3*dKappa0)/(3 + 3*rStar))` — algebraically identical to SymPy's `(-dr_c/3 + dκ0)/(r_c+1)` (multiply by 3/3).
- `dγ_W = -((deltaR - 9*dGamma0)/(9 + 9*rStar))` — matches SymPy's `(-dr_c/9 + dγ0)/(r_c+1)`.
- `exact compensated-branch slippage identity = 0` → PASS
- `pure-scale harmlessness = 0` → PASS

Engines still agree on the first-order pieces and both `expectZero` assertions; the agreement is now reached by genuinely different routes (chain rule vs. series-coefficient).

**Output freshness:** mathematica `.wl` mtime 1779984284 < output `.txt` mtime 1779989528 (delta +5244s); sympy `.py` mtime 1779945093 < output `.txt` mtime 1779989432. Both fresh post-edit.

## Material-change assessment

`material_change`: false. The edits change only the derivation mechanism inside the `.wl`; the asserted residuals, the boxed identity, the gate substitution, and all printed numerical/symbolic outputs are unchanged (same algebraic content, just a different path). No downstream unit depends on the internal symbol names or recipe shape of stage 160's Mathematica audit.

## Side observations (non-blocking)

- The tangential block still uses `rc` (not `rStar`) in its own `$Assumptions` (line 54) and `dgWTan` denominator (line 57). This is consistent with the directive's item 4 (do not modify the tangential block) and matches the SymPy script's local naming there, so it is not a regression — but it does mean the symbol `rc` reappears later in the file, just outside the load-bearing derivation block.
- Banner labels in both scripts still read "STAGE 160" correctly (no stale-stage-number artifact here — the original audit's mention of "STAGE 143" was elsewhere).

## Verdict justification

The single F1 finding is fully resolved: the load-bearing derivation block (lines 28-49) was restructured to an independent chain-rule route with all required symbol renames; forbidden substrings are confirmed absent; both `expectZero` calls still pass with non-tautological residuals; SymPy script and tangential block were left untouched per directive scope; output files are fresh and exit 0. No regressions, no new findings introduced. Verdict: `verified`.
