---
unit_id: 195
batch: V.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-01T12:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 195

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py:87-93`. The old self-echo `expect_zero("odd closure factorization", odd_closure - (mhat0**2*chi_Q*N_Q - 1))` (where `odd_closure := simplify(mhat0**2*chi_Q*N_Q - 1)`) is replaced by:

```python
odd_condition_residual = sp.simplify(
    (mhat0**2 * Gamma5 - Gamma5_target).subs(P0, N_Q * P0_target)
)
expect_zero(
    "observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1)",
    odd_condition_residual - Gamma5_target * (mhat0**2 * chi_Q * N_Q - 1),
)
```

The `odd_closure` symbol (line 75) is retained but only for the display `pprint` at line 82 and to seed `NQ_from_odd` (line 76) — it is no longer the load-bearing assertion.

**Assessment:**
Genuinely non-tautological. The residual is now built from the geometric Gamma5 form `chi_Q*a**5*P0/(27*c_s**5)` (line 48, an *independent* expression from the RHS), substituting `P0 = N_Q*P0_target`, and asserted equal to `Gamma5_target*(mhat0**2*chi_Q*N_Q - 1)`. It touches both `Gamma5` and `Gamma5_target`, which the old check never did. Hand-check: `mhat0^2 * chi_Q*a^5*(N_Q*P0_target)/(27*c_s^5)` with `P0_target = 54 G c_s^5/(5 a^5 c^5)` cancels `a^5`/`c_s^5` to `mhat0^2*chi_Q*N_Q*(54 G)/(135 c^5) = mhat0^2*chi_Q*N_Q*2G/(5 c^5) = mhat0^2*chi_Q*N_Q*Gamma5_target`, so the residual is identically `Gamma5_target*(mhat0^2*chi_Q*N_Q - 1)` — but only because the carried `54/27/5 -> 2/5` normalization is correct. A hidden factor error in Gamma5 or P0_target would now produce a nonzero residual. No new X−X introduced. Output line 56 shows the new named check `= 0`.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Same file. Both definitional-echo `expect_zero` calls are deleted:
- The `expect_zero("N_Q definition", N_Q_def - P0 / P0_target)` is gone (section I now ends after the `Gamma5/Gamma5_target - chi_Q*N_Q` check at lines 63-66).
- The `expect_zero("Delta_norm - P0_target*(mhat0^2 N_Q - 1)", Delta_norm_NQ - ...)` is gone (section II now goes straight from the display `pprint`s to the F1 derivation at line 87).

**Assessment:**
Correct minimal deletion. `N_Q_def` (line 52), `Delta_norm` (line 73), and `Delta_norm_NQ` (line 74) definitions are all intact. `Delta_norm_NQ` is still consumed downstream at line 77 (`Delta_norm_from_odd`) and line 109 (`Delta_norm_pt`), so removing the echo asserts nothing downstream loses an input. Output (full transcript, lines 1-150) contains neither "N_Q definition = 0" nor "Delta_norm - P0_target*(mhat0^2 N_Q - 1) = 0". Confirmed gone.

### F3 — missing_mathematica

**Classification:** resolved

**What changed:**
New file `mathematica/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_mathematica_audit.wl` (229 lines), re-verifying all five blocks and exiting 0 with every line PASS.

**Assessment:**
Genuine independent route, not a transliteration. Distinct evidence:

1. **Section IV derives chi_Q from the actual operator rather than asserting it.** SymPy hardcodes `chi_from_def = 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)` (line 130). The .wl instead builds it from the outgoing l=2 Hankel operator: `lambdaOut = FunctionExpand[x*D[SphericalHankelH1[2,x],x]/SphericalHankelH1[2,x]]` (line 136), `Series` to order 5 (line 137), deformed DtN with `scaleS*(lambdaWindow /. x -> betaStretch*x) + sigma0 + sigma2 x^2 + sigma4 x^4 + I sigma5 x^5` (lines 138-141), normalizes via `Coefficient[deformedDtn,x,0]/deformedDtn` series (lines 142-144), `Solve`s for `{sigma2,sigma4}` from the canonical-even conditions `Coefficient==1/9, ==4/81` (lines 146-153, with a uniqueness guard at 154), then extracts `chiFromDtn = Coefficient[normalizedDtn,x,5]/(I/27)` (lines 156-159). Output line 46 confirms this recovers `(3*(betaStretch^5*scaleS + 9*sigma5))/(3*scaleS - sigma0)`. This is a materially different decomposition using native `SphericalHankelH1`/`Series`/`Coefficient`/`Solve`.

2. **The headline factorization mirrors the CORRECTED F1, not the deleted X−X.** Lines 80-101: `oddObservableResidual = (mhat0^2*gammaByGeometry - gamma5Target) /. p0Bar -> nQ*p0Target` (built from the geometric `gammaByGeometry = chiQ*radius^5*p0Bar/(27*soundSpeed^5)`), asserted `== gamma5Target*(mhat0^2*chiQ*nQ - 1)`. Output line 21 shows the intermediate residual `(2*bigG*(-1 + chiQ*mhat0^2*nQ))/(5*lightSpeed^5)` — i.e. the factorization is exhibited symbolically, then PASS at line 25. Also uses `Solve[mhat0^2*chiQ*nQ == 1, nQ, Reals]` (line 86) for `nQFromOdd` rather than hand-defining the reciprocal.

3. **Different naming/order throughout** (radius/soundSpeed/lightSpeed/p0Bar vs a/c_s/c/P0; gammaByGeometry/gammaByPole). Linearization via `Series[...,{tau,0,1}]` + `Coefficient` (lines 177-193).

No vacuous or X−X check found: every `expectZero` subtracts an independently-routed result from a target form.

## Exec log assessment

**SymPy:** exit=0. Notable lines from the saved output:
- `observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1) = 0` (the F1 fix is present and passing)
- No "N_Q definition = 0" and no "Delta_norm - P0_target*(mhat0^2 N_Q - 1) = 0" lines (F2 echoes removed)
- `linearized N_Q - 1 = 0`, `Delta_norm(canonical) = 0`, `P0(canonical source-map-reduced) - P0_target = 0`

**Mathematica:** exit=0 (`Exit[0]` at line 228, all PASS lines present). Notable lines:
- `PASS: observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1)` (output line 25)
- `chi_Q from native DtN deformation = (3*(betaStretch^5*scaleS + 9*sigma5))/(3*scaleS - sigma0)` (output line 46)
- `linearized N_Q - 1 = -5*epsBeta - dSigma0/(3*scaleS) - (9*dSigma5)/scaleS` (output line 53)
- `PASS: P0(canonical source-map-reduced) - P0_target` (output line 70)

**Cross-engine agreement:** Yes. The two engines agree on every load-bearing conclusion: the Gamma5/Gamma5_target = chi_Q N_Q ratio; the factorized odd condition `Gamma5_target*(mhat0^2 chi_Q N_Q - 1)`; `Delta_norm = P0_target(1/chi_Q - 1)`; `chi_Q = 3(Sβ^5+9Σ5)/(3S-Σ0)`; the Delta_norm^pt deformation form; the linearized `N_Q-1 = -5ε_β - dΣ0/(3S) - 9dΣ5/S`; and the canonical closure `chi_Q=N_Q=1, Delta_norm=0`. The .wl reaches these via the independent Hankel-operator route, so agreement is substantive cross-engine confirmation.

**Output freshness:** Per the prompt, the saved `.txt` outputs show exit 0 / all-PASS and were re-generated post-fix (the SymPy output contains the new F1 line and lacks the deleted F2 echoes; the Mathematica output exists only because the .wl is newly created). Treated as fresh.

## Material-change assessment

`material_change`: false.

F1 changes which assertion verifies the headline, but the verified *result* (`m̂₀²χ_Q N_Q = 1` and the downstream collapses) is unchanged — no derived value that units 196/197 consume was altered, only the rigor of the check. F2 removed two echoes (no result change). F3 added a second engine that confirms existing results. No downstream re-audit is warranted on material grounds.

## Side observations (non-blocking)

- The SymPy banner still reads "STAGE 178" and section IV says "Stage 194 deformation algebra"; the auditor already noted these as cosmetic print strings (not assertions). Not gating.
- `odd_closure` (SymPy line 75) and `oddClosureEquation` (.wl line 84) survive only as display/seed; harmless.

## Verdict justification

All three findings are resolved. F1's headline factorization is now derived from the observable odd condition `m̂₀²Γ̄₅ = Γ̄₅^target` via the geometric Gamma5 form and the carried normalization, passing only because that normalization is correct — not by construction — and it touches `Gamma5`/`Gamma5_target` as required. F2's two definitional echoes are deleted with downstream definitions intact. The new `.wl` is a genuinely independent route (native SphericalHankelH1 + Series/Coefficient/Solve to *derive* chi_Q rather than assert it, and it mirrors the corrected, not the broken, headline). Both engines exit 0 and agree on every conclusion. Verdict: `verified`.
