---
unit_id: 070
batch: III.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 070

## Per-finding outcomes

### F1 — tautological_check (de-tautologize the sech-moment anchor)

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage070_gnls_wall_shell_sympy_audit.py:57` — added, right after `If_sym = sp.integrate(chi_phi**2, (xi, -oo, oo))` (L56):
  `expect_zero("sech-profile moment I_f = 2/3", sp.simplify(If_sym - sp.Rational(2, 3)))`.
  The pre-existing `I1/J1 - 4πa²ℓ` ratio line (L64) and the "expected 2/3" print (L65) are left in place as documentation, exactly as the directive prescribed.
- Mathematica `mathematica/moving_throat_pde_stage070_gnls_wall_shell_mathematica_audit.wl:86–90` — corrected the wrong I_g print annotation `8/15`→`14/15` (L86, `N[8/15,30]`→`N[14/15,30]`), and added two numeric assertions on the already-computed `IfNum`/`IgNum` against the in-scope `tol = 10^-10` (L87–90):
  `If[Abs[IfNum - 2/3] < tol, pass[...], fail[...]]` and `If[Abs[IgNum - 14/15] < tol, pass[...], fail[...]]`.
  The existing structural `I_1/J_1` `expectZero` (L96–98) and the symbolic `kappa`/`W_wall`/`Xi` assertions (L57/65/73) are untouched.

**Assessment:**
The new asserts are genuinely non-tautological. `If_sym` (SymPy) is the closed-form result of `sp.integrate((d/dξ sech ξ)², ξ, -∞, ∞)`; pinning it to `sp.Rational(2,3)` exercises the actual symbolic integration — a wrong integrand or a botched integration would yield a residual ≠ 0 and the script would raise. `IfNum`/`IgNum` (Mathematica) are NIntegrate results to 30-digit precision of `(f')²` and `(f'')²` for `f = sech`; comparing them to `2/3` and `14/15` makes the quadrature load-bearing. These are NOT the self-cancelling `(4πa²ℓ·If/Hw)/(If/Hw)` ratio the auditor flagged — that ratio cancels `If` identically and is profile-independent; the new asserts instead bind the profile moment to its analytic value. The `I_g = 14/15` constant is correct (`∫sech²=2, ∫sech⁴=4/3, ∫sech⁶=16/15 ⇒ ∫(f'')² = 2 − 16/3 + 64/15 = 14/15 ≈ 0.9333`), and the prior draft's `8/15` is fully purged (annotation corrected, assert uses `14/15`). No collateral edit beyond the directive: the inert `I_1/J_1` checks remain as documentation, the deliverable assertions are unchanged.

### F2 — stale self-labels (numbering)

**Classification:** resolved

**What changed:**
SymPy docstring, NUMBER-token-only:
- `:3` filename self-label `moving_throat_pde_stage53_...` → `moving_throat_pde_stage070_...` (3-digit, per directive).
- `:5` prose `SymPy audit for Stage 53:` → `SymPy audit for Stage 70:` (2-digit, per directive).

**Assessment:**
Correct and scope-respecting. The diff touches exactly the two docstring number tokens. The deferred cross-references and variable names are confirmed UNTOUCHED: py:57 (`Stage-47 integral`), py:58/59 (`Stage-47`/`Stage-48` normalization comments), wl:91 (`Stage-48`), wl:92 (`Stage-47`), and the variable names `J1_stage48` (py:63) / `J1Stage48` (wl:78) all remain at their pre-renumber labels, as the SCOPE GUARD required. The `STAGE 70` ledger banner was not padded. No collateral edits.

### F3 — stale_output (orchestrator-refreshed)

**Classification:** resolved (output refresh)

Both committed `.txt` transcripts were regenerated post-fix (mtimes below). They now carry the canonical `STAGE 70` (SymPy) / `STAGE 070` (Mathematica) banners, include the new sech-moment assertion lines, the corrected `14/15` annotation, and the full numeric-profile cross-check section. No stale `STAGE 53`/`STAGE 053` banner survives.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `kappa - expected = 0`
- `sech-profile moment I_f = 2/3 = 0`  ← new F1 assert PASSES (residual 0)
- `W_wall - expected = 0`, `Xi - W_wall = 0`

**Mathematica:** exit=0. Notable lines:
- `I_g (sech profile) = 0.9333…  (analytic 14/15 = 0.9333…)`  ← corrected annotation, value ≈0.9333 confirmed
- `PASS: sech-profile moment I_f = 2/3`
- `PASS: sech-profile moment I_g = 14/15`  ← **I_g=14/15 PASS confirmed** (≈0.9333, NOT 8/15)
- `PASS: kappa - expected`, `PASS: W_wall - expected`, `PASS: Xi - W_wall` (unchanged deliverables still pass)
- `PASS: kappa numeric profile check`, `PASS: W_wall numeric profile check`, `PASS: Xi = W_wall numeric profile check`

**Output freshness:** confirmed fresh. Script mtimes 2026-06-05 13:48:40 (.py) / 13:48:51 (.wl); committed `.txt` mtimes both 2026-06-05 13:58:32 — outputs are ~10 min NEWER than their scripts. Committed `.txt` content matches the exec logs (banners, sech asserts, `14/15` annotation, PASS lines all present in the committed transcripts).

## Material-change assessment

`material_change`: false.

The four paper deliverables (`kappa`, `W_wall`, `Xi`, `J_1` route) are algebraically unchanged — their assertions and printed closed forms are byte-identical to the pre-fix state (`kappa = 4(mc_{s,w}L/ħ)² + (I_g/I_f)(L/ℓ)²`, `W_wall = Xi = 4ρ_w²V0²L²/(ħ²c_{s,w}²ℓ²)`). F1 only ADDS internal anchor assertions on the sech profile moments (`I_f = 2/3`, `I_g = 14/15`), which are this stage's own scaffolding and are not consumed by any downstream unit (Stage 071 separately uses the tanh profile with `I_f = 1/3`). F2 is label-only. No downstream unit can depend on a changed derived result, so no `upstream_stale` propagation is materially warranted.

## Side observations (non-blocking)

- The corrected Mathematica `I_g` annotation is a genuine fix to a previously-misleading print (`8/15` was wrong); the original auditor report itself had repeated the wrong `8/15` in its F1 "Required change" and in the value-reconciliation table (report L99, L161), but the directive caught and corrected this to `14/15` (orchestrator-verified). The shipped code uses the correct `14/15` throughout. Non-blocking — the shipped state is correct.
- SymPy still has no `I_g` symbolic moment assertion (only `I_f`); the directive only asked for the `I_f` symbolic assert in SymPy and the `I_g` numeric assert in Mathematica, so this matches the directive. Not a gap.

## Verdict justification

Both findings are fully resolved. F1's de-tautologization is substantive: the new `I_f = 2/3` (SymPy symbolic) and `I_f = 2/3` / `I_g = 14/15` (Mathematica 30-digit NIntegrate) assertions genuinely exercise the sech-profile integration — they are not the self-cancelling `I_1/J_1` ratio and would fail on a wrong moment. The I_g value is the correct `14/15` (≈0.9333), with the prior `8/15` error purged from both the annotation and the assert, and the Mathematica log shows `PASS: sech-profile moment I_g = 14/15`. F2's label fix is number-token-only and respects the deferred-scope guard. All deliverable assertions (kappa/W_wall/Xi) are unchanged and still pass; both engines exit 0; committed outputs are refreshed with canonical banners and the new lines. `material_change: false`.
