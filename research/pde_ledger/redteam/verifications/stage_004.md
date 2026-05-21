---
unit_id: 004
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 004

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**

`scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py:32-44` (post-edit numbering).
The old block (formerly lines 37-39) constructed `by_parts_density = sp.diff(W*Q, w) - W*sp.diff(Q, w) - sp.diff(W, w)*Q` over abstract `W, Q = sp.Function(...)(w)` symbols and asserted it was zero — a product-rule identity that `sp.diff` enforces by construction. Codex replaced it with concrete decaying Gaussian densities:

```
lam_ibp = sp.Symbol("lam_ibp", positive=True)
W_ex = sp.exp(-w**2 / lam_ibp**2)
Q_ex = w * sp.exp(-w**2 / lam_ibp**2)
ibp_lhs = sp.integrate(W_ex * sp.diff(Q_ex, w), (w, -sp.oo, sp.oo))
ibp_rhs = -sp.integrate(sp.diff(W_ex, w) * Q_ex, (w, -sp.oo, sp.oo))
assert_zero("projection integration-by-parts (decaying test functions)",
            sp.simplify(ibp_lhs - ibp_rhs))
```

The previously unused `W = sp.Function("W")(w)` and `Q = sp.Function("Q")(w)` declarations were removed per the directive's mechanical-removal allowance (they had no other readers after this block was rewritten).

**Assessment:**

The replacement is exactly what the directive's "required change" block specified, character for character. Substantively, the new check is non-tautological:

- SymPy must actually evaluate two improper Gaussian integrals on `(-∞, ∞)`, not merely apply the product rule.
- The identity holds only because `W_ex · Q_ex` decays at both infinities, so the boundary term `[W·Q]_{-∞}^{∞}` vanishes. If the densities had been chosen non-decaying (e.g. polynomial probes, or one Gaussian and one constant), the integrals would either diverge or not satisfy the IBP identity.
- The label string `"projection integration-by-parts (decaying test functions)"` is honest about what is exercised.

The saved output transcript at `scripts/output/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.txt` shows `STATUS: PASS` with mtime 2026-05-21 11:26, post-dating the script mtime 01:09, so the assertion fires on the post-fix script without raising. No collateral edits beyond the directive-sanctioned removal of the unused abstract declarations.

### F2 — missing_verification_script

**Classification:** resolved

**What changed:**

New file `mathematica/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl` (124 lines) was created. It contains six independent blocks M1-M6, each terminating in an `If[FullSimplify[...] =!= 0, Print["FAIL: ..."]; Exit[1]]` guard, with `Print["M# residual = ...]` diagnostic above and `Print["PASS: ..."]` below. Final line prints `STATUS: PASS` and exits 0.

**Assessment:**

Per-claim substance:

- **M1** (`.wl:12-26`): Codex chose the combined-integrand formulation `Integrate[W·f' + W'·f, {w, -∞, ∞}]` rather than the directive's literally-stated two-separate-integrals form. The two are mathematically equivalent (and Mathematica's `Integrate` recognizes the combined integrand as `d/dw(W·f)` and evaluates to the boundary term, which vanishes by decay). This is the iter2 fix after an iter1 residual bug. The combined form is slightly less informative than two-side comparison, but it is not tautological: a non-decaying choice of densities would either fail to integrate or yield a nonzero boundary term. Acceptable as a substantive IBP-under-decay check.
- **M2** (`.wl:28-69`): Three independent cyclic permutations. Each builds the antisymmetric two-form `twoForm##` from explicit field substitutions (`twoForm02 = -E2[t,x,y,z]`, etc.), takes `D[...]` derivatives, and compares to the explicit Faraday vector form. Uses Mathematica's `D[]` and `FullSimplify` rather than SymPy's `sp.diff`. Each component has its own guard.
- **M3** (`.wl:71-82`): `Integrate[Exp[-w^2/lam^2], {w, -∞, ∞}, Assumptions -> lam > 0]` independently evaluates and compares to `Sqrt[Pi]*lam`.
- **M4** (`.wl:84-94`): `Integrate[Exp[-w^2/lam^2]^2, ...]` compared to `Sqrt[2*Pi]*lam/2`.
- **M5** (`.wl:96-107`): `matchedWeight = localizedProfile / profileArea`, then `Integrate[matchedWeight*localizedProfile, ...]` compared to `Sqrt[2]/2`. Derived from the integral definition above the assertion, not transliterated from `I_WZ`.
- **M6** (`.wl:109-121`): `pointSourceCoupling = mu0*(matchedWeight /. w -> 0)/overlapValue`, `volumeReducedCoupling = mu0/profileArea`, ratio compared to `Sqrt[2]` with `couplingAssumptions = lam > 0 && mu0 > 0`.

Naming conventions are Mathematica-idiomatic and do not mirror the SymPy intermediates: `localizedProfile/profileArea/matchedWeight/overlapValue/pointSourceCoupling/volumeReducedCoupling` instead of `Z/Z_int/W_match/I_WZ/mu_proj_delta/mu_red`. Each claim is derived from its own integral or algebraic definition immediately above the guard — no imports of the `.py` script's intermediates. The script uses `Assumptions ->` correctly (scaleAssumptions for `lam > 0`-only claims, couplingAssumptions for M6 which involves `mu0`).

The exec log shows all six residuals are literally `0` (not just `FullSimplify[...] === 0`-after-fiddling) and the script exits `0`. The new file exists at the directive-named path.

## Exec log assessment

**SymPy:** exit=n/a. The orchestrator did not capture a `stage_004_sympy.log` for this iteration. However, the saved transcript at `scripts/output/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.txt` (mtime 2026-05-21 11:26, after the script's 01:09 mtime) contains:

```
STEP 01 PROJECTED MAXWELL README AUDIT
Checked script inventory, projection identity, vector Bianchi signs, and Gaussian coupling mismatch.
STATUS: PASS
```

This is sufficient to confirm the post-fix script reached its final `print("STATUS: PASS")` without `AssertionError` (which would have aborted before reaching that line).

**Mathematica:** exit=0. Notable lines from `redteam/exec_logs/stage_004_mathematica.log`:

```
M1 residual = 0
M2 component 1 residual = 0
M2 component 2 residual = 0
M2 component 3 residual = 0
M3 residual = 0
M4 residual = 0
M5 residual = 0
M6 residual = 0
STATUS: PASS
```

All six claims pass with residuals literally simplifying to `0`.

**Output freshness:** Confirmed.
- `scripts/.../sympy_audit.py` mtime 2026-05-21 01:09 < `scripts/output/.../sympy_audit.txt` mtime 2026-05-21 11:26.
- `mathematica/.../mathematica_audit.wl` mtime 2026-05-21 11:22 < `mathematica/output/.../mathematica_audit.txt` mtime 2026-05-21 11:50.

Both saved transcripts post-date their source scripts, so they reflect the post-fix runs.

## Material-change assessment

`material_change`: false.

Neither edit changes a derived result that downstream units consume. F1 replaced a tautological assertion with a genuine integral check — the script's exported numerical/symbolic claims (Gaussian integrals A5-A8, Faraday sign checks A2-A4) are unchanged. F2 is purely additive: a new `.wl` file that re-derives the same claims in a second engine. No downstream unit reads variables out of this script.

## Side observations (non-blocking)

1. The orchestrator did not produce `stage_004_sympy.log` (the file is absent from `redteam/exec_logs/`), while it did produce `stage_004_mathematica.log`. The `.txt` saved transcript is fresh and confirms PASS, so this isn't a verification blocker, but the orchestrator may want to investigate why the sympy log capture step was skipped for this iteration.

2. M1 in the `.wl` deviates from the directive's literal "two-separate-integrals" specification, using a combined-integrand form instead. This is functionally equivalent and was reportedly the iter2 fix after an iter1 IBP residual bug. The combined form is slightly less informative as a diagnostic (it cannot localize which side of the IBP identity is wrong if a residual appears), but it is not tautological and matches the spirit of the claim. Not raising as a finding.

3. The combined-integrand M1 in Mathematica and the two-sided M1-equivalent in SymPy are not structurally identical formulations, which is consistent with the directive's "do NOT mirror SymPy variable choreography" requirement. The engine cross-check is now meaningful: both engines independently confirm the IBP identity under Gaussian decay, the Gaussian normalization `Sqrt[Pi]*lam`, the squared-norm `Sqrt[2*Pi]*lam/2`, the matched-kernel overlap `Sqrt[2]/2`, and the projection/reduction ratio `Sqrt[2]`.

## Verdict justification

Both findings are `resolved`. F1's new IBP check actually evaluates two improper Gaussian integrals and compares them, which is genuinely non-tautological — a misplaced sign or a non-decaying density choice would produce a nonzero residue. F2's new `.wl` re-derives each M1-M6 claim independently from integral or algebraic definitions, uses Mathematica idioms, employs distinct naming from the `.py`, and exits `0` with every residual literally `0`. Saved transcripts post-date their source scripts. No regressions visible in the diff; no collateral edits beyond directive-sanctioned ones. Verdict: `verified`.
