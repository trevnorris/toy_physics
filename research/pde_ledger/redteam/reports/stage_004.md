---
unit_id: 004
batch: I.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-20T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: missing
  engines_agree: n/a
  outputs_fresh: true
---

# Audit unit 004 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py`
- mathematica: `(missing)`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.txt`
- mathematica output: `(missing)`

## What the script claims to verify

The script is a bundle-index audit for the first "projected Maxwell" bundle. It performs four kinds of checks: (i) the three downstream bundle scripts (stages 005, 006, 007) exist on disk; (ii) an "integration-by-parts density" identity in `w`; (iii) the sign/index conventions that take the antisymmetric two-form's cyclic derivative (`∂_{[a} F_{bc]}`) to the vector Faraday form `B_t + curl(E)`; and (iv) closed-form Gaussian integrals plus a `√2` ratio between a "delta-source projection" coupling and a naive "reduction" coupling using a matched kernel `W = Z/∫Z`. The verdict applies to those four checks as asserted in the script.

## Assertion inventory

| #  | Script | Line   | Form                                                          | Anchored to claim? |
|----|--------|--------|---------------------------------------------------------------|--------------------|
| A1 | sympy  | 38–39  | `simplify(diff(W*Q,w) - W*diff(Q,w) - diff(W,w)*Q) == 0`      | no (tautology)     |
| A2 | sympy  | 53–55, 61–62 | Faraday component 1 cyclic-sum vs `B1_t + E3_y - E2_z`   | partial            |
| A3 | sympy  | 56–57, 61–62 | Faraday component 2 cyclic-sum vs `B2_t + E1_z - E3_x`   | partial            |
| A4 | sympy  | 58–59, 61–62 | Faraday component 3 cyclic-sum vs `B3_t + E2_x - E1_y`   | partial            |
| A5 | sympy  | 74     | `simplify(Z_int - sqrt(pi)*lam) == 0`                          | yes                |
| A6 | sympy  | 75     | `simplify(Z2_int - sqrt(2*pi)*lam/2) == 0`                     | yes                |
| A7 | sympy  | 76     | `simplify(I_WZ - sqrt(2)/2) == 0`                              | yes                |
| A8 | sympy  | 77     | `simplify(mu_proj_delta/mu_red - sqrt(2)) == 0`                | yes                |
| —  | sympy  | 28–30  | `FileNotFoundError` on missing bundle scripts                  | meta (file inventory) |

A1 is flagged as tautology below. A2–A4 are "partial" because the assertion only checks that the F-definitions defined immediately above (`F23=B1`, `F02=-E2`, etc.) are self-consistent with the hand-written vector form — it confirms sign conventions but does not exercise an independent identity. They are sufficient as sign/index checks per the comment "vector Bianchi signs," so no finding is raised on them. A5–A8 verify the Gaussian closed forms and the `√2` projection/reduction ratio non-trivially: the LHS is computed by SymPy from `sp.integrate` on the actual Gaussian densities and is compared to an independently expected closed form.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py:37-39`

**What's wrong:**

The comment above the first assertion claims to test the integration-by-parts identity:

```
# int W Q_w dw = [WQ] - int W_w Q dw.
by_parts_density = sp.diff(W * Q, w) - W * sp.diff(Q, w) - sp.diff(W, w) * Q
assert_zero("projection integration-by-parts density", by_parts_density)
```

But the residue tested is `d/dw(W·Q) − W·Q_w − W_w·Q`, which is the **product rule for differentiation**, not the integration-by-parts identity. SymPy's `sp.diff(W*Q, w)` is implemented by the product rule and returns exactly `W*Q'(w) + W'(w)*Q` for two unknown symbols-as-functions of `w`. The residue is therefore identically zero by construction of `sp.diff`, independent of any physical content. The assertion cannot fail no matter what `W` and `Q` are; it tests SymPy, not the script's claimed identity.

Worse, the comment misnames the check: the integration-by-parts identity `∫ W Q_w dw = [W Q] − ∫ W_w Q dw` involves an actual integral and a boundary term `[WQ]`, neither of which appears anywhere in the residue. The label "projection integration-by-parts density" in the assert message therefore overstates the check.

**Why this matters:**

A tautological check raises false confidence: future readers (or an automated coverage scan that counts assertions) will count this as a verified property of the projection step, when in fact it is just `True == True`. Either the projection integration-by-parts identity is verified by an actual integral test, or the assert+comment must be honest that all that's exercised is the product rule.

**Required change:**

Replace lines 37–39 so the assertion actually exercises integration-by-parts at the density level. One direct, in-script-derivable check: introduce concrete decaying densities and verify the integral identity holds numerically/symbolically with the boundary term explicit. Concretely, edit the block to:

```
# Verify integration-by-parts at density level: d(WQ)/dw = W Q_w + W_w Q
# is the product rule (built into sp.diff). The non-trivial content is the
# boundary-term cancellation under decay, exercised on concrete decaying W,Q.
lam_ibp = sp.Symbol("lam_ibp", positive=True)
W_ex = sp.exp(-w**2 / lam_ibp**2)
Q_ex = w * sp.exp(-w**2 / lam_ibp**2)
lhs = sp.integrate(W_ex * sp.diff(Q_ex, w), (w, -sp.oo, sp.oo))
rhs = -sp.integrate(sp.diff(W_ex, w) * Q_ex, (w, -sp.oo, sp.oo))
assert_zero("projection integration-by-parts (decay, boundary vanishes)",
            sp.simplify(lhs - rhs))
```

This actually exercises the integral identity under decay assumptions and cannot pass by tautology — if the sign convention or boundary handling were wrong, the residual would not simplify to zero. Cite: `moving_throat_pde_stage004_projected_maxwell_bundle_index_sympy_audit.py:37-39`.

**Verification:**

A re-run produces a transcript whose new line `projection integration-by-parts (decay, boundary vanishes)` is reached without `AssertionError`, and the script exits `0`. The verifier confirms the new label string appears in the script source at the named line range.

### F2 — missing_verification_script

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/` (expected `moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl`)

**What's wrong:**

The audit-prompt header explicitly lists the Mathematica script for unit 004 as `(missing)`, and the unit's manifest entry has `is_status_only_candidate: False`. Per the auditor charter, non-status, non-checkpoint units require both engines. Stage 004 currently has only a SymPy script; nothing independently re-derives, in a second engine, the three substantive claims A5–A8 (Gaussian closed forms and `√2` projection/reduction ratio) or A2–A4 (Bianchi-to-Faraday sign reduction).

**Why this matters:**

Without a second engine, the `√2` projection/reduction ratio and the Gaussian closed forms rest on a single CAS's implementation of `integrate` and `simplify`. A second-engine derivation is the policy guard against engine-specific quirks (assumption handling on improper integrals, branch cuts in `√`, simplification of `sqrt(2*pi)*lam/2` vs `lam/sqrt(2/pi)`, etc.). Until a `.wl` exists, the `engines_agree` field is `n/a` rather than `true`, and the unit's status as "verified" is incomplete.

**Required change:**

Create `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage004_projected_maxwell_bundle_index_mathematica_audit.wl` that independently derives — from the physical premises stated in the SymPy script's docstring and inline comments, NOT by line-by-line transliteration — the claims enumerated in the manifest below. Use Mathematica idioms (`Integrate`, `FullSimplify`, `Assuming`) and choose intermediate forms differently from the SymPy script. Each claim must terminate in a check of the form `If[FullSimplify[lhs - rhs] =!= 0, Print["FAIL: <label>"]; Exit[1]]`.

**Claim manifest:**

- **M1** — Density-level integration by parts under decay:
  For decaying test functions `W(w), Q(w)` on `(-∞, ∞)`, verify
  `∫_{-∞}^{∞} W(w) Q'(w) dw + ∫_{-∞}^{∞} W'(w) Q(w) dw = 0`
  using concrete decaying examples (e.g. `W = exp(-w^2/lam^2)`, `Q = w·exp(-w^2/lam^2)`). The boundary term `[WQ]_{-∞}^{∞}` vanishes by decay.

- **M2** — Cyclic Bianchi reduces to vector Faraday, component 1:
  With `F_{23} = B_1`, `F_{30} = E_3`, `F_{02} = -E_2`, verify
  `∂_t F_{23} + ∂_y F_{30} + ∂_z F_{02} = ∂_t B_1 + ∂_y E_3 − ∂_z E_2`.
  (Component 2 and component 3 analogous with cyclic permutation.)

- **M3** — Gaussian normalization integral:
  For `Z(w) = exp(-w²/λ²)` with `λ > 0`, verify
  `∫_{-∞}^{∞} Z(w) dw = √π · λ`.

- **M4** — Gaussian squared-norm integral:
  `∫_{-∞}^{∞} Z(w)² dw = (√(2π)/2)·λ`.

- **M5** — Matched-kernel coupling `I_{WZ}`:
  Define `W_match(w) = Z(w) / ∫Z`. Verify
  `∫_{-∞}^{∞} W_match(w)·Z(w) dw = √2 / 2`.

- **M6** — Delta-source projection/reduction ratio:
  Define `μ_{proj,δ} = μ_0 · W_match(0) / I_{WZ}` and `μ_{red} = μ_0 / ∫Z`. Verify
  `μ_{proj,δ} / μ_{red} = √2`.

The Mathematica script must NOT structurally mirror the SymPy script (no identical variable choreography, no identical intermediate symbol names). Derive each result from the integral definition above the assertion; do not import intermediates from the `.py` file.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 004` and confirm the new `.wl` script exists, contains assertions of the `If[FullSimplify[...] =!= 0, ... Exit[1]]` form for each of M1–M6, and exits `0`.

## Independent-derivation check (Mathematica)

No `.wl` exists for this unit; the independent-derivation check is not applicable beyond noting its absence (see F2).

## Engine cross-check

Single-engine only (SymPy). Engine cross-check is not applicable; `engines_agree: n/a` in the front-matter reflects this. Once F2 is resolved, this section should compare the final symbolic forms of `Z_int`, `Z2_int`, `I_WZ`, and `mu_proj_delta/mu_red` between engines.

## Verdict justification

The four Gaussian / projection-ratio assertions (A5–A8) hold up under attack: I traced each integral by hand (`∫exp(-w²/λ²) dw = √π·λ`, `∫exp(-2w²/λ²) dw = √(π/2)·λ`, `I_WZ = Z2_int/Z_int = √2/2`, and `μ_proj_δ/μ_red = (2/√(2π·λ²))/(1/(√π·λ)) = √2`), and they match the script's closed-form targets. The sign-convention checks A2–A4 also survive: substituting `F_{02}=-E_2` and cyclic permutations into `∂_t F_{23}+∂_y F_{30}+∂_z F_{02}` does produce `B_{1,t}+E_{3,y}-E_{2,z}`. Symbol assumptions (`lam>0`, `mu0>0`, `w` real) are consistent with the improper-integral closed forms used.

The two findings are: A1 is a product-rule tautology mislabeled as integration-by-parts (F1), and the unit lacks a Mathematica counterpart (F2). Verdict is `findings`, not `stop_cold` — neither finding is `UNFIXABLE` (both are mechanically patchable) and neither is `CRITICAL_DOWNSTREAM` (F1 just changes the test, not any quoted forward result; F2 is purely an additive artifact).
