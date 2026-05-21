---
unit_id: 018
batch: I.2
created_at: 2026-05-21T12:29:31-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-21T13:28:53-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 018

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:41-42`

**Issue:**
The assertion `assert_zero("compatibility equality", sp.simplify(K_from_norm - K_from_one_pole) - compatibility)` reduces to `expr - expr == 0` by construction: `K_from_norm - K_from_one_pole` is, by definition, `N0/Ptarget - 3*(MSigma+B2+Z2)**2/(B4+Z4)`, and `compatibility` is defined on the preceding line as exactly that expression. The assertion cannot fail regardless of physics.

**Required change:**
Replace lines 41-42 of `moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py` with a substantive cross-closure check.

Before:
```python
    compatibility = N0 / Ptarget - 3 * (MSigma + B2 + Z2) ** 2 / (B4 + Z4)
    assert_zero("compatibility equality", sp.simplify(K_from_norm - K_from_one_pole) - compatibility)
```

After:
```python
    compatibility = N0 / Ptarget - 3 * (MSigma + B2 + Z2) ** 2 / (B4 + Z4)
    N0_from_compat = sp.solve(compatibility, N0)[0]
    N0_from_equality = sp.solve(K_from_norm - K_from_one_pole, N0)[0]
    assert_zero("compatibility equality", N0_from_compat - N0_from_equality)
    assert_nonzero("mutated compatibility equality should fail", N0_from_compat - 2 * N0_from_equality)
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 018` and confirm the substantive check appears and the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py`
- summary: Replaced the compatibility identity assertion with independent N0 solves and a mutation check.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:70-71`

**Issue:**
The assertions `assert_zero("compensated K1", K1.subs(sol))` and `assert_zero("compensated H_even", H_even.subs(sol))` substitute `sol` — the `sp.solve` solution that *defines* `K1 = 0` and `H_even = 0` — back into the same expressions. Both assertions are algebraically guaranteed by the contract of `sp.solve` and cannot fail.

**Required change:**
Replace lines 70-71 of `moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py` with substantive cross-checks routed through `expected_dK` and `expected_dM` instead of `sol`.

Before:
```python
    assert_zero("compensated K1", K1.subs(sol))
    assert_zero("compensated H_even", H_even.subs(sol))
```

After:
```python
    assert_zero("expected slopes satisfy K1", K1.subs({dKSigma: expected_dK, dMSigma: expected_dM}))
    assert_zero("expected slopes satisfy H_even", H_even.subs({dKSigma: expected_dK, dMSigma: expected_dM}))
    assert_nonzero("mutated expected slopes fail K1", K1.subs({dKSigma: expected_dK + 1, dMSigma: expected_dM}))
```

Do not move or alter lines 69 (`D01_comp`) or 72-79 (`Xi1` block).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 018` and confirm the new substantive checks appear and the script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py`
- summary: Replaced solve back-substitution checks with direct closed-form slope checks and a mutation check.
- deviation: none

## F3 — missing_verification_script

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl` (create new file)

**Issue:**
No Mathematica audit script exists for stage 018. Manifest entry `is_checkpoint: false` and `is_status_only_candidate: false` require the two-engine policy. The new `.wl` file must independently re-derive (not transliterate the SymPy choreography) the claims enumerated in the manifest below.

**Required change:**
Create the file at the Target path. The script must:

1. Begin with a comment block stating the unit, that this is the Mathematica counterpart of the SymPy audit, and the claim list.
2. For each claim, write the symbolic setup *from physical premises* (definitions of `D0, D2, D4, K1, H_even, Xi1`, etc.) in Mathematica idiom — use `D[expr, var]`, `Solve[eqns, vars]`, `Integrate[expr, {w, -Infinity, Infinity}]`. Do NOT mirror the SymPy variable-by-variable choreography line-by-line; derive each claim independently.
3. Wrap each assertion in `If[FullSimplify[lhs - rhs] =!= 0, (Print["FAIL: <label>"]; Exit[1])]` for `assert_zero` claims, and `If[FullSimplify[expr] === 0, (Print["FAIL: <label> unexpectedly vanished"]; Exit[1])]` for `assert_nonzero` mutation claims.
4. End with `Print["STAGE 018 MATHEMATICA AUDIT PASS"]; Exit[0]`.

**Claim manifest:**

M1 — One-pole numerator identity:
$$ (u_4 - 4 u_2^2) - \frac{D_0 (B_4 + Z_4) - 3(M_\Sigma + B_2 + Z_2)^2}{D_0^2} = 0 $$
where $D_0 = K_\Sigma - B_0 - Z_0$, $D_2 = -(M_\Sigma + B_2 + Z_2)$, $D_4 = -(B_4 + Z_4)$, $u_2 = -D_2/D_0$, $u_4 = (D_2^2 - D_0 D_4)/D_0^2$.

M2 — One-pole `KSigma` closure: substituting $K_\Sigma \to B_0 + Z_0 + 3(M_\Sigma + B_2 + Z_2)^2/(B_4 + Z_4)$ into $u_4 - 4 u_2^2$ gives 0.

M3 — Normalization closure: substituting $K_\Sigma \to B_0 + Z_0 + N_0/P_{\rm target}$ into $N_0/D_0$ gives $P_{\rm target}$.

M3-mut — Mutation check: $(N_0/D_0)|_{K_\Sigma = K_{\rm from\_norm}} - 2 P_{\rm target} \neq 0$.

M4 — Compatibility cross-closure: with `compatibility = N0/Ptarget - 3(M+B2+Z2)^2/(B4+Z4)`, the `N0` solution of `compatibility = 0` equals the `N0` solution of `K_from_norm - K_from_one_pole = 0`. (Mirror F1's substantive form.)

M5 — Even-gate determinant:
$$ \det \begin{pmatrix} \partial K_1 / \partial \delta K_\Sigma & \partial K_1 / \partial \delta M_\Sigma \\ \partial H_{\rm even}/\partial \delta K_\Sigma & \partial H_{\rm even}/\partial \delta M_\Sigma \end{pmatrix} = \frac{1}{27} $$
where $K_1 = D_{21} + D_{01}/9$, $H_{\rm even} = D_{41} - (2/3) D_{21} - D_{01}/27$, $D_{01} = \delta K_\Sigma - B_{01} - Z_{01}$, $D_{21} = -(\delta M_\Sigma + B_{21} + Z_{21})$, $D_{41} = -(B_{41} + Z_{41})$.

M5-mut — `det + 1/27 != 0`.

M6 — Wall-stiffness slope: the closed-form $\delta K_\Sigma = B_{01} + Z_{01} + 27(B_{41} + Z_{41})$ together with $\delta M_\Sigma = -(B_{21} + Z_{21}) + 3(B_{41} + Z_{41})$ satisfies both $K_1 = 0$ and $H_{\rm even} = 0$. (Verify by direct substitution, not by re-running `Solve` and back-substituting — that path mirrors SymPy too closely and would be tautological.)

M7 — Residual amplitude: with the closed-form slopes substituted,
$$ \Xi_1 = \frac{N_{01}}{N_0} - \frac{27(B_{41}+Z_{41})}{K_\Sigma - B_0 - Z_0} $$
where $\Xi_1 = N_{01}/N_0 - D_{01}/D_0$.

M7-mut — Same with sign flipped on the second term, must not vanish.

M8 — Gaussian wall integrals: with $\beta(w) = \exp(-w^2/2)$,
$$ \int_{-\infty}^{\infty} \beta^2 \, dw = \sqrt{\pi}, \qquad \int_{-\infty}^{\infty} \left( \left(\frac{d\beta}{dw}\right)^2 + \beta^2 \right) dw = \frac{3\sqrt{\pi}}{2}. $$

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 018` and confirm the new `.wl` script is present, contains the `If[... Exit[1]]` assertions for each Mn above, and exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl`
- summary: Added the stage 018 Mathematica audit covering the manifest checks and mutation checks.
- deviation: none

## F4 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:79`

**Issue:**
The `residual Xi1` assertion at lines 72-75 follows by elementary substitution from A11 (`D01_comp = 27*(B41+Z41)`). The only independent path is to recompute `Xi1` from the *closed-form expected slopes* rather than from `sp.solve`'s output `sol`, then verify it equals the same right-hand side.

**Required change:**
Insert an additional check directly after line 79 (after the existing `assert_nonzero` mutation block on `Xi1`). Do not modify lines 72-79.

Insert after line 79:
```python
    Xi1_from_expected = Xi1.subs({dKSigma: expected_dK, dMSigma: expected_dM})
    assert_zero(
        "residual Xi1 from expected slopes",
        Xi1_from_expected - (N01 / N0 - 27 * (B41 + Z41) / (KSigma - B0 - Z0)),
    )
```

This routes the `Xi1` claim through `expected_dK`/`expected_dM` (the closed forms) instead of `sol`, providing an independent verification path that does not share intermediate simplification with A11.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 018` and confirm the new check appears after the existing mutation row and the script exits 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py`
- summary: Added an Xi1 residual check routed through the closed-form expected slopes.
- deviation: none
