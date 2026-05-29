---
unit_id: 119
batch: IV.3
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 119

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py:52`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:59`

**Issue:**
Section III only confirms the algebra of solving `4 L_W^2 / (pi^2 a^2) == (1 + rc)/3` for `L_W`. It does not exercise the family connection `r_c = rhat^2` that the notes use to derive the boxed `L_W` formula. The assertion is algebraically guaranteed by construction. We add an extra consistency check that ties the section's `rc` parameter back to the `rhat` symbol used in sections I-II, so the assertion no longer floats free of the rest of the parent-balance family.

**Required change (SymPy):**

After the existing `expect_zero("tube-length law", ...)` call at line 52 of `scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py`, insert two lines:

Before (around line 52):
```python
print("L_W =", sp.simplify(L_sel))
expect_zero("tube-length law", L_sel - sp.pi*a*sp.sqrt((1+rc)/3)/2)
```

After:
```python
print("L_W =", sp.simplify(L_sel))
expect_zero("tube-length law", L_sel - sp.pi*a*sp.sqrt((1+rc)/3)/2)
expect_zero(
    "tube-length law (rc -> rhat**2 link)",
    L_sel.subs(rc, rhat**2) - sp.pi*a*sp.sqrt((1+rhat**2)/3)/2
)
```

This new assertion confirms that with the notes-mandated identification `r_c = rhat^2`, the same `L_W` closed form holds with `rhat` (the parameter used in sections I-II). It is non-tautological because it ties the L_W section to the dimensionless-ratio section.

**Required change (Mathematica):**

After the existing `expectZero["tube-length law", ...]` call at line 59 of `mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl`, insert two lines. Note that Section III declares `$Assumptions = Element[{a, lW, rC}, Reals] && a > 0 && lW > 0 && rC > 0;` and does NOT declare `rHat`. To use `rHat` in the new assertion without breaking the section's assumptions, extend the `$Assumptions` for this assertion only, or extend the section-level assumption block. Use the section-level extension to keep the block self-contained:

Before (around line 53-59):
```mathematica
Clear[a, lW, rC];
$Assumptions = Element[{a, lW, rC}, Reals] && a > 0 && lW > 0 && rC > 0;

kappa0 = FullSimplify[4*lW^2/(Pi^2*a^2), Assumptions -> $Assumptions];
lSel = FullSimplify[lW /. First[Solve[kappa0 == (1 + rC)/3, lW, Reals]], Assumptions -> $Assumptions];

Print["L_W = ", fmt[lSel]];
expectZero["tube-length law", lSel - (Pi*a*Sqrt[(1 + rC)/3])/2];
```

After:
```mathematica
Clear[a, lW, rC, rHat];
$Assumptions = Element[{a, lW, rC, rHat}, Reals] && a > 0 && lW > 0 && rC > 0;

kappa0 = FullSimplify[4*lW^2/(Pi^2*a^2), Assumptions -> $Assumptions];
lSel = FullSimplify[lW /. First[Solve[kappa0 == (1 + rC)/3, lW, Reals]], Assumptions -> $Assumptions];

Print["L_W = ", fmt[lSel]];
expectZero["tube-length law", lSel - (Pi*a*Sqrt[(1 + rC)/3])/2];
expectZero[
  "tube-length law (rC -> rHat^2 link)",
  (lSel /. rC -> rHat^2) - (Pi*a*Sqrt[(1 + rHat^2)/3])/2
];
```

(The added `rHat` symbol in `Clear` and the assumptions list keeps the section self-contained; Section IV later resets `$Assumptions` with its own `Clear` block, so this change does not leak.)

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 119` and `redteam exec-mathematica 119` and confirm a new check line `tube-length law (rc -> rhat**2 link) = 0` appears in the SymPy output and `tube-length law (rC -> rHat^2 link)` appears as `PASS` in the Mathematica output, and both scripts exit 0.

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage119_parent_balance_sympy_audit.py:66-69`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:78-82`

**Issue:**
Section IV solves for `T_m` on both branches and only `print`s the result; no assertion is made. The notes (§5) give the explicit form `T_m = sqrt(2 Z_q K_s) / (J_s c_s sqrt(mu_0 L_W)) * 1 / (rhat +/- (1/2) sqrt(1 + rhat^2))`. Add assertions that the script's `Tm_sol_plus` and `Tm_sol_minus` match this notes-given form. (The factor-of-2 rearrangement makes the literal target form `2 sqrt(2 Z_q K_s) / (J_s c_s sqrt(mu_0 L_W) (2 rhat +/- sqrt(1 + rhat^2)))`, which matches the current printed SymPy output exactly.)

**Required change (SymPy):**

After the existing `print("T_m (- branch) =", Tm_sol_minus)` line (line 69) and before the `print("\nAll Stage 119 symbolic checks passed.")` line (line 71), insert two assertions:

```python
expect_zero(
    "T_m (+ branch) match",
    Tm_sol_plus - 2*sp.sqrt(2)*sp.sqrt(K_s)*sp.sqrt(Zq) /
        (J_s*sp.sqrt(L_W)*c_s*sp.sqrt(mu0)*(2*rhat + sp.sqrt(1+rhat**2)))
)
expect_zero(
    "T_m (- branch) match",
    Tm_sol_minus - 2*sp.sqrt(2)*sp.sqrt(K_s)*sp.sqrt(Zq) /
        (J_s*sp.sqrt(L_W)*c_s*sp.sqrt(mu0)*(2*rhat - sp.sqrt(1+rhat**2)))
)
```

The `expect_zero` helper uses `sp.simplify(sp.expand(...))` which should reduce the difference to 0 since the printed forms in the output already exhibit this exact closed form.

**Required change (Mathematica):**

The Mathematica `tMPlus` and `tMMinus` outputs are wrapped in `ConditionalExpression[..., condition]` (the saved transcript shows this at lines 39-40). The `expectZero` helper applies `FullSimplify[Together[Expand[expr]]]` and then checks `TrueQ[res === 0]`. `ConditionalExpression` will survive these passes, so a direct subtraction will not reduce to `0`. The safest fix is to strip `ConditionalExpression` heads before comparison.

After the existing `Print["T_m (- branch) = ", fmt[tMMinus]];` line (line 82) and before the `Print[""];` line (line 84), insert:

```mathematica
stripCE[expr_] := expr /. ConditionalExpression[v_, _] :> v;

expectZero[
  "T_m (+ branch) match",
  stripCE[tMPlus] - (2*Sqrt[2]*Sqrt[kS]*Sqrt[zQ]) /
    (jS*Sqrt[lW]*cSound*Sqrt[mu0]*(2*rHat + Sqrt[1 + rHat^2]))
];
expectZero[
  "T_m (- branch) match",
  stripCE[tMMinus] - (2*Sqrt[2]*Sqrt[kS]*Sqrt[zQ]) /
    (jS*Sqrt[lW]*cSound*Sqrt[mu0]*(2*rHat - Sqrt[1 + rHat^2]))
];
```

`stripCE` is a local helper (defined inline at insertion point); no other code depends on it. If the inline definition feels intrusive, hoist it to just below the existing `expectZero` definition at line ~24, but the inline placement is acceptable since it is only used here.

The section's current `$Assumptions` (line 64-66) includes `kS > 0 && zQ > 0 && mu0 > 0 && cSound > 0 && tM > 0 && jS > 0 && lW > 0`. It does NOT declare `rHat`. For the new assertions to simplify cleanly, add `rHat` to the `Element[...]` list in the existing `$Assumptions` line, matching what Section II uses. Specifically, change line 65 from:

```mathematica
  Element[{kS, rHat, zQ, mu0, cSound, tM, jS, lW}, Reals] &&
```

Wait — re-reading the actual file: line 65 already contains `Element[{kS, rHat, zQ, mu0, cSound, tM, jS, lW}, Reals] &&`. `rHat` is already there. No change to the assumption block is needed for F2.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 119` and `redteam exec-mathematica 119`. Expected output additions:
- SymPy: two new lines `T_m (+ branch) match = 0` and `T_m (- branch) match = 0`.
- Mathematica: two new `PASS: T_m (+ branch) match` and `PASS: T_m (- branch) match` lines.
Both scripts must exit 0.
