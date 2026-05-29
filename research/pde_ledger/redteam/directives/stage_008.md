---
unit_id: 008
batch: I.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-25T17:12:32-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive -- unit 008

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>`
block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append
`## Blocked: F<n>` with a question instead -- skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 -- tautological_check

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl` lines 31-46 (M2 Pair A), 49-67 (M2 Pair B), 88-131 (M4 block, specifically line 119), 152-171 (M7 block)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py` lines 184-201 (matched-Gaussian H=Z) and 237-241 (independent-Gaussian H=Z)

**Issue:**

On every concrete profile, the H=Z gauge-parameter check is computed as a ratio of two integrals
whose integrand strings are identical (`W * Z` for I_WZ and `W * Z` for I_WH with H=Z). The
algebra collapses `X / X = 1` regardless of what the integrator returns, so the check passes
even if `Integrate` fails or returns the wrong value. The paper card explicitly promises
"mutation guards distinguishing H=1 from H=Z" on concrete kernels; the H=1 side is genuine,
but the H=Z side is currently an algebraic identity, not a verification.

Fix structure: introduce a textually-distinct integrand for I_WH_HZ on each profile (re-typed
copy of the Z expression, not the Python/Mathematica variable carrying Z), then assert that
the two integrals collapse to the same closed form BEFORE the ratio check.

**Required change:**

Edit 1 -- SymPy script, matched-Gaussian H=Z (current lines 178-201, around line 186):

Before:
```python
Z = sp.exp(-w**2 / lam**2)
Z_int_gauss = sp.simplify(sp.integrate(Z, (w, -sp.oo, sp.oo)))
Z2_int_gauss = sp.simplify(sp.integrate(Z**2, (w, -sp.oo, sp.oo)))
W_match = sp.simplify(Z / Z_int_gauss)
I_WZ_match = sp.simplify(sp.integrate(W_match * Z, (w, -sp.oo, sp.oo)))
I_WH_H1_match = sp.simplify(sp.integrate(W_match, (w, -sp.oo, sp.oo)))
I_WH_HZ_match = sp.simplify(sp.integrate(W_match * Z, (w, -sp.oo, sp.oo)))
```

After:
```python
Z = sp.exp(-w**2 / lam**2)
Z_int_gauss = sp.simplify(sp.integrate(Z, (w, -sp.oo, sp.oo)))
Z2_int_gauss = sp.simplify(sp.integrate(Z**2, (w, -sp.oo, sp.oo)))
W_match = sp.simplify(Z / Z_int_gauss)
I_WZ_match = sp.simplify(sp.integrate(W_match * Z, (w, -sp.oo, sp.oo)))
I_WH_H1_match = sp.simplify(sp.integrate(W_match, (w, -sp.oo, sp.oo)))
# H=Z computed from a textually-distinct integrand so the integrator must independently
# reproduce the same closed form.
H_Z_match = sp.exp(-w**2 / lam**2)
I_WH_HZ_match = sp.simplify(sp.integrate(W_match * H_Z_match, (w, -sp.oo, sp.oo)))
```

Then insert a new assertion immediately before the existing `assert_zero("Gaussian H=Z gauge parameter", ...)`
line (currently line 201):
```python
assert_zero("Gaussian H=Z integrals match (independent computation)", I_WH_HZ_match - I_WZ_match)
```

Edit 2 -- SymPy script, independent-Gaussian H=Z (current lines 231-241, around line 239):

Before:
```python
I_WZ_indep = sp.simplify(sp.integrate(W_indep * Z, (w, -sp.oo, sp.oo)))
# H = Z:
I_WH_HZ_indep = sp.simplify(sp.integrate(W_indep * Z, (w, -sp.oo, sp.oo)))
xi_eff_HZ_indep = sp.simplify(xi * I_WZ_indep / I_WH_HZ_indep)
assert_zero("independent-profile H=Z gauge parameter", xi_eff_HZ_indep - xi)
```

After:
```python
I_WZ_indep = sp.simplify(sp.integrate(W_indep * Z, (w, -sp.oo, sp.oo)))
# H = Z: re-typed integrand so the integrator must independently reproduce the same closed form.
H_Z_indep = sp.exp(-w**2 / lam**2)
I_WH_HZ_indep = sp.simplify(sp.integrate(W_indep * H_Z_indep, (w, -sp.oo, sp.oo)))
assert_zero("independent-profile H=Z integrals match (independent computation)", I_WH_HZ_indep - I_WZ_indep)
xi_eff_HZ_indep = sp.simplify(xi * I_WZ_indep / I_WH_HZ_indep)
assert_zero("independent-profile H=Z gauge parameter", xi_eff_HZ_indep - xi)
```

Edit 3 -- Mathematica script, M2 Pair A (current lines 28-47):

Before (lines 34-36):
```
gaussGaugeWeight =
  Integrate[gaussianWeight*localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
```

After:
```
gaussGaugeWeight =
  Integrate[gaussianWeight*Exp[-w^2/lambda^2], {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
```

Insert a new equality check immediately before the existing `m2aGaugeResidual = ...` line
(currently line 38):
```
If[FullSimplify[gaussGaugeWeight - gaussOverlap] =!= 0,
  Print["FAIL: M2 Pair A H=Z integrals match"]; Exit[1]
];
Print["PASS: M2 Pair A H=Z integrals independently match"];
```

Edit 4 -- Mathematica script, M2 Pair B (current lines 49-68):

Before (lines 55-57):
```
lorentzGaugeWeight =
  Integrate[lorentzWeight*localizedGaussian, {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
```

After:
```
lorentzGaugeWeight =
  Integrate[lorentzWeight*Exp[-w^2/lambda^2], {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
```

Insert a new equality check immediately before `m2bGaugeResidual = ...` (currently line 59):
```
If[FullSimplify[lorentzGaugeWeight - lorentzOverlap] =!= 0,
  Print["FAIL: M2 Pair B H=Z integrals match"]; Exit[1]
];
Print["PASS: M2 Pair B H=Z integrals independently match"];
```

Edit 5 -- Mathematica script, M4 matched H=Z check (current lines 84-131):

After the existing `matchedDistributedSource = ...` definition block (currently around line 99,
before the `m4Residuals = {...}` list at line 100), add:
```
matchedGaugeOverlap =
  Integrate[matchedWeight*Exp[-w^2/lambda^2], {w, -Infinity, Infinity},
    Assumptions -> lambda > 0];
If[FullSimplify[matchedGaugeOverlap - matchedOverlap] =!= 0,
  Print["FAIL: M4 H=Z matched integrals match"]; Exit[1]
];
Print["PASS: M4 H=Z matched integrals independently match"];
```

Then in the `m4Residuals = {...}` list (currently line 100-108), replace the entry
`FullSimplify[xi*matchedOverlap/matchedOverlap - xi]` (currently line 104) with
`FullSimplify[xi*matchedOverlap/matchedGaugeOverlap - xi]`. Similarly replace the corresponding
`If[FullSimplify[xi*matchedOverlap/matchedOverlap - xi] =!= 0, ...]` block (currently lines 119-121)
with `If[FullSimplify[xi*matchedOverlap/matchedGaugeOverlap - xi] =!= 0, ...]`.

Edit 6 -- Mathematica script, M7 Lorentzian numeric H=Z (current lines 152-171):

Before line 152 (`lorentzMatchedSourceOverlap = ...`), add a fresh independent-integrand
gauge overlap:
```
lorentzGaugeOverlap =
  Integrate[lorentzWeight*Exp[-w^2/lambda^2], {w, -Infinity, Infinity},
    Assumptions -> lambda > 0 && sigma > 0];
```

Then in the `m7GaugeResiduals = ...` block (currently lines 155-159), change
`xi*lorentzOverlap/lorentzGaugeWeight - xi` to `xi*lorentzOverlap/lorentzGaugeOverlap - xi`.
Match the change in the corresponding `If[Chop[m7GaugeResiduals] =!= {0, 0}, ...]` check
(currently lines 165-167) -- the assertion text stays the same since the variable name is
already inside `m7GaugeResiduals`; only the residual definition needs the variable swap.

**Self-test traps checked:**
- Variable independence: no `D[...]`/`sp.diff` involved.
- Symmetry: integrands are products of positive even-w functions; integrals over (-inf, inf) are
  finite and nonzero.
- Trivial-case substitution: matched-Gaussian H=Z integrals are
  `(1/(sqrt(pi)*lambda)) * integrate(exp(-2*w^2/lambda^2), -inf, inf) = sqrt(2)/2` for both
  textually distinct integrands. Lorentzian H=Z integrals close to a Voigt-like form that
  Mathematica resolves the same way for both integrands. Both fix branches preserve the
  existing PASS values.
- Path specifications: both files exist; no new file paths to introduce.
- Paper round-trip: the fix neither introduces nor changes a constant or symbolic claim. The
  H=Z weight `H(w) = Z(w)` is preserved exactly; only the implementation textually duplicates
  the integrand so the integrator's output is independently verified to match.

**Verification command:**

After Codex applies, the verifier runs `redteam exec-sympy 008` and `redteam exec-mathematica 008`
and confirms (a) both scripts still exit 0 with `STATUS: PASS`, (b) new "H=Z integrals match"
print lines appear at the locations specified above, and (c) the H=Z ratio assertions still
pass (now substantively rather than by algebraic cancellation).

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage008_projected_maxwell_extension_mathematica_audit.wl`
  - `scripts/moving_throat_pde_stage008_projected_maxwell_extension_sympy_audit.py`
- summary: Added textually distinct H=Z gauge-overlap integrands and pre-ratio equality guards for the concrete Mathematica and SymPy profile checks.
- deviation: none
