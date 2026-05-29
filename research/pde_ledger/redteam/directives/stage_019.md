---
unit_id: 019
batch: I.2
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T04:07:07Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 019

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line range named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — symbol_assumption_error

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py:25-29`

**Issue:** Every symbol is declared with only `nonzero=True`. The physical setup is real (wall constants, projector moments, target-surface parameters) and the dimensionful prefactors `G, c_s, a, c, \widehat m_0` are physically positive. The Mathematica mirror at `mathematica/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_mathematica_audit.wl:62-69` correctly declares `Element[{...}, Reals]`, so the two engines verify the same identities under different domain assumptions. All current assertions pass because the residuals cancel without needing positivity, but the latent risk is that a future `simplify`/`radsimp`/`sqrt`-manipulation edit that does depend on sign would silently produce a wrong simplification and the assertion could pass for the wrong reason.

**Required change:**
Replace the symbol declarations at lines 25-29 (inside `main`). Before:

```python
    KSigma, MSigma = sp.symbols('KSigma MSigma', nonzero=True)
    B0, B2, B4, Z0, Z2, Z4 = sp.symbols('B0 B2 B4 Z0 Z2 Z4', nonzero=True)
    N0, N2, N4 = sp.symbols('N0 N2 N4', nonzero=True)
    mhat0, G, cs, a, c = sp.symbols('mhat0 G cs a c', nonzero=True)
    eps = sp.symbols('eps', nonzero=True)
```

After:

```python
    KSigma, MSigma = sp.symbols('KSigma MSigma', real=True, nonzero=True)
    B0, B2, B4, Z0, Z2, Z4 = sp.symbols('B0 B2 B4 Z0 Z2 Z4', real=True, nonzero=True)
    N0, N2, N4 = sp.symbols('N0 N2 N4', real=True, nonzero=True)
    mhat0, G, cs, a, c = sp.symbols('mhat0 G cs a c', positive=True)
    eps = sp.symbols('eps', real=True, nonzero=True)
```

Do not change KSigma/MSigma/B*/Z*/N* to `positive=True`. The M-root analysis genuinely sweeps both signs of `MSigma` (the positive and negative one-pole branches), and the response-sign criterion at lines 116-173 must not pre-commit `D_0>0` into the algebra — that positivity is what those numeric samples are testing.

Make no other edits. Do not modify any assertion text, any closed-form definition, any sample dictionary, or any `lines.append(...)` print line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 019` and confirm: (a) the script exits 0; (b) every existing assertion still passes (no `AssertionError` raised); (c) the saved-output line `positive-root numeric u2 = 0.6324555320336759` remains identical (numeric value unchanged); (d) every `... = PASS` line and the trailing `STATUS: PASS` remain present. The saved-output line `u2_on_positive_root = ...` is permitted to change form (the sign cosmetics may collapse to a simpler `sqrt`-of-positive expression) — this is the expected effect of the added `real=True`/`positive=True` assumptions.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage019_parent_throat_action_isotropic_bundle_sympy_audit.py`
- summary: Updated the stage 019 SymPy symbol assumptions to use real domains for algebraic symbols and positive domains for the physical prefactors.
- deviation: none
