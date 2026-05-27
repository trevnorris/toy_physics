---
unit_id: 066
batch: III.3
created_at: 2026-05-26T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T18:52:08Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 066

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py:66-74`

**Issue:** Under the `MONOTONICITY` banner the SymPy script only `print`s four partial derivatives (`dW/dV0^2, dW/da, dW/dL, dW/dell`) without asserting their signs, and entirely omits `dW/dJ_1` and `dW/dT_X`. Notes §3 lists six signed monotonicity directions as a stage deliverable; the Mathematica engine (`.wl` lines 73-78) asserts all six with `expectTrue`. The SymPy engine must independently certify the same six signs so the closure status `\StatusExactClosure{}` is honored by both engines.

**Required change:**
Replace the current block at lines 66-74 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py` with the following. This adds two missing derivative prints (`dW/dJ_1`, `dW/dT_X`) and six sign assertions. Do not change any line outside this range.

Before (lines 66-74):

```python
banner("MONOTONICITY")

# W_wall is manifestly monotone in V0^2, a^2, L^2, J1 and inverse ell.
Vp = sp.symbols("Vp", positive=True, real=True)
W_Vp = sp.simplify(W_wall.subs(V0**2, Vp))
print("dW/d(V0^2) =", sp.diff(W_Vp, Vp))
print("dW/da =", sp.diff(W_wall, a))
print("dW/dL =", sp.diff(W_wall, L))
print("dW/dell =", sp.diff(W_wall, ell))
```

After (replace lines 66-74 with):

```python
banner("MONOTONICITY")

# W_wall is manifestly monotone in V0^2, a^2, L^2, J1 and inverse ell, T_X
# (notes section 3: six signed directions).
Vp = sp.symbols("Vp", positive=True, real=True)
W_Vp = sp.simplify(W_wall.subs(V0**2, Vp))

dW_dV0sq = sp.simplify(sp.diff(W_Vp, Vp))
dW_da    = sp.simplify(sp.diff(W_wall, a))
dW_dL    = sp.simplify(sp.diff(W_wall, L))
dW_dell  = sp.simplify(sp.diff(W_wall, ell))
dW_dJ1   = sp.simplify(sp.diff(W_wall, J1))
dW_dTX   = sp.simplify(sp.diff(W_wall, TX))

print("dW/d(V0^2) =", dW_dV0sq)
print("dW/da =", dW_da)
print("dW/dL =", dW_dL)
print("dW/dell =", dW_dell)
print("dW/dJ1 =", dW_dJ1)
print("dW/dT_X =", dW_dTX)

assert sp.simplify(dW_dV0sq > 0) is sp.true, "dW/d(V0^2) should be positive"
assert sp.simplify(dW_da    > 0) is sp.true, "dW/da should be positive"
assert sp.simplify(dW_dL    > 0) is sp.true, "dW/dL should be positive"
assert sp.simplify(dW_dell  < 0) is sp.true, "dW/dell should be negative"
assert sp.simplify(dW_dJ1   > 0) is sp.true, "dW/dJ1 should be positive"
assert sp.simplify(dW_dTX   < 0) is sp.true, "dW/dT_X should be negative"
```

Notes for Codex:
- All symbols `V0, ell, a, L, J1, T_X` are declared `positive=True, real=True` at lines 47-49, and `Vp` is declared `positive=True, real=True` in this block. Under those assumptions SymPy resolves `simplify(monomial > 0)` to `sp.true` for manifestly-positive rational monomials, and `simplify(-monomial < 0)` similarly. The `is sp.true` comparison is the correct guard — a `Relational` that fails to simplify will *not* be `is sp.true`, and the assert will fire.
- Do not change the docstring header (the historical "Stage 49" label is a separate naming-drift issue outside this finding's scope).
- Do not modify the CONSTANT-COMPRESSIBILITY section (lines 76 onward).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 066` and confirm:
1. The script exits 0.
2. The output now contains lines `dW/dJ1 = 4*pi*L**2*V0**2*a**2/(T_X*ell)` (or sympy's canonical reordering thereof) and `dW/dT_X = -4*pi*J1*L**2*V0**2*a**2/(T_X**2*ell)`.
3. The four existing derivative print lines remain unchanged (modulo sympy's canonical ordering).
4. No new failures appear elsewhere in the transcript.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py`
- summary: Added the missing J1 and T_X monotonicity derivative prints and asserted all six required monotonicity signs.
- deviation: none
