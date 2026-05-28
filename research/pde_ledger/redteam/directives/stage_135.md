---
unit_id: 135
batch: IV.4
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 135

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what is specified. Do NOT touch paper.tex, notes/, or any prose documents. Do NOT run python or mathematica.

The two findings below are tightly coupled: F1 says "the only assertion is tautological"; F2 says "and the substantive deliverables are not tested." The required edit for both is a single edit batch on the SymPy script — apply F1's edit and F2's edit together as one revision of `scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py`. The Mathematica script does not require any edit.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:46-49`

**Issue:** The only SymPy assertion is the closure residual `Pi_* - Sigma_*(4 - S_*)`, where `Sigma_*` was just defined as `solve(Pi_* == sigma*(4 - S_*), sigma)`. The residual is identically zero by construction (modulo floating-point noise), so the `raise AssertionError` cannot fire. The check verifies arithmetic on `solve`, not the physics.

**Required change:**

Keep the residual print as a numerical sanity probe, but the load-bearing assertion is moved into F2's anchored checks. Concretely: leave lines 46-47 as-is (they print the residual for the transcript), and delete lines 48-49:

```python
if abs(residual) > 1e-12:
    raise AssertionError("Outlet-consistent threshold did not close.")
```

The replacement assertions are added per F2 below. Do not silently drop the residual print — only drop the conditional `raise`.

**Verification:** After applying, lines 48-49 of the prior revision (the `if abs(residual)...` block) are gone, and the residual is still printed at line 46-47. The new assertions from F2 are the only `assert`/`raise` statements in the script.

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:19-49`

**Issue:** The SymPy script does not assert any of the paper's six deliverables (the substitution `M_s + M_q*S_q -> Sigma_m*(4 - S_q)`; `0 < S_q(Pi_*) < 1`; numerical Sigma_m^*, M_s^*, M_q^*; mixed-lane correction). It only prints them. The Mathematica script has all six checks at lines 56, 67, 74-77; SymPy must reach parity.

**Required change:**

1. **Add a symbolic substitution check** immediately after the `S_q = sp.simplify(...)` line (currently line 28). Introduce two independent symbols `M_s_sym, M_q_sym` (positive=True is wrong for `M_q`, since `M_q` is negative — declare them as plain real symbols), then build the generic law and substitute:

   ```python
   M_s_sym, M_q_sym = sp.symbols("M_s_sym M_q_sym", real=True)
   generic_law = M_s_sym + M_q_sym * S_q
   reduced_law = generic_law.subs({M_s_sym: 4*Sigma_m, M_q_sym: -Sigma_m})
   expected_law = Sigma_m * (4 - S_q)
   residual_sub = sp.simplify(reduced_law - expected_law)
   print("M_s + M_q*S_q - Sigma_m*(4 - S_q) =", residual_sub)
   assert residual_sub == 0, f"Outlet substitution failed: residual = {residual_sub}"
   ```

   Insert this block between the existing `S_q = sp.simplify(...)` and the existing `Pi_eq = sp.simplify(Sigma_m * (4 - S_q))` (current line 29).

2. **Assert the inequality** `0 < S_q(Pi_*) < 1`. Replace the current line 41:

   ```python
   print("0 < S_q(Pi_star) < 1 ->", bool(0 < S_star < 1))
   ```

   with:

   ```python
   s_in_range = bool(0 < S_star < 1)
   print("0 < S_q(Pi_star) < 1 ->", s_in_range)
   assert s_in_range, f"S_q(Pi_*) out of range: S_star = {S_star}"
   ```

3. **Assert numerical anchors** against the notes-stated values. After the current `M_q_star = sp.N(-Sigma_star, 30)` (line 38) and the existing prints (lines 40-44), add:

   ```python
   Sigma_target = sp.Float("0.451485277739090", 30)
   M_s_target = sp.Float("1.80594111095636", 30)
   M_q_target = sp.Float("-0.451485277739090", 30)

   assert abs(Sigma_star - Sigma_target) < sp.Float("1e-12", 30), \
       f"Sigma_m^* mismatch: {Sigma_star} vs {Sigma_target}"
   assert abs(M_s_star - M_s_target) < sp.Float("1e-11", 30), \
       f"M_s^* mismatch: {M_s_star} vs {M_s_target}"
   assert abs(M_q_star - M_q_target) < sp.Float("1e-12", 30), \
       f"M_q^* mismatch: {M_q_star} vs {M_q_target}"
   print("Sigma_m^*, M_s^*, M_q^* anchored to notes values within tolerance.")
   ```

4. **Compute and assert the mixed-lane correction**. After step 3 above, add:

   ```python
   mixed_correction = sp.N(M_q_star * S_star, 30)
   mixed_target = sp.Float("-0.297111597463199", 30)
   print("M_q^* * S_q(Pi_*) =", mixed_correction)
   assert abs(mixed_correction - mixed_target) < sp.Float("1e-11", 30), \
       f"mixed-lane correction mismatch: {mixed_correction} vs {mixed_target}"
   ```

   The notes file (line 125-127) gives the reference value `-0.297111597463199151`; using 15 digits with tolerance 1e-11 is safely inside that precision.

5. Delete the conditional `raise` at lines 48-49 (per F1). The residual print at lines 46-47 stays.

After all five steps, the script will have five new substantive `assert` statements plus the existing residual print.

**Claim manifest** (so the verifier can confirm coverage):
- M1: `M_s + M_q*S_q` with `M_s=4Σ, M_q=-Σ` simplifies to `Σ*(4 - S_q)` — symbolic identity, step 1 above.
- M2: `0 < S_q(Pi_*) < 1` — boolean assertion, step 2 above.
- M3: `Sigma_m^* ≈ 0.451485277739090` — numeric anchor, step 3 above.
- M4: `M_s^* ≈ 1.80594111095636` — numeric anchor, step 3 above.
- M5: `M_q^* ≈ -0.451485277739090` — numeric anchor, step 3 above.
- M6: `M_q^* * S_q(Pi_*) ≈ -0.297111597463199` — numeric anchor, step 4 above.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 135` and confirm: (a) the new check at step 1 prints `M_s + M_q*S_q - Sigma_m*(4 - S_q) = 0` and the `assert` does not raise; (b) lines for M2-M6 each print and assert; (c) the script exits 0; (d) the residual print at lines 46-47 of the original remains visible in the new output.
