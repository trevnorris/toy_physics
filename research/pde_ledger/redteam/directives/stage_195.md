---
unit_id: 195
batch: V.3
created_at: 2026-06-01T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-01T11:49:55-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 195

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the script (`python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py`) and iterate until it exits 0 with all in-file checks passing.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py:89`

**Issue:** The stage's headline `\stagefield{Output}` is the boxed factorization `m̂₀²χ_Q N_Q = 1` (notes section 2, lines 130–141: derived from the observable odd condition `m̂₀²Γ̄₅ = Γ̄₅^target` plus the verified ratio `Γ̄₅/Γ̄₅^target = χ_Q N_Q`). The current load-bearing assertion at line 89 is `expect_zero("odd closure factorization", odd_closure - (mhat0**2 * chi_Q * N_Q - 1))`, where `odd_closure = sp.simplify(mhat0**2 * chi_Q * N_Q - 1)` (line 76). This is `X - X == 0`: it can never fail and never touches `Gamma5` or `Gamma5_target`. The named deliverable is asserted, not derived.

**Required change:**
Replace the tautological assertion at line 89 with a derivation of the factorization from the observable odd condition. Concretely, change line 89 from:

```python
expect_zero("odd closure factorization", odd_closure - (mhat0**2 * chi_Q * N_Q - 1))
```

to:

```python
odd_condition_residual = sp.simplify(
    (mhat0**2 * Gamma5 - Gamma5_target).subs(P0, N_Q * P0_target)
)
expect_zero(
    "observable odd condition factorizes as Gamma5_target*(mhat0^2 chi_Q N_Q - 1)",
    odd_condition_residual - Gamma5_target * (mhat0**2 * chi_Q * N_Q - 1),
)
```

`Gamma5`, `Gamma5_target`, `P0`, `N_Q`, `P0_target`, `chi_Q`, `mhat0` are all already defined above (lines 42–52). You may keep the `odd_closure` definition (line 76) and its `print`/`pprint` (lines 82–83) as display; just do not let the tautological subtraction be the verification. Do not add any new constant.

**Self-test (must hold before you trust the run):** `mhat0**2 * Gamma5` with `Gamma5 = chi_Q*a**5*P0/(27*c_s**5)` and `P0 = N_Q*P0_target = N_Q*54*G*c_s**5/(5*a**5*c**5)` equals `mhat0**2 * chi_Q * N_Q * 2*G/(5*c**5) = mhat0**2 * chi_Q * N_Q * Gamma5_target`, so `odd_condition_residual = Gamma5_target*(mhat0**2*chi_Q*N_Q - 1)` and the new `expect_zero` is identically zero — but only because the carried normalization is correct, not by construction.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 195` and confirms the new check line ("observable odd condition factorizes as ...") appears in the output and the script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py`
- summary: Replaced the tautological odd-closure subtraction with a residual derived from the observable odd condition after substituting `P0 = N_Q * P0_target`.
- deviation: none

## F2 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py:63` and `:88`

**Issue:** Two secondary `expect_zero` calls are definitional echoes that cannot fail.
- Line 63: `expect_zero("N_Q definition", N_Q_def - P0 / P0_target)` with `N_Q_def = sp.simplify(P0 / P0_target)` (line 52) → `simplify(P0/P0_target) - P0/P0_target == 0`.
- Line 88: `expect_zero("Delta_norm - P0_target*(mhat0^2 N_Q - 1)", Delta_norm_NQ - P0_target * (mhat0**2 * N_Q - 1))` with `Delta_norm_NQ = sp.simplify(Delta_norm.subs(P0, N_Q * P0_target))` (line 75) and `Delta_norm = mhat0**2*P0 - P0_target` (line 74) → `Delta_norm_NQ` IS `P0_target*(mhat0**2*N_Q - 1)`, so the subtraction is identically zero.

**Required change:**
Delete the two `expect_zero` calls:
- Remove line 63 entirely: `expect_zero("N_Q definition", N_Q_def - P0 / P0_target)`.
- Remove line 88 entirely: `expect_zero("Delta_norm - P0_target*(mhat0^2 N_Q - 1)", Delta_norm_NQ - P0_target * (mhat0**2 * N_Q - 1))`.

Leave the surrounding `print`/`pprint` display lines untouched (lines 54–61 in section I, and lines 80–87 in section II). Do not add any replacement assertion at these two spots; F1 already establishes the substantive section-II check. Do not alter the `N_Q_def`/`Delta_norm_NQ` definitions (they are still used downstream: `Delta_norm_NQ` is consumed at lines 78 and 105).

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 195` and confirms the output no longer contains the lines "N_Q definition = 0" and "Delta_norm - P0_target*(mhat0^2 N_Q - 1) = 0", and the script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py`
- summary: Removed the two definitional `expect_zero` checks while leaving the surrounding display and downstream definitions intact.
- deviation: none

## F3 — missing_mathematica

**Issue:** Stage 195 is dual-engine-capable (exact rational ratios, normalization-condition factorization, first-order Taylor series) but has no Mathematica `.wl`. Under the dual-engine rule, an independent second-engine verification is required wherever Mathematica can do the math.

**Required change (you design the route and write the script):**
Create `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_mathematica_audit.wl`.
- Independently re-verify EVERY load-bearing assertion in the CORRECTED SymPy script (after F1/F2 above) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_sympy_audit.py`. Read that script to enumerate the claims and their target conclusions; the paper card `paper/stages/stage_195.tex` and the notes file are the source of truth. In particular mirror the CORRECTED headline factorization (derive `m̂₀²χ_Q N_Q = 1` from the observable odd condition `m̂₀²Γ̄₅ = Γ̄₅^target`, NOT the deleted X−X self-echo), not the broken original.
- Use Mathematica-NATIVE primitives (`Series`+`Coefficient`, `Factor`/`Together`, `Solve`/`Reduce`, `D[...]`) via a DIFFERENT derivation route than the SymPy script — NOT a line-by-line port with the same variable names and step order. Reference an existing verified `.wl` (e.g. `mathematica/moving_throat_pde_stage187_*_mathematica_audit.wl`) ONLY for house idioms (`expectZero`, `$Assumptions` positivity, `stripCE`, the `math -script` convention).
- Assert cross-engine agreement: each conclusion must match the SymPy result.

**Anti-transliteration:** a `.wl` that merely re-types the SymPy closed forms and subtracts them is a transliteration and will be REJECTED at verification. Design a genuinely independent route. RUN it (`timeout 600 math -script <path>`) and iterate to exit 0; a timeout (124) is a failure — reformulate, don't raise the cap.

**Verification command:** the verifier runs `redteam exec-mathematica 195`, confirms exit 0 with all PASS lines, and reviews that the `.wl` is a genuinely independent route whose conclusions agree with the SymPy engine.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage195_source_map_reduction_of_canonical_outgoing_branch_mathematica_audit.wl`
- summary: Added an independent Mathematica audit that rederives the corrected source-map factorization, Packet-A collapse, deformation laws, linearization, and canonical closure.
- deviation: none
