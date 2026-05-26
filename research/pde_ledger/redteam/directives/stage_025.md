---
unit_id: 025
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-26T05:34:18Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 025

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:151`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:108`

**Issue:** `N0` is built at sympy line 80 (`N0 = sp.simplify(P**2 / Delta**2)`) and at math line 46 (`n0 = FullSimplify[p^2/delta^2, ...]`). The assertion at sympy line 151 / math line 108, `expect_zero("N0 - P^2/Delta^2", N0 - P^2/Delta^2)`, is structurally `(P^2/Delta^2) - P^2/Delta^2`, which is zero by construction and cannot fail. The check therefore exercises nothing.

**Required change:**

Replace the tautological check at each location with a substantive structural re-derivation that rebuilds `N_0` from the primitive raw symbols (OmegaU, OmegaW, GU, GW, R) without referring to the cached `P` or `Delta` shorthand. The replacement check then verifies that the cached pair `(P, Delta)` consistently reconstructs N_0 from the raw frequencies/couplings.

**Sympy edit** (file: `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`, line 151):

Replace the single line

```python
    expect_zero("N0 - P^2/Delta^2", N0 - P**2 / Delta**2)
```

with

```python
    P_raw = OmegaU**2 * GW + R * GU
    Delta_raw = OmegaU**2 * OmegaW**2 - R**2
    expect_zero("N0 reconstructed from raw symbols", N0 - P_raw**2 / Delta_raw**2)
```

This builds `P_raw` and `Delta_raw` locally from the primitive declared symbols `OmegaU, OmegaW, GU, GW, R` (no reference to the cached `P, Delta` returned by `zero_frequency_coefficients`), then checks that `N0` matches `P_raw**2 / Delta_raw**2`. If a future edit ever decoupled the cached `P, Delta` from their definitions, this check would fail.

**Mathematica edit** (file: `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`, line 108):

Replace the single line

```mathematica
  expectZero["N0 - P^2/Delta^2", n0 - p^2/delta^2];
```

with

```mathematica
  Module[{pRaw, deltaRaw},
    pRaw = omegaU^2*gW + r*gU;
    deltaRaw = omegaU^2*omegaW^2 - r^2;
    expectZero["N0 reconstructed from raw symbols", n0 - pRaw^2/deltaRaw^2];
  ];
```

Same logic: `pRaw` and `deltaRaw` are local rebuilds from the primitive declared symbols `omegaU, omegaW, gU, gW, r`; the residual `n0 - pRaw^2/deltaRaw^2` is non-trivially zero only if `n0` was correctly constructed from those primitives.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 025` and `redteam exec-mathematica 025`. Both must still exit 0 and the saved transcripts must contain the new check name `N0 reconstructed from raw symbols = 0` followed by `PASS: ...` (Mathematica) — the old `N0 - P^2/Delta^2` line should no longer appear.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Replaced the tautological cached N0 check with primitive-symbol reconstructions of P and Delta in both engines.
- deviation: none

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py` (insert at end of `normalization_formula()`, after the existing II.2 block — concretely after line 138)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl` (insert at end of `normalizationFormula[]` Module, before the `]` and `;` that close it — concretely insert after line 100, before line 101's `];`)

**Issue:** The paper card's `\stagefield{Checks}` list (`/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_025.tex:94`) enumerates four corollary checks, the fourth being "If P = 0, then N_0 = 0 and a positive quadrupole target cannot be reached." Neither engine verifies this corollary. The scripts must honor each paper-enumerated check.

**Required change:**

Add a new subsection II.3 (a block at the end of `normalization_formula` / `normalizationFormula`) that:
1. Substitutes a symbolic choice that forces `P = 0` — concretely, set `GW = -R*GU/OmegaU**2`, which gives `P = OmegaU**2 * (-R*GU/OmegaU**2) + R*GU = -R*GU + R*GU = 0`.
2. Asserts `N0` evaluates to zero under that substitution.
3. Asserts `(mhat**2 * P0_compact - target)` evaluates to `-target` under that substitution (so the LHS - RHS is strictly negative for any positive mhat — i.e., no positive-mhat solution).

**Sympy edit** (file: `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`, insert immediately after line 138, still inside `normalization_formula()`):

```python

    subbanner("II.3 — P = 0 corollary (paper Checks item 4)")
    # Forcing GW = -R*GU/OmegaU^2 gives P = OmegaU^2*GW + R*GU = 0 symbolically.
    P_zero_sub = {GW: -R * GU / OmegaU**2}
    N0_at_Pzero = sp.simplify(N0.subs(P_zero_sub))
    print(f"N0 at P=0 = {N0_at_Pzero}")
    expect_zero("N0 vanishes when P=0", N0_at_Pzero)
    residual_at_Pzero = sp.simplify((mhat**2 * P0_compact - target).subs(P_zero_sub))
    print(f"(mhat^2*P0 - target) at P=0 = {residual_at_Pzero}")
    expect_zero("(mhat^2*P0 - target) at P=0 equals -target", residual_at_Pzero + target)
```

Note: `expect_zero("(mhat^2*P0 - target) at P=0 equals -target", residual_at_Pzero + target)` checks that `residual = -target`, which is non-trivially true (it would fail if N_0 did not vanish or if `target` were not the constant the script claims it is). N0 and P0_compact are already in scope at this point (assigned at lines 101-103 within `normalization_formula()`); `target` is assigned at line 123.

**Mathematica edit** (file: `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`, insert just before the closing `];` of `normalizationFormula[]` at line 101):

```mathematica

  subbanner["II.3: P = 0 corollary (paper Checks item 4)"];
  Module[{pZeroSub, n0AtPzero, residualAtPzero},
    pZeroSub = {gW -> -r*gU/omegaU^2};
    n0AtPzero = FullSimplify[n0 /. pZeroSub, Assumptions -> $Assumptions];
    Print["N0 at P=0 = ", fmt[n0AtPzero]];
    expectZero["N0 vanishes when P=0", n0AtPzero];
    residualAtPzero = FullSimplify[(mhat^2*p0Compact - target) /. pZeroSub, Assumptions -> $Assumptions];
    Print["(mhat^2*P0 - target) at P=0 = ", fmt[residualAtPzero]];
    expectZero["(mhat^2*P0 - target) at P=0 equals -target", residualAtPzero + target];
  ];
```

The local `Module` guards the new symbol bindings so they do not leak. `n0`, `p0Compact`, and `target` are already in scope (assigned earlier in the outer `normalizationFormula[]` Module).

**Trivial-case pre-check (mental):** At `GW = -R*GU/OmegaU**2`:
- `P = OmegaU**2 * GW + R * GU = OmegaU**2 * (-R*GU/OmegaU**2) + R*GU = -R*GU + R*GU = 0`.
- `N0 = P**2 / Delta**2 = 0 / Delta**2 = 0` (Delta = OmegaU^2 OmegaW^2 - R^2 is nonzero generically).
- `P0_compact = P**2 / (Delta * (K*Delta - Delta*C^2/varpi^2 - Q)) = 0` (numerator vanishes).
- `mhat**2 * P0_compact - target = 0 - target = -target = -54 G c_s^5 / (5 a^5 c^5)`.
- `residual_at_Pzero + target = -target + target = 0`. Non-trivially zero — would fail if `target` were not its declared form or if `P0_compact` did not vanish when P=0.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 025` and `redteam exec-mathematica 025`. Both must still exit 0; the saved transcripts must contain the new lines `II.3 — P = 0 corollary` and `PASS: N0 vanishes when P=0` and `PASS: (mhat^2*P0 - target) at P=0 equals -target` (Mathematica), or the equivalent print/assert flow (SymPy).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Added the II.3 P=0 corollary checks showing N0 vanishes and the normalization residual equals -target.
- deviation: none
