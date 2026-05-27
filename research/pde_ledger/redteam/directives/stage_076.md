---
unit_id: 076
batch: III.4
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 076

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:61-65`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl:54-57`

**Issue:** The "Theta_w vs alternative-form derivation" check compares `Theta_w = 4 rho_w^2 mu_star^2 / (hbar^2 csw^2)` against `Theta_w_alt = (2 rho_w mu_star / (hbar csw))^2`. The two expressions differ only by the algebraic identity `(2 a)^2 = 4 a^2`; both use the same `mu_star_solved` symbol. The check does not exercise the enthalpy-lock factor `1/4` in `mu_* = lambda_mu m c_sw^2 / 4`. Replace it with a comparison against the closed-form target stated in notes §2.

**Required change (SymPy `.py`):**

Replace lines 61–65 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py`:

Before:
```python
Theta_w = sp.simplify(4 * rho_w**2 * mu_star_solved**2 / (hbar**2 * csw**2))
# Independent route: Theta_w as (2 rho_w mu_star / (hbar c_sw))^2
Theta_w_alt = sp.simplify((2 * rho_w * mu_star_solved / (hbar * csw))**2)
print("Theta_w (enthalpy lock) =", Theta_w)
expect_zero("Theta_w vs alternative-form derivation", Theta_w - Theta_w_alt)
```

After:
```python
Theta_w = sp.simplify(4 * rho_w**2 * mu_star_solved**2 / (hbar**2 * csw**2))
# Closed-form target from notes section 2: Theta_w = lambda_mu^2 m^2 rho_w^2 c_sw^2 / (4 hbar^2).
# This independently states the simplified form; the assertion below exercises the /4 factor
# in the enthalpy lock (mu_* = lambda_mu * m * c_sw^2 / 4).
Theta_target = sp.Rational(1, 4) * lambda_mu**2 * m**2 * rho_w**2 * csw**2 / hbar**2
print("Theta_w (enthalpy lock) =", Theta_w)
expect_zero("Theta_w under enthalpy lock", Theta_w - Theta_target)
```

**Required change (Mathematica `.wl`):**

Replace lines 54–57 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl`:

Before:
```
thetaW = FullSimplify[4*rhoW^2*muStarSolved^2/(hbar^2*cSw^2), Assumptions -> $Assumptions];
thetaCanonical = FullSimplify[(2*rhoW*muStarSolved/(hbar*cSw))^2, Assumptions -> $Assumptions];
Print["Theta_w (enthalpy lock) = ", fmt[thetaW]];
expectZero["Theta_w vs alternative-form derivation", thetaW - thetaCanonical];
```

After:
```
thetaW = FullSimplify[4*rhoW^2*muStarSolved^2/(hbar^2*cSw^2), Assumptions -> $Assumptions];
(* Closed-form target from notes section 2; the assertion exercises the 1/4 in the enthalpy lock. *)
thetaTarget = FullSimplify[(1/4)*lambdaMu^2*mpsi^2*rhoW^2*cSw^2/hbar^2, Assumptions -> $Assumptions];
Print["Theta_w (enthalpy lock) = ", fmt[thetaW]];
expectZero["Theta_w under enthalpy lock", thetaW - thetaTarget];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 076` and `redteam exec-mathematica 076`. The renamed check "Theta_w under enthalpy lock" should appear in both output files with residual `0`. Both scripts must still exit 0.

## F2 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:77-80`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl:68-71`

**Issue:** `ref_factor = 1/20` is hardcoded as the load-bearing input that produces the headline `25`. The current inline comment carries an undischarged `TODO(provenance)` token; the upstream anchor (notes §4) is unambiguous, so the TODO can be replaced with a direct citation.

**Required change (SymPy `.py`):**

Replace lines 77–80 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py`:

Before:
```python
# Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20.
# TODO(provenance): cite the upstream stage that fixes ref_factor. This factor is
# the load-bearing piece of the "25" in the normalized reference identity.
ref_factor = sp.Rational(1, 20)  # reference-branch convention: ell = a * ref_factor  (see F2 below for provenance)
```

After:
```python
# Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20.
# Source: Family-1 reference-branch description carried forward as input to this stage
# (notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md section 4).
# This factor is the load-bearing piece of the "25" in the normalized reference identity.
ref_factor = sp.Rational(1, 20)
```

**Required change (Mathematica `.wl`):**

Replace lines 68–71 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl`:

Before:
```
(* Reference-branch convention: ell = a * refFactor with refFactor = 1/20.
   TODO(provenance): cite the upstream stage that fixes refFactor. This factor is
   the load-bearing piece of the "25" in the normalized reference identity. *)
refFactor = 1/20;  (* reference-branch convention: ell = a * refFactor  (see F2 below for provenance) *)
```

After:
```
(* Reference-branch convention: ell = a * refFactor with refFactor = 1/20.
   Source: Family-1 reference-branch description carried forward as input to this stage
   (notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md section 4).
   This factor is the load-bearing piece of the "25" in the normalized reference identity. *)
refFactor = 1/20;
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 076` and `redteam exec-mathematica 076`. The numeric output should be byte-identical to the current transcripts; only the comment text changes. Both scripts must still exit 0.
