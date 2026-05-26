---
unit_id: 042
batch: III.1
created_at: 2026-05-26T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T07:39:05Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 042

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py:92`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl:70`

**Issue:**

The current `F_general` (sympy line 92) / `fGeneral` (mathematica line 70) is constructed from the *expected* closed-form overlap expressions (`Z_expected`, `S_expected` / `zExpected`, `sExpected`), then compared against `F_expected` (sympy line 96-100) / `fExpected` (mathematica line 71-75). Since `F_expected` is literally `Z_expected * S_expected / (1-xi)` after substitution, the residual `F_general - F_expected` is identically zero by construction. The assertion at sympy line 101 / mathematica line 82 cannot fail regardless of whether the underlying physical claim holds, so the F_(q,r,t) deliverable from notes Section 3 is not actually exercised.

The substantive checks for the Z and S overlaps remain at sympy lines 89-90 / mathematica lines 80-81 and are non-tautological; the F_general consequence check should be routed through the eigenvector-constructed forms `Z_overlap` / `zOverlap` and `S_overlap` / `sOverlap` so that the multiplication chains the simplifications established by A3/A4 (B2a/B2b) through the product.

**Required change:**

### SymPy script

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`, change line 92.

Before (line 92):
```python
F_general = sp.simplify(Z_expected * S_expected / (1 - xi))
```

After (line 92):
```python
F_general = sp.simplify(Z_overlap * S_overlap / (1 - xi))
```

Do not change lines 93-101 (the print, factor, `F_expected` definition, and `expect_zero` call are correct as-is — they will now form a non-tautological check because `F_general` is built from the constructed-from-eigenvector overlaps rather than the named closed forms).

### Mathematica script

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl`, change line 70.

Before (line 70):
```
fGeneral = FullSimplify[zExpected sExpected/(1 - xi), Assumptions -> $Assumptions];
```

After (line 70):
```
fGeneral = FullSimplify[zOverlap sOverlap/(1 - xi), Assumptions -> $Assumptions];
```

Do not change lines 71-82 (the `fExpected` definition, prints, and `expectZero[fGeneral - fExpected]` call are correct as-is).

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 042` and `redteam exec-mathematica 042` and confirm both scripts:

1. Exit 0.
2. Still print `F_general - expected = 0` (sympy) and `fGeneral - fExpected = 0 / PASS: F_general - expected` (mathematica). The mathematical content of the assertion is preserved — `Z_overlap = Z_expected` and `S_overlap = S_expected` are independently established by A3/A4 (B2a/B2b), so `Z_overlap * S_overlap / (1-xi) = Z_expected * S_expected / (1-xi) = F_expected` holds, but only after the simplifier chains through the eigenvector-form residuals. The check is now non-tautological because if Z_overlap or S_overlap did not reduce to their expected closed forms (e.g. if the eigenvector ratio were wrong or `D_qr` were misstated), the residual `F_general - F_expected` would not simplify to zero.

No other assertions, expressions, or saved outputs need to change.

## Applied: F1

- files_changed:
  - `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`
  - `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage042_rank2_selected_mode_mathematica_audit.wl`
- summary: Routed the F_general/fGeneral product through the eigenvector-derived overlap expressions instead of the expected closed-form overlaps.
- deviation: none
