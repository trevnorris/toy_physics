---
unit_id: 055
batch: III.2
created_at: 2026-05-22T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-22T17:40:42-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 055

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:56`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl:54`

**Issue:**
Both scripts assert a "KX/KW equivalence" check that is algebraically guaranteed by the prior `x_floor = 4 - pi^2/zeta_req` definition. The LHS form `1 - x/4` is typed in by hand rather than derived from `A_K`, so the assertion's residual is zero by pure substitution on the closed form of `x_floor`. Re-anchor the check to `A_K` (via `1/A_K` at `y = 0`) so a regression in either `A_K`'s definition or `x_floor`'s closed form will trip the assertion.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`, replace line 56 exactly.

Before (line 56):
```
expect_zero("KX/KW equivalence", (1 - x / 4).subs(x, x_floor) - pi**2 / (4 * zeta_req))
```

After (line 56):
```
expect_zero("KX/KW equivalence", (1 / AK).subs(y, 0).subs(x, x_floor) - pi**2 / (4 * zeta_req))
```

Do not modify any other line. In particular, leave line 55 (`KX_over_KW = sp.symbols(...)`) and the comment on line 54 (`# Equivalent stiffness-ratio form.`) unchanged.

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl`, replace line 54 exactly.

Before (line 54):
```
expectZero["KX/KW equivalence", ((1 - x/4) /. x -> xFloor) - Pi^2/(4 zetaReq)];
```

After (line 54):
```
expectZero["KX/KW equivalence", ((1/aK) /. y -> 0 /. x -> xFloor) - Pi^2/(4 zetaReq)];
```

Do not modify any other line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 055` and `redteam exec-mathematica 055`. Both must exit 0 and the captured outputs must still show `KX/KW equivalence = 0` (and `PASS: KX/KW equivalence` for Mathematica). The string `(1 / AK).subs(y, 0).subs(x, x_floor)` must appear at line 56 of the SymPy file, and `((1/aK) /. y -> 0 /. x -> xFloor)` must appear at line 54 of the Mathematica file.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl`
- summary: Re-anchored the KX/KW equivalence checks to the reciprocal A_K expression at y = 0.
- deviation: Mathematica's matching assertion was line 60 in the working tree due to pre-existing edits; only that matching assertion line was changed.
