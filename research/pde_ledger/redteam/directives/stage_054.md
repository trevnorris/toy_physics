---
unit_id: 054
batch: III.2
created_at: 2026-05-22T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-22T17:37:53-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 054

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:34`

**Issue:** The Mathematica script hardcodes the Neumann-derived coefficient as `bExpr = FullSimplify[a Tan[k ell], Assumptions -> $Assumptions];` at line 34, without ever solving the bottom Neumann condition. The subsequent assertion at line 40, `expectZero["Robin equation -> k tan(kL) - h", charEq/a + h - k Tan[k ell]];`, is therefore tautological by construction. The SymPy counterpart at line 33 uses `sp.solve(sp.Eq(sp.diff(psi, s).subs(s, L), 0), B)[0]` and is genuinely derivational; Mathematica must do the same.

**Required change:**
At line 34, replace the hardcoded form with an explicit `Solve` of the Neumann condition. Concretely:

Before (line 34):
```
bExpr = FullSimplify[a Tan[k ell], Assumptions -> $Assumptions];
```

After (line 34):
```
bExpr = FullSimplify[b /. First@Solve[(D[psi, s] /. s -> ell) == 0, b], Assumptions -> $Assumptions];
```

Do not change any other lines. The print on line 38 (`"B from Neumann bottom = ", fmt[bExpr]`) should continue to output `a*Tan[ell*k]` because the Solve yields exactly that expression after `FullSimplify`.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 054` and confirm the assertion `Robin equation -> k tan(kL) - h` still passes (output line 16: `PASS: Robin equation -> k tan(kL) - h`), the printed `B from Neumann bottom = a*Tan[ell*k]` is unchanged (output line 13), and the script exits 0. The verifier additionally inspects line 34 of the source to confirm it now contains `Solve[(D[psi, s] /. s -> ell) == 0, b]` (or an equivalent Solve of the same Neumann condition for `b`) rather than the hardcoded `a Tan[k ell]`.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl`
- summary: Replaced the hardcoded Neumann coefficient with an explicit Solve of the bottom Neumann condition for b.
- deviation: none

## F2 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl:78`

**Issue:** Line 78 defines `xFloor = FullSimplify[4 - 4/zetaReq, Assumptions -> zetaReq > 0];` — i.e. it writes in the claimed solution to `A_K,max = zeta_req` by hand. The assertion at line 83, `expectZero["x floor = 4 - 4/zeta_req", xFloor - (4 - 4/zetaReq)];`, reduces to `(4 - 4/zetaReq) - (4 - 4/zetaReq) == 0` and is tautological. The SymPy counterpart at line 83 uses `sp.solve(sp.Eq(AK_max, zeta_req), x)[0]` to actually invert `A_K,max = zeta_req`; Mathematica must do the same.

**Required change:**
At line 78, replace the hardcoded form with an explicit `Solve` that inverts `aKMax == zetaReq` for `x`. Concretely:

Before (line 78):
```
xFloor = FullSimplify[4 - 4/zetaReq, Assumptions -> zetaReq > 0];
```

After (line 78):
```
xFloor = FullSimplify[x /. First@Solve[aKMax == zetaReq, x], Assumptions -> zetaReq > 0];
```

Do not change any other lines. The print at line 82 (`"x floor at saturation = ", fmt[xFloor]`) should continue to output `4 - 4/zetaReq` because Solve returns exactly that. The downstream assertion at line 84 (`A_K,max(x_floor) - zeta_req`) continues to pass because the substitution `aKMax /. x -> xFloor` evaluates to `zetaReq` either way.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 054` and confirm assertions `x floor = 4 - 4/zeta_req` (output line 39) and `A_K,max(x_floor) - zeta_req` (output line 41) both still pass, and the script exits 0. The verifier additionally inspects line 78 of the source to confirm it now contains `Solve[aKMax == zetaReq, x]` (or an equivalent inversion of `aKMax == zetaReq` for `x`) rather than the hardcoded `4 - 4/zetaReq`.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl`
- summary: Replaced the hardcoded saturation floor with an explicit Solve of aKMax == zetaReq for x.
- deviation: none
