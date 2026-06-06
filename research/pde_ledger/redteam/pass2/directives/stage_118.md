---
unit_id: 118
batch: IV.3
created_at: 2026-06-06T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-06T17:08:16Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 118

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py:82`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl:90`

**Issue:** `K_q` is defined directly as its target closed form and then the "K_q closed form" check subtracts the identical literal — an `X − X` tautology that cannot fail. Unlike the sibling checks (`g_q`, `J_s`, `g_s`, `I_q`), it does not exercise any derived quantity. The load-bearing reduction (`∫(χ')²dz = π²/(4L_W²)`, already computed as `chi_grad`/`chiGrad` and asserted by the "D/N stiffness check") is what makes `K_q` non-arbitrary, but the K_q check ignores it. Tie `K_q` to that verified integral so the closed-form check becomes load-bearing. This must NOT change the printed `K_q` value or any downstream value.

**Required change:**

SymPy, line 82. Before:
```
K_q = sp.simplify((Zq/mu0) * (sp.pi**2*c_s**2/(4*L_W**2)))
```
After:
```
K_q = sp.simplify((Zq/mu0) * c_s**2 * chi_grad)
```
(`chi_grad` is the gradient integral computed at line 42 and asserted equal to `π²/(4*L_W²)` at line 50. Leave line 97's assertion `expect_zero("K_q closed form", K_q - (Zq/mu0) * sp.pi**2 * c_s**2 / (4*L_W**2))` UNCHANGED — it now non-tautologically confirms `chi_grad` reduced correctly.)

Mathematica, line 90. Before:
```
kQ = FullSimplify[(zQ/mu0)*(Pi^2*cSound^2/(4*lW^2)), Assumptions -> $Assumptions];
```
After:
```
kQ = FullSimplify[(zQ/mu0)*cSound^2*chiGrad, Assumptions -> $Assumptions];
```
(`chiGrad` is computed at line 45 and asserted equal to `Pi^2/(4*lW^2)` at line 54. Leave line 105's `expectZero["K_q closed form", kQ - (zQ/mu0)*Pi^2*cSound^2/(4*lW^2)]` UNCHANGED.)

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 118` and `redteam exec-mathematica 118` and confirms:
- both scripts exit 0 with all checks PASS;
- the printed `K_q` value is unchanged (`pi**2*Zq*c_s**2/(4*L_W**2*mu0)` in the SymPy transcript, `(cSound^2*Pi^2*zQ)/(4*lW^2*mu0)` in the Mathematica transcript);
- the "K_q closed form" residual still prints `0` / `PASS` — but now `K_q` is built from `chi_grad`/`chiGrad`, so the check depends on the gradient integral having reduced to `π²/(4L_W²)`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage118_parent_core_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage118_parent_core_mathematica_audit.wl`
- summary: Rebuilt `K_q`/`kQ` from the previously computed D/N gradient integral while preserving the existing closed-form checks.
- deviation: none
