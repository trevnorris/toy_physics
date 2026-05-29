---
unit_id: 055
batch: III.2
created_at: 2026-05-26
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-26T03:01:34-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 055

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For `paper_misalignment` findings, do nothing — the orchestrator holds for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" a paper_misalignment unless the user has explicitly chosen a direction in a follow-up directive. (None present in this directive.)

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl:51`

**Issue:**
Line 51 hardcodes `xFloor = FullSimplify[4 - Pi^2/zetaReq, Assumptions -> zetaReq > 0];`. Line 58's assertion `expectZero["x floor = 4 - Pi^2/zeta_req", xFloor - (4 - Pi^2/zetaReq)];` then reduces to `(4 - Pi^2/zetaReq) - (4 - Pi^2/zetaReq) == 0`, which is algebraically guaranteed by construction and cannot fail. By contrast the SymPy script (line 50) derives `x_floor` via `sp.solve(sp.Eq(zeta_max, zeta_req), x)[0]`. Re-anchor the Mathematica derivation to `Solve` so line 58 becomes a genuine cross-check between the Solve output and the paper-stated closed form.

Note: Line 59 (`(zetaMax /. x -> xFloor) - zetaReq`) is already a substantive cross-check and must NOT be modified — it remains the primary verification of D3. Only line 51 changes; line 58 is unchanged (it becomes non-tautological automatically once xFloor is Solve-derived).

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl`, replace line 51 exactly.

Before (line 51):
```
xFloor = FullSimplify[4 - Pi^2/zetaReq, Assumptions -> zetaReq > 0];
```

After (line 51):
```
xFloor = FullSimplify[x /. First[Solve[zetaMax == zetaReq, x]], Assumptions -> $Assumptions];
```

Do not modify any other line. In particular, lines 55 (`Print["x floor from zeta_max = zeta_req = ", fmt[xFloor]];`), 58 (`expectZero["x floor = 4 - Pi^2/zeta_req", xFloor - (4 - Pi^2/zetaReq)];`), 59 (`expectZero["zeta_max(x_floor) - zeta_req", (zetaMax /. x -> xFloor) - zetaReq];`), and 60 (KX/KW equivalence) must remain unchanged.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 055`. The script must exit 0 and the captured output must still show:
- `x floor from zeta_max = zeta_req = 4 - Pi^2/zetaReq` (Solve must yield this form),
- `x floor = 4 - Pi^2/zeta_req = 0` and `PASS: x floor = 4 - Pi^2/zeta_req`,
- `zeta_max(x_floor) - zeta_req = 0` and `PASS: zeta_max(x_floor) - zeta_req`,
- `KX/KW equivalence = 0` and `PASS: KX/KW equivalence`.

The string `Solve[zetaMax == zetaReq, x]` must appear at line 51 of the Mathematica file.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl`
- summary: Replaced the hardcoded x-floor expression with a Solve-derived expression using the shared assumptions.
- deviation: none

## F2 — symbol_assumption_error

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:27`

**Issue:**
Line 27 declares `alpha, x, y, zeta_req = sp.symbols("alpha x y zeta_req", positive=True, real=True)`. The paper requires `0 < x < 4` (closure max `pi^2/(4-x)` is only finite and positive on this open interval) and `0 < y < pi/2` (principal branch of `y tan y = eta`). The Mathematica script correctly declares `0 < x < 4 && y > 0` (line 36). The SymPy declarations silently relax the domain. No current assertion is invalidated, but the script does not document the paper-stated domain at the point of declaration.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`, replace line 27 exactly.

Before (line 27):
```
alpha, x, y, zeta_req = sp.symbols("alpha x y zeta_req", positive=True, real=True)
```

After (replace line 27 with the following 4 lines):
```
# Paper-stated domain: alpha > 0, 0 < x < 4 (Stage 054 softening), 0 < y < pi/2
# (principal branch of y tan y = eta), zeta_req > 0.  SymPy lacks compound
# symbol-level bounds, so positivity is declared here and the (0, 4) / (0, pi/2)
# constraints are exercised by the endpoint substitutions below.
alpha, x, y, zeta_req = sp.symbols("alpha x y zeta_req", positive=True, real=True)
```

Do not modify any other line. The symbol declaration itself stays identical; only a documentary comment block is inserted directly above it. This brings the SymPy script's domain documentation into alignment with the Mathematica script's `$Assumptions` (line 36 of the .wl file) and the paper/notes domain.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 055`. The script must exit 0 and the captured output must be byte-identical to the prior transcript (the change is a comment only; no symbolic algebra is altered). The comment block "`Paper-stated domain: alpha > 0, 0 < x < 4 ...`" must appear immediately above the `sp.symbols(...)` call at line 27 (now line 31 of the file).

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`
- summary: Added the requested paper-domain comment block immediately above the SymPy symbol declaration.
- deviation: none
