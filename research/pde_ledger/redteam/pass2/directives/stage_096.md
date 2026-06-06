---
unit_id: 096
batch: IV.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T23:17:05-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 096

Apply the finding below. After applying, append an `## Applied: F1` block under the finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If the required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F1` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py:88-116`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl:57-80`

**Issue:**
SECTION II hardcodes `eps_2 = sp.Integer(0)`, `eps_4 = sp.Integer(0)` (`.wl`: `eps2 = 0; eps4 = 0`) and evaluates the obstruction formula `c_pole = (1+eps_4)/(4*(1+eps_2)^2)` only at that degenerate point, where it collapses to the literal `1/4`. Consequently the five SECTION II assertions (`c_pole - 1/4`, `c_geom - 3/4`, `rho_alpha - 4/3`, `zeta_req - 1/3`, and the constant coefficients of `Yhat_cons - Yhat_expected`) are arithmetic-by-construction and cannot fail for any structurally wrong obstruction formula (wrong power on `(1+eps_2)`, wrong factor of 4, eps_2/eps_4 swapped). This is a checkpoint (higher bar): the obstruction formula's structure — the load-bearing object — must be exercised, not just evaluated at the point where it degenerates. The card's stated deliverable values (1/4, 3/4, 4/3, 1/3, Yhat closed form) are all correct and must remain unchanged; this fix only makes the eps→0 reduction load-bearing rather than tautological.

**Required change (SymPy, around lines 88–116):**

1. Replace the literal hardcoding with free symbols PLUS a general obstruction formula, keeping the static-limit evaluation:
   - Declare `eps_2, eps_4 = sp.symbols("eps_2 eps_4", real=True)`.
   - Define the general formula `c_pole_gen = (1 + eps_4) / (4 * (1 + eps_2)**2)`.
2. Add CAN-FAIL eps-sensitivity checks that pin the formula's structure to concrete distinct literals (these fail if a power/factor/symbol is wrong):
   - `expect_zero("c_pole|eps_4=1,eps_2=0 - 1/2", c_pole_gen.subs({eps_2: 0, eps_4: 1}) - sp.Rational(1, 2))`  — proves eps_4 enters linearly in the numerator.
   - `expect_zero("c_pole|eps_2=1,eps_4=0 - 1/16", c_pole_gen.subs({eps_2: 1, eps_4: 0}) - sp.Rational(1, 16))`  — proves eps_2 enters as a squared factor of 4 in the denominator.
   (Both target literals differ from 1/4, so a swapped/wrong-power formula fails them.)
3. Take the static limit FROM the general formula and keep the existing deliverable checks:
   - `eps_2_val = sp.Integer(0)`, `eps_4_val = sp.Integer(0)` (keep printing eps_2/eps_4 as before).
   - `c_pole = sp.simplify(c_pole_gen.subs({eps_2: eps_2_val, eps_4: eps_4_val}))`.
   - Leave `c_geom`, `rho_alpha`, `zeta_req`, `Yhat_cons`, `Yhat_expected` definitions (lines 94–98, 100–101) and the five existing `expect_zero` deliverable assertions (lines 109–116) UNCHANGED — they now consume the limit of the general formula, so they still print/assert 1/4, 3/4, 4/3, 1/3 and the Yhat closed form exactly.

**Required change (Mathematica, around lines 57–80):** mirror symmetrically:
   - `cPoleGen = (1 + eps4)/(4*(1 + eps2)^2)` with `eps2, eps4` left symbolic.
   - `expectZero["c_pole|eps4=1,eps2=0 - 1/2", (cPoleGen /. {eps2 -> 0, eps4 -> 1}) - 1/2];`
   - `expectZero["c_pole|eps2=1,eps4=0 - 1/16", (cPoleGen /. {eps2 -> 1, eps4 -> 0}) - 1/16];`
   - `cPole = FullSimplify[cPoleGen /. {eps2 -> 0, eps4 -> 0}];` then keep `cGeom`, `rhoAlpha`, `zetaReq`, `yhatCons`, `yhatExpected` (lines 60–64) and the five existing `expectZero` deliverable checks (lines 76–80) UNCHANGED.

Do NOT change SECTION I (orthogonality / Laplace eigenvalue) in either engine — it is already substantive. Do NOT change any deliverable value or the docstring's claimed values.

**Self-test (already performed by auditor):**
`(1+1)/(4·(1+0)^2) = 2/4 = 1/2` and `(1+0)/(4·(1+1)^2) = 1/16`; both differ from `1/4`, so the new checks are can-fail and non-tautological. At eps_2=eps_4=0 the formula gives `1/4`, so all existing downstream deliverable values are preserved. No derivatives are introduced, so no identically-zero-derivative trap.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 096` and `redteam exec-mathematica 096` and confirms: (a) the two new eps-sensitivity check lines appear in each output reporting `= 0` / `PASS` for the `1/2` and `1/16` literals; (b) the existing deliverable lines (`c_pole - 1/4 = 0`, `c_geom - 3/4 = 0`, `rho_alpha - 4/3 = 0`, `zeta_req - 1/3 = 0`, `Yhat_Q^cons - [...] = 0`) still report 0/PASS; (c) both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage096_geometry_lane_check_verdict_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage096_geometry_lane_check_verdict_mathematica_audit.wl`
- summary: Replaced the degenerate static obstruction setup with symbolic general formulas, added eps-sensitivity checks, and then evaluated the unchanged deliverables on the static limit.
- deviation: none
