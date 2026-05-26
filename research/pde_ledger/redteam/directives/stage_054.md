---
unit_id: 054
batch: III.2
created_at: 2026-05-26T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-26T02:59:06-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 054

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py` (insert new check between the current line 70 and line 72)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl` (insert new check between the current line 79 and line 81)

**Issue:** The paper's window claim `1 <= A_K <= 4/(4-x)` (eq. `app-stage054-AK-window` in `paper/stages/stage_054.tex` lines 36-40) is currently exercised in both scripts only at the two endpoints `y = pi/2` (DN limit, sympy line 69 / mma line 78) and `y -> 0+` (soft-mouth limit, sympy line 70 / mma line 79). Two endpoint values alone do not establish a window for a continuous function. The notes (`notes/stages/moving_throat_pde_stage054_robin_softening_support_lane.md` line 148) state explicitly: "the map `y -> eta = y tan y` is strictly increasing on (0,pi/2), and `A_K` is strictly decreasing in y" — i.e. monotonicity is the load-bearing piece that closes the window. Neither engine verifies this. The fix is to add one non-tautological symbolic identity for `dA_K/dy` in each engine. Combined with positivity of the prefactor on the declared domain `0 < x < 4, 0 < y < pi/2`, this certifies the window.

**Required change:**

**SymPy (`scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py`)**:

Between the current line 70 (`expect_zero("soft-mouth limit", AK_soft - 4 / (4 - x))`) and line 72 (`banner("PURE SOFTENING THRESHOLD")`), insert the following block (preserve the blank line above and below):

```python
# Monotonicity certificate: A_K is strictly decreasing in y on (0, pi/2)
# with 0 < x < 4. Verify the closed form of dA_K/dy; positivity of the
# prefactor 2*x*y/pi^2 on the assumed domain then implies dA_K/dy < 0,
# which together with the endpoint values brackets A_K in [1, 4/(4-x)].
dAK_dy = sp.diff(AK_sym, y)
dAK_dy_expected = -2 * x * y / (pi**2 * (1 - x / 4 + x * y**2 / pi**2) ** 2)
expect_zero("dA_K/dy closed form", dAK_dy - dAK_dy_expected)
print("Prefactor 2*x*y/pi^2 > 0 on 0<x<4, 0<y<pi/2 => dA_K/dy < 0 (monotone decreasing).")
```

Do not modify any existing assertion, print, or banner.

**Mathematica (`mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl`)**:

Between the current line 79 (`expectZero["soft-mouth limit", aKSoft - 4/(4 - x)];`) and line 81 (`ineqRhs = FullSimplify[1/zetaReq - 1 + x/4, ...];`), insert the following block (preserve a blank line above and below the insertion):

```mathematica
(* Monotonicity certificate: A_K is strictly decreasing in y on (0, Pi/2)
   with 0 < x < 4.  Verify the closed form of D[aKSym, y]; positivity of
   the prefactor 2 x y / Pi^2 on the declared domain then implies
   D[aKSym, y] < 0, closing the window 1 <= A_K <= 4/(4-x). *)
dAKdy = FullSimplify[D[aKSym, y], Assumptions -> $Assumptions];
dAKdyExpected = -2 x y / (Pi^2 (1 - x/4 + x y^2/Pi^2)^2);
expectZero["dA_K/dy closed form", dAKdy - dAKdyExpected];
Print["Prefactor 2*x*y/Pi^2 > 0 on 0<x<4, 0<y<Pi/2 => D[aKSym, y] < 0 (monotone decreasing)."];
```

The existing `$Assumptions` block at lines 53-55 already declares `0 < x < 4` and `y > 0`, which is sufficient for the FullSimplify call. Do not modify any other line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 054` and `redteam exec-mathematica 054`. Both scripts must exit 0. The SymPy transcript must contain a new line `dA_K/dy closed form = 0` (in the existing `EXACT SOFTENING FACTOR` section, just after the `soft-mouth limit = 0` line). The Mathematica transcript must contain new lines `dA_K/dy closed form = 0` and `PASS: dA_K/dy closed form` (just after `PASS: soft-mouth limit`). Existing PASS lines and outputs must be unchanged.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage054_robin_softening_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage054_robin_softening_mathematica_audit.wl`
- summary: Added the requested symbolic derivative checks certifying that A_K is monotone decreasing in y on the declared domain.
- deviation: none
