---
unit_id: 044
batch: III.1
created_at: 2026-05-22T00:00:00-06:00
findings_count: 5
stop_cold: null
applied: true
applied_at: 2026-05-22T12:43:45-06:00
findings_applied: 5
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 044

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:42-61`

**Issue:** The `.wl` script's quadratic-branch and `xi_phys` construction (WL lines 42-61) is a line-by-line port of the SymPy construction at `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:54-79`. Both engines manually write the same `B_cont`, `C_cont` and form `xi_phys = (-B + sqrt(disc))/2`. Mathematica must provide an independent algebraic route — concretely, it should derive `xi_phys` via `Solve` on `branchEq` and verify that the `Solve`-derived root agrees with the manually-written `(-bCont + Sqrt[deltaDisc])/2`. This makes the two engines do algebraically different things, so a typo in the manually-written coefficients would be caught.

**Required change:**

In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`, after the existing block ending at line 61, insert the following lines (do not delete the existing block — augment it):

```
(* Independent route: solve branchEq for xi via Mathematica's algebraic solver,
   then pick the root that vanishes at zero load and verify it matches xiPhys. *)
xiRoots = xi /. Solve[branchEq == 0, xi];
xiPhysSolve = SelectFirst[
  xiRoots,
  TrueQ[FullSimplify[(# /. {mMix -> 0, mSupp -> 0}) === 0,
                     Assumptions -> $Assumptions]] &
];
Print["xi_phys (from Solve) = ", fmt[xiPhysSolve]];
expectZero["xiPhys solve match", xiPhysSolve - xiPhys];
```

Place these lines immediately after the existing `expectZero["zero-load root", ...]` at WL line 61, before the `subbanner["2. Exact continuum-selected normalization function"];` at WL line 63.

If `SelectFirst` returns `Missing[...]` (no root vanishes at zero load with the default branch convention), fall back to:
```
xiPhysSolve = First[Select[xiRoots,
  TrueQ[Simplify[(# /. {mMix -> 0, mSupp -> 0}) == 0,
                 Assumptions -> $Assumptions]] &]];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-mathematica 044` and confirm the new `xi_phys (from Solve)` print and `PASS: xiPhys solve match` line appear in the saved output, and the script exits 0.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- summary: Added the independent Mathematica `Solve` route for `xiPhys`, including the directive fallback for the zero-load branch selection.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:139-144`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:115-120`

**Issue:** The "quadratic mismatch penalty" assertion (sympy line 144, WL line 120) is a pure variable rename — substituting `Rphi -> RU - mismatch` into `(RU - Rphi)^2` is algebraically guaranteed to yield `mismatch^2`. It anchors no physics. Replace with a check that extracts the `xi`-constant coefficient of `branch_eq` and compares to the claimed form `-delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*(RU-Rphi)^2`.

**Required change:**

(a) In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`, replace lines 139-144 with:

```python
subbanner("27.5 — Exact mismatch penalty")

# Extract the xi-constant coefficient of branch_eq (the numerator of n_req - M_supp).
# branch_eq = 9*(xi^2 + B_cont*xi + C_cont), so the constant-in-xi term is 9*C_cont.
C_from_branch_eq = sp.simplify(sp.Poly(branch_eq, xi).nth(0) / 9)
C_expected = sp.simplify(-delta * (Mmix + Msupp) + lambda0 * Mmix * Msupp * (RU - Rphi)**2)
expect_zero("mismatch penalty in C coefficient", C_from_branch_eq - C_expected)
```

Remove the `mismatch = sp.symbols("Delta_R", real=True)` line entirely (it was at the start of the deleted block, line 141).

(b) In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`, replace lines 115-120 with:

```
subbanner["5. Exact mismatch penalty"];

(* Extract the xi-constant coefficient of branchEq (= 9*(xi^2 + bCont*xi + cCont)). *)
cFromBranchEq = FullSimplify[CoefficientList[branchEq, xi][[1]]/9,
                             Assumptions -> $Assumptions];
cExpected = FullSimplify[-delta (mMix + mSupp) + lambda0 mMix mSupp (rU - rPhi)^2,
                         Assumptions -> $Assumptions];

expectZero["mismatch penalty in C coefficient", cFromBranchEq - cExpected];
```

Also remove `deltaR` from the `Clear[...]` at WL line 35 and from the `Element[..., Reals]` list at WL line 37 (since `deltaR` is no longer used anywhere).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 044` and `redteam exec-mathematica 044` and confirm:
- The saved output no longer contains "quadratic mismatch penalty" but contains "mismatch penalty in C coefficient" with residual `= 0` (sympy) and `PASS: mismatch penalty in C coefficient` (WL).
- Both engines exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- summary: Replaced the mismatch rename check with extraction of the `xi`-constant branch coefficient in both engines and removed the unused WL `deltaR` assumption.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:127-129`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:103, 112`

**Issue:** The "tracking total-loading law" assertion is algebraically guaranteed by the preceding "tracking collapse of n_req" assertion. In sympy, `(n_track - Msupp) - (G_q - (Mmix + Msupp))` expands to `n_track - G_q + Mmix`, which is exactly what A5 (`n_track - (G_q - Mmix)`) already asserts. Remove the redundant assertion in both engines.

**Required change:**

(a) In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`, delete lines 127-129 (the comment `# Setting actual support baseline ...`, the `tracking_equation = ...` assignment, and the `expect_zero("tracking total-loading law", tracking_equation)` call). Leave line 125 (`expect_zero("tracking collapse of n_req", ...)`) and line 131 onwards (`F_track = ...`) intact.

After deletion, the block from sympy line 121 should read:
```python
subbanner("27.4 — Interference-matched tracking surface")

n_track = sp.simplify(n_req.subs(Rphi, RU))
G_q = sp.simplify(xi * (delta + xi) / (delta + (1 + lambda0 * RU**2) * xi))
expect_zero("tracking collapse of n_req", n_track - (G_q - Mmix))

F_track = sp.simplify(F_cont.subs(Rphi, RU))
```

(b) In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`:
- Delete line 103 (the `trackingEquation = ...` assignment).
- Delete line 112 (the `expectZero["tracking total-loading law", trackingEquation];` call).
- Leave lines 101 (`nTrack = ...`), 102 (`gQ = ...`), 104 (`fTrack = ...`), 105-109 (`fTrackExpected = ...`), 111 (`expectZero["tracking collapse of n_req", ...]`), and 113 (`expectZero["tracking F collapse", ...]`) intact.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 044` and `redteam exec-mathematica 044` and confirm the line `tracking total-loading law = 0` (and `PASS: tracking total-loading law`) no longer appears in the saved outputs, while `tracking collapse of n_req = 0` and `tracking F collapse = 0` (with their PASS lines) still appear. Both engines exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- summary: Removed the redundant tracking total-loading assertion from both audit scripts while leaving the collapse checks intact.
- deviation: none

## F4 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:81-96`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:63-76`

**Issue:** `F_cont` is exercised only at two slices (`Rphi=1` and `Rphi=RU`). Add a third independent literal slice (`Rphi=2`) with an independently-transcribed expected form so a coefficient typo in `F_cont`'s general definition would surface.

**Required change:**

(a) In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`, insert the following block immediately after line 96 (`print("F_cont built successfully.")`) and before line 98 (`subbanner("27.3 — Minimal-kernel source-tied surface")`):

```python
# Third slice: independent literal R_phi = 2 to constrain bivariate dependence
# beyond the R_phi=1 and R_phi=R_U collapsed surfaces.
Rphi_lit = sp.Integer(2)
F_lit = sp.simplify(F_cont.subs(Rphi, Rphi_lit))
F_lit_expected = sp.simplify(
    (delta + (1 + lambda0 * RU * Rphi_lit) * xi)**2
    * (delta + (1 + lambda0 * Rphi_lit) * xi - Mmix * lambda0 * (RU - Rphi_lit) * (RU - 1))**2
    / ((1 - xi) * (
        (delta + xi - Mmix * lambda0 * RU * (RU - Rphi_lit))**2
        + lambda0 * (Mmix * (RU - Rphi_lit) + Rphi_lit * xi)**2
    )**2)
)
expect_zero("third-slice F at Rphi=2", F_lit - F_lit_expected)
```

(b) In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`, insert the following block immediately after the `Print["F_cont(xi) = ", fmt[fCont]];` at line 76 and before the `subbanner["3. Minimal-kernel source-tied surface"];` at line 78:

```
(* Third slice: independent literal rPhi = 2 to constrain bivariate dependence. *)
fLit = FullSimplify[fCont /. rPhi -> 2, Assumptions -> $Assumptions];
fLitExpected = FullSimplify[
  (delta + (1 + lambda0 rU 2) xi)^2
    (delta + (1 + lambda0 2) xi - mMix lambda0 (rU - 2) (rU - 1))^2/
    ((1 - xi) ((delta + xi - mMix lambda0 rU (rU - 2))^2
               + lambda0 (mMix (rU - 2) + 2 xi)^2)^2),
  Assumptions -> $Assumptions
];
expectZero["third-slice F at rPhi=2", fLit - fLitExpected];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 044` and `redteam exec-mathematica 044` and confirm:
- The saved sympy output contains `third-slice F at Rphi=2 = 0`.
- The saved Mathematica output contains `PASS: third-slice F at rPhi=2`.
- Both engines exit 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- summary: Added the independent literal `R_phi = 2` / `rPhi = 2` normalization slice checks to both engines.
- deviation: none

## F5 — symbol_assumption_error

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:50, 13-14`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:35, 37`

**Issue:** `sigma0` / `sigma_0` is declared in both scripts but never used anywhere. The docstring mentions a `sigma0=0` limit that is never tested (the actual check uses `Rphi=1`). Remove the unused symbol declaration and update the docstring to match what is actually verified.

**Required change:**

(a) In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`:
- Delete line 50 entirely: `sigma0 = sp.symbols("sigma_0", real=True)`.
- Edit line 13 from `4. The minimal-kernel limit sigma0=0 gives the exact source-tied closure.` to `4. The minimal-kernel limit R_phi=1 gives the exact source-tied closure.` (replace `sigma0=0` with `R_phi=1`).

(b) In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`:
- Edit line 35 from `Clear[xi, delta, mMix, mSupp, rU, rPhi, sigma0, deltaR];` to `Clear[xi, delta, mMix, mSupp, rU, rPhi];` (remove both `sigma0` and `deltaR`; the latter is removed under F2).
  - If F2 has not yet been applied at the time this edit is made, leave `deltaR` and only remove `sigma0`: `Clear[xi, delta, mMix, mSupp, rU, rPhi, deltaR];`.
- Edit line 37 from `Element[{xi, delta, mMix, mSupp, rU, rPhi, sigma0, deltaR}, Reals] &&` to `Element[{xi, delta, mMix, mSupp, rU, rPhi}, Reals] &&` (remove `sigma0` and `deltaR`).
  - Same conditional fallback as above if F2 not yet applied: leave `deltaR` if necessary.

**Verification command:**
After Codex applies, the verifier confirms no occurrence of `sigma0` or `sigma_0` remains in either script. Both engines still exit 0 since the symbol was unused.

## Applied: F5

- files_changed:
  - `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- summary: Removed the unused `sigma0` / `sigma_0` declarations and corrected the SymPy docstring to the verified `R_phi=1` limit.
- deviation: none
