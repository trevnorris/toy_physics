---
unit_id: 228
batch: VII.1
created_at: 2026-06-02T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-06-02T17:21:30-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 228

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

The `paper_misalignment` finding(s) (F2 & F3) have been RESOLVED by the user (2026-06-02) — see the `## RESOLVED` block at the END of this directive and apply the authorized notes-only edit there as part of this fix loop. (Codex applies notes/*.md; Claude reviews.)

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing.

Do NOT touch paper.tex or the appendix. The ONLY authorized prose/notes edit is the notes-only change specified in the `## RESOLVED` block at the END of this directive (user-authorized 2026-06-02); apply exactly that and make no other notes/prose edits.

## F1 — missing_verification_script (subtype: missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_mathematica_audit.wl`

**Issue:** No Mathematica audit exists for stage 228. The unit is not status-only and carries `ExactClosure` symbolic claims a second engine can independently verify. Create an independent `.wl` that re-derives every claim from the physical premises using NATIVE Mathematica primitives and a DIFFERENT decomposition than the SymPy script — this must NOT be a transliteration of the `.py`.

**Anti-transliteration requirements (mandatory):**
- Do NOT replicate the SymPy `np.roots`-on-a-numeric-companion-matrix route for the pole census. Use exact `Solve`/`Roots` (or `NSolve`/`Eigenvalues` with a genuinely different formulation) on the quartic `F(y)` in `y = omega^2`.
- Do NOT echo the SymPy variable choreography (`Kd = K*exp(eps*xK)`, etc. then `sp.diff(...,eps).subs(eps,0)`). Prefer `Series[..., {eps,0,1}]` + `Coefficient[..., eps, 1]` (or `D[...,eps]/.eps->0`) and `Coefficient`/`CoefficientArrays` for the row extraction, with `MatrixRank`/`NullSpace` for the counts and `Det`/`Factor` for the determinant.
- Use `FullSimplify` for symbolic-zero checks (strip any `ConditionalExpression[0,...]`), and a tolerance helper (`Abs[a-b] < tol`) for the numeric checks.

**Sample branch (identical physical premises — these are inputs, not derived here):**
`kappa = 2 Sqrt[2]/Pi`, `lamB = 1/2`, `lamU = 3/10`, `lamW = 2/5`, `lamR = 1/4`, `OmU = 1`, `OmW = 7/5`, `varpi = 2`, `M = 1`, `a = 1`, `c_s = 1`. Derive `K_compat` in-script as the SymPy does (`K_compat = B0 + Z0 + 3(1+B2+Z2)^2/(B4+Z4)` evaluated on the sample), do NOT hardcode its numeric value. Mixed slope basis order: `{xLU, xLW, xLR, xOU, xOW}`.

The carried upstream-static budgets are inputs (used as-is, NOT derived here): `threshold = 21.854566296358396`, `budget_both = 0.367930328492646`, `budget_nonempty = 0.737619063660757`.

**Claim manifest** (the new `.wl` must independently verify each):

- **M1** (split identity): `Xi_1 - 2 (pi_1 - delta_1) == 0` symbolically on the pure-transfer corridor, where `pi_1 = P_01/P`, `delta_1 = Delta_01/Delta`, `Xi_1 = N_01/N_0`. (`expectZero` via `FullSimplify`.)
- **M2** (exact `pi_1` row): coefficients of `pi_1` in `{xLU,xLW,xLR,xOU,xOW}` equal `{3/19, 16/19, 3/19, 32/19, 0}`.
- **M3** (exact `delta_1` row): coefficients of `delta_1` in `{xLU,xLW,xLR,xOU,xOW}` equal `{0, 0, 50/(25-98 Pi^2), 196 Pi^2/(98 Pi^2-25), 196 Pi^2/(98 Pi^2-25)}`. (NOTE: the `196` here is the SymPy-verified value; the notes' `247` is under user review per F2 — if the user later changes the script, this manifest entry must follow. Until then, use `196`.)
- **M4** (rank/nullity): on the `D01,D21,D41` mixed-coefficient matrix, `MatrixRank = 3`, nullity `= 2`; adding the `pi_1` row gives rank `4`, nullity `1`; adding the `delta_1` row gives rank `4`, nullity `1`; adding both gives rank `5`, nullity `0`.
- **M5** (reduced determinant != 0): on a basis of the pure-transfer nullspace `T`, `Det[{pi_coeff.T.T, delta_coeff.T.T}]` equals `196 (200 + 147 Pi^2)(80000 + 343225 Pi^2 + 43218 Pi^4)/(475 (8670000 + 14894275 Pi^2 + 2117682 Pi^4))` (the SymPy-verified form; notes' `247(251+215 Pi^2)...` is under user review per F3 — use `196(200+147 Pi^2)...` until resolved) and is nonzero.
- **M6** (unit generators): the oriented (`Xi_1 > 0`) Euclidean-unit nullspace generators satisfy `v_num ~= {-0.55551149, 0.31814576, -0.65766801, -0.04533730, -0.39447126}` and `v_den ~= {-0.26583993, 0.18448137, 0.94454459, 0.04984499, -0.02543112}` (tol 1e-8), with `pi_1(v_num)=0`, `delta_1(v_num) ~= -0.86805617`, `Xi_1(v_num) ~= 1.73611235`; `delta_1(v_den)=0`, `pi_1(v_den) ~= 0.34646608`, `Xi_1(v_den) ~= 0.69293215`; and `Xi_1(v_num) = -2 delta_1(v_num)`, `Xi_1(v_den) = 2 pi_1(v_den)`.
- **M7** (pole census): on the undeformed sample point, the two wall-like positive-omega roots of `F(y)` give `omega_- ~= 1.99753568`, `omega_+ ~= 4.94905432`, `R_Q,- ~= 30.19990756`, `R_Q,+ ~= 36.17118648`, and `P_0 ~= 0.00206979232` (tol ~2e-10 / 2e-15). `R_Q(omega) = 27 c_s^5/(a^5 omega^5 N(omega))` with `N(omega) = ((OmU^2-omega^2) GW + R GU)^2 / (((OmU^2-omega^2)(OmW^2-omega^2)-R^2)^2)`; `P_0 = N_0/(K - C^2/varpi^2 - Z_0)`.
- **M8** (first-order log-slopes): along `v_num`, `dln P_0 ~= 1.73611235`, `dln R_Q,+ ~= -0.52346582`, `dln R_Q,- ~= +0.71358484`, wall-frequency slopes negligible (`< 5e-5`); along `v_den`, `dln P_0 ~= 0.69293215`, `dln R_Q,+ ~= -0.35245541`, `dln R_Q,- ~= -0.23169484`, wall-frequency slopes negligible. Also `dln P_0 == Xi_1` to ~5e-7 on each branch. (Use a symmetric difference or symbolic `D`; if finite-difference, an independent step size is fine — the point is method independence, not matching SymPy's `1e-8`.)
- **M9** (ceilings + decisive inequality): dynamic ceilings `(num_both, num_nonempty) ~= (0.96253269, inf)` and `(den_both, den_nonempty) ~= (1.39592653, 1.42955095)`; transported static ceilings `num_both ~= 0.21192772`, `num_nonempty ~= 0.42486828`, `den_both ~= 0.53097598`, `den_nonempty ~= 1.06448959` (each `= budget / Xi_1`); and the decisive comparisons `num_both_dyn > num_both_stat`, `den_both_dyn > den_both_stat`, `den_nonempty_dyn > den_nonempty_stat` all hold. (`num_nonempty_dyn` is `Infinity` because one wall pole improves — exclude it from the finite-inequality check, as the SymPy script does.)

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 228` and confirms the new `.wl` appears at the Target path, exits 0, and all checks pass — and that it is an independent re-derivation, not a transliteration.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_mathematica_audit.wl`
- summary: Created the independent Mathematica audit for M1-M9 using series coefficient extraction, exact matrix ranks/nullspaces, quartic NSolve roots, and implicit dynamic-slope differentiation.
- deviation: none

## F2 — paper_misalignment (subtype: notes_contradicts_script)

**Paper side (notes):**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md:151` quote: `+\frac{247\pi^2}{98\pi^2-25}x_{\Omega_U}`
- `...:152` quote: `+\frac{247\pi^2}{98\pi^2-25}x_{\Omega_W}.`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage228_..._sympy_audit.py:194` quote: `196 * sp.pi**2 / (98 * sp.pi**2 - 25), 196 * sp.pi**2 / (98 * sp.pi**2 - 25)]`
- committed output `:14` quote: `delta_1= ... + 196*pi**2*xOU/(-25 + 98*pi**2) + 196*pi**2*xOW/(-25 + 98*pi**2)`

### Resolve before fix_loop

The notes state the `delta_1` row coefficient on `x_OmU` and `x_OmW` as `247 Pi^2/(98 Pi^2 - 25)`; the script (and its committed output) use `196 Pi^2/(98 Pi^2 - 25)`. By hand: `Delta_01 = 2 OmU^2 OmW^2 (x_OU + x_OW) - 2 R^2 x_LR`, and on the sample branch the `x_OU` coefficient of `delta_1 = Delta_01/Delta` is `2*(49/25)/((98 Pi^2 - 25)/(50 Pi^2)) = 196 Pi^2/(98 Pi^2 - 25)`. The script value is mathematically correct; the notes value appears to be a typo (the same `247`-vs-`196` slip seen in F3).

Possible directions (the user picks one):
- (a) Script is correct (recommended; matches hand derivation) → change notes lines 151-152 from `247\pi^2` to `196\pi^2`. No script change. This is a paper-side edit Codex applies only after explicit user authorization.
- (b) Notes are correct → re-derive; change the script's `expected_delta` entries and re-run sympy + the new mathematica audit.
- (c) Both trace to a third source that contradicts both → flag for deeper review.

F2 is RESOLVED by the user (see the ## RESOLVED block at the end) — apply the authorized notes edit.

## Applied: F2

- files_changed:
  - `notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md`
- summary: Replaced the two stale `247\pi^2/(98\pi^2-25)` delta-row note coefficients with `196\pi^2/(98\pi^2-25)`.
- deviation: none

## F3 — paper_misalignment (subtype: notes_contradicts_script)

**Paper side (notes):**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage228_..._sympy_audit.md:196` quote: `\frac{247(251+215\pi^2)(80000+343225\pi^2+43218\pi^4)}{475(8670000+14894275\pi^2+2117682\pi^4)}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage228_..._sympy_audit.py:250-252` quote: `196 * (200 + 147 * sp.pi**2) * (80000 + 343225 * sp.pi**2 + 43218 * sp.pi**4) / (475 * (8670000 + 14894275 * sp.pi**2 + 2117682 * sp.pi**4))`
- committed output `:21` quote: `det[...] = 196*(200 + 147*pi**2)*(80000 + 343225*pi**2 + 43218*pi**4)/(475*(8670000 + 14894275*pi**2 + 2117682*pi**4))`

### Resolve before fix_loop

The notes give the reduced determinant's leading factors as `247(251+215 Pi^2)`; the script (machine-derived by SymPy and asserted at line 254) gives `196(200+147 Pi^2)`. These are numerically distinct (`247*251 = 61997` vs `196*200 = 39200`), so it is a genuine mismatch, not a regrouping. This is the determinant companion of the F2 typo (same `247`-vs-`196` slip).

Possible directions (the user picks one):
- (a) Script is correct (recommended; SymPy-derived and consistent with F2's corrected `196`) → change notes line 196 to `\frac{196(200+147\pi^2)(80000+343225\pi^2+43218\pi^4)}{475(8670000+14894275\pi^2+2117682\pi^4)}`. No script change. Paper-side edit Codex applies only after explicit user authorization.
- (b) Notes are correct → re-derive; change the script's `expected_det` and re-run sympy + the new mathematica audit.
- (c) Both trace to a third source that contradicts both → flag for deeper review.

F3 is RESOLVED by the user (see the ## RESOLVED block at the end) — apply the authorized notes edit.

## Applied: F3

- files_changed:
  - `notes/stages/moving_throat_pde_stage228_numerator_denominator_split_of_the_pure_transfer_corridor_and_first_actual_dynamic_window_test_sympy_audit.md`
- summary: Replaced the stale reduced-determinant leading factors `247(251+215\pi^2)` with `196(200+147\pi^2)`.
- deviation: none

## RESOLVED — F2 & F3 paper_misalignment (USER-AUTHORIZED 2026-06-02)

Direction: correct the notes to match the verified SymPy script (authoritative; script:194/250-252, output:14/21). Notes-only; Codex applies, Claude reviews. Do NOT change the script, paper.tex, or appendix.
- F2: the δ_1 row coefficient `247 π²/(98 π² − 25)` (`247\pi^2/(98\pi^2-25)`) → `196 π²/(98 π² − 25)` (`196\pi^2/(98\pi^2-25)`).
- F3: the reduced-determinant leading factors `247 (251 + 215 π²)` → `196 (200 + 147 π²)`.
- Both are notes-only typos in the same `247`-vs-`196` / `215`-vs-`147` family. Acceptance: stale `247`/`215` in those two contexts gone; `196`/`200`/`147` present. Append `## Applied: F2` and `## Applied: F3`.
