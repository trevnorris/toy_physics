---
unit_id: 231
batch: VII.2
created_at: 2026-06-02T00:00:00-06:00
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-06-02T21:52:12-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 231

Apply each finding below in order (F1 — the authorized notes edit ONLY —, F2, F3, F4). After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

F1 is a `paper_misalignment` (notes-vs-script) that the USER HAS RESOLVED — see `## RESOLVED — F1` below. Apply F1 by making ONLY the one authorized notes edit specified there (correct the notes polynomial to match the script). Do NOT change the script's polynomial — it is correct.

If a non-`paper_misalignment` finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes beyond what each finding names. The ONLY authorized prose/notes edit is the F1 notes correction specified in `## RESOLVED — F1`; do not touch `paper.tex` or any other prose document.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md:98` quote: `\frac{(9\delta+11\xi)^3\bigl(81\delta^3+240\delta^2\xi+72\delta^2+297\delta\xi^2+189\xi^3\bigr)}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py:81-84` quote: `(9 * delta + 11 * xi) ** 3 * (81 * delta**3 + 189 * delta**2 * xi + 72 * delta**2 + 297 * delta * xi**2 + 121 * xi**3) / (81 * (1 - xi) ** 2 * (9 * delta**2 + 18 * delta * xi + 11 * xi**2) ** 3)`
- Output corroboration `/var/projects/toy_physics/research/pde_ledger/scripts/output/...sympy_audit.txt:14`: SymPy's own factorization prints `... (81*delta**3 + 189*delta**2*xi + 72*delta**2 + 297*delta*xi**2 + 121*xi**3) ...`

## RESOLVED — F1 (user direction 2026-06-02: correct notes to script)

Direction (a): the SymPy script polynomial is canonical (SymPy independently computed it via `Factor`; saved output line 14 prints the script's coefficients). The notes line 98 is a transcription typo. **Authorized notes edit (Codex APPLIES; Claude reviews):**

In `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md` line 98, inside the `\bigl( ... \bigr)` numerator, correct the two coefficients so the polynomial matches the script / SymPy `Factor`:
- `240\delta^2\xi` → `189\delta^2\xi`
- `189\xi^3` → `121\xi^3`

Resulting polynomial: `81\delta^3+189\delta^2\xi+72\delta^2+297\delta\xi^2+121\xi^3`. No script change. The new F2 `.wl` (M1) independently recomputes `Factor[D[F,x]]` and will cross-engine-corroborate the corrected coefficients.

## Applied: F1

- files_changed:
  - `notes/stages/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.md`
- summary: Corrected the two notes numerator coefficients to `189\delta^2\xi` and `121\xi^3` to match the script and SymPy factorization.
- deviation: none

## F2 — missing_verification_script (subtype missing_mathematica)

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl`

**Issue:** Unit 231 has only a SymPy script; the dual-engine policy requires an independent Mathematica verifier wherever Mathematica can independently verify the claims, which it can here (symbolic differentiation/factoring of rational functions, one-sided limits, threshold root-finding, and a symbolic product-law identity). The new `.wl` must re-derive the results from the physical premises using native Mathematica primitives and a DIFFERENT decomposition than the `.py` — not a line-by-line transliteration. You design the route; the requirements and acceptance criteria are below. The `_mathematica_audit.wl` filename suffix is mandatory and the file lives in `mathematica/`.

**Required change:**
Author the new `.wl` so that it independently establishes the claim manifest below. Use Mathematica-native machinery (e.g. `D`, `Together`, `Factor`, `Apart`, `Simplify`/`FullSimplify`, `Limit[..., Direction -> "FromBelow"]`, `Reduce`/`Resolve` or `FindRoot` for the threshold roots). Each claim must be checked with a hard failure on mismatch (e.g. `If[TrueQ[Simplify[lhs - rhs] == 0], Print["PASS ..."], Print["FAIL ..."]; Exit[1]]`), not a bare `Print`. Do NOT import or echo the `.py`'s intermediate `expected_*` polynomials verbatim as the derivation route; let Mathematica derive the derivative itself and compare to the published closed form.

**Claim manifest** (the new script must independently verify each):

- M1. Geometry/classifier derivatives and signs. With `F = (9d+11x)^4 / (81(1-x)(9d^2+18 d x+11 x^2)^2)`, `G = 9 x (x+d)/(9 d+11 x)`, `RND = 72 d^2 (1-x)/((9 d+11 x)(9 d^2+18 d x+11 x^2))`, verify `D[F,x] > 0`, `D[G,x] > 0`, `D[RND,x] < 0` for `0 <= x < 1, d > 0`. Confirm the published closed forms: `D[G,x] = 9(9 d^2+18 d x+11 x^2)/(9 d+11 x)^2`; and `D[F,x]` numerator factor equals `81 d^3 + 189 d^2 x + 72 d^2 + 297 d x^2 + 121 x^3` (i.e. Mathematica's own `Factor[D[F,x]]` should expose this polynomial — report what Mathematica gets so the F1 notes discrepancy is independently adjudicated).
- M2. Endpoints: `F(0,d) = 1`, `G(0,d) = 0`, `RND(0,d) = 8/(9 d)`, and `Limit[F, x -> 1, Direction -> "FromBelow"] = Infinity`.
- M3. Pullback compiler: `dRphys/dRtarget = D[RND,x] / D[F,x] < 0` on a representative grid of `(x,d)` points (you choose >= 5 distinct points with `0 < x < 1, d > 0`).
- M4. Threshold compiler: `Pc = c(9 d+11 x)(9 d^2+18 d x+11 x^2) - 72 d^2(1-x)`; verify `D[Pc,x] = 3(87 c d^2 + 198 c d x + 121 c x^2 + 24 d^2) > 0`; `Pc(0,d) = 9 d^2(9 c d - 8)`; and the onset law `delta_c = 8/(9 c)` makes `Pc(0, delta_c) = 0`.
- M5. Pulled-back thresholds: with `R* = 0.411024574532864 / 0.334368725711457 ~ 1.229255438463336` and `c = 1`, solve `RND(x_c, d) = c` for `x_c` and report `R = F(x_c, d)` on the three slices `d in {0.25, 0.50, 0.75}`. Acceptance: the `R_flip` (cap `R*`) and `R_den` (cap `1`) values must agree with the SymPy sample rows to ~1e-9 — `d=0.25: R_flip~1.330868539, R_den~1.393832566`; `d=0.50: R_flip~1.139956630, R_den~1.221087062`; `d=0.75: R_flip~1.0, R_den~1.071471867` — and `R_flip <= R_den`. Also verify collapse to onset above the onset deltas: `R_flip(0.80) = 1` (since `0.80 > 8/(9 R*) ~ 0.7231`) and `R_den(1.00) = 1` (since `1.00 > 8/9`).
- M6. Continuum placement product law: with `delta = delta0/(1-eps_eta)`, `M_mix = 8 Z_W (1+rho)^2/(pi^2(1-eps_eta)(1-eps_W))`, `R_target = Lam(1-eps_eta)(1-eps_W)^2/(Z_W(1+rho)^2)`, verify `Simplify[R_target * M_mix - 8 Lam (1-eps_W)/pi^2] == 0`.
- M7. Static-first survives pullback: `0.967282389363822 > 0.367930328492646` and `0.990581810705233 > 0.737619063660757`. (These are carried-forward upstream budgets; assert the strict inequalities. If you can cheaply express the surviving-inequality as a set-inclusion `RND(x,d) in [0, RND(0,d)]` over the M3 grid, add that as corroboration, but it is not required for acceptance.)

**Verification command:**
After Codex applies, the verifier runs `redteam exec-mathematica 231` and confirms each M-claim's PASS line appears AND the script exits 0, and the M5 numeric values match the SymPy sample rows.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl`
- summary: Added an independent Mathematica audit covering M1-M7 with symbolic derivative/sign checks, endpoint checks, sampled pullback monotonicity, threshold roots, product-law verification, and static-first inequalities.
- deviation: none

## F3 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py:180-181`

**Issue:** `s_minus_den = 0.411024574532864` and `s_minus_num = -0.334368725711457` are hardcoded carry-forward literals from stage 230 with no provenance citation. The numbers are not anchored in unit 231's card/notes (which state only `R_* ~ 1.2292554`). Risk is low because `R_*` is independently pinned by `assert_close(..., 1.229255438463336)` at line 186.

**Required change:**
Add a single provenance comment immediately above line 180 naming the specific stage-230 source (script filename and/or result anchor) from which `s_minus_den` and `s_minus_num` are carried. Do NOT change the numeric values. Example shape (Codex fills in the actual stage-230 source identifier — locate it; do not invent):
`    # Carried from Stage 230 selected-branch sign-flip data: s_minus from <stage230 source script / result anchor>`

If the exact stage-230 source cannot be located with confidence, append `## Blocked: F3` with the question rather than citing a guessed source.

**Verification command:**
`redteam exec-sympy 231`; a provenance comment naming the stage-230 source appears above line 180 and the script exits 0 with unchanged numeric output.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py`
- summary: Added a provenance comment tying `s_minus_den` and `s_minus_num` to the Stage 230 selected-branch dynamic-slope block.
- deviation: none

## F4 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py:270-278`

**Issue:** Block 7's static-first inequalities compare four hardcoded literals (`B_dyn_both_inf`, `B_dyn_nonempty_inf`, `B_stat_both`, `B_stat_nonempty`, lines 270-273). The inequality cannot fail unless a literal is retyped, and it does not connect the four budgets to the pullback the stage is supposed to prove. This is the headline "static-first survives pullback" deliverable, so a literal-vs-literal comparison is not load-bearing.

**Required change:**
Two parts, both minimal:
1. Add a provenance comment above line 270 naming the upstream stage (247) source script/result that defines these four budgets. Do NOT change the four values.
2. Connect the inequality to the pullback by adding, after the existing two `if not ...: raise` checks, a loop over the already-defined `probe_grid` asserting that the sampled physical classifier value `R_num(xi_val, delta_val)` lies in the physical subset `[0, R_num(0.0, delta_val)]` (i.e. `0 <= R_phys_sample <= R_ND(0,delta)`), which is the set-inclusion fact the paper's survival argument relies on. Use the existing `R_num` lambdify handle (defined at line 131). Do NOT introduce new symbols or change any existing assertion.

Self-test of part 2: `R_num` depends on both `xi` and `delta`, so the substitutions are non-trivial; `R_ND` is decreasing in `xi` (verified by A5) with `R_ND(0,delta) = 8/(9 delta)`, so for each grid point `0 < xi < 1` the value `R_num(xi,delta)` is strictly between `R_num(1,delta) >= 0` and `R_num(0,delta)`, making the bound a genuine, satisfiable, non-tautological check (it would fail if the wrong classifier were sampled).

**Verification command:**
`redteam exec-sympy 231`; a provenance comment appears above line 270, a new subset-bound assertion over `probe_grid` appears after line 278's block, and the script exits 0.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py`
- summary: Added static-first budget provenance and a `probe_grid` loop asserting sampled `R_num(xi, delta)` stays in `[0, R_num(0, delta)]`.
- deviation: Cited Stage 247's old notes numbering together with the canonical Stage 230 source script, because the current Stage 247 script/output do not define these four budget values.
