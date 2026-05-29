---
unit_id: 029
batch: II.1
created_at: 2026-05-25T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-26T05:38:32Z
findings_applied: 4
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 029

Apply each non-paper_misalignment finding below (F1, F2, F3) in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

For F4 (`paper_misalignment`), do nothing — the orchestrator is holding for user resolution. Do not edit `paper/stages/stage_029.tex`, the notes, or scripts to "fix" F4 unless the user has explicitly chosen a direction in a follow-up directive.

If a non-paper_misalignment finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding and continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch `paper/`, `notes/`, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — paper_misalignment (legacy stage-number drift in docstrings/banner)

**Subtype:** notes_contradicts_script

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:3,5`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:33`

**Note:** Although this is a `paper_misalignment` category, the resolution is entirely script-side (docstring/banner label only). No paper or notes edit is needed. Codex MAY apply this fix without waiting for user resolution — there is no ambiguity about the correct stage number (the file is named `stage_029`, the paper card is `stage_029.tex`, and the Mathematica script already prints "Stage 029 Mathematica audit passed." on its last line).

**Issue:** SymPy docstring header refers to "moving_throat_pde_stage12_dynamic_loading_sympy_audit.py" / "Stage 12"; Mathematica banner prints `STAGE 012 — DYNAMIC LOADING`. The file is in fact stage_029 in the current paper numbering.

**Required change:**
1. In `scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py`:
   - Change line 3 from `moving_throat_pde_stage12_dynamic_loading_sympy_audit.py` to `moving_throat_pde_stage029_dynamic_loading_sympy_audit.py`.
   - Change line 5 from `SymPy audit for Stage 12 of the moving-throat PDE program.` to `SymPy audit for Stage 029 of the moving-throat PDE program.`.
   - Leave the rest of the docstring (lines 7-25) untouched — its references to "Stage-11 loading parameter" describe a physical relationship to an earlier stage in the program and are not file-identification labels. Do NOT mass-replace "Stage 11" or "Stage-11" elsewhere.
2. In `mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl`:
   - Change line 33 from `banner["STAGE 012 — DYNAMIC LOADING"];` to `banner["STAGE 029 — DYNAMIC LOADING"];`.

**Verification command:**
`grep -nE 'Stage 12|STAGE 012|stage12_|stage_12' scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl` — should return no matches (matches on physical-content references like "Stage-11 ansatz" are NOT failures; only file-identification "Stage 12" / "STAGE 012" / "stage12_" / "stage_12" tokens are).

## F2 — insufficient_verification (selected odd coefficient not asserted as combined identity)

**Targets:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:246-283`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:191-229`

**Issue:** Paper eq `selected-odd` claims `delta D_-^odd(omega) = -i Gamma_port * beta_0 * (v.e_-)^2 * omega^5 + O(omega^7)`, where `Gamma_port = a^5/(27 c_s^5)`. Both scripts compute `beta_5 = Gamma_port * beta_0` and `kappa_sel^2 = (v.e_-)^2 |_{alpha=alpha_0}` but never assert that their printed `odd_projection`/`oddProjection` equals the paper's combined symbolic form. Add the assertion.

**Required change:**

(a) In sympy `selected_mode_projection`, insert the following BEFORE the final print at line 281 (i.e. after the existing `print("kappa_sel^2 = ...")` block at lines 274-275 and after the existing weak/strong limit checks at lines 278-279):

```python
    # Anchor the paper's eq selected-odd against the script-computed pieces.
    delta_D_paper = (
        -sp.I * Gamma_port
        * (Omega_U**2 * lambda_W + lambda_R * lambda_U) ** 2
        / (Omega_U**2 * Omega_W**2 - lambda_R**2 * sigma) ** 2
        * kappa_sel_sq
        * omega**5
    )
    delta_D_script = -sp.I * beta5 * kappa_sel_sq * omega**5
    expect_zero(
        "delta D_-^odd (script) - delta D_-^odd (paper formula)",
        delta_D_script - delta_D_paper,
    )
```

(b) In Mathematica, insert the following BEFORE `Print[""];` at line 231 (i.e. after the existing `Print["delta D_-^(odd)(omega) template = ", fmt[oddProjection]];` at line 229):

```mathematica
deltaDPaper = -I*GammaPort
              * (OmegaU^2*lambdaW + lambdaR*lambdaU)^2 / delta0^2
              * kappaSelSq
              * omega^5;
deltaDScript = -I*beta5*kappaSelSq*omega^5;
expectZero[
  "delta D_-^odd (script) - delta D_-^odd (paper formula)",
  deltaDScript - deltaDPaper
];
```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 029` and `redteam exec-mathematica 029` and confirm the new line `delta D_-^odd (script) - delta D_-^odd (paper formula) = 0` (sympy) and `PASS: delta D_-^odd (script) - delta D_-^odd (paper formula)` (mathematica) appears in each output and each script exits 0.

## Applied: F2

files_changed: scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:291-303; mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:222-231
summary: Added selected-odd combined-identity assertions tying the script-computed odd coefficient to the paper formula in both engines.
deviation: Parenthesized the Mathematica multi-line deltaDPaper RHS to preserve the requested formula safely in .wl script form.

## F3 — insufficient_verification (sympy lacks direct eigenvector projection cross-check)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:268-279`

**Issue:** Mathematica asserts `kappa_sel^2 closed-form vs eigenvector projection` (line 223) — Hellmann-Feynman vs direct `Eigensystem` projection of `v` onto `e_-`. SymPy only has the Hellmann-Feynman form plus weak/strong limit checks. Add the SymPy mirror of M8 so that sympy independently verifies the closed-form `kappa_sel^2` against a direct nullspace/eigenvector projection.

**Required change:**

In sympy `selected_mode_projection`, insert the following BEFORE the weak-loading limit check at line 278 (i.e. after `kappa_sel_sq = sp.simplify(kappa_sel_template.subs(al, alpha0))` at line 273 and its print at lines 274-275, but BEFORE the existing `expect_zero("weak-loading ...")` call):

```python
    # Direct eigenvector projection of v onto the lower eigenvector of K_eff(al).
    # Independent of the Hellmann-Feynman derivation above; the two must agree.
    K_eff_al = sp.Matrix(
        [
            [K0t - al * kappa0**2, -al * kappa0 * kappa1],
            [-al * kappa0 * kappa1, K1t - al * kappa1**2],
        ]
    )
    null = (K_eff_al - lam_minus_template * sp.eye(2)).nullspace()
    assert null, "lower-eigenvector nullspace is empty"
    vec_lo = null[0]
    norm_sq = sp.simplify((vec_lo.T * vec_lo)[0])
    kappa_sel_sq_direct_template = sp.simplify(((vec_lo.T * v)[0]) ** 2 / norm_sq)
    expect_zero(
        "kappa_sel^2 closed-form vs eigenvector projection",
        sp.simplify(kappa_sel_template - kappa_sel_sq_direct_template),
    )
```

The `nullspace` route (rather than `eigenvects()`) avoids SymPy's occasional `CRootOf` rendering on symbolic 2x2 matrices and is the more direct realization of "the eigenvector whose eigenvalue equals `lam_minus_template`".

**Verification command:**
After Codex applies, `redteam exec-sympy 029` outputs the new line `kappa_sel^2 closed-form vs eigenvector projection = 0` and the script exits 0.

**Self-test verification for F3 (auditor done):**
- At `al = 0`: `lam_minus_template = K0t`, `K_eff_al - K0t*I = diag(0, DeltaK_ax)`, nullspace = span((1, 0)). Then `(v . (1, 0))^2 / 1 = kappa_0^2`. Hellmann-Feynman at `al=0` also gives `kappa_0^2`. Difference = 0. PASS.
- At `al -> infty`: dominant eigenvalue `-> -al * sigma`, eigenvector `-> v / ||v||`. Then `(v . v)^2 / ||v||^2 = ||v||^2 = sigma`. Hellmann-Feynman limit also gives `sigma`. Difference -> 0. PASS.
- In the interior, both expressions are algebraic functions of `(al, DeltaK_ax)` (after the K0/Xi_0 cancellations); they're equal by the standard Hellmann-Feynman identity for a 2x2 symmetric matrix linear in a parameter (`K_eff_al = M_0 - al * v v^T`, with `M_0 = diag(K0t, K1t)`).

## Applied: F3

files_changed: scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:269-285
summary: Added the SymPy direct nullspace eigenvector projection check against the Hellmann-Feynman kappa_sel^2 template.
deviation: none

## F4 — paper_misalignment (alpha_crit verified in scripts, not in paper card)

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_029.tex:106` quote: *"Stage~029 outputs the Schur-complement split \eqref{eq:app-stage029-schur-split}, the static branch data \eqref{eq:app-stage029-static-data}, the outgoing transfer coefficient \eqref{eq:app-stage029-beta0}, and the selected odd coefficient \eqref{eq:app-stage029-selected-odd}."* (no mention of alpha_crit)
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_029.tex:108-113` quote: `\stagefield{Checks}` lists four items (I_2/vv^T separation, positivity of alpha_0, beta_0 as square/squared-determinant, projection-onto-e_- inserts (v.e_-)^2) — alpha_crit not present.
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage029_dynamic_loading.md:173-186` quote: *"### 3.1 Refined conservative softening threshold ... alpha_crit = (K_0 - Xi_0)(K_1 - Xi_0) / [ (K_1 - Xi_0) kappa_0^2 + (K_0 - Xi_0) kappa_1^2 ]."* — the notes do treat alpha_crit as a substantive stage result.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:189-194` quote:
  ```
  al = sp.symbols('alpha_load', real=True)
  det_template = sp.simplify(((K0t - al*kappa0**2)*(K1t - al*kappa1**2) - al**2*kappa0**2*kappa1**2))
  alpha_crit = sp.solve(sp.Eq(det_template, 0), al)[0]
  alpha_crit_expected = sp.simplify(K0t * K1t / (K1t * kappa0**2 + K0t * kappa1**2))
  print("alpha_crit =", alpha_crit)
  expect_zero("alpha_crit - expected", alpha_crit - alpha_crit_expected)
  ```
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:159-166` quote:
  ```
  al = Symbol["alphaLoad"];
  detTemplate = FullSimplify[Det[{{k0t - al*kappa0Sq, -al*kappa0*kappa1}, {-al*kappa0*kappa1, k1t - al*kappa1Sq}}], Assumptions -> $Assumptions];
  alphaCrit = FullSimplify[k0t*k1t/(k1t*kappa0Sq + k0t*kappa1Sq), Assumptions -> $Assumptions];
  Print["alpha_crit = ", fmt[alphaCrit]];
  expectZero["det(alpha_crit)", detTemplate /. al -> alphaCrit];
  ```

## Resolve before fix_loop

**RESOLVED** (2026-05-25, batch II.1 v2): User approved direction **(b)** — trim `alpha_crit` from stage 029 scripts. Stage 031 owns the refined threshold (boxed in `paper/stages/stage_031.tex:43` as `alpha_{\rm crit}=AB/(B\kappa_0^2+A\kappa_1^2)`, Output line 65; and verified in `scripts/moving_throat_pde_stage031_..._sympy_audit.py:87,94,116` and `mathematica/moving_throat_pde_stage031_..._mathematica_audit.wl:59,61,71`).

Applied via Codex apply session 2026-05-25 (see `redteam/resolutions/batch_II1_paper_alignment.md` § Apply log → Q2). Stage 029 scripts re-run post-trim, both exit 0.

F1 also resolved (relabel Stage 12 → Stage 029) per Q1 approval; see Apply log → Q1.

Remaining work for fix_loop: F2 and F3 only (both insufficient_verification, script-side).

## Applied: F1

files_changed: scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:3,5; mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:33
summary: Relabeled docstring "Stage 12" → "Stage 029" and Mathematica banner "STAGE 012" → "STAGE 029" per user-approved Q1 (a) in batch II.1 v2 paper_alignment resolution.
deviation: none

## Applied: F4

files_changed: scripts/moving_throat_pde_stage029_dynamic_loading_sympy_audit.py:17,189-194; mathematica/moving_throat_pde_stage029_dynamic_loading_mathematica_audit.wl:159-166
summary: Removed alpha_crit definition + assertion block from stage 029 audits per user-approved Q2 (b). Stage 031 owns alpha_crit derivation (destination-verified: paper/stages/stage_031.tex:43,65; scripts/moving_throat_pde_stage031_..._sympy_audit.py:87,94,116; mathematica/moving_throat_pde_stage031_..._mathematica_audit.wl:59,61,71). Both stage 029 scripts re-run post-trim, exit 0.
deviation: none
