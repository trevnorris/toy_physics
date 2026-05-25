---
unit_id: 076
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 076

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
Both scripts were restructured to derive the four claimed identities from independent routes rather than postulating them on both sides of `expect_zero` (see `redteam/exec_logs/stage_076_diff.patch`).

In `scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:34-55`:
- `P = K * rho**n_poly` introduces an exponent-indexed EOS (line 36).
- `U_per_mass = sp.integrate(P / rho**2, rho)` and `U = rho * U_per_mass` derive `U` from the polytrope integral (lines 41-42) rather than the original `U = K*rho**5/4`.
- `cs2 = sp.diff(P, rho) / m` and `h = sp.diff(U, rho)` are computed symbolically and then specialized to `n_poly = 5` for the assertion (lines 43-50).
- A `residual_n3` check raises `AssertionError` if the identity also holds at `n=3` (lines 52-55), preventing the n-independent failure mode the auditor flagged in the self-test notes.

In `scripts/...py:57-65` the `Theta_w` block now introduces a separate symbol `mu_star_sym`, solves the enthalpy-lock condition via `sp.solve`, and cross-checks against an independent algebraic form `(2 rho_w mu_star / (hbar c_sw))^2`.

In `scripts/...py:67-75` the healing block solves `csw - hbar/(2 m ell) == 0` for both `ell` and `csw` and substitutes the `csw`-solution into `Theta_w` to check the reduction. The substitution rule is no longer hand-typed.

In `scripts/...py:77-88` the `25` in the reference assertion is now derived from `(1/ref_factor)**2 / 16` rather than typed independently.

The Mathematica script (`mathematica/...wl:28-77`) mirrors these changes using `Integrate[]`, `Solve[]`, `ConditionalExpression` stripping, and the renamed intermediate `thetaCanonical`.

**Assessment:**
The edits substantively de-tautologize the four assertions:
- A1: `h_n5 - m*cs2_n5/4 = 0` now passes only because `U` was derived from `Integrate[P/rho^2]` at `n_poly=5`, AND because the parallel `n=3` residual is verified to be nonzero (output line 8: `n=3 residual (should be NONZERO) = 3*K*rho**2/4`). The pair (`expect_zero` at n=5, `raise` if zero at n=3) constitutes a genuine non-vacuous identity.
- A2: `Theta_w - Theta_w_alt = 0` now compares two algebraically distinct forms (`4 rho_w^2 mu_star^2/(hbar^2 c_sw^2)` vs `(2 rho_w mu_star/(hbar c_sw))^2`) of the same symbolic quantity, with `mu_star` obtained from `sp.solve`. This is still a fairly mild check (the two forms are equal by trivial squaring algebra), but at least the `mu_star` value flows through `solve` rather than being typed twice.
- A3: `Theta_w_in_ell - Theta_heal_target = 0` now substitutes the `solve`-derived `csw_from_ell` into `Theta_w`. The output line 11 confirms `ell from healing condition = hbar/(2*c_sw*m)` matches the prior hand-typed substitution.
- A4: `Theta_ref_norm - ref_target = 0` now derives `ref_target` from `(1/ref_factor)**2/16 * lambda_mu^2 rho_w^2`, so changing `ref_factor` automatically updates the target. This is still mostly a consistency check on `ref_factor`'s arithmetic, but the link is now algebraic rather than two independently-typed numbers.

The codex `## Applied: F1` block flags a deviation: the directive's literal `P = K*rho**(1 + 1/n_poly)` form (per-mass enthalpy convention) would have given `h = m (n+1)/n K rho^(1/n)` and `cs2 = (1+1/n) K rho^(1/n) / m`, so `m cs2 / 4 = (n+1)/(4n) K rho^(1/n)` and the identity `h = m cs2/4` would require `(n+1)/n = (n+1)/(4n)`, i.e. `1 = 1/4` — never true. The directive's literal proposal was internally inconsistent. Codex's substitution `P = K*rho**n_poly` (a power-law EOS with exponent `n`, not a polytrope index) repairs this: at `n=5`, `U_per_mass = K rho^(n-1)/(n-1)`, `U = K rho^n/(n-1)`, so `h = n K rho^(n-1)/(n-1)`, `cs2 = n K rho^(n-1)/m`, and `m cs2 / 4 = n K rho^(n-1)/4`. The identity `h = m cs2/4` then requires `n/(n-1) = n/4`, i.e. `n-1 = 4`, i.e. `n=5` — which is exactly what the script claims to verify. The deviation is therefore not just acceptable but necessary for the non-tautology to be substantive. (The auditor's own self-test notes had to retract a similar proposal for the same reason.)

Mathematica mirrors the same fix and additionally exercises `Solve` and `Integrate`, which are absent from the SymPy script's solution path. Both engines now have independent routes through their respective `solve`/`Solve` and `integrate`/`Integrate` calls.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:9` — docstring line 4 now reads `Theta_w = 25 lambda_mu^2 rho_w^2` (was `(25/4)`). This is the typo fix the user resolved (the assertion and the paper consistently use `25`, the docstring's `25/4` was stale).
- `scripts/...py:77-80` — three-line provenance comment above the `ref_factor = sp.Rational(1, 20)` line, plus the inline `(see F2 below for provenance)` note on the assignment line.
- `mathematica/...wl:68-71` — mirror provenance comment block above `refFactor = 1/20`.

**Assessment:**
The docstring fix is a doc-only change with no symbolic or numeric impact (output line 15: `Theta_w (reference branch, normalized wall units) = 25*lambda_mu**2*rho_w**2`, unchanged from the pre-fix saved output cited in the audit report). The provenance comments are non-executable annotations flagging the unresolved upstream citation. Both targets of F2 are addressed exactly as the directive instructed; nothing was silently changed in the assertion or the assertion's numeric target. Per the stage-specific context provided by the orchestrator, this is explicitly a docstring typo fix and NOT a material change.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/...wl:55` — `thetaCanonical` (new name, F1 step 3) replaces `thetaExpected`; semantics unchanged.
- `mathematica/...wl:64-66, 72` — `thetaHealReduced` (new name, F3 directive) replaces `thetaHealTarget` throughout the healing and reference blocks (4 occurrences, as required).
- `mathematica/...wl:33-34` — `Integrate[]` derivation of `uPerMass` and `uRho` (F1 step 2) has no SymPy analogue in choreography (SymPy also uses `integrate` but the Mathematica intermediate naming is now distinct).
- `mathematica/...wl:53, 60-61` — `Solve[]` calls for `muStarSolved`, `cSwFromEll`, `ellSolved` (F1 steps 3, 4) introduce derivation steps that SymPy reaches via `sp.solve` but with different intermediate names.

**Assessment:**
`grep -c thetaHealTarget` returns 0 (confirmed) and `grep -c thetaHealReduced` returns 4 (confirmed), meeting the directive's verification command exactly. The variable mapping between the two scripts is no longer 1-to-1 (`Theta_heal_expected (sympy)` ↔ `thetaHealReduced (math)` rather than the original `Theta_heal_expected` ↔ `thetaHealTarget`). The Mathematica script's `Solve`, `ConditionalExpression` stripping, and `Integrate` call chains exercise machinery absent from the SymPy script's flow even where the two engines reach equivalent symbolic conclusions. The cross-engine agreement on residuals (all four = 0) is now structurally non-trivial because the algebraic paths diverge in concrete ways.

## Exec log assessment

**SymPy:** exit=0 (codex_logs/076_iter1.txt: "`python3 scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py` exits 0"). The exec log file `redteam/exec_logs/stage_076_sympy.log` is not present on disk — only `stage_076_diff.patch` was retained in `exec_logs/`. The codex iteration log confirms a clean run, and the saved `.txt` output reflects post-fix content.

Notable lines from the saved `.txt` (`scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt`):
- line 5: `c_s^2(rho)  [n=5] = 5*K*rho**4/m` (matches pre-fix output, n=5 specialization preserves the original numeric)
- line 6: `h(rho)      [n=5] = 5*K*rho**4/4` (matches pre-fix)
- line 7: `n=5 enthalpy identity = 0` (PASS)
- line 8: `n=3 residual (should be NONZERO) = 3*K*rho**2/4` (NEW; explicit non-tautology evidence)
- line 11: `ell from healing condition = hbar/(2*c_sw*m)` (NEW; from `sp.solve`)
- line 15: `Theta_w (reference branch, normalized wall units) = 25*lambda_mu**2*rho_w**2` (matches pre-fix value exactly; F2 docstring-only fix did not change the numeric)

**Mathematica:** exit=0 (codex_logs/076_iter1.txt: "`math -script mathematica/...wl` exits 0"). Saved `.txt` (`mathematica/output/...txt`):
- line 9: `n=3 residual (should be NONZERO) = (3*kConst*rho^2)/4` (NEW; matches SymPy at the symbolic level)
- line 13: `ell from healing condition = hbar/(2*cSw*mpsi)` (NEW; from `Solve`)
- lines 8, 12, 15, 20: explicit `PASS:` lines for all four assertions
- line 18: `Theta_w (reference branch, normalized wall units) = 25*lambdaMu^2*rhoW^2` (unchanged from pre-fix)
- line 22: `Stage 076 Mathematica audit passed.`

**Output freshness:** Both `.txt` outputs (mtimes 1779513788 SymPy, 1779513799 Mathematica) are newer than their respective scripts (1779513577 SymPy, 1779513656 Mathematica), confirming post-fix regeneration. The diff patch was captured at 1779513743 — between the script edits and the output regeneration, consistent with the orchestrator's exec-then-capture sequence.

## Material-change assessment

`material_change`: false.

The only printed numeric/symbolic values that changed between the pre-fix and post-fix saved outputs are *additions* (new diagnostic lines: `n=3 residual = ...`, `ell from healing condition = ...`). Every value that existed in the pre-fix output is reproduced identically post-fix:
- `c_s^2(rho)` is now labeled `[n=5]` but the expression `5*K*rho**4/m` is unchanged.
- `h(rho)` is now labeled `[n=5]` but the expression `5*K*rho**4/4` is unchanged.
- `Theta_w (enthalpy lock)` = `c_sw^2*lambda_mu^2*m^2*rho_w^2/(4*hbar^2)` is unchanged.
- `Theta_w (healing lock)` = `lambda_mu^2*rho_w^2/(16*ell^2)` is unchanged.
- `Theta_w (reference branch, general a)` = `25*lambda_mu^2*rho_w**2/a^2` is unchanged.
- `Theta_w (reference branch, normalized wall units)` = `25*lambda_mu^2*rho_w**2` is unchanged.

No downstream unit can observe any change in the symbolic values exported by this stage. The provenance-route changes (sp.solve, sp.integrate, Solve, Integrate) and the F3 rename (thetaHealTarget → thetaHealReduced) are local to this stage's `.wl`. The F2 docstring fix is a comment-only change. Per the orchestrator's explicit note in the verifier prompt and per the project rule (provenance changes alone do not count as material), this stage stays non-material.

## Side observations (non-blocking)

1. The directive's F1 literal text (`P = K*rho^(1 + 1/n_poly)`) is internally inconsistent with the assertion `h = m cs^2 / 4`; codex's deviation to `P = K*rho^n_poly` was necessary and is documented in the `## Applied: F1` block. The auditor's self-test notes already retracted an analogous `expect_nonzero` instruction for the same reason, so this is a known issue in the audit drafting, not a verification defect.
2. A2's post-fix assertion `Theta_w - Theta_w_alt = 0` is still a relatively soft check — `(2 rho_w mu_star/(hbar c_sw))^2` and `4 rho_w^2 mu_star^2/(hbar^2 c_sw^2)` are equal by exponent algebra. Routing `mu_star` through `sp.solve` adds some structural separation, but a future reviewer may want a deeper independent route for `Theta_w` (e.g., from a primitive energy-per-healing-length definition). Flagged as observation, not a verification block.
3. The Mathematica banner still says `STAGE 059`. The SymPy banner says `STAGE 59`. The directive did not request a fix here, and the auditor mentioned it only in passing in F3. Non-blocking.

## Verdict justification

All three findings are `resolved`. F1's tautological assertions are now genuinely substantive: the n=5 identity is paired with an explicit n=3 fail-check that confirms `(P,U)` machinery is exponent-sensitive; `Theta_w`, the healing-length substitution, and the `25` reference target are all routed through `sp.solve` / `Solve` or derived algebraically from `ref_factor`. F2's docstring typo is fixed exactly as the human resolved (the orchestrator-provided context confirms the docstring was wrong, not the assertion), and provenance comments are added to both scripts. F3's `thetaHealReduced` rename is in place (4 occurrences, 0 of the old name). Both scripts exit 0; the saved `.txt` outputs are fresh and show the new diagnostic lines (`n=3 residual`, `ell from healing condition`) with no change to any pre-existing numeric or symbolic value. `material_change: false` because the only printed values that changed are additions; every previously-printed value is reproduced identically. No regressions in the diff.

stage 076: verified
