---
unit_id: 076
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md
  paper_appendix: present
---

# Audit unit 076 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_076.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 130 — the stage line is the only relevant entry)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.txt`

## What the paper claims

Stage 076 expresses the remaining wall-depth datum `Theta_w` in normalized Family-1 density variables. The `\stagefield{Output}` quotes the boxed identity

> \(\Theta_w = 25 \lambda_\mu^2 \rho_w^2\)

(stage_076.tex line 17). The appendix row at part03 line 130 echoes the same identity. The notes file derives the chain in four steps: (1) the exact n=5 EOS gives `h = m c_s^2 / 4`; (2) the local enthalpy lock `mu_* = lambda_mu h_w` reduces `Theta_w = lambda_mu^2 m^2 rho_w^2 c_(s,w)^2 / (4 hbar^2)`; (3) the Stage-57 healing-width lock `ell = hbar/(2 m c_(s,w))` reduces this to `Theta_w = lambda_mu^2 rho_w^2 / (16 ell^2)`; (4) the Family-1 reference-branch convention `ell/a = 1/20` with `a = 1` in normalized wall coordinates gives `Theta_w = 25 lambda_mu^2 rho_w^2`. The notes treat the `1/20` as a reference-branch input quoted from the explicit Family-1 wall description, not as something re-derived here.

## What the script claims to verify

The sympy docstring lists four checks: (1) the n=5 enthalpy identity `h = m c_s^2 / 4`, (2) the exact algebraic reduction of `Theta_w` under `mu_* = lambda_mu h_w`, (3) the healing-lock reduction to `lambda_mu^2 rho_w^2 / (16 ell^2)`, and (4) the reference-branch form `Theta_w = 25 lambda_mu^2 rho_w^2`. Internally the polytropic exponent is parameterized as `n_poly`, with `U = rho * integrate(P/rho^2, rho)` (sympy line 41–42 / mathematica line 34), `cs2 = diff(P, rho)/m`, `h = diff(U, rho)`, then specialized to `n_poly = 5` for the assertion. The reference factor `25` is *computed* in-script from `(1/ref_factor)^2 / 16` with `ref_factor = 1/20`, with a `TODO(provenance)` inline comment acknowledging the input is not derived here. A non-tautology guard verifies the `n=3` residual is nonzero. The Mathematica script mirrors all four checks using `Integrate`, `Solve`, and `FullSimplify`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `h = m c_s^2 / 4` for n=5 (notes §1) | A1/B1 — `h_n5 - m*cs2_n5/4 == 0` derived from polytrope `U` integral | match |
| `Theta_w = lambda_mu^2 m^2 rho_w^2 c_sw^2 / (4 hbar^2)` (notes §2) | A2/B2 — `Theta_w - Theta_w_alt == 0` where alt is the same expression in `(2x)^2` form | partial (algebraic identity only — see F1) |
| `Theta_w = lambda_mu^2 rho_w^2 / (16 ell^2)` (notes §3) | A3/B3 — substitute `csw_from_ell` derived via `solve`, check against hand-typed target | match |
| `Theta_w = 25 lambda_mu^2 rho_w^2` (stage_076.tex line 17, notes §4) | A4/B4 — `Theta_ref_norm - (1/ref_factor)^2/16 * lambda_mu^2 rho_w^2 == 0` with `ref_factor = 1/20` | match |
| Reference-branch input `ell/a = 1/20` (notes §4) | hardcoded literal `ref_factor = 1/20` with `TODO(provenance)` comment | partial (acknowledged in-script; not derived here) |

Paper-side numeric `25` is reproduced exactly. No `target_mismatch`, no `value_mismatch`. The headline output aligns; the underlying check structure has one residual tautology (A2) and one acknowledged-but-not-cited provenance input (`1/20`).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50 | `expect_zero("n=5 enthalpy identity", h_n5 - m*cs2_n5/4)` | claim 1 (n=5 enthalpy identity) | yes (`U` from `integrate(P/rho^2, rho)` parameterized in `n_poly`; n=3 non-tautology guard at line 52–55) |
| A1b | sympy | 53 | non-tautology guard: `n=3 residual` nonzero | claim 1 negative control | yes (residual `3 K rho^2 / 4` in output line 8) |
| A2 | sympy | 65 | `expect_zero("Theta_w vs alternative-form derivation", Theta_w - Theta_w_alt)` | claim 2 | no — pure algebraic identity `(2x)^2 == 4 x^2`; both expressions are `mu_star_solved`-based with the same numerator |
| A3 | sympy | 74 | `expect_zero("healing-lock reduction", Theta_w_in_ell - Theta_heal_target)` | claim 3 | yes (`csw_from_ell` is `sp.solve(healing_condition, csw)[0]` on line 71) |
| A4 | sympy | 88 | `expect_zero("normalized reference factor", Theta_ref_norm - ref_target)` | claim 4 (Output line 17) | partial — `ref_target` is computed from `ref_factor = 1/20`; the value `1/20` is hardcoded (see F2) |
| B1 | math | 43 | `expectZero["n=5 enthalpy identity", hRhoN5 - mpsi*csSqN5/4]` | claim 1 | yes (mirror of A1; `uPerMass = Integrate[press/rho^2, rho]` on line 34) |
| B1b | math | 48 | guard: `n=3` residual nonzero | claim 1 negative control | yes |
| B2 | math | 57 | `expectZero["Theta_w vs alternative-form derivation", thetaW - thetaCanonical]` | claim 2 | no — mirror of A2; same `(2x)^2 vs 4x^2` identity |
| B3 | math | 65 | `expectZero["healing-lock reduction", thetaWInEll - thetaHealReduced]` | claim 3 | yes (`cSwFromEll` via `Solve` on line 60) |
| B4 | math | 77 | `expectZero["normalized reference factor", thetaRefNorm - refTarget]` | claim 4 | partial — mirror of A4 |

A1 / A3 / A4 are now genuine in v2 (the v1 critique that `U := P/4` was hand-typed has been addressed by deriving `U` from the polytrope integral; the v1 critique that `cs_sub` was hand-typed has been addressed by `sp.solve`; the v1 docstring/assertion 25 vs 25/4 mismatch is resolved — the docstring now matches the paper's `25`). The remaining script-side issues are A2 (which became purely algebraic after the rewrite — both sides come from the same `mu_star_solved` symbol via `(2*x)^2 = 4*x^2`) and the `1/20` provenance, which the script self-flags via inline `TODO(provenance)`.

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:61-65`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl:54-57`

**What's wrong:**
The "Theta_w vs alternative-form derivation" check compares two expressions that are *trivially the same* up to a `(2a)^2 = 4 a^2` rewrite. From the sympy script:

```
Theta_w     = sp.simplify(4 * rho_w**2 * mu_star_solved**2 / (hbar**2 * csw**2))
Theta_w_alt = sp.simplify((2 * rho_w * mu_star_solved / (hbar * csw))**2)
...
expect_zero("Theta_w vs alternative-form derivation", Theta_w - Theta_w_alt)
```

`(2 a / b)^2 = 4 a^2 / b^2` identically for any `a, b` — there is no physical content in this check. Both sides use the same `mu_star_solved` and `hbar` and `csw`. Mathematica reproduces the same pattern on lines 54–57 (`thetaW` versus `thetaCanonical`).

The paper's claim 2 (notes §2: `Theta_w = lambda_mu^2 m^2 rho_w^2 c_sw^2 / (4 hbar^2)`) is supposed to be verified by **deriving** `Theta_w` from the definition `Theta_w := 4 rho_w^2 mu_*^2 / (hbar^2 c_sw^2)` *and substituting* the enthalpy-lock value `mu_* = lambda_mu m c_sw^2 / 4`, then checking the simplified form against the notes target. The check that is actually present skips that step: it just verifies that `(2x)^2 = 4 x^2`, which is true for any `x`.

**Why this matters:**
If the enthalpy-lock substitution `mu_* = lambda_mu m c_sw^2 / 4` were written as `mu_* = lambda_mu m c_sw^2 / 5` (a different factor), A2 would still pass (both sides would simply use the wrong `mu_star_solved`). The check does not exercise the load-bearing `/4` factor in the enthalpy lock — only the algebraic identity. Claim 2 from the notes is *not* actually checked by A2; the substantive check would compare `Theta_w` to the closed-form target `lambda_mu^2 m^2 rho_w^2 c_sw^2 / (4 hbar^2)` written *independently* from the construction of `mu_star_solved`.

**Required change:**
Replace `Theta_w_alt` (sympy line 63) with a hand-typed target form expressed in the simplified variables — i.e., the closed form from notes §2:

```python
Theta_target = sp.Rational(1, 4) * lambda_mu**2 * m**2 * rho_w**2 * csw**2 / hbar**2
expect_zero("Theta_w under enthalpy lock", Theta_w - Theta_target)
```

This way the check verifies that the `Theta_w` derived from `mu_star_solved` via the enthalpy-lock equation simplifies to the closed-form expression the notes claim, rather than verifying a trivial algebraic identity. Apply the analogous change in the Mathematica script: replace `thetaCanonical = FullSimplify[(2*rhoW*muStarSolved/(hbar*cSw))^2, ...]` (line 55) with `thetaTarget = FullSimplify[(1/4)*lambdaMu^2*mpsi^2*rhoW^2*cSw^2/hbar^2, Assumptions -> $Assumptions]` and update the `expectZero` argument on line 57 accordingly.

**Verification:**
After Codex applies, the sympy output line for "Theta_w vs alternative-form derivation" should be renamed (e.g., "Theta_w under enthalpy lock") with residual `0`. The mathematica output should similarly show the renamed check passing. If someone modifies the enthalpy-lock factor in line 59 (`enthalpy_lock = mu_star_sym - lambda_mu * m * csw**2 / 4`) to `/5`, the new A2 check must fail. Both scripts must still exit 0.

### F2 — hardcoded_result

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py:78-80`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl:68-71`

**What's wrong:**
`ref_factor = sp.Rational(1, 20)` (sympy line 80; `refFactor = 1/20` mathematica line 71) is the load-bearing input that produces the headline `25` (since `(1/ref_factor)^2 / 16 = 400/16 = 25`). The inline comment acknowledges the missing provenance:

```python
# Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20.
# TODO(provenance): cite the upstream stage that fixes ref_factor. This factor is
# the load-bearing piece of the "25" in the normalized reference identity.
```

The notes file §4 introduces `ell/a = 1/20` as a Family-1 reference-branch convention (notes line 101). That is the upstream anchor; the script's TODO comment can be discharged by citing it. This is not a paper-misalignment (the paper and the script agree on `25` and on `1/20`), but the TODO is currently undischarged and the citation belongs in the script.

**Why this matters:**
The `1/20` is the sole non-derived input feeding the headline `25`. Without a citation, a downstream reader cannot verify in-script which upstream stage / which notes section the value comes from. The math itself is correct.

**Required change:**
Replace the `TODO(provenance)` comment block on sympy lines 77–80 with a comment naming the upstream source. Suggested text:

```python
# Reference-branch convention: ell = a * ref_factor with ref_factor = 1/20,
# from the Family-1 reference-branch description carried forward as input to this stage
# (notes/stages/moving_throat_pde_stage076_n5_wall_depth_lock.md §4).
ref_factor = sp.Rational(1, 20)
```

Apply the analogous comment change in the Mathematica `.wl` file on lines 68–71. No assertion changes; the numeric output should be byte-identical after the edit.

**Verification:**
After Codex applies, the inline comment in both scripts cites the notes §4 anchor. The `TODO(provenance)` token disappears. All assertion residuals unchanged; both scripts still exit 0.

## Independent-derivation check (Mathematica)

In v2, both engines now route the load-bearing steps through CAS solvers:
- `U` is derived via `Integrate[press/rho^2, rho]` (mathematica line 34) and `sp.integrate(P/rho**2, rho)` (sympy line 41).
- `csw_from_ell` is derived via `Solve[healingCondition == 0, cSw]` (mathematica line 60) and `sp.solve(healing_condition, csw)` (sympy line 71).
- The reference factor `25` is computed from `(1/ref_factor)^2 / 16` rather than hand-typed.

Both scripts use the same algorithmic skeleton (declare positive-real symbols → derive `U` via polytrope integral → specialize to n=5 → assert against `m cs^2 / 4` → introduce `mu_star_sym` → solve enthalpy lock → form `Theta_w` → solve healing condition → substitute → check against target → substitute `ref_factor` → check `25` reduction). The Mathematica script is structurally a port — same variable choreography, same intermediate sequence, same assertion order — but in v2 each engine *does* perform genuine derivations (Integrate, Solve) rather than copy-paste algebra. The remaining transliteration concern is reduced from "echoes algebra" to "follows the same recipe", which the audit policy permits when both engines reach the result via real symbolic operations.

The one remaining tautology (A2) is mirrored exactly in B2 (lines 54–57): both engines compare `(2x)^2` against `4 x^2` from the same `mu_star_solved`. Fixing A2 per F1 should also fix B2 (the directive specifies parallel edits). No additional `mathematica_transliteration` finding is raised on top of F1.

## Engine cross-check

Both engines produce zero residuals on all four assertions and identical symbolic intermediates:

| Assertion | SymPy `.txt` | Mathematica `.txt` |
|-----------|-------------|---------------------|
| n=5 enthalpy identity | `0` (line 7) | `0` (lines 7–8 PASS) |
| n=3 residual (nonzero guard) | `3*K*rho**2/4` (line 8) | `(3*kConst*rho^2)/4` (line 9) |
| Theta_w vs alternative-form | `0` (line 10) | `0` (lines 11–12 PASS) |
| healing-lock reduction | `0` (line 12) | `0` (lines 14–15 PASS) |
| normalized reference factor | `0` (line 16) | `0` (lines 19–20 PASS) |

Symbolic intermediates match:

| Quantity | SymPy | Mathematica |
|----------|-------|-------------|
| `c_s^2(rho) [n=5]` | `5*K*rho**4/m` | `(5*kConst*rho^4)/mpsi` |
| `h(rho) [n=5]` | `5*K*rho**4/4` | `(5*kConst*rho^4)/4` |
| `Theta_w` (enthalpy lock) | `c_sw**2*lambda_mu**2*m**2*rho_w**2/(4*hbar**2)` | `(cSw^2*lambdaMu^2*mpsi^2*rhoW^2)/(4*hbar^2)` |
| `ell` from healing condition | `hbar/(2*c_sw*m)` | `hbar/(2*cSw*mpsi)` |
| `Theta_w` (healing lock) | `lambda_mu**2*rho_w**2/(16*ell**2)` | `(lambdaMu^2*rhoW^2)/(16*ell^2)` |
| `Theta_w` (ref branch, general a) | `25*lambda_mu**2*rho_w**2/a**2` | `(25*lambdaMu^2*rhoW^2)/a^2` |
| `Theta_w` (ref branch, normalized) | `25*lambda_mu**2*rho_w**2` | `25*lambdaMu^2*rhoW^2` |

Outputs are fresh: sympy script May 22 23:19, sympy output May 22 23:23; mathematica script May 22 23:20, mathematica output May 22 23:23. No `stale_output` finding.

Cosmetic note (not a finding): both banner strings still read "STAGE 59" / "STAGE 059" inside the body (sympy line 30, mathematica line 26), while the final mathematica line says "Stage 076 Mathematica audit passed." The stage was renumbered to 076; the banner strings were not updated. Not blocking; the file/path/notes all align on 076.

## Verdict justification

Verdict: `findings`, count 2, paper alignment: aligned. The headline identity `Theta_w = 25 lambda_mu^2 rho_w^2` is reproduced exactly in both engines and matches the paper card, notes §4, and the appendix row. v1's three findings have been mostly addressed: A1 now derives `U` via a polytrope integral with an `n=3` non-tautology guard, A3 routes the healing-lock substitution through `sp.solve` / `Solve`, A4 computes `25` from `ref_factor = 1/20` rather than hand-typing it, and the v1 docstring/assertion `25/4` vs `25` mismatch is resolved (docstring now reads `25`). What remains: A2 is still a trivial algebraic identity `(2x)^2 = 4 x^2` that does not exercise the enthalpy-lock factor of `1/4` (F1); and `ref_factor = 1/20` is hardcoded with an undischarged `TODO(provenance)` comment although the upstream anchor (notes §4) is unambiguous (F2). Neither finding propagates downstream — the *output values* are correct and would remain correct after the fixes — so neither is `CRITICAL_DOWNSTREAM`. No `paper_misalignment`: the paper, notes, appendix row, docstring, and assertion all now agree on `25`. No `UNFIXABLE`: both findings are local edits.

## Self-test notes

(1) **Variable independence / derivatives**: For F1's proposed `Theta_target = (1/4) * lambda_mu^2 * m^2 * rho_w^2 * csw^2 / hbar^2`, this is `csw`-dependent and `m`-dependent. The `Theta_w` expression already has the same symbolic content after `mu_star_solved` is substituted (since `Theta_w = lambda_mu^2 m^2 rho_w^2 csw^2 / (4 hbar^2)` per the v1 output line 9 and the current output line 9 — verified equal). The residual `Theta_w - Theta_target` is symbolically zero only because the `1/4` is correctly placed; if the enthalpy lock factor were changed to `/5`, the residual would not cancel. Verified non-tautological. (2) **Symmetry/parity**: No symmetric-domain integrals introduced; `Integrate[K rho^(n-2), rho]` is a polynomial primitive, no parity trap. (3) **Trivial-case pre-check**: For F1's new check, mentally substitute `lambda_mu = 1, m = 1, rho_w = 1, csw = 2, hbar = 1`: `mu_star_solved = 1 * 1 * 4 / 4 = 1`; `Theta_w = 4 * 1 * 1 / (1 * 4) = 1`; `Theta_target = (1/4) * 1 * 1 * 1 * 4 / 1 = 1`; residual = 0. With `enthalpy_lock` factor changed to `/5`: `mu_star_solved = 4/5`; `Theta_w = 4 * 1 * (4/5)^2 / (1 * 4) = 16/25`; `Theta_target = (1/4) * 4 = 1`; residual = `16/25 - 1 != 0`. Correctly distinguishes. (4) **Path specifications**: No `missing_verification_script` findings; both files exist at the paths confirmed in the prompt. (5) **Paper round-trip**: F1's fix replaces a tautological check with a hand-typed target form that *exactly* equals the notes §2 expression `lambda_mu^2 m^2 rho_w^2 c_sw^2 / (4 hbar^2)`. F2's fix is comment-only and changes no numeric content. Neither introduces a new paper_misalignment; the `25` headline remains the verified output.
