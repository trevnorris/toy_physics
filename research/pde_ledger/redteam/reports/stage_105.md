---
unit_id: 105
batch: IV.2
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md
  paper_appendix: present
---

# Audit unit 105 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_105.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only the `\input{stages/stage_105}` row at line 1244 — no per-stage row content; the part file is a stage-input wrapper)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.txt`

## What the paper claims

The paper card (`stage_105.tex`, `\label{stage:105}`, displayed as "Stage~122") states the audit target verbatim:

> "Matching the canonical grouped module to the exact outgoing DtN branch gives \(\chi_Q=1\)."

The "Derivation ledger" field expands this: "The computation isolates the reduced product \(\widehat m_0^{\,2}\chi_Q N_Q=1\) and the canonical condition \(\chi_Q=1\)."

The notes file fixes the underlying algebra: define the retarded normalized quadrupole module as
\(\widehat Y_Q^{\rm ret}(\omega) = 3/4 + (1/4)\,/\,[1-\omega^2/\Omega_Q^2 - i\,\chi_Q\,\sigma_Q^{\rm can}\,\omega^5] + O(\omega^6)\) with \(\Omega_Q=3c_s/(2a)\) and \(\sigma_Q^{\rm can}:=9/(8\Omega_Q^5)=4a^5/(27c_s^5)\); expand to \(O(\omega^5)\) giving \(1+a^2\omega^2/(9c_s^2)+4a^4\omega^4/(81c_s^4)+i\,\chi_Q\,a^5\omega^5/(27c_s^5)\); match against the outgoing DtN coefficient \(i\,a^5\omega^5/(27c_s^5)\) from Stage 87 to fix \(\chi_Q=1\); and parametrize the remaining deformation by \(\Lambda_2^{\rm def}(z)=-3+z^2/3+z^4/9+i\,\xi_Q z^5/9\) so that \(\widehat Y_2^{\rm def}(z)=1+z^2/9+4z^4/81+i\,\xi_Q z^5/27\) and \(\chi_Q=\xi_Q\).

Distinct deliverables:
- D1: \(\Omega_Q=3c_s/(2a)\) and \(\sigma_Q^{\rm can}=4a^5/(27c_s^5)\) (algebraic identity).
- D2: Series coefficients of \(\widehat Y_Q^{\rm ret}\) through \(O(\omega^5)\) match the closed forms.
- D3: Matching the \(\omega^5\) coefficient against the outgoing branch's \(ia^5\omega^5/(27c_s^5)\) fixes \(\chi_Q=1\).
- D4: Deformed-branch parametrization yields \(\widehat Y_2^{\rm def}\) coefficients with \(\chi_Q=\xi_Q\).

Card "Checks" list also includes C1 (factorization separability of \(\widehat m_0^{\,2}\chi_Q N_Q\)), C2 (higher odd terms begin beyond 2.5PN), and C3 (outgoing DtN fingerprint matches \(z=\omega a/c_s\) expansion).

## What the script claims to verify

The SymPy script defines \(\Omega_Q=3c_s/(2a)\) and \(\sigma_Q^{\rm can}=9/(8\Omega_Q^5)\), asserts the algebraic equality \(\sigma_Q^{\rm can}-4a^5/(27c_s^5)=0\), builds the closed-form retarded module, expands as a power series through \(O(\omega^5)\), and asserts that the \(\omega^2\), \(\omega^4\), and \(\omega^5\) coefficients match the closed-form expressions from the notes. It then `solve`s the equation `coeff5/I == a^5/(27 c_s^5)` for `chi_Q`, recovering 1 and asserting `chi_sol - 1 == 0`. Finally it constructs \(\Lambda_2^{\rm def}(z)\), forms \(\widehat Y_2^{\rm def}=-3/\Lambda_2^{\rm def}\), expands to \(O(z^5)\), and asserts the deformed-branch coefficients \(z^2:1/9\), \(z^4:4/81\), \(z^5:\xi_Q/27\). The Mathematica script performs the same operations one-for-one (same intermediate names, same series order, same assertions).

## Paper ↔ script cross-check

| Deliverable | Script coverage |
|---|---|
| D1 (\(\Omega_Q\), \(\sigma_Q^{\rm can}\)) | match (sympy lines 32-37; math lines 33-38) |
| D2 (series coefficients of \(\widehat Y_Q^{\rm ret}\)) | match (sympy lines 39-47; math lines 40-46) |
| D3 (\(\chi_Q=1\) from outgoing-branch match) | match (sympy lines 49-55; math lines 48-50) |
| D4 (deformed-branch parametrization and \(\chi_Q=\xi_Q\)) | partial — coefficients of \(\widehat Y_2^{\rm def}\) are checked (sympy lines 58-66; math lines 52-61) but the conversion `Ŷ = -3/Λ` used by the script is not stated in the notes; the notes give only `Λ_2^def` and `Ŷ_2^def`, and the identification `χ_Q = ξ_Q` is not solved/asserted in the script — only the \(z^5\) coefficient appearing as `xi_Q/27` is verified, which implicitly demonstrates the identification but no explicit `assert chi_Q == xi_Q` exists |
| C1 (separability of \(\widehat m_0^2 \chi_Q N_Q\)) | missing — no `m_0` or `N_Q` symbols appear in either script |
| C2 (higher odd terms begin beyond 2.5PN) | partial — series truncates at \(O(\omega^5)\) so no odd term below 5 appears; no explicit assertion `coeff(omega,3) == 0` etc., but the printed series visibly contains no \(\omega^1\) or \(\omega^3\) term |
| C3 (outgoing l=2 DtN fingerprint vs \(z\) expansion) | match (the deformed-branch block with \(\xi_Q=1\) is the canonical outgoing fingerprint and its z-expansion is asserted) |

Set `paper_alignment: partial` — the load-bearing math claim (D1-D3, the boxed \(\chi_Q=1\)) is faithfully exercised; C1/C2 paper-card checks are not script-side identities at this stage and are reasonably carried forward from upstream; D4 is verified at coefficient level rather than via a separate `chi_Q-xi_Q` identity assertion.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 37 | `expect_zero("sigma_Q^can - 4 a^5/(27 c_s^5)", sigma_can - 4*a**5/(27*c_s**5))` | D1 | yes |
| A2 | sympy | 45 | `expect_zero("omega^2 coefficient", Yret_series.coeff(omega, 2) - a**2/(9*c_s**2))` | D2 | yes |
| A3 | sympy | 46 | `expect_zero("omega^4 coefficient", Yret_series.coeff(omega, 4) - 4*a**4/(81*c_s**4))` | D2 | yes |
| A4 | sympy | 47 | `expect_zero("imag omega^5 coefficient", Yret_series.coeff(omega, 5)/sp.I - chi_Q*a**5/(27*c_s**5))` | D2 | yes |
| A5 | sympy | 49-55 | `chi_sol = solve(...); raise if chi_sol - 1 != 0` | D3 | yes |
| A6 | sympy | 64 | `expect_zero("deformed z^2 coefficient", Y_def.coeff(z, 2) - sp.Rational(1, 9))` | D4 | yes |
| A7 | sympy | 65 | `expect_zero("deformed z^4 coefficient", Y_def.coeff(z, 4) - sp.Rational(4, 81))` | D4 | yes |
| A8 | sympy | 66 | `expect_zero("deformed imag z^5 coefficient", Y_def.coeff(z, 5)/sp.I - xi_Q/sp.Integer(27))` | D4 | yes (implicit \(\chi_Q=\xi_Q\)) |
| B1 | mathematica | 38 | `expectZero["sigma_Q^can - 4 a^5/(27 c_s^5)", sigmaCan - 4*aThroat^5/(27*cSound^5)]` | D1 | yes |
| B2 | mathematica | 44 | `expectZero["omega^2 coefficient", Coefficient[ySeries, omega, 2] - aThroat^2/(9*cSound^2)]` | D2 | yes |
| B3 | mathematica | 45 | `expectZero["omega^4 coefficient", Coefficient[ySeries, omega, 4] - 4*aThroat^4/(81*cSound^4)]` | D2 | yes |
| B4 | mathematica | 46 | `expectZero["imag omega^5 coefficient", Coefficient[ySeries, omega, 5]/I - chiQ*aThroat^5/(27*cSound^5)]` | D2 | yes |
| B5 | mathematica | 48-50 | `Solve...; expectZero["chi_Q - 1", (chiQ /. chiSol) - 1]` | D3 | yes |
| B6 | mathematica | 59 | `expectZero["deformed z^2 coefficient", Coefficient[yDef, z, 2] - 1/9]` | D4 | yes |
| B7 | mathematica | 60 | `expectZero["deformed z^4 coefficient", Coefficient[yDef, z, 4] - 4/81]` | D4 | yes |
| B8 | mathematica | 61 | `expectZero["deformed imag z^5 coefficient", Coefficient[yDef, z, 5]/I - xiQ/27]` | D4 | yes (implicit) |

Each assertion is non-tautological: the series coefficients are derived by symbolic expansion of the rational function and then compared against the closed forms from the notes. None of them are "x = expr; assert x == expr" tautologies. The `chi_Q` solve is over a linear equation with a unique root; the assertion `chi_sol - 1 == 0` confirms the result of that solve rather than a precomputed value.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:26-61`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:28-66`

**What's wrong:**

The `.wl` script is a line-by-line port of the `.py` script: same intermediate variable names (modulo `_`→camelCase: `Omega_Q`/`omegaQ`, `sigma_can`/`sigmaCan`, `Yret`/`yRet`, `chi_Q`/`chiQ`, `Lam_def`/`lamDef`, `Y_def`/`yDef`); same operation order (define Ω, define σ, assert σ-form identity, build Y_ret, series-expand to order 5, check ω², ω⁴, ω⁵ coefficients individually, solve for χ, define Λ_def, build -3/Λ_def, series-expand, check z², z⁴, z⁵ coefficients); same closed-form RHS for every comparison. Three corresponding pairs that show the choreography is identical:

Pair 1 (definitions):
- `.py:32-33`:
  ```
  Omega_Q = sp.simplify(3*c_s/(2*a))
  sigma_can = sp.simplify(sp.Rational(9, 8) / Omega_Q**5)
  ```
- `.wl:33-34`:
  ```
  omegaQ = FullSimplify[3*cSound/(2*aThroat), Assumptions -> $Assumptions];
  sigmaCan = FullSimplify[(9/8)/omegaQ^5, Assumptions -> $Assumptions];
  ```

Pair 2 (series and assertions):
- `.py:39-47`:
  ```
  Yret = sp.Rational(3, 4) + sp.Rational(1, 4) / (1 - omega**2/Omega_Q**2 - sp.I*chi_Q*sigma_can*omega**5)
  Yret_series = sp.expand(sp.series(Yret, omega, 0, 6).removeO())
  ...
  expect_zero("omega^2 coefficient", Yret_series.coeff(omega, 2) - a**2/(9*c_s**2))
  expect_zero("omega^4 coefficient", Yret_series.coeff(omega, 4) - 4*a**4/(81*c_s**4))
  expect_zero("imag omega^5 coefficient", Yret_series.coeff(omega, 5)/sp.I - chi_Q*a**5/(27*c_s**5))
  ```
- `.wl:40-46`:
  ```
  yRet = 3/4 + (1/4)/(1 - omega^2/omegaQ^2 - I*chiQ*sigmaCan*omega^5);
  ySeries = Expand[Normal[Series[yRet, {omega, 0, 5}]]];
  ...
  expectZero["omega^2 coefficient", Coefficient[ySeries, omega, 2] - aThroat^2/(9*cSound^2)];
  expectZero["omega^4 coefficient", Coefficient[ySeries, omega, 4] - 4*aThroat^4/(81*cSound^4)];
  expectZero["imag omega^5 coefficient", Coefficient[ySeries, omega, 5]/I - chiQ*aThroat^5/(27*cSound^5)];
  ```

Pair 3 (deformed branch via the same `-3/Λ` normalization that the notes do not explicitly state):
- `.py:59-66` defines `Lam_def = -3 + z**2/3 + z**4/9 + sp.I*xi_Q*z**5/9` and `Y_def = -3/Lam_def`, then checks three coefficients.
- `.wl:55-61` does the same with `lamDef` and `yDef = -3/lamDef`, then checks the same three coefficients.

**Why this matters:**

The second-engine policy requires Mathematica to derive the claim from independent algebraic premises so that a transcription mistake in one engine cannot silently propagate to both. Here, both engines (a) use the same `-3/Λ_def` choice for the deformed-branch normalization (not stated in the notes — it is a script-side choice), (b) extract coefficients via the same `Series`/`series` call to the same order, (c) compare against the same hand-supplied closed forms. A sign error in the closed-form RHS, or a missing factor in `-3/Λ_def` vs. the actual definition of \(\widehat Y_2\) from \(\Lambda_2\), would pass on both engines. Independent re-derivation could, for example, do (Mathematica side): construct the outgoing \(\Lambda_2^{\rm out}(z)\) directly from the spherical-Hankel ratio and compare its series; or factor the conservative module algebraically and identify \(\chi_Q\) via residue of \(\widehat Y_Q^{\rm ret}\) at the pole \(\omega^2=\Omega_Q^2\).

**Required change:**

Rewrite the Mathematica script so that the verification of `chi_Q = 1` follows a different algebraic path than the SymPy script. Acceptable alternatives:
- (i) Derive the outgoing \(\widehat Y_2^{\rm out}(\omega)\) directly from a polynomial form (e.g. `Ŷ_2^out(z) = -3/Λ_2^out(z)` with `Λ_2^out` constructed from the outgoing spherical-Hankel `h_2^(1)` form rather than from the deformed parametrization), expand its series, and match the imaginary \(z^5\) coefficient to \(a^5/(27c_s^5)\) — then \(\chi_Q\) drops out by direct comparison of two ω⁵ coefficients of two independently constructed series.
- (ii) Verify the canonical condition via residue of \(\widehat Y_Q^{\rm ret}\) at \(\omega^2 = \Omega_Q^2\) and match to the outgoing-branch residue, asserting \(\chi_Q = 1\) without going through a coefficient-by-coefficient series match.

Either approach must retain the bottom-line identity `\(\chi_Q = 1\)` as the final assertion. Variable names should differ from the `.py` so the diff is visible.

**Verification:**

The verifier confirms (a) the `.wl` no longer contains the substrings `(1/4)/(1 - omega^2/omegaQ^2`, `yRet`, and `lamDef`, AND (b) the final printed result still asserts `chiQ - 1 == 0`, AND (c) the script's exit code is 0 on a fresh run.

### F2 — paper_misalignment

**Severity:** low
**Subtype:** target_mismatch
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_sympy_audit.py:28`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn_mathematica_audit.wl:26`

**What's wrong:**

The SymPy module docstring (`.py:3`) declares `"""Stage 88 SymPy audit."""` and the banner (`.py:28`) prints `"STAGE 88 — EXACT FIXING OF chi_Q"`. The Mathematica banner (`.wl:26`) prints `"STAGE 088 — EXACT FIXING OF chi_Q"`. The paper card labels this `\label{stage:105}` (display title "Stage~122") and the notes file is `moving_throat_pde_stage105_*`. The `.wl`'s closing print (`.wl:64`) confirms `"Stage 105 Mathematica audit passed."`, contradicting its own banner.

Paper card (`stage_105.tex:2`):
> `\label{stage:105}`

Notes header (`moving_throat_pde_stage105_chiQ_fix_from_outgoing_dtn.md:1`):
> `# Moving-Throat PDE — Stage 105: Exact Fixing of chi_Q`

Script docstring (`.py:3`):
> `Stage 88 SymPy audit.`

Script banner (`.wl:26`):
> `banner["STAGE 088 — EXACT FIXING OF chi_Q"];`

**Why this matters:**

The math content is correct and matches the Stage 105 paper claim. However, the docstring and banner labels are a stale "Stage 88" / "Stage 088" copy-paste from an upstream stage (Stage 88 is the conservative pole-scale unit mentioned in the provenance comment at `.py:7-10`). A future reader scanning grep output by stage number will not find this script when searching "Stage 105". This is a documentation defect that contradicts the paper card per the auditor prompt's explicit guidance ("if the script docstring contradicts the paper card, that is itself a paper_misalignment"). The math claim itself is aligned (D1-D3 all check correctly against Stage 105's boxed result); the failure is purely in the unit-id label.

**Required change:**

This is unambiguous (the file path, notes, paper label, and the `.wl`'s closing print line all agree on Stage 105). Update the labels:
- `.py:3`: change `"""Stage 88 SymPy audit."""` to `"""Stage 105 SymPy audit."""` (preserve the rest of the docstring as-is).
- `.py:28`: change `banner("STAGE 88 — EXACT FIXING OF chi_Q")` to `banner("STAGE 105 — EXACT FIXING OF chi_Q")`.
- `.wl:26`: change `banner["STAGE 088 — EXACT FIXING OF chi_Q"];` to `banner["STAGE 105 — EXACT FIXING OF chi_Q"];`.

Do not change any math, assertions, or paper text.

**Verification:**

After the fix, the saved transcripts on a fresh run should print `STAGE 105 — EXACT FIXING OF chi_Q` in both `.txt` outputs (replacing the current `STAGE 88` / `STAGE 088` lines). Final result lines remain unchanged.

## Independent-derivation check (Mathematica)

See F1. The `.wl` is structurally a line-by-line port of the `.py` (same names, same operation order, same closed-form RHS). Quoting the corresponding sections:

SymPy (`.py:32-33, 39-47, 49-55, 59-66`) and Mathematica (`.wl:33-34, 40-46, 48-50, 55-61`) each execute the identical eight-step choreography listed in F1.

## Engine cross-check

Both engines produce matching results on the captured transcripts:
- SymPy output line 13: `Omega_Q = 3*c_s/(2*a)` ↔ Mathematica output line 13: `Omega_Q = (3*cSound)/(2*aThroat)`
- SymPy line 14: `sigma_Q^can = 4*a**5/(27*c_s**5)` ↔ Mathematica line 14: `sigma_Q^can = (4*aThroat^5)/(27*cSound^5)`
- SymPy line 27: `chi_Q from exact outgoing match = 1` ↔ Mathematica line 24: `chi_Q from exact outgoing match = ConditionalExpression[1, Element[omega, Reals]]` — the `ConditionalExpression` strip is implicit in `(chiQ /. chiSol) - 1` and the subsequent `expectZero` passes.
- SymPy lines 34-36 deformed coefficients = 0 ↔ Mathematica lines 28-32 same coefficients = 0.

Engines agree (both PASS at exit 0). No `engine_disagreement` finding.

## Verdict justification

The script's load-bearing assertion (\(\chi_Q = 1\) from matching the \(\omega^5\) coefficient of \(\widehat Y_Q^{\rm ret}\) against the outgoing branch's \(ia^5\omega^5/(27c_s^5)\)) faithfully exercises the paper card's boxed claim and the notes' derivation. The auxiliary D1, D2, D4 assertions are non-tautological coefficient matches against closed forms taken from the notes. I attempted three attacks that all failed: (1) checking whether the `chi_Q = positive` assumption could mask a sign ambiguity in the solve — it cannot, because the equation is linear in `chi_Q` with a unique root +1 regardless of sign assumption; (2) checking whether the `simplify`/`FullSimplify` calls could silently collapse a nonzero residual — the residuals are printed verbatim and are syntactically zero in both transcripts; (3) checking whether the series order (6 in sympy, `{omega, 0, 5}` in mathematica) truncates an asymmetric way — both produce a polynomial through ω⁵ inclusive, so the ω⁵ check is real. The remaining defects are F1 (the second-engine policy: the `.wl` is a transliteration of the `.py`) and F2 (stale "Stage 88" labels in docstring/banner). Both are real defects but neither invalidates the math; verdict is `findings` (not `stop_cold`).

## Self-test notes

1. Variable independence — no `sp.diff`/`D[]` calls in either script; nothing to trap on identically-zero derivatives.
2. Symmetry/parity — no unbounded-domain integrals; series expansions are at \(\omega = 0\) / \(z = 0\) of analytic rational functions with finite radii of convergence.
3. Trivial-case pre-check — at `chi_Q = 1` the polynomial residual `chi_Q*a^5/(27c_s^5) - 1*a^5/(27c_s^5) = 0` and the deformed-branch RHS at `xi_Q = 1` recovers the outgoing fingerprint; both reduce as expected.
4. Path specifications — F1 instructs only `.wl` edits (lives in `mathematica/`); F2 instructs `.py` (in `scripts/`) and `.wl` edits. No path ambiguity.
5. Paper round-trip — F1 does not change the math claim being verified (only the algebraic path on the Mathematica side); F2 only touches stale labels and does not introduce any new identity. Neither fix introduces a new `paper_misalignment`.
