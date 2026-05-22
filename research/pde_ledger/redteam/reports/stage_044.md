---
unit_id: 044
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00-06:00
verdict: findings
stop_cold: null
findings_count: 5
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 044 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.txt`

File timestamps (mtimes):
- .py = Apr 1 12:39
- .py.txt = May 11 12:42 (newer than script — fresh)
- .wl = May 11 11:56
- .wl.txt = May 11 12:50 (newer than script — fresh)

## What the script claims to verify

The script claims to verify a "continuum-selected rank-2 closure" for the moving-throat PDE. Specifically: (1) substituting baseline support load `M_supp` and support direction `R_phi` into a Stage-24 support theorem produces a quadratic in `xi` with explicit coefficients `B_cont`, `C_cont`; (2) the physical root of that quadratic vanishes when both loads vanish; (3) a normalization function `F_cont(xi)` arises from inserting `q = sqrt(lambda0) R_U`, `r = sqrt(lambda0) R_phi` into a Stage-25 normalization law; (4) `R_phi = 1` ("sigma0 = 0 minimal-kernel limit" per docstring) gives a source-tied closure; (5) `R_phi = R_U` collapses the rank-2 branch to a total-loading tracking law; (6) the mismatch penalty is `+lambda0 M_mix M_supp (R_U - R_phi)^2`. Both engines run, both report PASS.

## Assertion inventory

| #  | Script      | Line  | Form                                                                     | Anchored to claim? |
|----|-------------|-------|--------------------------------------------------------------------------|--------------------|
| A1 | sympy       | 73    | `expect_zero("quadratic branch equation", branch_eq - 9*branch_expected)` | yes                |
| A2 | sympy       | 79    | `expect_zero("zero-load root", xi_phys.subs({Mmix:0,Msupp:0}))`           | yes                |
| A3 | sympy       | 118   | `expect_zero("source-tied n", n_source - n_source_expected)`              | partial            |
| A4 | sympy       | 119   | `expect_zero("source-tied F", F_source - F_source_expected)`              | partial            |
| A5 | sympy       | 125   | `expect_zero("tracking collapse of n_req", n_track - (G_q - Mmix))`       | yes                |
| A6 | sympy       | 129   | `expect_zero("tracking total-loading law", ... )`                         | NO (redundant w/ A5) |
| A7 | sympy       | 137   | `expect_zero("tracking F collapse", F_track - F_track_expected)`          | partial            |
| A8 | sympy       | 144   | `expect_zero("quadratic mismatch penalty", C_rewritten - C_expected)`     | NO (tautological)  |
| B1 | mathematica | 56    | `expectZero["quadratic branch equation", branchEq - 9 branchExpected]`    | yes                |
| B2 | mathematica | 61    | `expectZero["zero-load root", xiPhys /. {mMix->0, mSupp->0}]`             | yes                |
| B3 | mathematica | 96    | `expectZero["source-tied n", nSource - nSourceExpected]`                  | partial            |
| B4 | mathematica | 97    | `expectZero["source-tied F", fSource - fSourceExpected]`                  | partial            |
| B5 | mathematica | 111   | `expectZero["tracking collapse of n_req", nTrack - (gQ - mMix)]`          | yes                |
| B6 | mathematica | 112   | `expectZero["tracking total-loading law", trackingEquation]`              | NO (redundant w/ B5) |
| B7 | mathematica | 113   | `expectZero["tracking F collapse", fTrack - fTrackExpected]`              | partial            |
| B8 | mathematica | 120   | `expectZero["quadratic mismatch penalty", cRewritten - cExpected]`        | NO (tautological)  |

## Findings

### F1 — mathematica_transliteration

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:42-120`

**What's wrong:**
The `.wl` script is a structural line-by-line port of the `.py` script. Every conceptual step in SymPy has a direct one-to-one Mathematica counterpart with identical variable choreography and identical intermediate quantities:

- SymPy lines 60-62 define `n_req` as `(xi*(delta+xi) - Mmix*(delta + (1 + lambda0*RU^2)*xi)) / (delta + (1 + lambda0*Rphi^2)*xi - Mmix*lambda0*(RU - Rphi)^2)`. WL lines 45-47 define `nReq` with exactly the same algebraic form.
- SymPy line 66 forms `branch_eq = numerator of (n_req - Msupp)`. WL line 49 forms `branchEq = Numerator[Together[nReq - mSupp]]` — identical operation, identical name.
- SymPy lines 70-73 define `B_cont`, `C_cont`, `branch_expected = xi^2 + B*xi + C` and check `branch_eq - 9*branch_expected == 0`. WL lines 50-56 do the same construction with `bCont`, `cCont`, `branchExpected` and check `branchEq - 9*branchExpected`.
- SymPy lines 75-79 compute `xi_phys = (-B + sqrt(disc))/2` and substitute Mmix=Msupp=0. WL lines 58-61 do the same.
- SymPy lines 84-92 define `D_cont`, `F_cont` from the same algebraic templates; WL lines 65-74 mirror this. Same substitution sequence Rphi=1 and Rphi=RU follows in both.

This violates the second-engine independence policy. The Mathematica script does not re-derive any result from the underlying premises by an algebraically distinct route (e.g., `Solve` on the quadratic, `Reduce` on the support equation, alternate parametrization). It echoes the SymPy algebra and confirms FullSimplify of the same expressions equals zero.

**Why this matters:**
If a coefficient sign or factor error were transcribed from a common source into the SymPy script, copying the same construction into Mathematica would propagate the same error and both engines would "agree." The point of two engines is to catch such errors.

**Required change:**
Restructure the `.wl` script so at least one of `xi_phys`, `F_cont`, or the quadratic coefficients is obtained by an algebraically distinct route. Concretely:
- Replace WL lines 49-59 with a derivation that uses `Solve[branchEq == 0, xi]` (Mathematica's algebraic solver) to obtain the two roots, then select the physical root by checking which root vanishes when `mMix -> 0, mSupp -> 0`, and verify that this `Solve`-derived `xiPhys` matches the manually-formed `(-bCont + Sqrt[deltaDisc])/2`. This is an independent algebraic path; if the manual coefficients are wrong, the `Solve` result will not match.

**Verification:**
The verifier should see a new `Solve[branchEq == 0, xi]` (or `Reduce`) call in the .wl, plus a new `expectZero["xiPhys solve match", xiPhysFromSolve - xiPhys]` assertion. Output should print both forms and the residual.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:139-144`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:115-120`

**What's wrong:**
The "quadratic mismatch penalty" assertion is purely a variable-renaming check.

`C_cont` is defined at SymPy line 71 as
```
C_cont = -delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*(RU - Rphi)^2
```
Then at lines 141-144:
```
mismatch = sp.symbols("Delta_R", real=True)
C_rewritten = C_cont.subs(Rphi, RU - mismatch)
C_expected = -delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*mismatch^2
expect_zero("quadratic mismatch penalty", C_rewritten - C_expected)
```
Substituting `Rphi -> RU - mismatch` into `(RU - Rphi)^2` gives `(RU - (RU - mismatch))^2 = mismatch^2` by elementary algebra. The "expected" form is just the same expression with `(RU-Rphi)^2` replaced by `mismatch^2`. This assertion cannot fail by construction — it is algebraically guaranteed by the SymPy `subs`. It tests no physics; it tests SymPy's substitution rule.

The Mathematica copy at WL line 117-120 does the same thing (`cRewritten = cCont /. rPhi -> rU - deltaR`) and is equally tautological.

**Why this matters:**
Docstring claim 6 — that the rank-2 mismatch penalty is `+lambda0 M_mix M_supp (R_U - R_phi)^2` — is presented as a result, but the script's check only verifies that "when you rename `(R_U - R_phi)` to `Delta_R`, you get `Delta_R^2`." It does not anchor the form of `C_cont` to any independent derivation, so the claim's substantive content is unverified.

**Required change:**
Replace the trivial rename with a substantive check. The substantive content of claim 6 is that the `C` coefficient of the branch quadratic (the constant term in `xi^2 + B*xi + C = 0`) has the form `-delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*(RU-Rphi)^2`. Anchor this by extracting `C` from the already-derived `branch_eq` (the numerator of `n_req - Msupp`, line 66) and comparing.

Concrete edit at sympy lines 139-144:
```python
subbanner("27.5 — Exact mismatch penalty")

# Extract the constant-in-xi part of branch_eq (which equals 9*(xi^2 + B*xi + C)).
C_from_branch_eq = sp.simplify(sp.Poly(branch_eq, xi).nth(0) / 9)
C_expected = sp.simplify(-delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*(RU-Rphi)**2)
expect_zero("mismatch penalty in C coefficient", C_from_branch_eq - C_expected)
```
Make the analogous change in the .wl at lines 115-120 using `CoefficientList[branchEq, xi]` to pull off the constant term and divide by 9, then compare.

**Verification:**
The verifier should see (a) the substitution `Rphi -> RU - deltaR` removed, (b) a new line extracting the constant-in-xi coefficient from `branch_eq` (or `branchEq`), (c) a comparison against the manually-written `(RU - Rphi)^2` form. The assertion can now fail if `C_cont` was defined wrong — it actually tests claim 6.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:127-129`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:103, 112`

**What's wrong:**
The "tracking total-loading law" assertion (A6/B6) is algebraically guaranteed once A5/B5 hold.

At sympy line 125 (A5) the script asserts `n_track - (G_q - Mmix) == 0`, i.e. `n_track == G_q - Mmix`. Then at line 128-129:
```python
tracking_equation = (n_track - Msupp) - (G_q - (Mmix + Msupp))
                  = n_track - Msupp - G_q + Mmix + Msupp
                  = n_track - G_q + Mmix
                  = (n_track - (G_q - Mmix))
```
This is bit-for-bit the same expression as A5 — `Msupp` cancels in `(n_track - Msupp) - (G_q - (Mmix + Msupp))`. There is no new content. WL lines 103, 112 (B6) are identical.

The inline comment at sympy line 127 ("Setting actual support baseline equal to the exact required support load yields total-loading law") suggests a substantive intent, but the algebra performed does not enact any such constraint — `Msupp` appears symbolically in both factors with opposite signs and cancels.

**Why this matters:**
Docstring claim 5 promises that the rank-2 branch collapses to a "one-direction tracking law with total loading M_mix + M_supp." A5/B5 actually verify the collapse (`n_track = G_q - Mmix`), but the total-loading content (`M_mix + M_supp`) is asserted only by a tautology. The substantive claim is unverified.

**Required change:**
Replace the redundant `tracking_equation` check with one that genuinely uses both A5 and the support equation `n_req = M_supp`. The substantive claim is: when `R_phi = R_U` and we impose `n_req = M_supp`, the resulting equation in `xi` is `G_q(xi) = M_mix + M_supp`, i.e. the single-direction Stage-24 law with the loads added.

Concrete edit at sympy lines 127-129:
```python
# When R_phi = R_U, the support equation n_req = M_supp combined with
# n_track = G_q - Mmix gives G_q = Mmix + Msupp (total-loading law).
tracking_total = sp.simplify(G_q.subs(xi, xi) - (Mmix + Msupp))
# Anchor it via the support-equation residual at R_phi=R_U:
support_residual_track = sp.simplify((n_req.subs(Rphi, RU) - Msupp) - (G_q - (Mmix + Msupp)))
expect_zero("tracking total-loading law", support_residual_track)
```
Wait — re-examination shows `support_residual_track` algebraically equals the existing `tracking_equation` and is still tautological once `n_track = G_q - Mmix` is established. A genuinely non-tautological replacement is to verify the structural prediction that `G_q` is the SAME function in the rank-1 (single-direction) and rank-2 tracking cases, i.e. that the rank-2 surface with `Rphi=Ru` recovers the Stage-24 single-direction support theorem with `Mmix -> Mmix + Msupp`. To do that without reading other units, the script must first define the rank-1 single-direction form independently. Since this is not currently in-script, the safer correction is to DELETE the redundant assertion (and its WL counterpart) rather than dress it up.

Final required change: at sympy lines 127-129 and WL lines 103, 112, delete the redundant `tracking_equation` / `trackingEquation` and its assertion. Keep A5/B5 (the substantive collapse) and add a one-line comment noting that the total-loading interpretation follows by inspection: substituting `n_track = G_q - Mmix` into `n_track = Msupp` (the support condition) gives `G_q = Mmix + Msupp`.

**Verification:**
The verifier should see the `tracking_equation` block and the corresponding WL `trackingEquation` line removed (or replaced by a non-tautological check that actually uses the support condition `n_req = M_supp`). The PASS count drops by one in each engine, both still exit 0.

### F4 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:81-96`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:63-76`

**What's wrong:**
Docstring claim 3 states that `F_cont` arises from "inserting `q = t R_U` and `r = t R_phi` into the Stage-25 normalization law." However, in this script `F_cont` is *defined* directly at lines 88-92 (sympy) / 69-74 (wl):

```python
F_cont = (delta + (1 + lambda0*RU*Rphi)*xi)^2
       * (delta + (1 + lambda0*Rphi)*xi - Mmix*lambda0*(RU-Rphi)*(RU-1))^2
       / ((1 - xi) * D_cont^2)
```
No Stage-25 expression appears in this script, and no derivation is performed. The block ends with `print("F_cont built successfully.")` — there is no assertion at all between lines 88 and 96.

The general form of `F_cont` is then exercised only at two slices: `Rphi = 1` (A4/B4, "source-tied") and `Rphi = R_U` (A7/B7, "tracking"). Two slice checks against manually-written expected expressions are consistent but cannot constrain the full bivariate dependence on `(Rphi, Mmix, Msupp, RU)`. Many alternative forms for `F_cont` would pass both slice checks.

**Why this matters:**
A reader of the script cannot conclude from the assertions alone that the displayed `F_cont` is correct as written; the assertions only verify it at two collapsed surfaces. If a coefficient in the cross-term `lambda0*RU*Rphi` or `lambda0*(RU-Rphi)*(RU-1)` were mistyped, both slice checks could still pass (e.g., if the error vanishes at `Rphi=1` and `Rphi=RU` simultaneously, which is possible for certain coefficient mistakes).

**Required change:**
Add a third slice check at a parameter value that is algebraically distinct from `Rphi=1` and `Rphi=R_U`. The simplest choice is `Rphi = 2` (a literal value that is neither `1` nor `RU`), with the expected form written out independently. Concrete addition between sympy lines 96 and 98:

```python
# Third slice: independent literal R_phi to constrain bivariate dependence.
Rphi_lit = sp.Integer(2)
F_lit = sp.simplify(F_cont.subs(Rphi, Rphi_lit))
F_lit_expected = sp.simplify(
    (delta + (1 + lambda0 * RU * Rphi_lit) * xi)**2
    * (delta + (1 + lambda0 * Rphi_lit) * xi - Mmix * lambda0 * (RU - Rphi_lit) * (RU - 1))**2
    / ((1 - xi) * ((delta + xi - Mmix * lambda0 * RU * (RU - Rphi_lit))**2
       + lambda0 * (Mmix * (RU - Rphi_lit) + Rphi_lit * xi)**2)**2)
)
expect_zero("third-slice F at Rphi=2", F_lit - F_lit_expected)
```
And the analogous block in the .wl with `rPhi -> 2`.

Note: this is a sanity check that the slice substitution is internally consistent, not an independent derivation. But because `F_lit_expected` is constructed by an independent textual transcription of the same generic template at a *literal* value of `Rphi`, any typo that mis-distributes `Rphi` (e.g., dropping a `Rphi` in a sum or replacing `RU*Rphi` with `RU+Rphi`) will produce a nonzero residual at the literal value.

**Verification:**
The verifier should see new `Rphi_lit` / `rPhi -> 2` blocks in both engines, each ending in an `expect_zero` / `expectZero` of a third-slice residual, with the printed residual equal to 0 in the saved output.

### F5 — symbol_assumption_error

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:50`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl:35-38`

**What's wrong:**
SymPy line 50 declares `sigma0 = sp.symbols("sigma_0", real=True)`. WL lines 35-38 include `sigma0` in `Clear[...]` and the `$Assumptions`. The symbol `sigma0` never appears anywhere else in either script.

Docstring claim 4 in sympy lines 13-14 says the minimal-kernel "sigma0=0 limit" gives the source-tied closure. However, the actual check (sympy line 102-119, WL line 80-97) substitutes `Rphi -> 1`, not `sigma0 -> 0`. The script does not derive or assert the equivalence `sigma0 = 0 <=> Rphi = 1`, so claim 4 as stated is unverified by the assertions; only an alternative formulation (`Rphi = 1`) is verified.

**Why this matters:**
A declared-but-unused symbol indicates one of two things: (a) leftover scaffolding from an earlier draft that should be removed, or (b) a verification step that was intended but not implemented. Leaving `sigma0` in place mis-suggests that the script tests a `sigma0 -> 0` limit when it does not.

**Required change:**
Either remove the declaration of `sigma0` if no `sigma0`-based check is intended, or add a real assertion that uses `sigma0`. Given that adding a new physical relation would expand scope (forbidden), the conservative fix is to remove the dead declaration:

- Delete sympy line 50: `sigma0 = sp.symbols("sigma_0", real=True)`.
- In the .wl, remove `sigma0` from `Clear[xi, delta, mMix, mSupp, rU, rPhi, sigma0, deltaR]` at line 35 (yielding `Clear[xi, delta, mMix, mSupp, rU, rPhi, deltaR]`) and from `Element[{xi, delta, mMix, mSupp, rU, rPhi, sigma0, deltaR}, Reals]` at line 37 (yielding `Element[{xi, delta, mMix, mSupp, rU, rPhi, deltaR}, Reals]`).
- Edit the sympy docstring line 13 to read "The minimal-kernel limit R_phi=1 gives the exact source-tied closure." (replacing `sigma0=0` with the variable actually substituted).

**Verification:**
The verifier confirms `sigma0` / `sigma_0` no longer appears anywhere in either script (`grep -i 'sigma' moving_throat_pde_stage044*` returns nothing other than possibly the docstring update).

## Independent-derivation check (Mathematica)

Strict transliteration. Side-by-side correspondence:

- SymPy 57-62 `n_req` <-> WL 44-48 `nReq` — identical algebraic form.
- SymPy 66 `branch_eq = ... .as_numer_denom()[0]` <-> WL 49 `branchEq = Numerator[Together[nReq - mSupp]]` — same operation.
- SymPy 70-72 `B_cont, C_cont, branch_expected` <-> WL 50-52 `bCont, cCont, branchExpected` — same coefficients, same template.
- SymPy 73 `expect_zero("quadratic branch equation", branch_eq - 9*branch_expected)` <-> WL 56 `expectZero["quadratic branch equation", branchEq - 9 branchExpected]` — same assertion with same factor of 9.
- SymPy 75-79 quadratic root construction <-> WL 58-61 — same `(-B + sqrt(disc))/2` and same Mmix=Msupp=0 substitution.
- SymPy 84-92 `D_cont, F_cont` <-> WL 65-74 `dCont, fCont` — same algebraic templates.
- SymPy 100-119 source-tied (Rphi=1) <-> WL 80-97 — same substitution sequence and same expected forms.
- SymPy 123-137 tracking (Rphi=RU) <-> WL 101-113 — same substitutions and same expected forms.
- SymPy 141-144 mismatch (Rphi -> RU - mismatch) <-> WL 117-120 — same substitution and same expected form.

There is no independent algebraic route in the .wl. This is the basis of F1.

## Engine cross-check

Both engines exit 0 with all assertions passing. Printed `n_req^(cont)`, `xi_phys`, `n_source`, `F_source`, and `F_cont(xi)` agree in algebraic form:

- SymPy `n_req^(cont)` (output line 17-22): `(2*Mmix*RU^2*xi + 9*Mmix*delta + 9*Mmix*xi - 9*delta*xi - 9*xi^2) / (2*Mmix*RU^2 - 4*Mmix*RU*Rphi + 2*Mmix*Rphi^2 - 2*Rphi^2*xi - 9*delta - 9*xi)`.
- Mathematica `n_req^(cont)` (output line 17): `(xi*(delta + xi) - mMix*(delta + xi + (2*rU^2*xi)/9))/(delta - (2*mMix*(rPhi - rU)^2)/9 + xi + (2*rPhi^2*xi)/9)`. Multiplying numerator and denominator by 9 reproduces the SymPy form modulo overall sign convention (SymPy carries an extra factor of `-1/-1` from its `factor`-style rearrangement). Algebraically equivalent.

`xi_phys`, `n_source`, `F_source`, `F_cont` outputs likewise agree at the level of algebraic equality after clearing the factor-of-9 conventions.

Engine cross-check passes. However, because the .wl mirrors the .py construction step-for-step (F1), agreement here is weaker evidence than two independent derivations would provide.

## Verdict justification

Both engines run and both report PASS on a self-consistent set of algebraic identities, but the audit reveals five real issues:

1. The Mathematica script is a transliteration of the SymPy script (F1, high). Both engines verifying the same construction agree by construction; a typo in a shared template would not be caught.
2. The "quadratic mismatch penalty" assertion (F2, high) is a pure variable rename and does not anchor docstring claim 6.
3. The "tracking total-loading law" assertion (F3, medium) is algebraically redundant with the preceding "tracking collapse" assertion.
4. `F_cont`'s general form (F4, medium) is exercised only at two collapsed slices; a third independent slice is needed to constrain the bivariate dependence.
5. `sigma0` is declared in both engines but never used (F5, low); the docstring promises a `sigma0=0` check that is never performed.

None of these findings rises to `UNFIXABLE` (the underlying math is internally consistent) or `CRITICAL_DOWNSTREAM` (no sign flip, no constant change). Verdict: `findings`.

Attacks tried that the script survived:
- Zero-load substitution in `xi_phys`: passes correctly given `delta > 0`.
- Substituting `Rphi = R_U` actually collapses the rank-2 numerator and denominator as advertised (A5 is substantive).
- The `quadratic branch equation` check (A1/B1) genuinely verifies that the numerator of `n_req - M_supp` factors as `9*(xi^2 + B_cont*xi + C_cont)` with the manually-written `B_cont`, `C_cont` — independent reconstruction.
- The two source-tied (Rphi=1) and tracking (Rphi=RU) F-slice expressions are independently transcribed, so a typo in either branch would surface.

## Self-test notes

- **Variable independence:** No `sp.diff` / `D[...]` calls appear in this audit unit; the script is purely algebraic identity checking. No derivative-trap risk.
- **Symmetry/parity:** No integrals over unbounded domains; no parity-based vanishing claims. N/A.
- **Trivial-case pre-check for F4's proposed addition:** Substituting `Rphi -> 2` into the SymPy `F_cont` formula (lines 88-92) yields `(delta + (1 + lambda0*RU*2)*xi)^2 * (delta + (1 + lambda0*2)*xi - Mmix*lambda0*(RU-2)*(RU-1))^2 / ((1-xi) * ((delta + xi - Mmix*lambda0*RU*(RU-2))^2 + lambda0*(Mmix*(RU-2) + 2*xi)^2)^2)`. The proposed `F_lit_expected` in F4 uses exactly this template with `Rphi_lit = 2`, so the residual `F_lit - F_lit_expected` is identically 0 by construction — the assertion passes trivially. This is the correct outcome for a sanity-check on the substitution; it would fail only if `F_cont`'s definition were edited to introduce an `Rphi`-dependent typo. Confirmed safe.
- **Trivial-case pre-check for F2's proposed replacement:** `branch_eq` equals `9*(xi^2 + B_cont*xi + C_cont)` (this is asserted by A1). Extracting the constant-in-xi coefficient gives `9*C_cont`; dividing by 9 yields `C_cont = -delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*(RU-Rphi)^2`. Comparing this to `C_expected = -delta*(Mmix+Msupp) + lambda0*Mmix*Msupp*(RU-Rphi)^2` gives zero. Passes correctly. This check is non-tautological because `branch_eq` was constructed by an algebraically distinct path (numerator of the support equation) rather than by writing `C_cont` directly.
- **Path specifications:** Both scripts are present at canonical paths; no `missing_verification_script` finding, so no path-spec issue for new files. F1 modifies the existing `.wl`; path verified at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage044_continuum_selected_rank2_mathematica_audit.wl`.
