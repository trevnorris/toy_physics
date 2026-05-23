---
unit_id: 063
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 063 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage063_parent_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.txt`

## What the script claims to verify

According to the docstring and inline comments, the script verifies five claims about parent-overlap fail/succeed thresholds derived from a microscopic gain `G_micro = rho_star * g_phi^2 * O_sp^2 / (m * cs_star_sq * K_X * N_ss)`: (1) the `g_phi^2` thresholds `g_(phi,fail/suff)^2 = m cs* K_X N_ss G_(fail/suff)/(rho* O_sp^2)` solve `G_micro = G_(fail/suff)`; (2) substituting `O_sp^2 = C^2 N_ss N_pp` yields equivalent coherence-form thresholds; (3) the coherence ratio `C_suff^2/C_fail^2` reduces to `G_suff/G_fail`; (4) the best-case (perfect-alignment) gain is `G_max = rho* g_phi^2 N_pp/(m cs* K_X)` so that `G_micro = G_max C^2`; (5) inserting the Stage-44 closures `G_(fail) = Pe_req/(kappa Delta_inf)`, `G_(suff) = Pe_req/(kappa Delta_0)`, and `kappa = K_X L^2/T_X` cancels `K_X` from the threshold prefactor, yielding the explicit forms in terms of `T_X, Pe_req, Delta_0, Delta_inf`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|--------|------|------|--------------------|
| A1 | sympy | 51-55 | `g_fail_sq.subs(Osp^2, C2*Nss*Npp) - m cs* KX G_fail/(rho* Npp C2) == 0` | no (pure substitution+cancellation; both sides built from same primitive) |
| A2 | sympy | 56-60 | `g_suff_sq.subs(Osp^2, C2*Nss*Npp) - m cs* KX G_suff/(rho* Npp C2) == 0` | no (mirror of A1) |
| A3 | sympy | 67 | `C_suff_sq/C_fail_sq - G_suff/G_fail == 0` | no (C_fail_sq and C_suff_sq differ only by G_fail vs G_suff; ratio is identically G_suff/G_fail) |
| A4 | sympy | 70-74 | `G_micro.subs(Osp^2, C2*Nss*Npp) - G_max*C2 == 0` | no (G_max defined to make this zero by construction) |
| A5 | sympy | 80-84 | `g_fail_sq.subs(G_fail->Pe/(kappa Delta_inf), kappa->KX L^2/TX) - m cs* TX Nss Pe/(rho* Osp^2 L^2 Delta_inf) == 0` | no (K_X explicitly cancels by simple algebra; RHS is hand-rearranged copy of LHS) |
| A6 | sympy | 85-89 | mirror of A5 with G_suff, Delta_0 | no (mirror of A5) |
| A7 | sympy | 91-95 | A5 with O_sp^2 -> C^2 Nss Npp substituted | no (same identity, additional substitution) |
| A8 | sympy | 96-100 | A6 with O_sp^2 -> C^2 Nss Npp substituted | no (same identity, additional substitution) |
| B1 | math | 42-45 | port of A1 | no |
| B2 | math | 46-49 | port of A2 | no |
| B3 | math | 55 | port of A3 (`cSuffSq/cFailSq - gSuff/gFail`) | no |
| B4 | math | 58 | port of A4 | no |
| B5 | math | 64-68 | port of A5 | no |
| B6 | math | 69-73 | port of A6 | no |
| B7 | math | 74-78 | port of A7 | no |
| B8 | math | 79-83 | port of A8 | no |

Every assertion in both engines is algebraically guaranteed by the script's own definitions. Walk-through for A3: `C_fail_sq = m*cs_star_sq*KX*G_fail/(rho_star*g_phi^2*N_pp)` (sympy line 63) and `C_suff_sq = m*cs_star_sq*KX*G_suff/(rho_star*g_phi^2*N_pp)` (line 64) differ only in the G factor, so `C_suff_sq/C_fail_sq = G_suff/G_fail` regardless of any physical input. Walk-through for A4: `G_max = rho_star*g_phi^2*N_pp/(m*cs_star_sq*K_X)` (line 70) is literally `(rho_star*g_phi^2*Osp^2)/(m*cs_star_sq*K_X*N_ss)` with `Osp^2 -> C2*Nss*Npp` and `C2` factored out — the assertion just unfactors what was factored. Walk-through for A5: after substituting `kappa -> KX*L^2/TX`, the explicit `K_X` in the `g_fail_sq` numerator cancels the `K_X` in the substituted `1/kappa`; the RHS on line 83 is identical to the LHS after that cancellation. None of these checks can fail unless SymPy's simplifier is broken.

## Findings

### F1 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:51-100`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl:42-83`

**What's wrong:**
Every assertion in both scripts is an algebraic identity guaranteed by the script's own definitions, not a verification of any threshold theorem. Concretely:

(a) `g_fail_sq` is **defined** on sympy line 45 as `m*cs_star_sq*KX*Nss*G_fail/(rho_star*Osp**2)`. The script never derives this from solving `G_micro = G_fail`; it simply writes it down. The coherence substitutions on lines 51-60 then check that `g_fail_sq.subs(Osp^2, C2*Nss*Npp)` equals the very same expression with `Osp^2` replaced and `N_ss` cancelled — i.e. `m*cs_star_sq*KX*Nss*G_fail/(rho_star*C2*Nss*Npp) - m*cs_star_sq*KX*G_fail/(rho_star*Npp*C2)`, which is `0` by elementary cancellation independent of any physics.

(b) `C_fail_sq` and `C_suff_sq` (lines 63-64) are written down explicitly differing only in `G_fail` vs `G_suff`; the ratio check on line 67 is then `G_suff/G_fail - G_suff/G_fail == 0` after the common factor cancels.

(c) `G_max` (line 70) is the rearrangement of `G_micro|_{Osp^2 -> C^2 N_ss N_pp}` after factoring out `C^2`; the check on lines 71-74 just multiplies back the factor that was pulled out.

(d) For lines 80-100, substituting `kappa = K_X L^2/T_X` into `G_fail = Pe_req/(kappa Delta_inf)` gives `T_X Pe_req/(K_X L^2 Delta_inf)`, so `g_fail_sq` becomes `m cs* K_X N_ss * T_X Pe_req / (K_X L^2 Delta_inf rho* Osp^2)`; `K_X` cancels trivially, leaving exactly the RHS on line 83. The script does not verify any derivation step (e.g. that `G_micro >= G_(fail/suff)` is the actual fail/succeed criterion, or that `kappa = K_X L^2/T_X` is the correct identification with Stage-44 closures, or that `G_max` is the supremum over alignments); it just substitutes and rearranges, then asserts the rearrangement.

The Mathematica file mirrors this structure assertion-for-assertion (lines 42-83), so the same critique applies.

**Why this matters:**
The docstring announces five thresholds-theorem claims, but none of them are exercised. If the underlying physics derivation were wrong (e.g. a missing factor of 2 in `G_micro`, or `Pe_req` evaluated at the wrong reference state, or `G_max` defined with `N_ss` instead of `N_pp`), every assertion in both scripts would still PASS because the assertions only test that SymPy/Mathematica can rearrange expressions. The script provides no evidence for or against the threshold theorem; it is a SymPy/Mathematica rearrangement self-test.

**Required change:**
Add at least one substantive check in each script that exercises the underlying derivation rather than restating it. Two are recommended; either is sufficient:

1. **Threshold-from-G_micro derivation.** Define `G_micro` symbolically (as already done on sympy line 41 / mathematica line 35), set up the equation `G_micro == G_fail` (with `g_phi` unknown), and use `sympy.solve(G_micro - G_fail, g_phi**2)` (`Solve[gMicro == gFail, gPhi^2]` in Mathematica) to *derive* `g_phi^2_fail`. Then assert that the solver's result equals the manually written `g_fail_sq` from line 45. This anchors the threshold formula to the gain definition; if `G_micro` is ever changed, the solver result will differ from the hand-rearranged form and the check will fail.

2. **G_max as the alignment supremum.** Since the script comment on line 69 says "Best-case gain at perfect alignment", add a check that under the Cauchy bound `O_sp^2 <= N_ss N_pp` (i.e. `C^2 <= 1`), `G_micro` attains `G_max` exactly at `C^2 = 1`. This can be done in sympy via `expect_zero("G_max - G_micro at C^2=1", G_max - G_micro.subs(Osp**2, Nss*Npp))` (perfect alignment ⇒ C=1 ⇒ Osp^2 = Nss Npp). If `G_max` were defined with a wrong factor (e.g. `N_ss` instead of `N_pp`), this check would fail.

Each script must add both checks (or at minimum check 1, which is the deeper derivation). The exact edits are spelled out in the directive.

**Verification:**
- `scripts/.../stage063...sympy_audit.py`: after the existing checks, two new `expect_zero` calls appear ("g_fail^2 from solve(G_micro=G_fail)" and "G_max == G_micro|_{C^2=1}"), and the output file gains corresponding lines.
- `mathematica/.../stage063...mathematica_audit.wl`: two new `expectZero` blocks with the same content.
- Both scripts continue to exit 0.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl:35-83`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:41-100`

**What's wrong:**
The `.wl` file is a line-by-line port of the `.py` file. Variable-by-variable, assertion-by-assertion:

| sympy line | mathematica line | identical content |
|---|---|---|
| 41 `G_micro = rho_star*g_phi^2*Osp^2/(m*cs_star_sq*KX*Nss)` | 35 `gMicro = rhoStar*gPhi^2*oSP^2/(m*csStarSq*kX*nSS)` | yes |
| 45 `g_fail_sq = m*cs_star_sq*KX*Nss*G_fail/(rho_star*Osp^2)` | 36 `gFailSq = m*csStarSq*kX*nSS*gFail/(rhoStar*oSP^2)` | yes |
| 46 `g_suff_sq = m*cs_star_sq*KX*Nss*G_suff/(rho_star*Osp^2)` | 37 `gSuffSq = m*csStarSq*kX*nSS*gSuff/(rhoStar*oSP^2)` | yes |
| 51-55 coherence substitution g_fail^2 | 42-45 same | yes |
| 56-60 coherence substitution g_suff^2 | 46-49 same | yes |
| 63 `C_fail_sq = m*cs_star_sq*KX*G_fail/(rho_star*g_phi^2*Npp)` | 51 `cFailSq = m*csStarSq*kX*gFail/(rhoStar*gPhi^2*nPP)` | yes |
| 64 `C_suff_sq` | 52 `cSuffSq` | yes |
| 67 `C_suff^2/C_fail^2 - G_suff/G_fail` | 55 same | yes |
| 70 `G_max = rho_star*g_phi^2*Npp/(m*cs_star_sq*KX)` | 57 `gMax = rhoStar*gPhi^2*nPP/(m*csStarSq*kX)` | yes |
| 71-74 `G_micro.subs - G_max*C^2` | 58 same | yes |
| 76-78 `G_fail_sub`, `G_suff_sub` | 60-61 `gFailSub`, `gSuffSub` | yes |
| 80-100 four threshold-with-kappa checks | 64-83 four same checks | yes |

Every intermediate object (`gMicro`, `gFailSq`, `gSuffSq`, `cFailSq`, `cSuffSq`, `gMax`, `gFailSub`, `gSuffSub`) and every RHS target in the eight `expectZero` calls is a camelCased rename of the corresponding sympy expression. The Mathematica engine performs no independent derivation — it just renames the variables and prints the same identities. This violates the second-engine policy: if the sympy script has an error in its setup (e.g. wrong placement of `N_ss` in `g_fail_sq`), the Mathematica script reproduces exactly the same wrong setup and confirms it.

**Why this matters:**
A second engine exists to catch errors the first engine cannot. A literal transliteration provides no second-engine check; it just runs the same algebra in a different syntax. The substantive bugs that the policy is designed to catch (wrong factor, wrong substitution, wrong simplification rule) propagate identically through the port.

**Required change:**
After F1's substantive checks land, the Mathematica script's new blocks should derive the same conclusion through a path the sympy script does not use. Specifically:

1. For the "g_fail^2 from solve" check, use `Solve[gMicro == gFail, gPhi^2]` in Mathematica (vs `sp.solve(G_micro - G_fail, g_phi**2)` in sympy) — but additionally, in Mathematica, perform an independent algebraic check via `Reduce[gMicro == gFail && gPhi > 0, gPhi]` and assert that the resulting `gPhi^2` matches `gFailSq`. The `Reduce` path solves the implicit equation rather than just re-arranging.
2. For the "G_max at C^2=1" check, in Mathematica derive `gMax` via `Maximize[{gMicro /. oSP^2 -> c2*nSS*nPP, 0 < c2 <= 1, nSS > 0, nPP > 0, rhoStar > 0, gPhi > 0, m > 0, csStarSq > 0, kX > 0}, c2]` (the assumption-bounded Cauchy maximum), and assert the Maximize output's first component equals `gMax`. The sympy version uses the algebraic substitution `Osp^2 -> Nss*Npp` directly.

Both checks then derive the same physical result via two genuinely different solver paths, satisfying the independent-derivation requirement.

**Verification:**
The Mathematica script gains at least one `Solve`/`Reduce`/`Maximize` invocation that does not appear in the sympy script, and both engines' output transcripts show the same final symbolic equality through different intermediate steps.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:33-39`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage063_parent_thresholds_mathematica_audit.wl:28-33`

**What's wrong:**
All symbols including `O_sp` and `C2` are declared `positive=True, real=True` (sympy lines 33-39) / `> 0` (Mathematica lines 30-33). But the Cauchy-Schwarz bound that motivates `G_max` (script comment line 69 says "Best-case gain at perfect alignment") requires `O_sp^2 <= N_ss N_pp` (i.e. `C^2 <= 1`). The script never asserts this constraint, and the symbol `C_sp_sq` (the squared coherence) is allowed to be arbitrarily large positive. The "best-case" interpretation of `G_max` is therefore not anchored: with `C^2` unbounded, `G_micro = G_max C^2` can exceed `G_max`, contradicting the comment.

Because the assumption `C^2 <= 1` is missing, even a substantive maximum-attainment check would have to be added by hand. The current script's positive-only assumption on `C_sp_sq` is consistent with positivity but does not implement the physical bound; this leaves the "best case" claim unverified.

**Why this matters:**
The docstring's claim (3) "Cauchy upper-gain factorization" rests on `C^2 <= 1`. Without that bound active in `$Assumptions` / sympy assumptions, no simplification can use it, and no check can fail when `G_max` is replaced by a larger value. Together with F1 this means the "best-case" semantics is not exercised at all.

**Required change:**
The new "G_max as alignment supremum" check from F1 must be performed under the Cauchy constraint explicitly:
- sympy: in the new `expect_zero` call, evaluate at `Osp**2 = Nss*Npp` (i.e. `C^2 = 1`) — this is the saturating case. No assumption change to existing symbols is required; the substitution itself implements the Cauchy saturation. Add a comment that `C^2 <= 1` is the Cauchy-Schwarz bound and `C^2 = 1` is perfect alignment.
- Mathematica: same — perform the check at the saturating substitution, with a comment naming `C^2 <= 1` as the Cauchy bound.

No change to global `$Assumptions` is needed (and would be risky); the saturating substitution is sufficient.

**Verification:**
The new `expect_zero` / `expectZero` call from F1 includes a comment naming the Cauchy bound, and the saturating substitution `Osp^2 -> Nss*Npp` (sympy) / `oSP^2 -> nSS*nPP` (Mathematica) is used. Output gains a line like `G_max - G_micro at C^2=1 = 0`.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration. Three corresponding-section quotes from the script files:

1. Setup:
   - sympy line 41: `G_micro = sp.simplify(rho_star * g_phi**2 * Osp**2 / (m := sp.symbols('m', positive=True, real=True)) / cs_star_sq / KX / Nss)`
   - mathematica line 35: `gMicro = FullSimplify[rhoStar*gPhi^2*oSP^2/(m*csStarSq*kX*nSS), Assumptions -> $Assumptions];`
   
   Same expression, same factor order, same variable identity.

2. Ratio check:
   - sympy line 67: `expect_zero("C_suff^2/C_fail^2 - G_suff/G_fail", sp.simplify(C_suff_sq / C_fail_sq - G_suff / G_fail))`
   - mathematica line 55: `expectZero["C_suff^2/C_fail^2 - G_suff/G_fail", FullSimplify[cSuffSq/cFailSq, Assumptions -> $Assumptions] - gSuff/gFail];`
   
   Identical assertion, same name, same residual form.

3. Final kappa substitution:
   - sympy lines 80-84: `expect_zero("KX*g_fail threshold with kappa inserted", g_fail_sq.subs({G_fail: G_fail_sub, kappa: KX * L**2 / TX}) - sp.simplify(m * cs_star_sq * TX * Nss * Pe_req / (rho_star * Osp**2 * L**2 * Deltainf)))`
   - mathematica lines 64-68: `expectZero["KX*g_fail threshold with kappa inserted", FullSimplify[gFailSq /. gFail -> gFailSubKappa, Assumptions -> $Assumptions] - m*csStarSq*tX*nSS*peReq/(rhoStar*oSP^2*ell^2*deltaInf)];`
   
   Same substitution choreography (substitute `G_fail`, then `kappa`), same RHS target rearranged identically.

The Mathematica script's only structural departure is to pre-compute `gFailSubKappa` once on line 62 rather than nesting two `subs` calls; this is a cosmetic difference, not an independent derivation path. Verdict: transliteration confirmed (F2).

## Engine cross-check

Both engines produce zero residual for every assertion and exit 0. They agree at the level claimed. The agreement is uninformative because both engines run the same algebra (see F2), but there is no contradiction between them. Outputs are fresh: sympy script mtime 2026-04-01, sympy output mtime 2026-05-11; mathematica script mtime 2026-05-11 11:56, mathematica output mtime 2026-05-11 12:55. No `stale_output` finding.

## Verdict justification

Both scripts run to completion and assert what they claim, but every assertion is algebraically guaranteed by the script's own definitions — substitute-and-rearrange identities, not anchored verifications of the threshold theorem. The Mathematica script is a literal port of the SymPy script and provides no independent algebraic path. The Cauchy-Schwarz bound underpinning the "best-case gain" interpretation is never used. The verdict is `findings` (not `clean`); not `stop_cold` because the fixes (derive `g_fail^2` via `solve` from `G_micro = G_fail`, and check `G_max == G_micro` at `C^2 = 1`) are local, mechanical, and do not propagate to downstream units' inputs — downstream units consume the *threshold formulas* (`g_fail_sq`, `G_max`, etc.), which are written down explicitly in this script's source; the new checks merely confirm their construction is consistent. Attacks tried that failed to escalate to UNFIXABLE: searched for sign errors in the `G_micro -> g_fail_sq` rearrangement (none — algebraic identity); searched for symbol-domain inconsistencies (`m, K_X, rho_star, etc.` all positive, consistent with their physical setup); searched for `simplify` calls hiding branch ambiguity (none — all expressions are rational in positive symbols).

## Self-test notes

Walked through each proposed new `expect_zero` mentally. For "g_fail^2 from solve": `sp.solve(G_micro - G_fail, g_phi**2)` returns the unique positive root `g_phi^2 = m cs* K_X N_ss G_fail/(rho* O_sp^2)` since `G_micro = G_fail` is linear in `g_phi^2`; the residual against `g_fail_sq` on sympy line 45 simplifies to 0 — confirmed substantive (would fail if `G_micro` were defined with a wrong factor). For "G_max at C^2=1": substituting `Osp**2 -> Nss*Npp` into `G_micro` gives `rho* g_phi^2 Nss Npp/(m cs* K_X Nss) = rho* g_phi^2 Npp/(m cs* K_X) = G_max` — residual is 0; would fail if `G_max` were declared with `N_ss` in place of `N_pp`. Variable independence: `solve(G_micro - G_fail, g_phi**2)` treats all other symbols as parameters — no derivative-of-constant trap. Symmetry/parity: no integrals are introduced. Trivial-case pre-check: setting `m = K_X = rho_star = cs_star_sq = N_ss = N_pp = O_sp = g_phi = 1` and `G_fail = 1` gives `g_fail_sq = 1` and `solve` returns `[1]`, residual 0; gives `G_max = 1` and at `C^2 = 1` returns `G_micro = 1`, residual 0. Paths: `.py` lives under `scripts/`, `.wl` lives under `mathematica/`, both confirmed in the existing tree.
