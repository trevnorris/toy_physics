---
unit_id: 036
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage036_support_feasibility_frontier.md]
  paper_appendix: present
---

# Audit unit 036 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_036.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage036_support_feasibility_frontier.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.txt`

(Note: a prior v1 directive at `/var/projects/toy_physics/research/pde_ledger/redteam/directives/stage_036.md` was already applied — its F1, F3, F4, F5 removed the previous round of tautological checks, and F2/F6 added the symbolic kappa derivation that is now present at sympy:123-139 / mathematica:124-138. The two findings below are v2 catches against the **current** state of the scripts, after those v1 fixes landed.)

## What the paper claims

Stage 036 isolates the dimensionless support-feasibility frontier left over once the Stage-035 normalization locus has selected `xi_req`. The card defines (i) the dimensionless mixed baseline `M_mix = 8 alpha_mix / (pi^2 A) = 8 Chi^2 / (pi^2 A Omega_U^2 Delta_0)`, (ii) the support-feasibility function `G(xi,delta) = 8 alpha_req / (pi^2 A) = 9 xi (xi+delta) / (9 delta + 11 xi)`, (iii) the required support loading `g_{B,req}^2/varpi^2 = (pi^2 A / 8)[G(xi_req,delta) - M_mix]`, (iv) the manifestly-positive monotonicity derivative `partial_xi G = 9(9 delta^2 + 18 delta xi + 11 xi^2)/(9 delta + 11 xi)^2`, (v) the endpoints `G(0,delta)=0` and `G_max(delta) = 9(1+delta)/(9 delta + 11)`, (vi) the near-onset expansion `G = xi - 2 xi^2/(9 delta) + O(xi^3)`, and (vii) the final admissibility test `R_target >= 1`, `F(xi_req,delta)=R_target`, `M_mix <= G(xi_req,delta)`. The Output line lists items (i), (ii), (iii), (iv), and the final-test triple as the load-bearing deliverables; the notes file restates the same content with more derivation context (`alpha_req`-split, `G_max` as the dimensionless form of the refined stability bound, parametric frontier interpretation).

## What the script claims to verify

The SymPy and Mathematica audits each: (1) print the closed forms of `G`, `F`, `M_mix`, `R_target`; (2) assert `(a_req - alpha_mix) - (pi^2 A/8)(G - M_mix) == 0` as the `g_{B,req}` factorization; (3) compute `dG/dxi` and check it matches the closed-form numerator `9(9 delta^2 + 18 delta xi + 11 xi^2)/(9 delta + 11 xi)^2`; (4) check `G(0,delta)=0` and `G_max = 9(1+delta)/(9 delta + 11)` via direct substitution / limit; (5) at a sample point `(delta=1, xi=1/2)`, check `R_target >= 1`, check `F(xi_req,delta) - R_target_host == 0` from a host-numeric microscopic kappa expansion, AND symbolically derive `R_target_sym` from the `kappa_0^2 = 8/pi^2`, `kappa_1^2 = 16/(9 pi^2)` expansion in `(A_sym, beta0_sym, xi, delta)` and check it equals `F`; (6) compare two sample `M_mix` values to `G(xi_req,delta)`; (7) check the near-onset series of `G` matches `xi - 2 xi^2/(9 delta)`. The Mathematica script adds an independent positivity argument via the discriminant of the numerator polynomial.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check(s) | Status |
|---|---|---|
| `M_mix` closed form | `Mmix_expr = 8 alpha_mix/(pi^2 A)` printed; used downstream | match |
| `G(xi,delta)` closed form | `G` printed; verified `G == 8 a_req/(pi^2 A)` (Mathematica line 57); definitional in sympy | match |
| `g_{B,req}^2/varpi^2 = (pi^2 A/8)(G - M_mix)` | sympy lines 68-71; Mathematica lines 59-62 | match (but tautological — see F2) |
| `dG/dxi = 9(9 d^2 + 18 d xi + 11 xi^2)/(9 d + 11 xi)^2 > 0` | sympy lines 74-77; Mathematica lines 68-81 (+ discriminant `-72 delta^2`) | match |
| `G(0,delta)=0`, `G_max = 9(1+delta)/(9 delta+11)` | sympy lines 78-82; Mathematica lines 82-92 | match |
| Near-onset `G = xi - 2 xi^2/(9 d) + O(xi^3)` | sympy lines 152-155; Mathematica lines 152-160 (coefficient-by-coefficient) | match |
| Final test `R_target >= 1` | sympy line 98; Mathematica line 110 | match (single sample point) |
| Final test `F(xi_req,delta) = R_target` | sympy lines 116-119 + symbolic 123-139; Mathematica lines 120-138 | match (host-sample + symbolic kappa re-derivation) |
| Final test `M_mix <= G(xi_req,delta)` | sympy lines 140-149; Mathematica lines 139-148 | **mismatch by tautology** — sample inequalities are constructed as `G_sample ± 1/10`, so the inequalities `G_sample - 1/10 <= G_sample` and `G_sample + 1/10 > G_sample` are arithmetic identities, not feasibility tests (F1) |

paper_alignment: aligned (all paper deliverables have script-side checks; the feasibility-witness check is structurally weak rather than misaligned with paper content).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 68-71 | `expect_zero("g_B,req^2/varpi^2 - (pi^2 A/8)(G - M_mix)", ...)` | g_{B,req} factorization | partial — identity follows from hardcoded `a_req` and `G` having matching closed forms |
| A2 | sympy | 77 | `expect_zero("dG/dxi - manifestly positive form", dG - dG_target)` | monotonicity | yes |
| A3 | sympy | 78 | `expect_zero("G(0,delta)", G.subs(xi,0))` | endpoint | yes |
| A4 | sympy | 82 | `expect_zero("G_max - closed form", Gmax - Gmax_target)` | G_max | yes |
| A5 | sympy | 87-90 | `expect_zero("final-test support inequality <-> nonnegative required support loading", ...)` | feasibility ↔ g_{B,req}^2 ≥ 0 | partial — same identity as A1 with xi→xi_req |
| A6 | sympy | 98 | `expect_true("admissible sample: R_target >= 1", ...)` | R_target ≥ 1 (final test) | yes (one sample) |
| A7 | sympy | 116-119 | `expect_zero("F(xi_req,delta) - R_target(host)", ...)` | F(xi_req,delta) = R_target | yes (numeric host) |
| A8 | sympy | 136-139 | `expect_zero("symbolic kappa derivation: F(xi,delta) - R_target_sym", ...)` | F(xi,delta) = R_target | yes (genuine: derives F from kappa expansion symbolically) |
| A9 | sympy | 140-144 | `expect_true("M_mix <= G(xi_req,delta)", Mmix_good <= G_sample, ...)` with `Mmix_good = G_sample - 1/10` | feasibility (M_mix ≤ G) | **no — tautological: (x − 1/10) ≤ x** |
| A10 | sympy | 145-149 | `expect_true("inadmissible sample: support deficit blocks the branch", Mmix_bad > G_sample, ...)` with `Mmix_bad = G_sample + 1/10` | feasibility (deficit) | **no — tautological: (x + 1/10) > x** |
| A11 | sympy | 155 | `expect_zero("G near-onset series through O(xi^2)", ...)` | near-onset expansion | yes |
| B1 | math | 57 | `expectZero["G - 8 alpha_req/(Pi^2 A)", ...]` | G definition | yes (a single algebraic simplification) |
| B2 | math | 59-62 | `expectZero["g_B,req^2/varpi^2 - (Pi^2 A/8)(G - M_mix)", ...]` | g_{B,req} factorization | partial — same tautology as A1 |
| B3 | math | 72-75 | `expectZero["dG/dxi positivity polynomial: ... == 11 xi^2 + 18 delta xi + 9 delta^2", ...]` | monotonicity | yes |
| B4 | math | 81 | `expectZero["dG/dxi numerator discriminant equals -72 delta^2", ...]` | manifest positivity (independent argument) | yes |
| B5 | math | 82 | `expectZero["G(0,delta)", gTarget /. xi -> 0]` | endpoint | yes |
| B6 | math | 92 | `expectZero["G_max - 9(1+delta)/(9delta+11)", ...]` | G_max | yes |
| B7 | math | 99-102 | `expectZero["final-test support inequality <-> nonnegative required support loading", ...]` | feasibility ↔ g_{B,req}^2 ≥ 0 | partial — same tautology as A5 |
| B8 | math | 110 | `expectTrue["R_target >= 1", fSample >= 1, ...]` | R_target ≥ 1 | yes (one sample) |
| B9 | math | 120-123 | `expectZero["F(xi_req,delta) - R_target(host)", ...]` | F(xi_req,delta) = R_target | yes (numeric host) |
| B10 | math | 134-138 | `expectZero["symbolic kappa derivation: F(xi,delta) - R_target_sym", ...]` | F(xi,delta) = R_target | yes (symbolic) |
| B11 | math | 139-143 | `expectTrue["M_mix <= G(xi_req,delta)", mMixGood <= gSample, ...]` with `mMixGood = gSample - 1/10` | feasibility | **no — tautological** |
| B12 | math | 144-148 | `expectTrue["inadmissible sample: support deficit blocks the branch", mMixBad > gSample, ...]` with `mMixBad = gSample + 1/10` | feasibility (deficit) | **no — tautological** |
| B13 | math | 152-160 | `expectZero["near-onset c0 = 0" / "c1 = 1" / "c2 = -2/(9 delta)", ...]` | near-onset expansion | yes (independently reads coefficients via `Coefficient[Series[...], ...]`) |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:96-99,140-149`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:108-109,139-148`

**What's wrong:**
The "admissible / inadmissible M_mix witness" checks for the paper-card feasibility condition `M_mix <= G(xi_req,delta)` are constructed so that they cannot fail regardless of what `G_sample` is. In sympy:

```python
Mmix_good = sp.simplify(G_sample - sp.Rational(1, 10))   # line 96
Mmix_bad  = sp.simplify(G_sample + sp.Rational(1, 10))   # line 97
...
expect_true(
    "admissible sample: M_mix <= G(xi_req,delta)",
    bool(Mmix_good <= G_sample),     # (x - 1/10) <= x  is always True
    f"M_mix={Mmix_good}, G={G_sample}",
)
expect_true(
    "inadmissible sample: support deficit blocks the branch",
    bool(Mmix_bad > G_sample),       # (x + 1/10) > x   is always True
    f"M_mix={Mmix_bad}, G={G_sample}",
)
```

The Mathematica script mirrors this exactly (`mMixGood = gSample - 1/10`, `mMixBad = gSample + 1/10`, then tests `mMixGood <= gSample` and `mMixBad > gSample` at lines 108-109, 139-148). Both inequalities are arithmetic identities on the rationals: `x - 1/10 <= x` and `x + 1/10 > x` are unconditionally true. They never test the physical feasibility condition `M_mix <= G(xi_req,delta)` because `M_mix` is being *defined* relative to `G_sample` rather than constructed independently (e.g. from a chosen `Chi, Omega_U, Delta_0, A` quartet) and then compared.

The output transcripts dutifully report `M_mix=53/145, G=27/58` and `M_mix=82/145, G=27/58` — and yes, those numerical inequalities happen to be correct, but the test code did not have to verify them: it pre-shifted them by ±1/10 from `G_sample` itself.

**Why this matters:**
This is the final-admissibility leg of the paper's Output (eq. `app-stage036-final-test`). The paper card lists the inequality `M_mix <= G(xi_req,delta)` as one of three load-bearing components of the admissibility test, and a reader would reasonably believe these two `expect_true` calls demonstrate (a) that the inequality holds for a physically admissible parameter pick and (b) that it fails when the support gap is overdrawn. They do neither. If someone later flipped the sign in the definition of `G` (or of `M_mix`), these checks would still pass — they don't even touch `M_mix_expr`. The v1 audit missed this entirely; the v2 paper-grounded read caught it because the paper card frames the inequality as the third leg of the admissibility test, which forced re-examination of the witness.

**Required change:**
Construct the sample `M_mix` from independent microscopic data, not from `G_sample`. At `(delta=1, xi=1/2)`, `G_sample = 27/58 ≈ 0.4655`. Provide a concrete admissible witness by choosing `(Chi, Omega_U, Delta_0, A) = (1, 1, 1, 29)`: then `M_mix_expr.subs(...) = 8/(29 pi^2) ≈ 0.0280`, well below `G_sample`. Provide a concrete inadmissible witness by choosing `(Chi, Omega_U, Delta_0, A) = (1, 1, 1, 1)`: then `M_mix_expr.subs(...) = 8/pi^2 ≈ 0.8106`, above `G_sample`. Numerical-evaluate both with `sp.N(...)` / `N[...]` and assert the inequality.

Concretely (sympy), replace lines 96-99 and lines 140-149 with:

```python
Mmix_admissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 29}))
Mmix_inadmissible = sp.N(Mmix_expr.subs({Chi: 1, OmegaU: 1, Delta0: 1, A: 1}))
G_sample_n = sp.N(G_sample)
expect_true(
    "admissible sample: M_mix < G(xi_req,delta)",
    bool(Mmix_admissible < G_sample_n),
    f"M_mix={Mmix_admissible}, G={G_sample_n}",
)
expect_true(
    "inadmissible sample: support deficit blocks the branch",
    bool(Mmix_inadmissible > G_sample_n),
    f"M_mix={Mmix_inadmissible}, G={G_sample_n}",
)
```

Mirror the same construction in Mathematica (lines 108-109, 139-148): replace `mMixGood = gSample - 1/10` and `mMixBad = gSample + 1/10` with `mMixAdmissible = N[mMix /. {Chi -> 1, OmegaU -> 1, Delta0 -> 1, A -> 29}]` and `mMixInadmissible = N[mMix /. {Chi -> 1, OmegaU -> 1, Delta0 -> 1, A -> 1}]`, and the `expectTrue` calls accordingly with `gSampleN = N[gSample]`. Both witnesses must be derived from `M_mix_expr / mMix` evaluated at independent parameter tuples — not by shifting `G_sample`.

**Verification:**
After the fix, the script's printed `M_mix=` and `G=` values for the two cases must no longer be `G_sample ± 1/10` (i.e. `53/145` and `82/145`) but instead be independently computed numbers approximately equal to `0.0280` (admissible) and `0.8106` (inadmissible). The verifier should confirm both lines now show `M_mix` numerators that are not arithmetically derived from `G_sample`.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:60-71`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage036_support_feasibility_frontier_sympy_audit.py:87-90`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:41-62`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage036_support_feasibility_frontier_mathematica_audit.wl:99-102`

**What's wrong:**
The "g_B,req^2/varpi^2 - (pi^2 A/8)(G - M_mix) == 0" check and its `xi -> xi_req` echo are identically zero by construction. In sympy:

```python
a_req = sp.simplify(9 * sp.pi**2 * A * xi * (xi + delta) / (8 * (9 * delta + 11 * xi)))   # line 60
gBreq_sq_over_varpi2 = sp.simplify(a_req - alpha_mix)                                      # line 61
...
expect_zero(
    "g_B,req^2/varpi^2 - (pi^2 A / 8) (G - M_mix)",
    gBreq_sq_over_varpi2 - (sp.pi**2 * A / 8) * (G - Mmix_expr),                          # line 68-71
)
```

Substituting the literals: `gBreq_sq_over_varpi2 = a_req - alpha_mix`, `Mmix_expr = 8 alpha_mix/(pi^2 A)`, so the residual is
`(a_req - alpha_mix) - (pi^2 A/8) G + alpha_mix = a_req - (pi^2 A/8) G`.
With `G = 9 xi(xi+delta)/(9 delta+11 xi)` and `a_req = 9 pi^2 A xi(xi+delta)/(8(9 delta+11 xi))`, the residual is identically zero — both quantities are independent hardcoded literals constructed to match. The check verifies the algebraic identity `G == 8 a_req/(pi^2 A)`, which is the same form a reader would write down by inspection. The line-87 `expect_zero` repeats the same identity with `xi -> xi_req`. Same issue in the Mathematica script (lines 41-62 set `alphaReq`, `g`, `gBReqSqOverVarpi2` from the same literals; lines 99-102 mirror the substitution).

**Why this matters:**
This is the only check that exercises the paper's defining equation `g_{B,req}^2/varpi^2 = (pi^2 A/8)(G - M_mix)`. Because both sides are hardcoded closed forms, the check confirms only that the author copied the same formula twice consistently. A genuine check would derive one side from the Stage-035 source (e.g. start from the split `alpha_req = alpha_mix + g_{B,req}^2/varpi^2` and verify the resulting `g_{B,req}^2/varpi^2` reduces to the boxed form when `alpha_req` is the Stage-035 closed form), not just assert `LHS == LHS`. This is less critical than F1 because the algebra is short and visually checkable, and because the kappa-derived F→R_target_sym chain (sympy:123-139, math:124-138) provides an independent anchor for `F` (and hence indirectly for `G = 8 alpha_req/(pi^2 A)` when `alpha_req` is read out of `F`'s structure). But the local lines 68-71 / 60-63 and 87-90 / 99-102 do not, by themselves, prove the factorization.

**Required change:**
Minimal version (preferred): add an inline comment above each of the four assertions clarifying that the check is a definitional self-consistency on the hardcoded forms, and pointing the reader to A8 / B10 (the symbolic kappa derivation) for the genuine anchor of the closed form. Specifically:

- sympy, before line 68: add comment `# Definitional identity: a_req and G are both hardcoded closed forms,`
  `# so this just confirms a_req = (pi^2 A / 8) G. The genuine anchor`
  `# for F (and hence for the closed form of G) is the symbolic kappa derivation`
  `# below at "symbolic kappa derivation: F(xi,delta) - R_target_sym".`
- sympy, before line 87: add comment `# Same definitional identity, with xi -> xi_req.`
- mathematica, before line 59: analogous comment.
- mathematica, before line 99: analogous comment.

Stronger version (optional, not required for this audit pass): express `g_{B,req}^2/varpi^2` symbolically as a free `alpha_req_sym - alpha_mix` (without baking the Stage-035 closed form into `alpha_req_sym`), substitute the Stage-035 closed form at the end, and confirm the boxed factorization emerges. This converts the check from a definitional identity into a real derivation step, but it is a structural rewrite rather than a one-line patch.

**Verification:**
For the minimal version: comments are present above the four assertions explicitly flagging them as definitional self-consistency. The printed residuals continue to be 0 in both saved outputs. The verifier confirms the comments landed at the named lines and the scripts still exit 0.

## Independent-derivation check (Mathematica)

The Mathematica script is **partially independent, partially a transliteration**. The arithmetic core (`alphaReq`, `g`, `fTarget`, `alphaMix`, `mMix`, `rTarget`, the host-sample point `(delta=1, xi=1/2, A=3, beta0=5)`, the same `kappa_0^2 = 8/Pi^2`, `kappa_1^2 = 16/(9 Pi^2)`, the symbolic kappa derivation at lines 124-138) is a direct line-by-line port of the SymPy script with renamed identifiers — same variable choreography, same expressions, same chosen sample. The two engines clearly share an author template.

However, Mathematica also performs three genuinely independent moves that the SymPy script does not:
1. Lines 68-75 derive `dG/dxi` and multiply through by `(9 delta + 11 xi)^2/9` to extract the polynomial `11 xi^2 + 18 delta xi + 9 delta^2`, then verify it matches term by term — a different algebraic route than sympy's direct `dG - dG_target` simplify.
2. Lines 78-81 compute `Discriminant[11 xi^2 + 18 delta xi + 9 delta^2, xi]` and confirm `disc = -72 delta^2`, giving an algebraic-discriminant argument for positivity that sympy does not perform.
3. Lines 152-160 use `Coefficient[Series[gTarget, {xi, 0, 2}], xi, k]` to *extract* the near-onset coefficients `c0, c1, c2` and independently check each against `0, 1, -2/(9 delta)` — sympy instead does a single `series` minus target comparison.

These three independent contributions clear the transliteration bar at the borderline (this restructuring is the legacy of the v1 F5 fix that was already applied — see the v1 directive). I do not file a `mathematica_transliteration` finding, but I note the structural overlap on the core algebra so downstream readers do not over-credit the Mathematica engine as a fully independent re-derivation. The two together are best read as a single algebraic derivation cross-checked by two algebra engines, which is acceptable for this stage's content (no nontrivial branch choices, no integration with multiple sign conventions).

## Engine cross-check

Both engines produce identical results at every shared check:

| Check | sympy output | Mathematica output |
|---|---|---|
| `G(xi,delta)` | `9*xi*(delta + xi)/(9*delta + 11*xi)` | `(9*xi*(delta + xi))/(9*delta + 11*xi)` |
| `F(xi,delta)` | `-(9*delta + 11*xi)**4/(81*(xi - 1)*(...)**2)` | `-1/81*(9*delta + 11*xi)^4/((-1 + xi)*(...)^2)` |
| `M_mix` | `8*Chi**2/(pi**2*A*Delta_0*Omega_U**2)` | `(8*Chi^2)/(A*Delta0*OmegaU^2*Pi^2)` |
| `R_target` | `pi**2*A*NQ/(8*beta0)` | `(A*NQ*Pi^2)/(8*beta0)` |
| `dG/dxi` residual | `0` | `0` |
| `dG/dxi` polynomial (math only) | n/a | `9*delta^2 + 18*delta*xi + 11*xi^2` |
| `Discriminant` (math only) | n/a | `-72*delta^2` |
| `G(0,delta)` | `0` | `0` |
| `G_max(delta)` | `9*(delta + 1)/(9*delta + 11)` | `(9*(1 + delta))/(11 + 9*delta)` |
| `R_target` sample | `1414562/558009` | `1414562/558009` |
| `F - R_target_host` (sample) | `0` | `0` |
| `F - R_target_sym` (symbolic kappa) | `0` | `0` |
| `M_mix sample` (admissible) | `M_mix=53/145, G=27/58` | `M_mix=53/145, G=27/58` |
| `M_mix sample` (inadmissible) | `M_mix=82/145, G=27/58` | `M_mix=82/145, G=27/58` |
| Near-onset series / coeffs | `xi - 2*xi**2/(9*delta)` | `c0=0, c1=1, c2=-2/(9*delta)` |

`engines_agree: true`. The agreement on the M_mix sample is also why F1 is not auto-caught by the engine cross-check: both engines compute the *same tautological inequality* and both report PASS, but neither tests the physical inequality.

Output freshness: sympy `.py` mtime is 2026-05-21 17:37, output `.txt` mtime is 17:40 — fresh. Mathematica `.wl` mtime 17:37, output `.txt` mtime 17:40 — fresh. No `stale_output` finding.

## Verdict justification

`findings`, not `clean`. The paper-side content is faithfully covered: all six numbered deliverables of the Output line (M_mix, G, g_{B,req} factorization, dG/dxi, G_max, near-onset expansion, final test) have script-side checks, and most of them are genuine — the dG/dxi check (especially Mathematica's polynomial and discriminant route), the endpoint and limit, the near-onset series (Mathematica reads coefficients independently), the symbolic kappa-based re-derivation of F (`F == R_target_sym` symbolically in `(xi, delta, A_sym, beta0_sym)`), and the host-sample numeric `F - R_target_host = 0` all hold up under attack. The `R_target >= 1` check is real (single sample at `F = 1414562/558009 ≈ 2.535`). I verified by hand that the symbolic kappa derivation reduces to `F` when `kappa_0^2 = 8/pi^2`, `kappa_1^2 = (2/9) kappa_0^2`, so that assertion is a genuine non-trivial identity.

What does not hold up is the third leg of the final admissibility test (`M_mix <= G(xi_req,delta)`): both engines test arithmetic identities of the form `(G_sample - 1/10) <= G_sample` and `(G_sample + 1/10) > G_sample`, which can never fail (F1). The `g_{B,req}^2/varpi^2` factorization check is also definitional in both engines (both sides hardcoded closed forms), recorded as F2 at low severity — the kappa-derived chain indirectly anchors the closed form, and the algebra is one line and visible, so a comment is sufficient.

Paper alignment is `aligned`. No `paper_misalignment`. The verdict is `findings` with `stop_cold: null` — Codex can fix both findings mechanically by (F1) replacing the witness construction with parameter-derived `M_mix` evaluations and (F2) adding clarifying comments above the four definitional assertions.

## Self-test notes

- Variable independence: in the F1 fix, `M_mix` for the admissible/inadmissible witnesses must be constructed from independent symbols `(Chi, Omega_U, Delta_0, A)` so the resulting numerical value is not algebraically tied to `G_sample`. I checked the proposed numerical witnesses: `M_mix_adm = 8/(29 pi^2) ≈ 0.0280`, `G_sample = 27/58 ≈ 0.4655` — admissible. `M_mix_inadm = 8/pi^2 ≈ 0.8106` — exceeds `G_sample`. The two witnesses straddle `G_sample` non-trivially, and the choice `(Chi=1, OmegaU=1, Delta0=1)` is consistent with the script's `Chi, OmegaU, Delta0 > 0` assumptions.
- Symmetry/parity: no integrals on unbounded domains; n/a.
- Trivial-case pre-check: F1's proposed witnesses give concrete inequalities (0.0280 < 0.4655 < 0.8106) that are non-vacuous; F2's required change is comment-only, so no trivial case applies.
- Path specifications: no missing scripts; both `.py` (scripts/) and `.wl` (mathematica/) already exist at the correct directories.
- Paper round-trip: F1 fix replaces the witness arithmetic with independent `M_mix` evaluations; the paper card does not state numeric witness values, so there is no risk of introducing a new paper_misalignment. F2 fix is comment-only; no risk.
