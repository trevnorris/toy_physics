---
unit_id: 146
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage146_positive_deformation_expansion.md]
  paper_appendix: present
---

# Audit unit 146 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_146.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage146_positive_deformation_expansion.md`
- part appendix: `/var/projects/toy_projects/.../paper/appendices/stage_appendix_part04.tex` → actual: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.txt`

## What the paper claims

Stage 146 is a finite mouth-profile corrections ledger step (stage_146.tex:7). The notes (the authoritative derivation source) state the question: given the unique regular Family-1 canonical exponential mouth branch with fixed constants `Pi_* ≈ 1.50882951349316`, `g_* ≈ 0.758035078944663`, `S_* ≈ 0.658075937605429` (notes lines 13-17), how must the canonical point move under a positive, normalized, *non-exponential* mouth-source deformation. The stage's deliverables are: (1) the exact convex deformation family `Σ_ε = (1-ε)Σ_* + ε·ς` preserving positivity and normalization (notes 41-49); (2) the exact affine moment laws `ḡ_ε = g_* + ε(ḡ_ς - g_*)`, `S̄_ε = S_* + ε(S̄_ς - S_*)` showing the first-order problem is two-dimensional in `(ḡ_ς, S̄_ς)` (notes 70-89); (3) the compensation law `δΠ = -ε(ḡ_ς - g_*)/g'_*` keeping the overlap pinned to `g_*` (notes 114-120); (4) the numeric slopes `g'_* ≈ 0.0714453558083195`, `S'_* ≈ 0.0483709542125041` (notes 124-129); and (5) the first non-exponential correction `δS_q = ε[(S̄_ς - S_*) - (S'_*/g'_*)(ḡ_ς - g_*)]` (notes 133-142). The .tex `\stagefield{Verification}` names this stage's audit target as exactly the SymPy/Mathematica transcripts. The appendix row (part04:31) classifies stages 146-153 as "first-order mouth-profile rigidity and finite mouth-only correction."

## What the script claims to verify

The SymPy script verifies the closed-form moment formulas `g(Π)` (py:24) and `S_q(Π)` (py:25-28) equal their defining integrals `∫Σ·cos(πx/2)` and `∫Σ·K_q` (py:33-53, with a numeric fallback at 4 samples for the S integral), prints both moments, runs three numeric kernel cross-checks (py:62-68), locates `Π_*` as the root of `g(Π) = g_-^{F1}` (py:71-74), reports `g_*, S_*, g'_*, S'_*` (py:76-83), forms the symbolic `δΠ` and `δS` retuning laws (py:87-94), and finally tests an "affine law (integral form)" residual for `ḡ_ε` and `S̄_ε` against a concrete bump profile `ς = 6x(1-x)` at `ε ∈ {1/10, 1/2}` (py:100-127). The Mathematica `.wl` performs the same sequence with the same formulas, same root condition, and the same affine-residual construction.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| (1) convex family Σ_ε preserving normalization/positivity | Σ_ε built py:101; normalization of ς=6x(1-x) implicit (∫=1) | match (structurally used) |
| (2a) affine law ḡ_ε = g_* + ε(ḡ_ς-g_*) | py:114 residual; but subtracts `gminus`, so measures root residual (1-ε)(g(Π_*)-gminus), not the affine structure | partial — label over-claims (F1) |
| (2b) affine law S̄_ε = S_* + ε(S̄_ς-S_*) | py:115 residual; subtracts `Sformula(Π_*)` = same source as Sbar_phys → identically zero | mismatch — tautological (F2) |
| (3) δΠ = -ε(ḡ_ς-g_*)/g'_* | py:87 builds dPi = -eps*(gbar-gminus)/g'(Π_*) | match (symbolic form printed) |
| (4) g'_* ≈ 0.07144…, S'_* ≈ 0.04837… | py:78-83 derive by sp.diff; outputs match notes to all digits | match |
| (5) δS_q = ε[(S̄_ς-S_*) - (S'_*/g'_*)(ḡ_ς-g_*)] | py:88 builds dS = eps*(Sbar - S(Π_*)) + S'(Π_*)*dPi | match (symbolic form printed) |
| g_* ≈ 0.75803…, S_* ≈ 0.65807…, Π_* ≈ 1.50882… | py:74,76-77 | match |

`paper_alignment: aligned` — every paper-side deliverable value reconciles with the script outputs; the two findings are *within* the script's own verification scope (an over-claimed label and a tautological check), not paper↔script value disagreements.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 42/53 | `expect_zero(gDirect - gFormula)` / numeric-sample fallback for S | moment formulas (feeds 2-5) | yes |
| A2 | sympy | 67-68 | `abs(integral - formula) < 1e-12` at 3 Π | S_q formula numeric anchor | yes |
| A3 | sympy | 73 | `Pi_star = nsolve(gPi - gminus, 1.5)` | Π_* root (deliverable) | yes (defines Π_*) |
| A4 | sympy | 122-125 | `abs(g_res) < 1e-25`, `abs(S_res) < 1e-25` | labeled "affine law"; g_res = root residual, S_res = identically 0 | g: partial (mislabeled); S: NO (tautological) |
| A5 | mathematica | 48/51 | `expectZero[gDirect-gFormula]`, `expectZero[sDirect-sFormula]` | moment formulas | yes |
| A6 | mathematica | 59 | `expectApprox[..., 10^-12]` at 3 Π | S_q numeric anchor | yes |
| A7 | mathematica | 79 | `expectApprox["Pi_* compensation point", gStar, gMinus, 10^-20]` | Π_* closes g = g_- | yes |
| A8 | mathematica | 108/118 | `Abs[gEpsSampleN] < 10^-25` / `Abs[sEpsSampleN] < 10^-25` | same "affine law" mislabel + S tautology as A4 | g: partial; S: NO |

## Findings

### F1 — insufficient_verification (over-claimed label)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:114, 120, 126`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:103, 106, 109`

**What's wrong:**
The g-residual is `gbar_phys - (gminus + eps*(gbar_v - gminus))` (py:114). By linearity of the integral, `gbar_phys = (1-eps)*g(Π_*) + eps*gbar_v`, so the residual reduces **exactly** to `(1-eps)*(g(Π_*) - gminus)` — i.e. `(1-eps)` times the *root residual* of the `Π_*` solve `g(Π_*) = gminus`. The inline comment (py:106-108) is honest about this mechanism ("the residual collapses by linearity ... to (1-eps)*(gPi(Pi_*) - gminus)"). But the printed label and PASS line call it "`g_eps affine law (integral form)`" (py:120, 126). What is actually measured is **root-solver accuracy / g-minus closure**, NOT the affine (convexity-in-ε) structure the paper deliverable (2a) describes. The affine structure itself is true by linearity of the integral for *any* root accuracy and is never the quantity being bounded. The `.wl` reproduces the identical construction (wl:103) and label (wl:106). This is the CAVEAT the audit brief named: the residual is not an independent `g_minus` guard, and it must not be read (in label) as verifying the affine law.

**Why this matters:**
A reader trusting the label believes Stage 146's affine moment law has been independently exercised against a concrete deformation. It has not — the bound is a root-accuracy check that is already covered by A7 (`Pi_* compensation point`, wl:79). The over-claim inflates the apparent coverage of deliverable (2a).

**Required change:**
Rename the label string at py:120 and py:126 (and wl:106, wl:109) so it states what is measured, e.g. `g_eps root-closure residual (1-eps)*(g(Pi_*)-g_minus)` rather than `g_eps affine law`. Do NOT change the residual expression or the tolerance — the numeric check is sound, only its description is wrong.

**Verification:**
After the fix, py:120/py:126 print a label that no longer says "affine law" for the g-residual; the residual value (~1e-51) is unchanged. The verifier confirms the script still exits 0 and the printed label matches the measured quantity.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:115, 124, 127`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:113, 116, 118`

**What's wrong:**
The S-residual is `Sbar_phys - (Sformula.subs(Pi, Pi_star) + eps*(Sbar_v - Sformula.subs(Pi, Pi_star)))` (py:115). Unlike the g-residual, this subtracts `Sformula(Π_*)` — which is **the same quantity** as the Σ-part of `Sbar_phys`: by linearity `Sbar_phys = (1-eps)*[∫Σ(Π_*)·K_q] + eps*Sbar_v`, and `∫Σ(Π_*)·K_q = Sformula(Π_*)` (verified exactly upstream at py:44-53). Substituting, the residual is `(1-eps)*Sformula(Π_*) + eps*Sbar_v - Sformula(Π_*) - eps*Sbar_v + eps*Sformula(Π_*) = 0` **identically**, regardless of any physics. The only thing the ~1e-50 output residual reflects is integration-engine round-off between `sp.integrate` and `Sformula.subs`. So `assert abs(S_res) < 1e-25` (py:124) cannot fail on the affine law; it is a tautology (`x - x = 0`). The asymmetry with the g-residual is exactly why g produces a meaningful (1-eps)(g(Π_*)-gminus) while S produces literal zero: g subtracts the *independent* algebraic target `gminus`, S subtracts its *own* derived value. The `.wl` (wl:113) copies the same tautological construction.

**Why this matters:**
Deliverable (2b) — the affine law for `S̄_ε` — has **no non-tautological script-side check**. The `S_eps affine law` PASS is guaranteed by construction and would still pass if `Sformula` were wrong, if `K_q` were wrong, or if the convex family were broken, because both sides are built from the same `Sformula(Π_*)`. The stage's S-side affine claim is effectively unverified.

**Required change:**
Make the S-residual subtract an *independent* anchor, mirroring how the g-residual subtracts `gminus`. Concretely, compute the affine prediction from the **separately-evaluated** moments rather than re-using `Sformula.subs(Pi, Pi_star)` for both `Sbar_phys`'s base and the prediction. One faithful route: form `S_at_star_indep = sp.N(sp.integrate((Sigma.subs(Pi, Pi_star))*Kq, (x,0,1)), 50)` and `S_eps_pred = (1-eps)*S_at_star_indep + eps*Sbar_v`, then test `Sbar_phys_numeric - S_eps_pred` — but since both still trace to the same exact integral this stays ~round-off. The genuinely discriminating check is to assert the affine law against the *closed-form* `Sformula(Π_*)` while computing `Sbar_phys` by an *independent* numeric quadrature of `Σ_ε·K_q` (e.g. high-precision `mpmath`/`NIntegrate` of the integrand, not symbolic linearity), so a wrong `Sformula` or wrong `K_q` would surface. Note: this is a `## Resolve before fix_loop` design question — the correct independent anchor is not mechanically obvious, so route to the user rather than auto-applying.

**Verification:**
After the fix, py:115/124 (and wl:113/118) test `S̄_ε` against an anchor not algebraically identical to `Sbar_phys`'s own base term, so perturbing `Sformula` or `K_q` would make the assertion fail; the script still exits 0 on the correct forms.

## Independent-derivation check (Mathematica)

The `.wl` is a **transliteration**, not an independent derivation.

1. **Closed-form moments — identical, not re-derived.** `.wl:41` `gFormula = 2*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1))` is the character-for-character image of `.py:24` `gPi = 2*Pi*(2*Pi*sp.exp(Pi)+sp.pi)/((4*Pi**2+sp.pi**2)*(sp.exp(Pi)-1))`. Likewise `.wl:42` `sFormula = p*(kap*Tanh[kap] + p*(Exp[-p]*Sech[kap] - 1))/((1 - Exp[-p])*(kap^2 - p^2))` mirrors `.py:26-28`. Both engines then *verify the same hand-supplied formula* against `Integrate[sigma*...]` rather than each deriving the closed form from scratch.

2. **Root condition — same choreography.** `.wl:64-66` `rF1 = Sqrt[(12*(37/20)^2)/Pi^2 - 1]; gMinus = ...; pStar = p /. FindRoot[gFormula == gMinus, {p, 1.5}, ...]` is the same variable choreography as `.py:71-73` `rF1 = sp.sqrt(12*sp.Rational(37,20)**2/sp.pi**2 - 1); gminus = ...; Pi_star = sp.nsolve(gPi - gminus, 1.5)`.

3. **Affine residual — same construction including the same g/S asymmetry.** `.wl:103` `gEpsRes = ((1 - eps)*(gDirect /. p -> pStar) + eps*gBarV) - (gMinus + eps*(gBarV - gMinus))` and `.wl:113` `sEpsRes = ((1 - eps)*(sDirect /. p -> pStar) + eps*sBarV) - (sStar + eps*(sBarV - sStar))` reproduce the *exact* py:114/115 structure — g anchored to `gMinus`/`gminus`, S anchored to its own `sStar`/`Sformula(Π_*)`. The `.wl` therefore inherits the identical F1 over-claim and F2 tautology, which is itself strong evidence of line-by-line porting: an independent re-derivation would not be expected to reproduce the same subtle S-side tautology. This is a `mathematica_transliteration` characteristic; I record it here rather than as a separate finding because both engines are present, agree, and the transliteration does not by itself change a verified value — but it means the second engine provides no independent cross-check of the S affine law (it shares F2).

## Engine cross-check

Both engines agree to full reported precision:

| value | sympy (.txt) | mathematica (.txt) |
|---|---|---|
| Π_* | 1.5088295134931555830055507559542749287786371931531 (l.27) | 1.5088295134931555830055507559542749287786371931530784314262 (l.20) |
| g_* | 0.758035078944662826919680890414 (l.28) | 0.7580350789446628269196808904141104577505… (l.21) |
| S_* | 0.658075937605429274616601849160 (l.29) | 0.6580759376054292746166018491599488… (l.22) |
| g'_* | 0.0714453558083195211894603019881 (l.30) | 0.071445355808319521189460301988135… (l.23) |
| S'_* | 0.0483709542125040992653572761725 (l.31) | 0.0483709542125040992653572761724701… (l.24) |

Affine residuals: sympy ~1e-50/1e-51 (l.199-202), mathematica ~1e-58 (g) and `0`/round-off (S) (l.33-37). Both pass their 1e-25 bounds. No `engine_disagreement`.

## Verdict justification

`findings` (2). Every load-bearing *value* in the stage reconciles exactly between the scripts and the notes (Π_*, g_*, S_*, g'_*, S'_* all match to all printed digits; the δΠ and δS_q symbolic laws are built in the form the notes state), so `paper_alignment: aligned` and there is no `paper_misalignment`. The findings are *internal* to the script's verification: F1 (low) — the g-side "affine law" PASS actually measures `(1-eps)·(g(Π_*)-gminus)` root closure, an over-claimed label the inline comment itself contradicts; F2 (medium) — the S-side "affine law" residual is algebraically identically zero (subtracts its own derived `Sformula(Π_*)`), a tautology that leaves deliverable (2b) without a discriminating check. I attacked: the `168π²`→`100π²` stale-constant trap (the output's `4107 - 100π²` is correct: `12·(37/20)²·100 = 4107`, denominator `100`, not `168`); symbol domains (`Pi positive,real`, `0≤x≤1` consistent with the mouth interval); the kernel formula numeric anchors (genuine, 1e-12); engine agreement (exact). The two surviving findings are the affine-law checks; both .wl and .py share them because the .wl is a transliteration.

## Self-test notes

Checked: (1) Variable independence — `sp.diff(gPi, Pi)`, `sp.diff(Sformula, Pi)` (py:78-79) both depend on `Pi`; nonzero, fine. (2) Linearity collapse — manually expanded both affine residuals: g → `(1-eps)*(g(Π_*)-gminus)` (root residual, nonzero), S → `0` identically (tautology); this is the basis of F1/F2. (3) Stale-constant trap — confirmed `100π²` is correct, not stale `168π²`. (4) Outputs fresh — both .txt mtimes (21:46) newer than scripts (21:37/21:39); no `stale_output`. F2's required change is flagged as a design question (independent anchor not mechanically obvious) and is routed via the directive rather than a blind mechanical edit.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 5 deliverable values checked, 0 misaligned.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| Π_* = 1.50882951349316 | py:74, .txt:27 / wl:68, .txt:20 | notes l.13 `\Pi_* \approx 1.50882951349316` | MATCH |
| g_* = 0.758035078944663 | py:80, .txt:28 / wl:75, .txt:21 | notes l.15 `\mathfrak g_*\approx 0.758035078944663` | MATCH |
| S_* = 0.658075937605429 | py:81, .txt:29 / wl:76, .txt:22 | notes l.16 `\mathcal S_*\approx 0.658075937605429` | MATCH |
| g'_* = 0.0714453558083195 | py:82, .txt:30 / wl:77, .txt:23 | notes l.125 `\mathfrak g'_* \approx 0.0714453558083195` | MATCH |
| S'_* = 0.0483709542125041 | py:83, .txt:31 / wl:78, .txt:24 | notes l.126 `\mathcal S'_* \approx 0.0483709542125041` | MATCH |
| δΠ = -ε(ḡ_ς-g_*)/g'_* (symbolic) | py:87, .txt:36-63 / wl:81, .txt:31 | notes l.114-120 boxed `\delta\Pi = -\epsilon\frac{\bar g_\varsigma-\mathfrak g_*}{\mathfrak g'_*}` | MATCH (form) |
| δS_q = ε[(S̄_ς-S_*)-(S'_*/g'_*)(ḡ_ς-g_*)] (symbolic) | py:88, .txt:64+ / wl:82, .txt:32 | notes l.133-142 boxed | MATCH (form) |

INTERNAL (scaffolding, no finding expected in prose): `rF1 = sqrt(12·(37/20)²/π² - 1)` and `gminus = rF1 - sqrt(1+rF1²)/2` (the Family-1 algebraic target that `Π_*` solves to; its numeric value *is* g_* = 0.758…, which is reported); the three kernel-check sample integrals at Π∈{1, 3/2, 5/2}; the four g/S numeric-sample diffs; the affine-law residuals (~1e-50/1e-58); `varsigma_test = 6x(1-x)` test profile; `gbar_v, Sbar_v` test-profile moments. None are stage deliverables.

All seven emitted deliverable values (5 numeric constants + 2 symbolic laws) reconcile with the notes; the `.tex` card is terse-by-design and carries none numerically, which is legitimate since the notes carry them. No MISMATCH, no MISSING-DELIVERABLE → no value-reconciliation `paper_misalignment` finding.
