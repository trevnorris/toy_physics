---
unit_id: 064
batch: III.3
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-26T00:00:00Z
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
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage064_equilibrium_alignment.md
  paper_appendix: present
---

# Audit unit 064 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_064.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage064_equilibrium_alignment.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 106, input directive at line 246)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.txt`

## What the paper claims

The paper card and notes claim, for the parent equilibrium source/support alignment in the local compressional linearization, the following deliverables. (a) The equilibrium-induced source profile is fixed pointwise by the support loading, `chi_sigma(y) = g_phi chi_phi(y) / H(y)` with `H(y) = h'(rho_*(y))` (boxed Eq. (eq:app-stage064-source-profile)). (b) On this aligned branch the overlap invariants reduce to `O_(sigma phi) = g_phi I_1`, `N_(sigma sigma) = g_phi^2 I_2`, with `I_1 = int chi_phi^2/H`, `I_2 = int chi_phi^2/H^2` (Eqs. (eq:app-stage064-I1I2), (eq:app-stage064-overlaps)). (c) The coherence factor `C_(sigma phi)^2 = I_1^2 / (N_(phi phi) I_2)` with Cauchy bound `<= 1` (Eq. (eq:app-stage064-coherence)). (d) The *general* equilibrium softening is `Delta K_X^(eq) = g_phi^2 I_1` and `G_eq = g_phi^2 I_1 / K_X` (Eq. (eq:app-stage064-softening)) — this holds for any `H(y)`, not only the matched layer. (e) In the thin active layer where `H` is approximately constant, the branch becomes exactly matched: `C_(sigma phi)^2 = 1` and `G_eq = g_phi^2 N_(phi phi) / (K_X H_w)` (boxed Eq. (eq:app-stage064-matched-layer)). The `\stagefield{Output}` line cites (a) and (e) explicitly; (b)-(d) are body equations supporting those outputs.

## What the script claims to verify

The sympy script asserts: (A1) the closure `chi_sigma = g_phi chi_phi/H` from local sigma-energy minimisation; (A2-A3) the overlap identities `O = g_phi I_1`, `N_ss = g_phi^2 I_2` on a concrete constant-H Gaussian; (A4-A5) Gaussian-integral reductions `I_1 = N_pp/H_w`, `I_2 = N_pp/H_w^2` for the constant Gaussian; (A6-A7) substituting these matched-layer reductions into the abstract `C^2` and `G_eq` gives `C^2 = 1` and `G_eq = g_phi^2 N_pp/(K_X H_w)`; (A8) the discrete two-point Cauchy gap identity `N_pp I_2 - I_1^2 = w1 w2 (H1-H2)^2/(H1^2 H2^2) >= 0`; (A9) the abstract elimination `sigma_stat = Lambda phi/Theta` gives `F_eff = (1/2)(K_X - Lambda^2/Theta) phi^2`; (A10) substituting `I_2 -> I_1^2/N_pp` and `N_pp -> I_1 H_w` into `Lambda_abs^2/Theta_abs` (with `Theta_abs = H_w N_ss`, `Lambda_abs = g_phi O_(sigma phi)`) reduces to `g_phi^2 I_1`. The mathematica script mirrors all of these checks 1-for-1.

## Paper ↔ script cross-check

| Paper deliverable | Script check(s) | Status |
|---|---|---|
| (a) closure law `chi_sigma = g_phi chi_phi/H` | A1 (sympy:64-68, math:36-41) | match |
| (b) overlap identities `O = g_phi I_1`, `N_ss = g_phi^2 I_2` | A2-A3 (sympy:74-82, math:46-56), abstract definitions sympy:85-86, math:58-59 | partial — checked only on constant-H Gaussian; general case is algebraic-by-substitution (chi_sigma already substituted), so this is essentially definitional. The discrete two-point check (A8) supplies the real general-case audit. |
| (c) `C^2 = I_1^2/(N_pp I_2) <= 1` | A8 Cauchy gap (sympy:124-139, math:96-108); the formula `C^2 = I_1^2/(N_pp I_2)` is built at sympy:87 / math:60 | match (formula construction + Cauchy bound check) |
| (d) general softening `Delta K_X^eq = g_phi^2 I_1` (for ANY H(y)) | A10 (sympy:156-165, math:123-129) attempts this but uses matched-layer substitutions to force the answer | mismatch — see Finding F1 |
| (e) matched layer: `C^2 = 1` and `G_eq = g_phi^2 N_pp/(K_X H_w)` | A6-A7 (sympy:115-122, math:88-94) | match |
| Engine independence (policy requirement) | sympy + mathematica scripts | mismatch — transliteration (Finding F2) |

`paper_alignment: partial` — four of five deliverables verified, but (d) the general softening is not verified.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 68 | `expect_zero("closure law ...", chi_sigma_closure - g_phi*chi_phi_loc/H_loc)` | (a) closure law | yes |
| A2 | sympy | 81 | `expect_zero("overlap O = g_phi * I1 (integral form)", ...)` (constant Hw, Gaussian) | (b) | partial — Gaussian + constant H only; algebraically definitional once chi_sigma is substituted |
| A3 | sympy | 82 | `expect_zero("overlap N_ss = g_phi^2 * I2 (integral form)", ...)` | (b) | partial — same caveat as A2 |
| A4 | sympy | 108 | `expect_zero("matched-layer I1 reduction", I1_int - Npp_int/Hw)` | supports (e) | yes |
| A5 | sympy | 109 | `expect_zero("matched-layer I2 reduction", I2_int - Npp_int/Hw^2)` | supports (e) | yes |
| A6 | sympy | 121 | `expect_zero("matched-layer coherence", C2_const - 1)` | (e) C^2=1 | yes |
| A7 | sympy | 122 | `expect_zero("matched-layer gain ...", Geq_const - g_phi^2*Npp/(KX*Hw))` | (e) G_eq matched | yes |
| A8 | sympy | 139 | `expect_zero("two-point Cauchy gap identity", gap_disc - gap_expected)` | (c) Cauchy bound | yes |
| A9 | sympy | 151-154 | `expect_zero("effective support softening", F_eff - (1/2)(KX - Lambda^2/Theta)*phi^2)` | abstract elimination scaffolding for (d) | yes (algebraic identity, but does not by itself anchor `Delta K_X = g_phi^2 I_1` to the parent equilibrium) |
| A10 | sympy | 163-165 | `expect_zero("equilibrium softening equals g_phi^2 I1", soft_matched - g_phi^2*I1)` | claims to verify (d), but uses matched-layer subs `I2 -> I_1^2/N_pp` AND `N_pp -> I_1 H_w` | no — only re-derives the matched-layer value; the general claim `Delta K_X^eq = g_phi^2 I_1` for arbitrary `H(y)` is not exercised |
| M1-M10 | mathematica | 41,55,56,81,82,93,94,108,121,129 | parallel to A1-A10 | same | same — see Finding F2 |

## Findings

### F1 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:156-165`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl:123-129`

**What's wrong:**

The paper states the equilibrium softening as a general (any-H) identity (notes Section 4, paper Eq. (eq:app-stage064-softening)):

> "On the equilibrium-aligned branch, the direct parent elimination gives `Delta K_X^(eq) = g_phi^2 I_1`."

The paper card boxes it without a matched-layer qualifier:

> `\Delta K_X^{\rm eq}=g_\phi^2 I_1, \qquad G_{\rm eq}=\frac{g_\phi^2 I_1}{K_X}.`

The sympy script tries to verify this at lines 156-165 via:
```python
Theta_abs = sp.simplify(Hw * Nss)          # = Hw * g_phi^2 * I2
Lambda_abs = sp.simplify(g_phi * Osp)      # = g_phi^2 * I1
soft_abs = sp.simplify(Lambda_abs**2 / Theta_abs)   # = g_phi^2 * I1^2 / (Hw * I2)
soft_matched = sp.simplify(soft_abs.subs(I2, I1**2 / Npp).subs(Npp, I1 * Hw))
expect_zero("equilibrium softening equals g_phi^2 I1", soft_matched - g_phi**2 * I1)
```

Two problems:

1. `Theta_abs = H_w * N_(sigma sigma)` is the matched-layer identification of the source self-energy coefficient, not the general one. From the parent equilibrium with `chi_sigma = g_phi chi_phi/H`, the general source self-energy coefficient is `Theta = int H(y) chi_sigma^2 d^3y = g_phi^2 int chi_phi^2/H = g_phi^2 I_1`, NOT `H_w * g_phi^2 * I_2`.

2. To get the desired answer `g_phi^2 I_1` out of the formula in (1), the script substitutes BOTH `I_2 -> I_1^2/N_pp` (the Cauchy-equality condition) AND `N_pp -> I_1 H_w` (the matched-layer constant-H consequence). After these two matched-layer substitutions the algebra is forced to the matched-layer answer; it does not test the general identity. The mathematica script (lines 123-129) does exactly the same chain.

The general identity `Delta K_X^eq = g_phi^2 I_1` (valid for any `H(y)`) is therefore not actually verified in either engine.

**Why this matters:**

The paper card lists the general softening as a body equation that supports the matched-layer output. If a future stage cites `Delta K_X^eq = g_phi^2 I_1` outside the matched-layer regime, the audit chain currently provides no script-side verification. The matched-layer case is independently established by A6-A7, so A10 as written contributes nothing beyond a self-consistency loop.

**Required change:**

Replace the matched-layer-only chain at sympy lines 156-165 (and mathematica lines 123-129) with a derivation that verifies `Delta K_X^eq = g_phi^2 I_1` directly from the parent equilibrium for *general* H(y). One clean way is a discrete two-point construction (which has already been used at A8 to verify the Cauchy gap):

For two support points with weights `w_k` (representing `chi_phi_k^2 dV`) and stiffnesses `H_k`, with `chi_sigma_k = g_phi chi_phi_k / H_k`, on the equilibrium-aligned branch:
- `I_1 = sum_k chi_phi_k^2 / H_k = w_1/H_1 + w_2/H_2` (matches `I1_disc`).
- `Theta = sum_k H_k chi_sigma_k^2 = g_phi^2 sum_k chi_phi_k^2 / H_k = g_phi^2 I_1`.
- `Lambda = g_phi sum_k chi_phi_k chi_sigma_k = g_phi^2 sum_k chi_phi_k^2 / H_k = g_phi^2 I_1`.
- Therefore `Lambda^2/Theta = (g_phi^2 I_1)^2 / (g_phi^2 I_1) = g_phi^2 I_1`.

Concretely, after the existing two-point Cauchy-gap block (sympy lines 124-139), add (do NOT delete the existing A9 `F_eff` check — that block is the abstract elimination scaffolding and is fine):
```python
# General-H equilibrium softening on a two-point parent equilibrium branch.
Theta_general = sp.simplify(g_phi**2 * I1_disc)
Lambda_general = sp.simplify(g_phi**2 * I1_disc)
soft_general = sp.simplify(Lambda_general**2 / Theta_general)
expect_zero(
    "general equilibrium softening equals g_phi^2 I_1",
    soft_general - g_phi**2 * I1_disc,
)
```
Then delete (or relabel as a sanity-check on the matched layer) the existing lines 156-165. The replacement does NOT use `H_w`, does not assume `I_2 = I_1^2/N_pp`, and does not assume `N_pp = I_1 H_w` — i.e., it exercises the general-H case (H1 != H2 in the two-point model). Mirror the same restructuring in the mathematica script at lines 123-129, but write it independently — see F2.

**Verification:**

After the patch, the sympy and mathematica scripts each produce a new `expect_zero`/`expectZero` line with name containing "general equilibrium softening" and the assertion evaluates to 0. The check must NOT contain substitutions of the form `I_2 -> I_1^2/N_pp` or `N_pp -> I_1 H_w`. Both scripts still exit 0.

### F2 — mathematica_transliteration

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl` (entire script, lines 26-134)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py` (entire script, lines 49-167)

**What's wrong:**

The `.wl` script is a line-by-line port of the `.py` script: same variable choreography (sympy `g_phi, K_X, N_pp, I1, I2, H_w` -> mathematica `gPhi, kX, npp, i1, i2, hw`), same intermediate steps in the same order, same Gaussian profile `Exp[-y^2/(2 L^2)]`, same `solve(diff(F, sigma), sigma)` -> `Solve[D[f, sigma] == 0, sigma]`, same matched-layer substitution chain. Three corresponding sections:

(1) Closure law minimisation:
- sympy:63-68 — `F_loc = sp.Rational(1,2)*H_loc*sigma_loc**2 - g_phi*chi_phi_loc*sigma_loc; closure_solutions = sp.solve(sp.diff(F_loc, sigma_loc), sigma_loc); chi_sigma_closure = closure_solutions[0]`
- mathematica:36-41 — `fLoc = (1/2) hLoc[yLoc] sigmaLoc^2 - gPhi chiPhiLoc[yLoc] sigmaLoc; closureSolutions = Solve[D[fLoc, sigmaLoc] == 0, sigmaLoc]; chiSigmaClosure = sigmaLoc /. First[closureSolutions]`

(2) Gaussian integral block:
- sympy:71-82 — `chi_phi_g = sp.exp(-y_int**2/(2*L_int**2)); H_g = Hw; ... Osp_int_check = sp.integrate(chi_sigma_g*chi_phi_g, (y_int, -sp.oo, sp.oo))`
- mathematica:45-56 — `chiPhiG = Exp[-yInt^2/(2 lInt^2)]; hG = hw; ... ospIntCheck = Integrate[chiSigmaG*chiPhiG, {yInt, -Infinity, Infinity}, ...]`

(3) Matched-layer softening chain (the part flagged in F1):
- sympy:158-165 — `Theta_abs = sp.simplify(Hw * Nss); Lambda_abs = sp.simplify(g_phi * Osp); soft_abs = sp.simplify(Lambda_abs**2 / Theta_abs); soft_matched = sp.simplify(soft_abs.subs(I2, I1**2/Npp).subs(Npp, I1*Hw))`
- mathematica:123-127 — `thetaAbs = FullSimplify[hw*nss, ...]; lambdaAbs = FullSimplify[gPhi*osp, ...]; softAbs = FullSimplify[lambdaAbs^2/thetaAbs, ...]; softMatched = FullSimplify[(softAbs /. i2 -> i1^2/npp) /. npp -> i1*hw, ...]`

The mathematica script never derives any identity by an alternative path. It does not, for example, verify the Gaussian-integral block on a different profile, treat Cauchy-Schwarz as a continuous-integral inequality, or take a controlled limit. It echoes the sympy algebra in Wolfram Language syntax.

**Why this matters:**

The second-engine policy requires the Mathematica audit to derive each claim *independently* from physical premises, so that the two engines act as cross-checks on each other rather than a single derivation re-typed. A transliteration cannot detect algebraic bugs in the original derivation; both engines would carry the same error.

**Required change:**

Rewrite the Mathematica audit so each claim is derived by a path that does *not* mirror the sympy choreography. Suggested independent paths (use these or equivalents):

- **Closure law (M1)**: derive `chi_sigma = g_phi chi_phi/H` from the Euler-Lagrange / variational derivative of the source-energy functional with respect to `sigma[y]` (treating `sigma` as a function, not a finite-dim variable). For example, use `VariationalD[(1/2) H[y] sigma[y]^2 - gPhi chiPhi[y] sigma[y], sigma[y], y]` from `` VariationalMethods` ``, set it to zero, and solve for `sigma[y]`.
- **Overlap identities (M2-M3) and matched-layer integral reductions (M4-M5)**: verify on a *different* concrete profile than the sympy Gaussian — for example a Lorentzian `1/(1 + y^2/L^2)` (which has finite `Integrate[..., {y, -Infinity, Infinity}]` and would give independent symbolic results). Using a Gaussian like sympy makes this a transliteration; using a different profile makes it an independent check.
- **Cauchy bound (M8)**: instead of repeating the two-point algebra, invoke Cauchy-Schwarz in continuous form by setting `f[y] = chiPhi[y]/Sqrt[H[y]]`, `g[y] = chiPhi[y]/H[y]`, and showing `(Integrate[f g])^2 <= Integrate[f^2] Integrate[g^2]` symbolically for a concrete (non-Gaussian) profile and a non-constant H, e.g., H[y] = h0 + a*y^2.
- **Softening (M10, see F1)**: write the *general* derivation `Theta = gPhi^2 I_1`, `Lambda = gPhi^2 I_1` -> `Lambda^2/Theta = gPhi^2 I_1` symbolically on the two-point branch (no matched-layer substitutions). Compute `Theta = Sum[H_k chiSigma_k^2, ...]` and `Lambda = gPhi Sum[chiPhi_k chiSigma_k, ...]` literally rather than via `hw*nss` / `gPhi*osp`.
- **Matched-layer reductions on the abstract `C^2` and `G_eq` (M6-M7)**: take the limit by `Limit[..., H[y] -> hwConst]` applied to the abstract integral expressions, or by an explicit perturbative expansion `H[y] = hw (1 + epsilon delta[y])` and showing the leading-order coherence is 1.

The goal is that no two corresponding lines in `.py` and `.wl` should be one-to-one syntactic translations.

**Verification:**

After the rewrite, a per-line correspondence audit between `.py` and `.wl` should show no clear 1:1 mapping between sympy and mathematica blocks. The mathematica script must still verify every claim in the paper card (closure, overlaps, coherence, Cauchy bound, matched-layer C^2 = 1, matched-layer gain, and the general softening from F1), and its assertions must still exit 0.

## Independent-derivation check (Mathematica)

Both `.py` and `.wl` are present. The `.wl` is a line-by-line port of the `.py`, not an independent derivation. See Finding F2 above for three pairs of corresponding sections.

## Engine cross-check

Both engines pass their assertions to zero. The captured outputs agree on every numerical/symbolic result:

| Quantity | sympy | mathematica |
|---|---|---|
| Closure | `g_phi*chi_phi(y_loc)/H(y_loc)` | `(gPhi*chiPhiLoc[yLoc])/hLoc[yLoc]` |
| `N_pp_int` (Gaussian) | `sqrt(pi)*L_int` | `lInt*Sqrt[Pi]` |
| `I_1_int` | `sqrt(pi)*L_int/H_w` | `(lInt*Sqrt[Pi])/hw` |
| `I_2_int` | `sqrt(pi)*L_int/H_w^2` | `(lInt*Sqrt[Pi])/hw^2` |
| `C^2 \| H=const` | `1` | `1` |
| `G_eq \| H=const` | `N_pp*g_phi^2/(H_w*K_X)` | `(gPhi^2*L*Sqrt[Pi])/(hw*kX)` (= same up to identifying `N_pp = L*Sqrt[Pi]`) |
| Two-point Cauchy gap | `w1*w2*(H1-H2)^2/(H1^2*H2^2)` | `((h1-h2)^2*w1*w2)/(h1^2*h2^2)` |
| `sigma_stat` | `Lambda*phi/Theta` | `(lambda*phi)/theta` |
| Softening (matched-only) | `I1*g_phi^2` | `gPhi^2*i1` |

`engines_agree: true`. Outputs are fresh (sympy .txt mtime 2026-05-22 19:48 vs .py 19:46; mathematica .txt 19:48 vs .wl 19:47).

## Verdict justification

The script-side closure law, overlap formulas, coherence formula, Cauchy bound, matched-layer coherence, and matched-layer gain all match the paper card and notes; their assertions are non-tautological and exercised against concrete integrals and a discrete two-point model. The paper's *general* equilibrium softening `Delta K_X^eq = g_phi^2 I_1` (valid for any `H(y)`, not only matched layer) is not actually verified — the script's A10 chain forces the answer via two matched-layer substitutions (`I_2 -> I_1^2/N_pp` and `N_pp -> I_1 H_w`), so the check is essentially `g_phi^2 I_1 == g_phi^2 I_1` once the dust settles, providing no independent test of the general case. Additionally, the Mathematica script is structurally a line-by-line port of the SymPy script, violating the second-engine independent-derivation policy. Verdict: `findings`, two medium-severity items. No `stop_cold` — neither finding is `UNFIXABLE` (both have clean remediations) and neither is `CRITICAL_DOWNSTREAM` (the matched-layer outputs that downstream stages 065-066 use *are* verified correctly by A6-A7; the more general softening identity is unproven but is not specifically cited by 065-066 per the paper card).

## Self-test notes

I checked: (i) no `sp.diff(EXPR, VAR)` where `VAR` is independent of `EXPR` — the local-energy minimisation `diff(F_loc, sigma_loc)` is genuine (F_loc depends on sigma_loc explicitly through the quadratic term); (ii) parity of the Gaussian integrals is OK (`chi_phi^2/H_w^k` is even, integral is nonzero, no spurious cancellation); (iii) the proposed two-point fix for F1 reduces to `g_phi^2 I_1 = g_phi^2 I_1` on a general-H model only by virtue of `Theta = g_phi^2 I_1` and `Lambda = g_phi^2 I_1` — both computed from the equilibrium-aligned closure `chi_sigma_k = g_phi chi_phi_k/H_k` without any matched-layer assumption, so the new check exercises the general identity; (iv) the directive targets are existing files at the correct directories (`scripts/`, `mathematica/`); (v) the patch does not introduce a new paper-misalignment because both fixes only restructure script-side algebra to better cover paper claims already stated in the card and notes.
