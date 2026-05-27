---
unit_id: 077
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage077_family1_theta_extraction.md
  paper_appendix: present
---

# Audit unit 077 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_077.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage077_family1_theta_extraction.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (line 272: `\input{stages/stage_077}` — no inline appendix prose for this stage; the part-level file is a stage-include table of contents)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.txt`

## What the paper claims

The stage card states the Purpose as "Stage~077 evaluates the natural wall-density weighting and a conservative lower envelope", with Inputs "The explicit radial wall profile and the canonical support weighting", and Output "Natural and conservative Family-1 wall-depth data \eqref{eq:app-stage077-theta-chi}--\eqref{eq:app-stage077-theta-J}". The two boxed deliverables are
- `\Theta_w^{(\chi)} = 25 \lambda_\mu^2 \langle\rho_r^2\rangle_\chi \simeq 4.06863235008162 \lambda_\mu^2` (natural shell-weighted datum)
- `\Theta_w^{(J)}  = 25 \lambda_\mu^2 \langle\rho_r\rangle_\chi^2 \simeq 0.927552032539308 \lambda_\mu^2` (conservative Jensen lower envelope)

together with the moments `\langle\rho_r\rangle_\chi \simeq 0.192619005556493` and `\langle\rho_r^2\rangle_\chi \simeq 0.162745294003265`. The notes carry the additional binding inputs: the Family-1 branch uses `alpha_r = 10`, the radial wall profile `rho_r(xi) = [1 - alpha_r S(xi)^2]_+^{1/4}` with `S(xi) = (1+tanh xi)/2`, the canonical support weight `chi_phi(xi) = S'(xi) = sech^2(xi)/2`, the exact cut point `xi_* = artanh(2/sqrt(alpha_r) - 1) ≈ -0.3855810692` (for `alpha_r=10`), and the support normalization `I_f = int chi^2 dxi = 1/3`. The coefficient `25` is carried from Stage 76 via `\eqref{eq:app-stage076-theta}`.

## What the script claims to verify

The SymPy script verifies (i) the symbolic identity `I_f = 1/3` for the canonical support normalization, (ii) the symbolic cut-point identity `1 - alpha_r * S(xi_*)^2 = 0`, (iii) the numerical shell-weighted moments `<rho>_chi` and `<rho^2>_chi` at `alpha_r = 10` against 50-digit literals, (iv) the numerical wall-depth data `Theta_w^(chi) = 25 <rho^2>_chi` and `Theta_w^(J) = 25 <rho>_chi^2` against 50-digit literals, and (v) the Jensen ordering `Theta_w^(chi) >= Theta_w^(J) > 0`. The Mathematica script independently re-derives the same identities using a different integration choreography (`NIntegrate` with `WorkingPrecision -> 60` and a distinct `rhoSqNum = Sqrt[1 - alphaR S^2]` integrand rather than squaring a fourth root) and enforces the same numeric targets to 28+ digits via `expectApprox` plus the same Jensen ordering. The docstrings and assertions are anchored to the same physics the paper card describes.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `<rho_r>_chi ≈ 0.192619005556493` | sympy `expect_close("<rho>_chi", R1, 0.19261900555649309777..., 1e-28)` and mathematica `expectApprox["<rho>_chi numeric check", ...]` | match |
| `<rho_r^2>_chi ≈ 0.162745294003265` | sympy `expect_close("<rho^2>_chi", R2, 0.16274529400326462037..., 1e-28)` and mathematica `expectApprox["<rho^2>_chi numeric check", ...]` | match |
| `Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2` (boxed) | sympy `expect_close("Theta_w^(chi)", Theta_chi, 4.0686323500816155092..., 1e-26)` and mathematica analogue | match |
| `Theta_w^(J) ≈ 0.927552032539308 lambda_mu^2` (boxed) | sympy `expect_close("Theta_w^(J)", Theta_J, 0.92755203253930797183..., 1e-27)` and mathematica analogue | match |
| Inequality `<rho_r^2>_chi >= <rho_r>_chi^2` ⇒ `Theta_w^(chi) >= Theta_w^(J)` (notes §5) | sympy `if not (Theta_chi >= Theta_J > 0)` and mathematica `expectTrue["Theta_w^(chi) >= Theta_w^(J) > 0", ...]` | match |
| Input: `alpha_r = 10` (notes §1) | both scripts use `alpha_num = mp.mpf('10')` / `alphaNum = 10` | match |
| Input: `S(xi) = (1+tanh xi)/2` (notes §1) | both scripts define `S = (1+tanh(xi))/2` and use it consistently | match |
| Input: `chi_phi = sech^2 xi / 2` (notes §2) | sympy `chi = sp.diff(S, xi)`; mathematica `D[sXi, xi]` | match |
| Input: `xi_* = artanh(2/sqrt(alpha_r) - 1)` (notes §1) | sympy `xi_star = sp.simplify(sp.atanh(2/sp.sqrt(alpha_r) - 1))`; mathematica `xiStar = FullSimplify[ArcTanh[2/Sqrt[alphaR] - 1], ...]` with symbolic cut-point identity then enforced | match |
| Input: `I_f = 1/3` (notes §2) | sympy `expect_zero("I_f - 1/3", If - 1/3)` and mathematica `expectZero["I_f - 1/3", ifMom - 1/3]` | match |
| Coefficient `25` (notes §3, paper card via Stage 76) | both scripts use literal `25` in `Theta_chi = 25*R2` and `Theta_J = 25*R1^2` | match (carried Input from Stage 076; not derived in 077 by design) |

Every paper-side deliverable has a corresponding script-side check, and every script-side load-bearing assertion exercises a paper-side claim. `paper_alignment: aligned`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 39 | `expect_zero("I_f - 1/3", If - sp.Rational(1,3))` | Input: `I_f = 1/3` | yes |
| A2 | sympy | 47 | `expect_zero("1 - alpha_r * S(xi_*)**2", rho_quartic_at_star)` | Input: cut point `xi_*` formula | yes |
| A3 | sympy | 88-93 | `expect_close("<rho>_chi", R1, 0.19261900555649309777..., 1e-28)` | `<rho_r>_chi ≈ 0.192619005556493` | yes |
| A4 | sympy | 94-99 | `expect_close("<rho^2>_chi", R2, 0.16274529400326462037..., 1e-28)` | `<rho_r^2>_chi ≈ 0.162745294003265` | yes |
| A5 | sympy | 100-105 | `expect_close("Theta_w^(chi)", Theta_chi, 4.0686323500816155092..., 1e-26)` | boxed `Theta_w^(chi)` | yes |
| A6 | sympy | 106-111 | `expect_close("Theta_w^(J)", Theta_J, 0.92755203253930797183..., 1e-27)` | boxed `Theta_w^(J)` | yes |
| A7 | sympy | 112-113 | `if not (Theta_chi >= Theta_J > 0): raise` | Jensen ordering (notes §5) | yes |
| A8 | mathematica | 50 | `expectZero["I_f - 1/3", ifMom - 1/3]` | Input: `I_f = 1/3` | yes |
| A9 | mathematica | 55 | `expectZero["1 - alphaR*S[xi_*]^2", rhoQuarticAtStar]` | Input: cut point `xi_*` formula | yes |
| A10 | mathematica | 92 | `expectApprox["<rho>_chi numeric check", r1, 0.19261900555649309777..., 10^-28]` | `<rho_r>_chi` value | yes |
| A11 | mathematica | 93 | `expectApprox["<rho^2>_chi numeric check", r2, 0.16274529400326462037..., 10^-28]` | `<rho_r^2>_chi` value | yes |
| A12 | mathematica | 94 | `expectApprox["Theta_w^(chi) numeric check", thetaChi, 4.0686323500816155092..., 10^-26]` | boxed `Theta_w^(chi)` | yes |
| A13 | mathematica | 95 | `expectApprox["Theta_w^(J) numeric check", thetaJ, 0.92755203253930797183..., 10^-27]` | boxed `Theta_w^(J)` | yes |
| A14 | mathematica | 96 | `expectTrue["Theta_w^(chi) >= Theta_w^(J) > 0", thetaChi >= thetaJ && thetaJ > 0]` | Jensen ordering (notes §5) | yes |

Every assertion is anchored to a paper-side claim. No tautological constructions remain. No assertions test claims absent from the paper.

## Findings

None.

## Independent-derivation check (Mathematica)

The two scripts are not transliterations. Three structural differences confirm independent derivation:

1. **Integration choreography differs.** SymPy uses `mp.quad` with explicit breakpoints `[-mp.inf, -4, xi_cut]` to assist adaptive quadrature (lines 66-68). Mathematica uses `NIntegrate` from `-Infinity` to `xiCut` with `WorkingPrecision -> 60` and `Quiet[..., NIntegrate::precw]` (lines 67-78) — no manual breakpoints.

2. **Integrand for the second moment differs.** SymPy builds `rho_num(x)` as the fourth root and then squares it inside the integrand for `<rho^2>` (lines 59-61, 68). Mathematica avoids the square-then-root by defining `rhoSqNum[x] := Sqrt[1 - alphaNum*sNum[x]^2]` directly (line 65) — a numerically distinct evaluation path.

3. **Positivity / clip handling.** SymPy explicitly clips `val**0.25 if val > 0 else 0` (line 61); Mathematica does not clip but relies on the integration bound being `xiCut = xi_*`, where the bracket is non-negative on the integration domain.

The fact that all four numerical results agree to 50 digits across these two genuinely different code paths is a strong cross-engine check, not a transliteration artifact.

## Engine cross-check

Numerical agreement at the printed precision (extracted verbatim from the saved output transcripts):

```
                              SymPy mp.quad (dps=50)                                    Mathematica NIntegrate (WP=60)
<rho>_chi    = 0.19261900555649309777068139356018510792903510747507     0.19261900555649309777068139356018510792903510747506717457...
<rho^2>_chi  = 0.16274529400326462037087418498629868328210821103971     0.16274529400326462037087418498629868328210821103971427483...
denominator  = 0.33333333333333333333333333333333333333333333333333     0.33333333333333333333333333333333333333333333333333268403...
Theta_chi    = 4.0686323500816155092718546246574670820527052759928      4.06863235008161550927185462465746708205270527599285687082...
Theta_J      = 0.92755203253930797183993260663904217023332624032789     0.92755203253930797183993260663904217023332624032789843141...
xi_*(10)     = -0.38558106921542562403635498846713378847348301441599    -0.38558106921542562403635498846713378847348301441599100022...
```

Symbolic results also agree: `chi_phi(xi) = sech^2(xi)/2`, `I_f = 1/3`, `xi_* = -atanh(1 - 2/sqrt(alpha_r))` (SymPy form) = `-(1/2)*Log[-1 + Sqrt[alphaR]]` (Mathematica form). Both engines verify `1 - alpha_r * S(xi_*)^2 = 0` symbolically. `engines_agree: true`.

## Output freshness

Script and output mtimes (epoch seconds):
- sympy script: 1779513885; sympy output: 1779513990 — output fresher than script by 105 s.
- mathematica script: 1779513885; mathematica output: 1779514003 — output fresher than script by 118 s.

`outputs_fresh: true`.

## Verdict justification

The v2 paper-grounded audit confirms full alignment between the paper card, the source notes, and the two engines' scripts. Every paper-side deliverable (both boxed Theta_w results, both moments, the Jensen ordering, and the carried inputs from the notes — alpha_r=10, S, chi, xi_*, I_f) has a corresponding non-tautological script-side assertion in both engines. The two engines derive the results via genuinely different integration paths and agree to 50 digits. All three v1 findings (missing symbolic cut-point identity in both engines, the tautological xi_* numeric self-check in Mathematica, and missing SymPy `expect_close` assertions on the numerical results) have been fully resolved by the prior Codex pass and are visible in the current scripts and output transcripts. Adversarial attacks tried and failed: (a) checked whether `alpha_r=10` is actually exercised — yes, in both numerical extraction blocks; (b) checked whether the coefficient `25` matches the paper (carried from Stage 76 by design and consistent); (c) checked whether the cut-point identity holds symbolically under the declared assumptions (`alpha_r positive` is sufficient because `tanh(atanh(z))=z` is formal); (d) checked whether the moment integrand handles the positive-part clip correctly (SymPy clips explicitly; Mathematica relies on the integration upper bound = xi_* keeping the bracket non-negative — both correct); (e) checked whether the assertion tolerances (1e-28, 1e-27, 1e-26) are tight enough to catch real numerical drift (residuals are ~1e-50, well below tolerance, so a real drift would be caught). Verdict: `clean`.

## Self-test notes

Variable-independence: no `sp.diff` or `D` traps — `chi = sp.diff(S, xi)` and `D[sXi, xi]` both involve real `xi`-dependence in `S`. Parity: the `I_f` integrand `sech^4(xi)/4` is even on a symmetric domain, integral is nonzero (1/3) — assertion has real content. Trivial-case: at `alpha_r → ∞`, `xi_* → atanh(-1) → -∞`, so the upper integration bound retreats and `<rho>_chi → 0` — limiting behavior consistent with the wall vanishing. Cut-point identity reduces to `1 - alpha_r * (1/sqrt(alpha_r))^2 = 0` algebraically, confirmed in the output. Numeric targets in the four `expect_close` calls match the paper's quoted decimals to all printed digits. No new paper_misalignment introduced. No directive written (zero findings).
