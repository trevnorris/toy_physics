---
unit_id: 033
batch: II.1
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-04T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage033_microscopic_normalization_equation.md]
  paper_appendix: present
---

# Audit unit 033 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_033.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage033_microscopic_normalization_equation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part02.tex` (row 56 references this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.txt`

## What the paper claims

Stage 033 rewrites the selected-branch 2.5PN normalization test in the microscopic couplings of the first finite-throat kernel model and separates three questions: existence of a stable selected branch, entry into the normalization window, and an exact hit of the universal target. The `\stagefield{Output}` (`/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_033.tex:111`) states: "Stage~033 outputs the coupling-level definitions \eqref{eq:app-stage033-defs}, the stability gate \eqref{eq:app-stage033-stability}, the microscopic target \eqref{eq:app-stage033-micro-target}, the onset condition \eqref{eq:app-stage033-onset}, and the weak-loading expansion \eqref{eq:app-stage033-weak-loading}." Distinct deliverables: (1) the coupling definitions `A`, `Delta_0`, `X`, `beta_0`, `alpha_0`; (2) the three-part stability gate `Delta_0>0`, `A>0`, `alpha_0<alpha_crit` with `alpha_crit=9 pi^2 A(A+ΔK_ax)/(8(11A+9ΔK_ax))`; (3) the microscopic selected target `X^2/Delta_0^2 · s_-^2/(kappa0^2 lambda_-)=N_Q^target`; (4) the onset condition `N_-(0)=beta_0 kappa0^2/A=X^2 kappa0^2/(A Delta_0^2) ≤ N_Q^target`; (5) the weak-loading expansion of `N_-(alpha_0)` (both generic and exact-D/N-constant closed forms). The notes additionally carry the exact monotonicity formula `dN_-/dalpha_0 = beta_0(2 s_- s_-' lambda_- + s_-^3)/(kappa0^2 lambda_-^2)` (`...stage033...md:157-161`) and the onset-stiffness rearrangement `K_0 ≥ g_U^2/Omega_U^2 + kappa0^2 X^2/(N_Q^target Delta_0^2)` (`...stage033...md:185`).

## What the script claims to verify

The SymPy script defines `s_-(alpha0)`, `lambda_-(alpha0)`, `N_- = beta0 s_-^2/(kappa0^2 lambda_-)` from finite-throat constants `kappa0^2=8/pi^2`, `kappa1^2=16/(9 pi^2)`, then asserts: (33.1) the monotonicity formula `dN/dalpha0 = beta0(2 s_- s_-' lambda_- + s_-^3)/(kappa0^2 lambda_-^2)` matches autodiff `diff(N_-,alpha0)`; (33.2) the two `alpha_crit` forms are equal; (33.3) `N_-(0)=beta0 kappa0^2/A`; (33.4) the first weak-loading coefficient equals both the generic and the closed-form expression; (33.5) the microscopic rewrite, the solved `K0_onset`, and its closed-form `gU^2/OmegaU^2 + kappa0^2 X^2/(NQ Delta0^2)`; (33.6) a non-tautological gate check that the substituted `alpha_crit_mic - alpha0_mic` denominator matches the claimed denominator up to a parameter-free constant (`den_ratio.is_number`). The Mathematica script mirrors these and adds (a) a `N_-(0)|_{K0_onset}=NQ` round-trip and (b) a two-point rational numeric cross-check of the monotonicity identity to 30 digits.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Defs `A`,`Delta_0`,`X`,`beta_0`,`alpha_0` (.tex:17-33) | py:81-85 microscopic rewrite, printed | match |
| Stability gate `alpha_crit` closed form (.tex:47-48) | py:57-60 `alpha_crit - alpha_crit_target` | match |
| Substituted coupling gate (.tex:53-56) | py:102-131 / wl:113-131 den_ratio.is_number | match |
| Microscopic target `N_-=...=N_Q^target` (.tex:60-66) | py:44 `Nminus` assembled + 33.5 beta0 sub | match (definition assembled) |
| Onset `N_-(0)=beta0 kappa0^2/A` (.tex:78-80) | py:63-65 `N0 - beta0 kappa0^2/A` | match |
| Onset stiffness `K_0 ≥ ...` (notes:185) | py:87,97-100 `K0_onset` closed form | match |
| Weak-loading generic coef (.tex:99) | py:69-70,73 `coef1 - coef1_target` | match |
| Weak-loading closed coef (.tex:107) | py:71,74 `coef1 - coef1_target_closed` | match |
| Monotonicity (notes:157-161) | py:51-54 `dN - dN_formula`; wl:146-154 numeric | match |

Dominant pattern: aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `expect_zero(dN - dN_formula)` | monotonicity (notes:157-161) | yes |
| A2 | sympy | 60 | `expect_zero(alpha_crit - alpha_crit_target)` | stability gate (.tex:47) | yes |
| A3 | sympy | 65 | `expect_zero(N0 - beta0 kappa0^2/A)` | onset value (.tex:78) | yes |
| A4 | sympy | 73 | `expect_zero(coef1 - coef1_target)` | weak-loading generic (.tex:99) | yes |
| A5 | sympy | 74 | `expect_zero(coef1 - coef1_target_closed)` | weak-loading closed (.tex:107) | yes |
| A6 | sympy | 97-100 | `expect_zero(K0_onset - closed form)` | onset stiffness (notes:185) | yes |
| A7 | sympy | 120-122 | `assert den_ratio.is_number` | substituted gate (.tex:53) | yes |
| A8 | sympy | 125-129 | `expect_zero(... gate_num reconstruction)` | (admittedly tautological-by-reconstruction; labeled) | partial (honest) |
| A9 | math | 60 | `expectZero(dN - dNFormula)` | monotonicity | yes |
| A10 | math | 65 | `expectZero(alphaCrit - alphaCritClosed)` | stability gate | yes |
| A11 | math | 69 | `expectZero(n0 - beta0 kappa0Sq/A)` | onset value | yes |
| A12 | math | 81-82 | `expectZero(coef1 - coef1Target/Closed)` | weak-loading | yes |
| A13 | math | 107 | `expectZero(n0Mic|_{K0_onset} - NQ)` | onset target round-trip | yes |
| A14 | math | 108-111 | `expectZero(k0Onset - closed form)` | onset stiffness | yes |
| A15 | math | 122-125 | `If[!NumericQ[denRatio], fail]` | substituted gate | yes |
| A16 | math | 146-154 | numeric `Abs[...]>10^-20 → fail` (×2) | monotonicity (independent route) | yes |

A8 is honestly flagged in-script as "tautological by reconstruction; substantive check is den_ratio.is_number above" — it is not a hidden tautology; the load-bearing gate check is A7/A15.

## Findings

### F1 — stale_output

**Severity:** low (informational, non-blocking)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage033_microscopic_normalization_equation_sympy_audit.txt` (mtime 2026-05-26 00:44:33)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage033_microscopic_normalization_equation_mathematica_audit.txt` (mtime 2026-05-26 00:44:46)
- vs scripts `.py`/`.wl` (both mtime 2026-06-03 15:59:11)

**What's wrong:**
Both saved-output transcripts are older than their scripts. I confirmed via `git diff c6ee7d7 e2a4780` that the only change in the Jun-03 commit (e2a4780, "numbering reconciliation Phase 1") to these two scripts was the banner strings: `STAGE 16.x` → `STAGE 33.x` (sympy) and `STAGE 016` → `STAGE 033` (mathematica). No mathematical line changed. The committed outputs accordingly still print the stale banners — sympy output line 3 `STAGE 016 — ...`, lines 3/30/36/etc.; mathematica output lines 3/30/etc. `STAGE 16.1 ... 16.6`. Every numeric/symbolic RESULT in both transcripts matches what the current scripts would produce.

**Why this matters:**
Cosmetic only: the banner labels in the captured `.txt` no longer match the scripts. No result is affected. The orchestrator's standard fresh re-run regenerates the transcripts with the corrected banners.

**Required change:**
None to script source. Re-run both scripts to refresh the committed `.txt` (orchestrator already does this on verify). No Codex directive needed (no script-source defect).

**Verification:**
After a fresh run, sympy output line 3 reads `STAGE 33.1 ...` and mathematica output line 3 reads `STAGE 033 — ...`; all PASS lines remain and residuals stay 0.

## Independent-derivation check (Mathematica)

The analytic portion of the `.wl` is a close port of the `.py` (same `lambdaMinus`/`sMinus`/`nMinus` construction, same two-form `alphaCrit` check, same `Solve[n0Mic==NQ,K0]`, same `denRatio` gate logic). Examples of 1:1 correspondence: py:42-44 `lambda_minus/s_minus/Nminus` ↔ wl:43-48; py:57-60 `alpha_crit` two-form ↔ wl:62-65; py:104-122 gate `den_ratio` ↔ wl:114-125. This is borderline-transliteration for the analytic checks. However the `.wl` carries two genuinely independent verification routes the `.py` lacks: (1) `n0Mic|_{K0->k0Onset} - NQ == 0` round-trip (wl:107) confirming the solved onset stiffness actually drives `N_-(0)` to the target; and (2) a two-point rational numeric cross-check (wl:138-154) verifying `dN - dNFormula = 0` to 30 digits at `numericRule1`/`numericRule2`, which is structurally distinct from the analytic `FullSimplify` path and would catch a shared algebra error. The script's own comment (wl:133-137) states exactly this intent. Because of the independent numeric arbiter, I do not raise `mathematica_transliteration` (it would be a likely false positive for this stage), but I note the analytic checks are closely parallel.

## Engine cross-check

Both engines agree on every closed form (cross-checked against the committed transcripts):
- `lambda_-`: sympy out l5 ↔ math out l5 — identical.
- `s_-(alpha0)`: sympy out l6 ↔ math out l6 — identical.
- `N_-(alpha0)`: sympy out l7 ↔ math out l5-26 (pretty-printed, same expression).
- `alpha_crit = (9 A(A+DeltaK) Pi^2)/(88A+72DeltaK)`: math out l10; equals sympy `9 pi^2 A(A+DK)/(8(11A+9DK))` since `8(11A+9DK)=88A+72DK`.
- `N_-(0)=8 beta0/(pi^2 A)`: math out l13, sympy out l38 — identical.
- weak coef `64 beta0(8A+9DK)/(9 A^2 DK pi^4)`: math out l16, sympy out l44 — identical.
- `K0_onset`: math out l27 (expanded), sympy out l57; both pass the closed-form identity assertion (`= 0`).
- gate `den_ratio`: math out l34 `-9 Pi^2`, sympy out l132 `9 pi^2`. The sign difference is purely a `fraction()` (SymPy) vs `Numerator/Denominator` (WL) sign-convention artifact on the rational `gate_diff`; both are `is_number`/`NumericQ` and both PASS. Not a math disagreement — `gate_diff = num/den` is the same rational either way.
- monotonicity: both engines `= 0` (analytic) and the `.wl` numeric residuals are `0` (`0\`\`78.83...` is WL high-precision zero) at both test points.

Engines agree.

## Verdict justification

Attacks attempted and failed: (1) hand-derived `alpha_crit` closed form from generic — matches exactly (`8(11A+9ΔK)/(9pi^2)` denominator). (2) Hand-derived the weak-loading closed coefficient and `N_-(0)` — match. (3) Checked the monotonicity identity for hidden tautology — it is NON-tautological: `dN` is autodiff, `dN_formula` is the hand form, and equality requires the spectral relation `dlambda_-/dalpha_0 = -s_-`, which I verified holds exactly for the script's `s_-`/`lambda_-` definitions; the check genuinely exercises it. (4) Scrutinized the gate check (py:113-129 / wl:113-131) — the load-bearing assertion `den_ratio.is_number` is built from `together/cancel(alpha_crit_mic - alpha0_mic)` WITHOUT reference to the claimed denominator, so it can fail if the claimed gate form is wrong; the subsequent `expect_zero` is honestly labeled tautological-by-reconstruction. (5) Verified `Delta_0` uses the same `sigma=88/(9pi^2)` the card names. (6) Confirmed paper↔script symbol mapping: script `DeltaK` is the card's `ΔK_ax`; `NQ` is `N_Q^target`. I read the card, the notes, and the appendix row; every stage-033 deliverable maps to a substantive, non-tautological script check, and both engines agree. The only finding is a cosmetic `stale_output` (banner-label-only drift from the Jun-03 numbering commit; no math changed), which is informational and self-heals on re-run. Verdict: `findings` (one informational stale_output), substance is clean; no `paper_misalignment`, no `stop_cold`.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 16 deliverable values checked, 0 misaligned.

I enumerated every RESULT/deliverable value the scripts emit (from `.py`/`.wl` source + committed `.txt` outputs; nothing executed) and located each in the `.tex` card and `.md` notes. All reconcile. (Outputs are stale only in banner text — see F1; every result value in the transcripts matches the current scripts, so the reconciliation is based on script source confirmed against the committed numeric/symbolic output.)

| value | source (py / wl + output line) | .tex / .md location | status |
|---|---|---|---|
| `A = K0 - gU^2/OmegaU^2` | py:81; sympy.txt:21 / wl.txt:51 | .tex:21 (eq-defs); .md:44 | MATCH |
| `Delta0 = OmegaU^2 OmegaW^2 - 88 gR^2/(9pi^2)` (=`Omega_U^2 Omega_W^2 - g_R^2 sigma`) | py:82; sympy.txt:22 / wl.txt:52 | .tex:21; .md:46 | MATCH |
| `Chi = OmegaU^2 gW + gR gU` | py:83; sympy.txt:23 / wl.txt:53 | .tex:24; .md:48 | MATCH |
| `beta0 = Chi^2/Delta0^2` | py:84; sympy.txt:24 / wl.txt:54 | .tex:30; .md:54 | MATCH |
| `alpha0 = gB^2/varpi^2 + Chi^2/(OmegaU^2 Delta0)` | py:85; sympy.txt:25 / wl.txt:55 | .tex:32; .md:56 | MATCH |
| `alpha_crit = 9 pi^2 A(A+DK)/(8(11A+9DK))` (=`(9A(A+DK)Pi^2)/(88A+72DK)`) | py:58; sympy.txt:10 / wl.txt:32 | .tex:47-48; .md:141 | MATCH |
| Microscopic target `N_- = X^2/Delta0^2 · s_-^2/(kappa0^2 lambda_-) = N_Q^target` | py:44+33.5; sympy.txt:7 | .tex:60-66; .md:77-95 | MATCH (definition assembled) |
| `N_-(0) = 8 beta0/(pi^2 A)` (=`beta0 kappa0^2/A`) | py:63-65; sympy.txt:13 / wl.txt:38 | .tex:78-80; .md:173 | MATCH |
| `N_-(0)` micro = `8 Chi^2/((K0-gU^2/OmegaU^2)(...)^2)` (=`X^2 kappa0^2/(A Delta0^2)`) | py:86; sympy.txt:26 / wl.txt:56 | .tex:79-80; .md:173 | MATCH |
| weak-loading generic coef `beta0 kappa0^2(4A kappa1^2 + DK kappa0^2)/(A^2 DK)` | py:70; (asserted, sympy.txt:17) | .tex:99; .md:202 | MATCH |
| weak-loading closed coef `64 beta0(8A+9DK)/(9 pi^4 A^2 DK)` | py:71; sympy.txt:16 / wl.txt:44 | .tex:107; .md:209-211 | MATCH |
| `N_-(alpha0)` to O(alpha0) series | py:68; sympy.txt:44 (math) / sympy.txt series | .tex:103-109; .md:207-211 | MATCH |
| `K0_onset = gU^2/OmegaU^2 + kappa0^2 Chi^2/(NQ Delta0^2)` | py:87,97-100; sympy.txt:27 / wl.txt:57 | .md:185 (notes carries it; .tex onset is the `X^2 ≤ ...` form eq:84-89) | MATCH |
| substituted stability gate inequality (`gB^2/varpi^2 + ... < alpha_crit(mic)`) | py:102-131; wl:113-131 | .tex:51-56; .md:145-146 | MATCH |
| onset condition `N_-(0) ≤ N_Q^target` ⇔ `X^2 ≤ N_Q^target A Delta0^2/kappa0^2` | implied by A6 / wl:107 round-trip | .tex:84-89; .md:177-181 | MATCH |
| monotonicity formula `dN/dalpha0 = beta0(2 s_- s_-' lambda_- + s_-^3)/(kappa0^2 lambda_-^2)` | py:51-54; wl:54-60,146-154 | .md:157-161 | MATCH |

INTERNAL scaffolding (accounted for, no finding): `kappa0_sq=8/pi^2`, `kappa1_sq=16/(9pi^2)`, `sigma=88/(9pi^2)`, `delta_kappa`, `Kprod`, `R`, `lambda_-(alpha0)` general closed form, `s_-(alpha0)` general closed form, `gate_den_claim`, `gate_num_target`, computed/claimed denominator expansions, `den_ratio (=±9 pi^2)`, numeric residuals (`=0`), monotonicity-formula residual (`=0`). (`lambda_-` and `s_-` general forms are upstream spectral data carried in from Stages 031/032; the card only reports their `alpha0=0` values `s_-(0)=kappa0^2`, `lambda_-(0)=A` at .tex:71-74, so the general forms are legitimately intermediate, not stage-033 prose deliverables.)

## Self-test notes

Checked: (1) variable-independence — `diff(Nminus, alpha0)` is well-posed (`Nminus` genuinely depends on `alpha0` via `s_-`/`lambda_-`/`R`); the `subs(alpha0,0)` checks are evaluations, not zero-derivatives. (2) The monotonicity `assert_zero` is non-trivial because it requires `lambda_-' = -s_-`, which I verified holds exactly. (3) Trivial-case: `N_-(0)` reduces to `8 beta0/(pi^2 A)` as asserted; `alpha_crit` algebra reduces as hand-derived. (4) No missing-script findings (both engines present), so no path-spec trap. (5) The lone finding (stale_output) prescribes no script-source edit and introduces no new paper_misalignment. Conclusion: substance clean; only a cosmetic banner-label staleness that self-heals on re-run.
