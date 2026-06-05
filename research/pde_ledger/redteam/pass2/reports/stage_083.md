---
unit_id: 083
batch: III.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage083_family1_direct_operator_window.md]
  paper_appendix: present
---

# Audit unit 083 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_083.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage083_family1_direct_operator_window.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 144)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.txt`

## What the paper claims

Stage 083 evaluates the Family-1 support window *directly* through the operator-selected fixed-point branch (rather than only through endpoint bounds). The card's `\stagefield{Output}` is: **"Direct operator-selected support window for Family--1."** The body equation (eq:app-stage083-zeta-direct) defines `zeta_F1^op = zeta_phys(Pe_*(Xi_F1,37,12321/5),37;12321/5)`, and the card explicitly stresses "This calculation is a numerical placement of exact functions, not an independent fit." The notes expand this into four stated deliverables: (1) the operator-selected transport bias `Pe_*` is monotone in wall/source strength `Xi` (proved via implicit differentiation `dPe_*/dXi = Delta/[1 - Xi ∂_Pe Delta] > 0`); (2) the explicit Family-1 support ratio and branch-product thresholds are therefore monotone in `Xi`; (3) inserting the natural shell-weighted (`chi`) and Jensen-floor (`J`) wall data reproduces the Stage-078/080/081 transport, zeta, and Pi/C_mix windows directly from the coupled operator; (4) the natural branch already sits essentially at the hard Family-1 ceiling `zeta_max^(F1) ≈ 2.46752922945601`. The notes quote concrete numbers for every window (`Delta_0`, `Delta_inf`, `Xi^(chi)`, `Xi^(J)`, the four `Pe`, four `zeta`, four `Pi/C_mix`, and the ceiling).

## What the script claims to verify

The docstring (py lines 4-11) lists three checks: (1) the exact implicit-differentiation formula for `dPe_*/dXi`; (2) the exact Family-1 support/source constants `Delta_0`, `Delta_inf`; (3) the direct operator-selected `Pe`, `zeta`, and `Pi/C_mix` windows for the `chi` and `J` branches. The load-bearing *assertions* are: the implicit-derivative residual (py:47 / wl:49), two "defining-equation residuals" for `Delta_0`/`Delta_inf` (py:75,81 / wl:81,91), the `y_F1` eigenvalue root residual (py:87 / wl:101), an `Omega(Pe)` "identity residual" (wl:154 only), a monotonicity spot-check on `dzeta/dPe` at four `Pe` values (py:144-150 / wl:186-194), and — Mathematica only — a battery of numeric `expectApprox` checks pinning `Delta_0`, `Delta_inf`, the four `Pe`, four `zeta`, and the ceiling to external literals (wl:169-179). The SymPy script asserts only the implicit-derivative, the two `Delta` residuals, the `y_F1` residual, and monotonicity; every window value (`zeta`, `Pi/C_mix`, the four `Pe`) is **printed but never asserted** in SymPy.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (1) `dPe_*/dXi = Delta/(1 - Xi ∂_Pe Delta) > 0` (monotone in Xi) | implicit-derivative residual py:47 / wl:49; monotonicity spot-check py:144 / wl:186 | match (closed-form identity proven; positivity is a 4-point spot-check, honestly labelled) |
| (2) support-ratio / branch-product thresholds monotone in Xi | `dzeta/dPe>0` spot-check py:144 / wl:186 | partial (4 sample points, not a global proof; comment is honest) |
| `Delta_0(F1) ≈ 1.733e-4`, `Delta_inf(F1) ≈ 2.014e-2` | residual py:75,81 / wl:81,91 (tautological); numeric `expectApprox` wl:169-170 (real anchor, Mathematica only) | partial — see F1 |
| `Xi^(chi) ≈ 556995.77`, `Xi^(J) ≈ 126981.87` | printed py:111-112 / wl:124-125 | match (value), not asserted |
| four `Pe` windows (96.53 / 11220.54 / 22.01 / 2558.02) | printed py:119-122; numeric `expectApprox` wl:171-174 | match; SymPy print-only — see F2 |
| four `zeta` windows + ceiling 2.46752922945601 | printed py:152-156; numeric `expectApprox` wl:175-179 | match; SymPy print-only — see F2 |
| four `Pi/C_mix` windows + ceiling 3.46752922945601 | printed py:158-162 / wl:196-200 (= 1 + zeta) | match; never asserted in either engine |

All numeric deliverables agree exactly with the notes (see Value Reconciliation). The paper alignment is **aligned** on values; the findings below are about *verification quality*, not value mismatch.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 47 | `expect_zero(dPe_formula - dPe_expected)` | claim 1 (monotonicity formula) | yes — solves the implicit eq independently, compares to closed form |
| A2 | sympy | 75 | `expect_zero(denom·Delta0 - numer)` | `Delta_0` constant | **no — tautological** (Delta0 := numer/denom) |
| A3 | sympy | 81 | `expect_zero(denom·DeltaInf - numer)` | `Delta_inf` constant | **no — tautological** |
| A4 | sympy | 87 | `raise if |y·tan y - eta| > 1e-25` | `y_F1` eigenvalue | yes — genuine root residual of nsolve |
| A5 | sympy | 144-150 | `raise if dzeta/dPe ≤ 0` at 4 Pe | claim 2 monotonicity | partial — 4-point spot-check |
| (zeta/Pi/Pe windows) | sympy | 111-162 | `print` only | claims 3,4 (the Output) | **no assertion** — see F2 |
| B1 | mathematica | 49 | `expectZero(dpeFormula - dpeExpected)` | claim 1 | yes |
| B2 | mathematica | 81 | `expectZero(denom·delta0 - numer)` | `Delta_0` | **no — tautological** (same defect as A2) |
| B3 | mathematica | 91 | `expectZero(denom·deltaInf - numer)` | `Delta_inf` | **no — tautological** |
| B4 | mathematica | 101 | `expectApprox(y·Tan y - eta, 0, 1e-20)` | `y_F1` | yes |
| B5 | mathematica | 103 | `expectApprox(aF1 - aF1Indep, 0, 1e-30)` | `A_F1` | **no — tautological** (aF1Indep = (kappa+(Pi/2)^2)/(kappa+yRoot^2) = aF1; `Pi^2/4 ≡ (Pi/2)^2`) |
| B6 | mathematica | 154 | `expectZero(omega·denom - numer)` | `Omega(Pe)` def | **no — tautological** (omega := numer/denom; comment admits "this is the definition we typed in") |
| B7 | mathematica | 169-179 | `expectApprox` vs external literals | `Delta_0`,`Delta_inf`,4 Pe,4 zeta,ceiling | yes — real numeric anchor (only such anchor in either engine) |
| B8 | mathematica | 186-194 | `Sign[dzeta/dPe]===1` at 4 Pe | claim 2 monotonicity | partial — 4-point spot-check |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:57-81`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:55-103`

**What's wrong:**
The "defining-equation residual" checks for `Delta_0` and `Delta_inf` are tautological by construction, and the inline comments falsely advertise them as independent verifications.

SymPy defines (py:57-60):
```
Delta0_F1 = sp.simplify(eta_F1*(cosh(alpha)-1) / (alpha**2*(alpha*sinh(alpha)+eta*cosh(alpha))))
```
then checks (py:71-75):
```
delta0_residual = (alpha**2*(alpha*sinh(alpha)+eta*cosh(alpha))) * Delta0_F1 - eta*(cosh(alpha)-1)
expect_zero("Delta_0(F1) defining-equation residual", delta0_residual)
```
Since `Delta0_F1 = numer/denom`, the residual is `denom·(numer/denom) - numer ≡ numer - numer ≡ 0`, identically, for *any* closed form typed into `Delta0_F1`. A Cosh↔Sinh paste error in the closed form would appear identically on both sides and still cancel to 0. The comment at py:68-70 asserts the opposite: *"These are non-tautological because they express Delta_* as the unique solution of a linear equation, not as a chosen closed form."* That rationale is incorrect as implemented — the script never solves a linear equation for `Delta_*`; it pins `Delta_*` to a chosen closed form and then multiplies it back by its own denominator.

The Mathematica `.wl` is worse on the comment side: lines 62-71 advertise an *"Independent BVP derivation"* of the two-point boundary value problem `u'' - alpha² u = 0` with a Robin wall condition, but the code at wl:77-91 performs the **identical** tautological residual (`denom·delta0F1 - numer`) and never sets up or solves any BVP (`DSolve`/`NDSolve` appear nowhere). Two further `.wl` checks are tautological for the same reason: `A_F1 independent vs closed-form` (wl:102-103) compares `aF1 = (kappa+Pi^2/4)/(kappa+yRoot^2)` against `aF1Indep = (kappa+(Pi/2)^2)/(kappa+yRoot^2)`, but `Pi^2/4 ≡ (Pi/2)^2`, so the diff is identically 0 — it is the same expression typed twice; and `Omega(Pe) identity residual` (wl:137-154) is `omega·denom - numer` where `omega := numer/denom` (the comment at wl:142-146 at least admits "It is mathematically the definition we typed in").

The one check that genuinely anchors the closed forms is the Mathematica numeric battery (wl:169-170, `expectApprox` vs external literals `1.73302079021525e-4` and `2.01447565540522e-2`). The **SymPy** script has no such anchor — its only checks on `Delta_0`/`Delta_inf` are the tautological residuals (py:75,81); it merely `print`s the numeric values (py:92-93) without asserting them.

**Why this matters:**
The card's `\stagefield{Inputs}` carries `Delta` and `Pe_*` forward, and these constants set every downstream window. If a future edit perturbs the closed-form Cosh/Sinh expressions in the SymPy script, *no SymPy assertion would catch it* — the tautological residual passes regardless, and there is no numeric anchor. The comments actively mislead a future maintainer into believing these closed forms are independently verified.

**Required change:**
- SymPy (py:71-93): replace the two tautological `expect_zero` residuals with non-tautological checks. Either (a) add numeric anchors `expect_close(sp.N(Delta0_F1,30), 1.73302079021525e-4, atol=1e-16)` and the analogous `2.01447565540522e-2` for `Delta_inf` (matching the literals the notes report and the `.wl` already pins), AND/OR (b) re-derive `Delta_0`/`Delta_inf` independently as boundary values of the linear ODE `u'' - alpha² u = 0` with the stated Robin/edge conditions via `sp.dsolve`, then `expect_zero(Delta0_indep - Delta0_F1)`. Correct the comment at py:68-70 — it currently claims non-tautology that does not hold for the residual form.
- Mathematica (wl:62-103,137-154): EITHER actually perform the advertised BVP derivation (solve `u'' - alphaF1² u == 0` with the Robin wall + unit-edge conditions via `DSolve`, extract the boundary values, and `expectZero[deltaIndep - delta0F1]`), OR delete the misleading "Independent BVP derivation" comment block (wl:62-71) and the tautological residual checks, relying on the real numeric anchors at wl:169-170. Likewise either remove or correct the false "Independent" framing on `A_F1` (wl:96-103) and `Omega` (wl:139-146): the `A_F1` "independent" form is byte-identical (`Pi^2/4 = (Pi/2)^2`) and provides zero coverage.

**Verification:**
After the fix, the SymPy script must contain at least one non-tautological assertion that would fail if the `Delta_0`/`Delta_inf` closed form were perturbed (a numeric `expect_close` against the external literal, or a `dsolve`-based independent derivation). The `.wl` either gains a real `DSolve` BVP check or sheds the misleading comment block. Output transcripts should still show all checks PASS, and a deliberate sign-flip test (mentally: swap `cosh`→`sinh` in `Delta0_F1`) should now flip a check from PASS to FAIL.

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:111-162`

**What's wrong:**
The stage's `\stagefield{Output}` is "Direct operator-selected support window for Family--1" — i.e. the four `Pe`, four `zeta`, four `Pi/C_mix` windows and the ceiling. In the SymPy script every one of these load-bearing deliverables is computed and **`print`ed but never asserted** (py:111-122 for `Xi`/`Pe`, py:152-162 for `zeta`/`Pi`). The only SymPy assertions are the implicit-derivative (A1, genuine), the two tautological `Delta` residuals (F1), the `y_F1` root residual (A4, genuine), and monotonicity (A5, a 4-point spot-check). Consequently the SymPy engine does **not** verify the stage's actual Output at all — it verifies two upstream ingredients (`dPe/dXi`, `y_F1`) plus a monotonicity sample, and prints the rest. The Mathematica `.wl` does pin the windows numerically (wl:171-179), so cross-engine the deliverables are anchored once; but the SymPy "audit" advertised in `\stagefield{Verification}` is, on its own, an unasserted print of the answer. Combined with F1, the SymPy script contains no assertion that exercises the stage's headline claim.

**Why this matters:**
The card lists the SymPy script as a verification artefact for this stage. A reader who runs only the SymPy script (or trusts its exit-0 as confirmation of the windows) is misled: the SymPy exit-0 confirms only the derivative identity, the eigenvalue root, and a monotonicity sample — not the four windows the stage exists to produce. The reconciliation depends entirely on the `.wl` for those values.

**Required change:**
In the SymPy script (py:111-162), add `expect_close`-style numeric assertions for the deliverable windows so the SymPy engine actually verifies the Output, mirroring the `.wl` battery (wl:171-179): assert `Pe_minus_chi ≈ 96.5285247264385`, `Pe_plus_chi ≈ 11220.5441626259`, `Pe_minus_J ≈ 22.0062226330754`, `Pe_plus_J ≈ 2558.01892349205`, `zeta_minus_chi ≈ 2.46622291347846`, `zeta_plus_chi ≈ 2.46752913273870`, `zeta_minus_J ≈ 2.44257571477179`, `zeta_plus_J ≈ 2.46752736855058`, `zeta_max_F1 ≈ 2.46752922945601`, each against the literal the notes report, with tolerances matching the `.wl` (1e-10 for Pe, 1e-12 for zeta). Add a small `expect_close` helper (the script currently has only `expect_zero`).

**Verification:**
After the fix, `grep -c "expect_close\|assert" sympy_audit.py` increases; the SymPy transcript shows the new window checks passing; perturbing any `Theta_*_coeff` literal (py:106-107) by 1e-3 now makes the SymPy script exit nonzero (previously it would still exit 0).

## Independent-derivation check (Mathematica)

The `.wl` is **not** a clean line-by-line transliteration in surface syntax (it adds the `expectApprox` numeric battery the SymPy lacks, and uses `FindRoot`/`DSolve`-style idioms), so I do not raise `mathematica_transliteration` as a standalone finding. However, the *structural* choreography of the load-bearing symbolic checks is identical: same closed-form `Delta_0`/`Delta_inf` (wl:55-60 ≡ py:57-64), same tautological residual construction (wl:77-91 ≡ py:71-81), same `Omega` closed form, same `zeta = A·Omega²`, same monotonicity sample at `{10,100,1000,10000}` (wl:186 ≡ py:144). The `.wl` advertises an "Independent BVP derivation" (wl:62-71) it never executes — so the *claimed* independence is absent precisely where it is advertised. The genuine added value of the `.wl` is the numeric anchoring (wl:169-179), which is independent of SymPy's algebra. Net: the engines are not echoing transliterations, but the symbolic verification is shared scaffolding; the `.wl`'s independence is carried entirely by its numeric pins, not by the BVP it claims. Folded into F1.

## Engine cross-check

Both engines agree to high precision on every value:

| value | sympy out | mathematica out |
|---|---|---|
| `Delta_0(F1)` | 0.000173302079021525149057156196550 | 0.0001733020790215251490571561965499… |
| `Delta_inf(F1)` | 0.0201447565540521594271032956099 | 0.0201447565540521594271032956099… |
| `y_F1` | 1.52948248371469964992710762240 | 1.5294824837146996499271076224… |
| `A_F1` | 1.00005192880219532865933408371 | 1.0000519288021953286593340837… |
| `Pe_-^(chi)` | 96.5285247264385161753086456051 | 96.5285247264385 |
| `zeta_-^(chi)` | 2.46622291347846457552701694323 | 2.466222913478465 |
| `zeta_max^(F1)` | 2.46752922945601223332958450157 | 2.467529229456012233329584501570… |

`engines_agree: true`. No `engine_disagreement`.

## Value Reconciliation (pass-2 augmentation)

Every deliverable value the scripts emit was located in the notes `.md` (the natural carrier; the terse `.tex` card legitimately omits intermediate constants). The `.tex` body equation (eq:app-stage083-zeta-direct) carries the symbolic Output form; the numeric windows live in the notes.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Delta_0(F1) = 1.73302079021525e-4` | py:92 / wl:105; out py:15 / wl:16 | notes:81 | MATCH |
| `Delta_inf(F1) = 2.01447565540522e-2` | py:93 / wl:106; out py:16 / wl:17 | notes:83 | MATCH |
| `Xi_chi coeff = 556995.768726174` | py:111; out py:19 | notes:105 | MATCH |
| `Xi_J coeff = 126981.873254631` | py:112; out py:20 | notes:107 | MATCH |
| `Pe_-^(chi) = 96.5285247264385` | py:119; out py:21 | notes:111 | MATCH |
| `Pe_+^(chi) = 11220.5441626259` | py:120; out py:22 | notes:113 | MATCH |
| `Pe_-^(J) = 22.0062226330754` | py:121; out py:23 | notes:115 | MATCH |
| `Pe_+^(J) = 2558.01892349205` | py:122; out py:24 | notes:117 | MATCH |
| `zeta_-^(chi) = 2.46622291347846` | py:152; out py:29 | notes:127 | MATCH |
| `zeta_+^(chi) = 2.46752913273870` | py:153; out py:30 | notes:129 | MATCH |
| `zeta_-^(J) = 2.44257571477179` | py:154; out py:31 | notes:131 | MATCH |
| `zeta_+^(J) = 2.46752736855058` | py:155; out py:32 | notes:133 | MATCH |
| `zeta_max^(F1) = 2.46752922945601` | py:156; out py:33 | notes:137 | MATCH |
| `Pi_suff^(chi)/C_mix = 3.46622291347846` | py:158; out py:34 | notes:169 | MATCH |
| `Pi_fail^(chi)/C_mix = 3.46752913273870` | py:159; out py:35 | notes:171 | MATCH |
| `Pi_suff^(J)/C_mix = 3.44257571477179` | py:160; out py:36 | notes:173 | MATCH |
| `Pi_fail^(J)/C_mix = 3.46752736855058` | py:161; out py:37 | notes:175 | MATCH |
| `Pi_max^(F1)/C_mix = 3.46752922945601` | py:162; out py:38 | notes:179 | MATCH |
| `Theta_chi_coeff = 4.06863235008162` (input literal) | py:106 / wl:119 | notes:99 | MATCH |
| `Theta_J_coeff = 0.927552032539308` (input literal) | py:107 / wl:120 | notes:101 | MATCH |
| prefactor `136900` | py:108-109 / wl:121-122 | notes:71 (`136900 Theta_w`) | MATCH |

INTERNAL (scaffolding, not expected in prose, no finding): `y_F1 ≈ 1.5294824837…` (eigenvalue intermediate; notes give only the defining eq `y tan y = eta`), `A_F1 ≈ 1.0000519288…` (support amplitude intermediate; notes give only the formula `(kappa+pi²/4)/(kappa+y²)`), the four `dzeta/dPe` monotonicity sample values, and all residual/diff scaffolding.

reconciliation: complete; 21 deliverable values checked, 0 misaligned.

## Verdict justification

`findings: 2`, both script-side `verification-quality` defects — no value mismatch, no paper misalignment, no stop-cold. Every numeric deliverable the scripts emit reconciles exactly with the notes (21/21), the body symbolic form matches the card, and both engines agree to ~30 digits. What does **not** hold up: (F1) the `Delta_0`/`Delta_inf` "defining-equation residual" checks are tautological by construction in both engines — `denom·(numer/denom) − numer ≡ 0` — yet the SymPy comment claims they are "non-tautological" and the Mathematica comment claims an "Independent BVP derivation" that the code never performs (the `A_F1` and `Omega` "independent" checks are likewise expression-typed-twice tautologies); and (F2) the SymPy engine never asserts a single one of the stage's headline window deliverables (`Pe`/`zeta`/`Pi` windows are printed, not checked), so the SymPy "audit" alone verifies only the derivative identity, the eigenvalue root, and a monotonicity sample. Cross-engine the windows are saved by the Mathematica `expectApprox` battery (wl:169-179), which is why this is `findings` (medium) rather than stop-cold: the values are correct and anchored once, but the SymPy script's claimed verification is hollow and the comments overstate the independence. I attacked the implicit-derivative identity (genuine, A1), the `y_F1` root residual (genuine, A4), the monotonicity spot-check (honest, partial), and the value reconciliation (all match) — those hold. I confirmed I read the card, notes, and appendix row; the script values match the paper, and the findings concern verification rigor, not correctness of the answer.

## Self-test notes

Checked: (1) variable independence — the only derivative is `dzeta/dPe` (py:143/wl:185), and `zeta_F1` genuinely depends on `Pe_sym`/`ppSym`, so it is non-trivially nonzero (confirmed by the positive sample outputs); no zero-derivative trap. (2) parity/integral — no unbounded-domain integrals in this unit; n/a. (3) trivial-case — the tautological residuals (F1) reduce to literal `0` for *any* input by construction, which is exactly the defect; the proposed numeric-anchor and `dsolve`-BVP replacements would instead reduce to a literal nonzero residual if the closed form were perturbed, so they are non-tautological as required. (4) path specs — F2 edits live in `scripts/...sympy_audit.py`; F1 touches both `scripts/...py` and `mathematica/...wl`; no new-script paths needed. (5) paper round-trip — the prescribed numeric anchors use exactly the literals the notes already report (`1.73302079021525e-4`, `2.01447565540522e-2`, and the window values at notes:111-179), so no new `paper_misalignment` is introduced.
