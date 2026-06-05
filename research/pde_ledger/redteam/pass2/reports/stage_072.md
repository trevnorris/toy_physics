---
unit_id: 072
batch: III.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-05T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage072_explicit_branch_thresholds.md]
  paper_appendix: present
---

# Audit unit 072 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_072.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage072_explicit_branch_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 122 references this stage)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.txt`

## What the paper claims

Stage 072 takes the universal matched-branch support/source theorem (the matched window `W_wall ≶ Pe_req/Delta_{inf,0}` from Stages 066/069) and re-expresses it directly in the three canonical branch controls `(chi_s, Lambda_ell, Upsilon_w)`. The card's `\stagefield{Output}` is: *"A threshold-surface map in the canonical branch variables."* The body equation (`eq:app-stage072-surfaces`) states `Upsilon_w Lambda_ell^2 ≶ Pe_req / Delta_{inf,0}(kappa(chi_s,Lambda_ell), eta=Lambda_ell)` with `kappa = 4 chi_s^2 + (4/5)Lambda_ell^2`, `eta = Lambda_ell`, `W_wall = Upsilon_w Lambda_ell^2`. The notes enumerate the explicit deliverables: the two threshold surfaces `Upsilon_fail = Pe_req/[Lambda_ell^2 Delta_inf]`, `Upsilon_suff = Pe_req/[Lambda_ell^2 Delta_0]`; the physical wall-amplitude versions `V0_fail^2`, `V0_suff^2`; and two asymptotic regimes — shell-gradient-dominated (`Upsilon_fail ~ 2 Pe_req/(sqrt(5) Lambda_ell)`, `Upsilon_suff ~ (4/5)(1+2/sqrt(5)) Pe_req`) and compression-dominated (`Upsilon_fail ~ 2 Pe_req chi_s/Lambda_ell^2`, `Upsilon_suff ~ 4 Pe_req chi_s^2(Lambda_ell+2 chi_s)/Lambda_ell^3`).

## What the script claims to verify

Both engines (1) substitute the explicit branch formulas `kappa = 4 chi_s^2 + (4/5)Lambda_ell^2`, `eta = Lambda_ell` into the canonical closed forms of `Delta_0` and `Delta_inf` (the Stage-058 forms `Delta_0 = eta(cosh α−1)/(α^2 W)`, `Delta_inf = (cosh α + (eta/α)sinh α − 1)/W`, `W = α sinh α + eta cosh α`, `α = sqrt(kappa)`), then form `Upsilon_fail`, `Upsilon_suff`, `V0_fail^2`, `V0_suff^2`. They then build hand-written shell-dominated (`α → c Lambda_ell`, `c = 2/sqrt(5)`) and compression-dominated (`α → 2 chi_s`) leading-order forms and confirm via residual-zero asserts that (a) the ratio of the full closed form to the hand-built form tends to 1 in the appropriate limit, and (b) the hand-built form reproduces the closed-form asymptotic thresholds quoted in the notes. All asserts are residual-zero (`expect_zero` / `expectZero` with `Exit[1]` on nonzero).

## Paper ↔ script cross-check

| Paper/notes deliverable | Script-side check | Status |
|---|---|---|
| `kappa = 4 chi_s^2 + (4/5)Lambda_ell^2`, `eta = Lambda_ell` | py:29-30, wl:33-34 (definitions) | match |
| Surface eq `W_wall ≶ Pe_req/Delta` with `W_wall = Upsilon_w Lambda_ell^2` | `Upsilon_fail/suff = Pe_req/(Lambda_ell^2 Delta_inf/0)` py:42-43, wl:46-47 | match |
| `Upsilon_fail`, `Upsilon_suff` closed surfaces | printed py:49-50, wl:53-54; exercised by asymptotic asserts | match |
| `V0_fail^2`, `V0_suff^2` | py:52-54, wl:56-58 | match |
| Shell `Upsilon_fail ~ 2 Pe_req/(sqrt5 Lambda_ell)` | `shell fail asymptotic` py:79-80, wl:85 | match |
| Shell `Upsilon_suff ~ (4/5)(1+2/sqrt5)Pe_req` | `shell suff asymptotic` py:81-82, wl:86 | match |
| Comp `Upsilon_fail ~ 2 Pe_req chi_s/Lambda_ell^2` | `compression fail asymptotic` py:106, wl:110 | match |
| Comp `Upsilon_suff ~ 4 Pe_req chi_s^2(Lambda_ell+2 chi_s)/Lambda_ell^3` | `compression suff asymptotic` py:107-108, wl:111-114 | match |
| Leading-order validity of shell/comp forms | ratio→1 limit asserts py:75-78,97-100, wl:78-81,103-106 | match |

`paper_alignment: aligned`. The Stage-058 canonical `Delta_0`/`Delta_inf` forms used here are identical to those derived and verified upstream (stage 058 sympy lines 84,90), so the carried-forward inputs are correct.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 75-76 | `expect_zero(Delta0/Delta0_shell limit − 1)` | shell leading-order validity | yes |
| A2 | sympy | 77-78 | `expect_zero(DeltaInf/DeltaInf_shell limit − 1)` | shell leading-order validity | yes |
| A3 | sympy | 79-80 | `expect_zero(Pe/(L^2 DeltaInf_shell) − Upsilon_fail_shell)` | shell fail threshold | yes |
| A4 | sympy | 81-82 | `expect_zero(Pe/(L^2 Delta0_shell) − Upsilon_suff_shell)` | shell suff threshold | yes |
| A5 | sympy | 97-98 | `expect_zero(Delta0/Delta0_comp limit − 1)` | comp leading-order validity | yes |
| A6 | sympy | 99-100 | `expect_zero(DeltaInf/DeltaInf_comp limit − 1)` | comp leading-order validity | yes |
| A7 | sympy | 106 | `expect_zero(Upsilon_fail_comp − 2 Pe chi_s/L^2)` | comp fail threshold | yes |
| A8 | sympy | 107-108 | `expect_zero(Upsilon_suff_comp − 4 Pe chi_s^2(L+2chi_s)/L^3)` | comp suff threshold | yes |
| B1-B8 | mathematica | 78-81,85-86,103-106,110-114 | `expectZero[...]` mirror of A1-A8 | same deliverables | yes |

All eight residual checks compare two **independently constructed** expressions: a hand-built leading-order form (built directly from the regime substitution `α → c Lambda_ell` or `α → 2 chi_s`) against either the limit of the full canonical closed form (A1/A2/A5/A6) or the notes-quoted asymptotic constant (A3/A4/A7/A8). None is tautological: if the hand-built form were wrong, the ratio would not be 1 and the residual would not vanish. The printed `Upsilon_fail/suff`, `V0_fail/suff^2` (py:49-55, wl:53-59) are not independently asserted but are pure substitutions whose correctness is implied by the canonical Delta forms; the load-bearing checks are the eight asymptotic asserts plus the four ratio→1 limits.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.txt:3,31` (banner reads `STAGE 55`)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.txt:3` (banner reads `STAGE 055`)

**What's wrong:**
Both committed transcripts predate the current scripts and were produced by the pre-renumber (055) revision. mtimes: sympy `.py` `2026-05-26 13:08`, sympy `.txt` `2026-05-22 21:38`; mathematica `.wl` `2026-05-26 13:08`, mathematica `.txt` `2026-05-22 21:39`. Both `.txt` are ~4 days older than their scripts. The current scripts emit `STAGE 072` (py:24, wl:26), but the captured transcripts emit `STAGE 55` (sympy out line 3, line 31 "STAGE 55 THEOREM LEDGER") and `STAGE 055` (mathematica out line 3). The numeric/symbolic content of the transcripts otherwise matches what the current scripts produce (all closed forms, all `= 0` residuals, all PASS lines are correct), so this is a freshness/label drift, not a content disagreement.

**Why this matters:**
The committed transcript is the auditable record of the stage; a stale banner means the captured run does not correspond to the current source. Low risk because content agrees, but the orchestrator's independent re-run should refresh both `.txt` so banners read `STAGE 072`.

**Required change:**
Re-run both scripts and overwrite the committed transcripts so the banners read `STAGE 072` (sympy) / `STAGE 072` (mathematica) and the ledger header reads `STAGE 072 THEOREM LEDGER`. No source edit needed beyond refreshing the captured output.

**Verification:**
After re-run, `scripts/output/...stage072..._sympy_audit.txt` line 3 reads `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES` and line 31 `STAGE 072 THEOREM LEDGER`; `mathematica/output/...stage072..._mathematica_audit.txt` line 3 reads `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES`. All residual/PASS lines unchanged.

## Independent-derivation check (Mathematica)

The `.wl` mirrors the `.py` choreography closely: identical `kappa`, identical `delta0`/`deltaInf` closed forms (wl:37-44 ≡ py:33-40), identical hand-built shell forms `delta0Shell`/`deltaInfShell` with `cShell = 2/Sqrt[5]` (wl:63-67 ≡ py:61-66), identical comp forms (wl:90-91 ≡ py:88-89), and identical residual targets (wl:85-86,110-114 ≡ py:79-82,106-108). This is structurally a transliteration. However, the stage's claim is a closed-form substitution + asymptotic-limit identity, where the only "derivation" freedom is which canonical Delta form to plug in (fixed upstream) and which limit to take. Both engines independently execute `FullSimplify`/`simplify` and `Limit`/`limit` rather than echoing a cached numeric result, so the algebra is genuinely re-run on each side. Per the established batch policy for pure identity/limit stages, the mirror is acceptable and not raised as a hard `mathematica_transliteration` finding — but it is noted (`independent_wl: transliteration`).

## Engine cross-check

The two engines agree at the level they claim. Both print the same canonical `Delta_0`/`Delta_inf` (algebraically identical, differing only in cosmetic simplification: e.g. sympy writes `cosh(2*sqrt5*sqrt(...)/5) − 1` while Mathematica writes `Sinh[...]^2` via the half-angle identity — same quantity). Both reduce all eight asymptotic residuals to `0` (sympy out lines 15-18,21-22,27-28; mathematica out lines 23-32,43-52, all PASS). The one cosmetically alarming line — mathematica out line 22 `DeltaInf shell leading-order ratio = 2/Sqrt[5] + (5 + 2*Sqrt[5])^(-1)` — is exactly 1 (`2/Sqrt[5] + 1/(5+2 Sqrt[5]) = 1`), which is why the subsequent `expectZero[ratio − 1]` at wl:80-81 yields 0 and PASS (mathematica out lines 25-26). The `Limit::alimv` warnings (assumptions involving the limit variable ignored) are benign — the limits resolve correctly. No engine disagreement.

## Verdict justification

The scripts faithfully and non-tautologically verify exactly what the paper card and notes claim: the branch-control re-expression of the matched support/source window, its two explicit threshold surfaces, the physical wall-amplitude versions, and both asymptotic regimes. The canonical `Delta_0`/`Delta_inf` inputs match the upstream Stage-058 forms verbatim. Attacks tried and failed: (1) checked whether the asymptotic asserts are tautological — they compare independently constructed forms, not `x == x`; (2) checked whether the shell/comp hand-built forms are the true leading order — confirmed by the ratio→1 limits and by hand (large-α reduction); (3) checked the `2/Sqrt[5] + 1/(5+2 Sqrt[5])` ratio that the Mathematica transcript prints unsimplified — it is exactly 1; (4) checked every deliverable value against the notes — all reconcile; (5) checked positivity assumptions — all symbols `positive=True`/`> 0`, consistent with the physical setup (all are amplitudes/ratios). The only defect is `stale_output`: both committed transcripts carry the pre-renumber `STAGE 55/055` banner and predate the current scripts; content otherwise agrees. Verdict `findings` (one low-severity, non-blocking stale-output), `paper_alignment: aligned`.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 13 values checked, 0 misaligned.

The paper card (`stage_072.tex`) is terse — it states only the surface inequality `eq:app-stage072-surfaces` and the abstract `Delta_{inf,0}` — but the per-stage notes (`...stage072_explicit_branch_thresholds.md`) carry every explicit closed form. Per the reconciliation guards, a value living correctly in the `.md` notes is a MATCH.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `kappa = 4 chi_s^2 + (4/5)Lambda_ell^2` | py:29, wl:33; sympy out L5 | tex:18 (`kappa(chi_s,Lambda_ell)`), md:16/37 | MATCH |
| `eta = Lambda_ell` | py:30, wl:34; sympy out L6 | tex:18 (`eta=Lambda_ell`), md:18/39 | MATCH |
| `W_wall = Upsilon_w Lambda_ell^2` | (LHS of surface eq) py:42-43 | tex:16, md:19/41 | MATCH |
| `Upsilon_fail = Pe_req/[Lambda_ell^2 Delta_inf]` | py:42, wl:46; sympy out L9 | md:51-52 | MATCH |
| `Upsilon_suff = Pe_req/[Lambda_ell^2 Delta_0]` | py:43, wl:47; sympy out L10 | md:54-55 | MATCH |
| `V0_fail^2 = hbar^2 c_sw^2 Upsilon_fail/(4 rho_w^2)` | py:52, wl:56; sympy out L11 | md:69-70 | MATCH |
| `V0_suff^2 = hbar^2 c_sw^2 Upsilon_suff/(4 rho_w^2)` | py:53, wl:57; sympy out L12 | md:72-73 | MATCH |
| `Upsilon_fail_shell = 2 Pe_req/(sqrt5 Lambda_ell)` | py:62, wl:64; sympy out L23 | md:99 | MATCH |
| `Upsilon_suff_shell = (4/5)(1+2/sqrt5) Pe_req` | py:63, wl:65; sympy out L24 | md:101 | MATCH |
| `Upsilon_fail_comp = 2 Pe_req chi_s/Lambda_ell^2` | py:90→106, wl:92→110; sympy out L25 | md:119 | MATCH |
| `Upsilon_suff_comp = 4 Pe_req chi_s^2(Lambda_ell+2 chi_s)/Lambda_ell^3` | py:91→108, wl:93→111; sympy out L26 | md:121 | MATCH |
| `Delta_0` closed form (eta(cosh α−1)/(α^2 W)) | py:33-36, wl:37-40; sympy out L7 | tex:18 abstract; canonical from Stage 058 | MATCH (carried) |
| `Delta_inf` closed form ((cosh α+(eta/α)sinh α−1)/W) | py:37-40, wl:41-44; sympy out L8 | tex:18 abstract; canonical from Stage 058 | MATCH (carried) |

INTERNAL (scaffolding, no prose expected): `alpha = sqrt(kappa)`, `W = alpha sinh alpha + eta cosh alpha` (intermediate), `c = 2/sqrt(5)` (shell substitution constant), `Delta0_shell`/`DeltaInf_shell`/`Delta0_comp`/`DeltaInf_comp` (hand-built leading-order intermediates), the four ratio→1 limit values, and all `= 0` residual / PASS flags.

All 13 deliverable values reconcile (to notes for the explicit forms, to the abstract surface eq + upstream Stage 058 for the Delta closed forms). Zero MISMATCH, zero MISSING-DELIVERABLE.

## Self-test notes

Trap 1 (variable independence / zero derivatives): N/A — no `sp.diff`/`D` calls; the stage is algebraic identities plus `Limit`. Trap 2 (integral parity): N/A — no integrals. Trap 3 (trivial-case): each `expect_zero` compares two independently built expressions (hand-built leading form vs. limit of the full canonical form, or vs. the notes constant), so a vanishing residual is a real match, not a constructed tautology; I hand-verified the large-α reductions for both shell and comp regimes and they reproduce the asserted constants. Trap 4/5 (paths / paper round-trip): the only fix is a transcript refresh — no source edit, no risk of introducing a new paper_misalignment. I also confirmed the Mathematica `2/Sqrt[5]+1/(5+2 Sqrt[5])` ratio equals exactly 1, so the PASS is legitimate.
