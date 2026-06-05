---
unit_id: 064
batch: III.3
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
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage064_equilibrium_alignment.md]
  paper_appendix: present
---

# Audit unit 064 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_064.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage064_equilibrium_alignment.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row at line 106; `\input` at line 246)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.txt`

## What the paper claims

Stage 064 proves that the parent equilibrium equations tie the support-induced source profile to the support profile rather than leaving the coherence `C_(sigma phi)` arbitrary. The bottom-line `\stagefield{Output}` (line 56) reads: *"Equilibrium source/support alignment \eqref{eq:app-stage064-source-profile} and matched-layer gain \eqref{eq:app-stage064-matched-layer}."* The distinct deliverables are: (1) the boxed aligned source profile `chi_sigma(y) = g_phi chi_phi(y)/H(y)` with `H = h'(rho_*)` (eq:app-stage064-source-profile, lines 16-19); (2) the overlap invariants `O_(sigma phi)=g_phi I_1`, `N_(sigma sigma)=g_phi^2 I_2` with `I_1=∫chi_phi^2/H`, `I_2=∫chi_phi^2/H^2` (lines 24-34); (3) the boxed coherence `C^2 = I_1^2/(N_(phi phi) I_2) ≤ 1` (eq:app-stage064-coherence, lines 38-39); (4) the eliminated-source softening `Delta K_X^eq = g_phi^2 I_1`, `G_eq = g_phi^2 I_1/K_X` (eq:app-stage064-softening, lines 44-46); and (5) the matched-layer (`H≈H_w`) reductions `C^2=1`, `G_eq = g_phi^2 N_(phi phi)/(K_X H_w)` (eq:app-stage064-matched-layer, lines 51-53). The notes add the derivation source (local static linear-response closure, §1) and the Cauchy–Schwarz argument with equality iff `1/H` is constant on `supp(chi_phi)` (§3, line 132).

## What the script claims to verify

Both scripts assert: the variational/linear-response minimiser of the local source-energy gives the closure `chi_sigma = g_phi chi_phi/H`; the overlap identities `O = g_phi I_1` and `N_ss = g_phi^2 I_2` (checked on a concrete profile — Gaussian in SymPy, Lorentzian in Mathematica); the coherence form `C^2 = I_1^2/(N_pp I_2)`; the matched-layer reductions `I_1 = N_pp/H_w`, `I_2 = N_pp/H_w^2` yielding `C^2 = 1` and `G_eq = g_phi^2 N_pp/(K_X H_w)`; a strict-positivity Cauchy gap `N_pp I_2 − I_1^2 ≥ 0` (discrete two-point in SymPy, continuous variable-H Lorentzian in Mathematica); the general-H softening `Lambda^2/Theta = g_phi^2 I_1` for arbitrary `H`; and the eliminated-source softening `F_eff = (1/2)(K_X − Lambda^2/Theta)phi^2`. Every check is a `expect_zero`/`expectZero` residual that raises/exits on nonzero.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| (1) `chi_sigma = g_phi chi_phi/H` | py L64-68 (`solve` of dF/dsigma, residual vs `g_phi chi_phi/H`); wl L37-44 (`VariationalD` + `Solve`) | match |
| (2) `O = g_phi I_1`, `N_ss = g_phi^2 I_2` | py L74-82 (Gaussian integrals); wl L52-82 (Lorentzian integrals) | match |
| (3) `C^2 = I_1^2/(N_pp I_2) ≤ 1` | py L87,121 + L129-139 (discrete gap ≥0); wl L86,93 + L106-153 (continuous gap, `>=0` asserts) | match |
| (4) `Delta K_X^eq = g_phi^2 I_1`, `G_eq = g_phi^2 I_1/K_X` | py L88,148-159,166-174; wl L87,160-180,182-201 | match |
| (5) matched-layer `C^2=1`, `G_eq = g_phi^2 N_pp/(K_X H_w)` | py L108-122; wl L97-104 | match |

`paper_alignment: aligned`. Every paper-side deliverable has a faithful, non-tautological script-side counterpart in both engines. The appendix row (part03 line 106) matches the card Output verbatim.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 68 | `expect_zero(closure − g_phi chi_phi/H)` | claim 1 | yes |
| A2 | sympy | 81 | `expect_zero(O_int − g_phi I1_int)` | claim 2 | yes |
| A3 | sympy | 82 | `expect_zero(Nss_int − g_phi^2 I2_int)` | claim 2 | yes |
| A4 | sympy | 108-109 | `expect_zero(I1_int − Npp/Hw)`, `(I2_int − Npp/Hw^2)` | claim 5 | yes |
| A5 | sympy | 121-122 | `expect_zero(C2_const − 1)`, gain vs `g_phi^2 Npp/(KX Hw)` | claim 5 | yes |
| A6 | sympy | 139 | `expect_zero(gap_disc − w1 w2 (H1−H2)^2/(H1^2 H2^2))` | claim 3 (gap≥0) | yes |
| A7 | sympy | 156-159 | `expect_zero(soft_general − g_phi^2 I1_disc)` | claim 4 (general H) | yes |
| A8 | sympy | 171-174 | `expect_zero(F_eff − (1/2)(KX−Lambda^2/Theta)phi^2)` | claim 4 | yes |
| A9 | math | 43-44 | `expectZero(closure − gPhi chiPhi/hLoc)` | claim 1 | yes |
| A10 | math | 81-82 | `expectZero(osp − gPhi i1)`, `(nss − gPhi^2 i2)` (Lorentzian) | claim 2 | yes |
| A11 | math | 93 | `expectZero(c2 − i1^2/(npp i2))` | claim 3 | yes |
| A12 | math | 97-98,103-104 | matched-layer I1/I2 reductions, `c2Const−1`, gain | claim 5 | yes |
| A13 | math | 144,150 | `expectZero(pairGap − pairExpected)`, `(cauchyGap − cauchyExpected)` + `>=0` | claim 3 | yes |
| A14 | math | 169-172 | `expectZero` Theta/Lambda integrand `= gPhi^2 chiPhi^2/hFun` | claim 4 (general H) | yes |
| A15 | math | 179-180 | `expectZero(softGeneral − gPhi^2 i1Integral)` | claim 4 | partial (bookkeeping; load carried by A14) |
| A16 | math | 198-201 | `expectZero(effEnergy − (1/2)(kX−mixCoeff^2/sourceCoeff)supportAmp^2)` | claim 4 | yes |

A15 is near-tautological in isolation (`thetaGeneral` and `lambdaGeneral` are assigned the identical expression at wl L175-176, so `Lambda^2/Theta` collapses to that expression by construction), but the substantive content — that both Theta and Lambda equal `gPhi^2 I_1` — is exercised non-tautologically by A14's integrand-equality checks. Not a standalone finding.

## Findings

### F1 — stale_output

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.txt` (whole file)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage064_equilibrium_alignment_mathematica_audit.txt` (whole file)

**What's wrong:**
Both saved `.txt` transcripts predate the current scripts and no longer reflect them — this is both an mtime and a content divergence.

mtimes:
- SymPy `.py` = 2026-06-03 15:59:11; SymPy `.txt` = 2026-05-22 19:48:42 (≈12 days stale).
- Mathematica `.wl` = 2026-06-03 15:59:11; Mathematica `.txt` = 2026-05-22 19:48:50 (≈12 days stale).

Content divergence (SymPy): the current `.py` has a `GENERAL-H EQUILIBRIUM SOFTENING CHECK` section (L141-159) and the matched-layer integrals use a Gaussian. The committed `.txt` instead shows an obsolete `ELIMINATED-SOURCE SOFTENING CHECK` block with `Lambda^2/Theta (closure form) = I1**2*g_phi**2/(H_w*I2)` and `Lambda^2/Theta (matched layer) = I1*g_phi**2` (txt L45-46) — lines that no longer exist in the current script. The banner reads `STAGE 47 — ...` / `STAGE 47 AUDIT PASSED` (txt L3, L50) where the current `.py` prints `STAGE 64` (py L49, L176).

Content divergence (Mathematica): the current `.wl` has a `CONTINUOUS CAUCHY BOUND CHECK` section (L106-153) using a continuous Lorentzian with variable `H` and hardcoded closed forms `Pi^2*L^2/(64*hw^2)` and `L^2*(15*Pi^2/128 − 256/225)/hw^3`. The committed `.txt` instead shows a `discrete two-point` gap `((h1−h2)^2 w1 w2)/(h1^2 h2^2)` (txt L34-37) that the current `.wl` does not produce, and a `Lambda^2/Theta` block (txt L42-43) absent from the current `.wl`. The banner reads `STAGE 047` (txt L3) where the current `.wl` prints `STAGE 064` (wl L26).

**Why this matters:**
The committed transcripts are not evidence for the current scripts; a reader trusting `output/*.txt` would be reading a previous revision's results. The orchestrator's independent re-run will regenerate them.

**Required change:**
Refresh both `.txt` outputs by re-running the current scripts (orchestrator exec / sed-refresh of the committed transcripts). No script-logic change is required for this finding.

**Verification:**
After refresh, SymPy `.txt` banner reads `STAGE 064`, contains a `GENERAL-H EQUILIBRIUM SOFTENING CHECK` section with `Theta (general, two-point)`/`Lambda (general, two-point)` lines, and no longer contains `Lambda^2/Theta (closure form)`. Mathematica `.txt` banner reads `STAGE 064`, contains a `CONTINUOUS CAUCHY BOUND CHECK` section with `continuous pair Cauchy gap` and `expected continuous gap` lines, and no longer contains the discrete `((h1−h2)^2 w1 w2)/(h1^2 h2^2)` gap. Both `.txt` mtimes newer than their script mtimes.

### F2 — paper_misalignment (subtype: target_mismatch — stale self/cross numbering labels in SymPy source)

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:2` (`"""Moving-Throat PDE — Stage 47 SymPy audit`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:25` (`reproducing the Stage-45/46 best-alignment formulas.`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:122` (`"matched-layer gain vs Stage-45 best-alignment formula"`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:180` (`In the matched-layer limit this reproduces the Stage-45/46 best-alignment formulas.`)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:49` (`banner("STAGE 64 — ...")` — 2-digit, sibling banners use 3-digit `064`)

**What's wrong:**
The canonical stage is **064** (paper card `\label{stage:064}` line 2; appendix row part03 line 106). The SymPy source carries stale pre-renumber self-labels: the docstring says `Stage 47` (47 + 17 = 64) and the matched-layer comparison comments/labels say `Stage-45/46` (45→062, 46→063). The notes are explicit that the best-alignment formulas come from **Stages 062/063** ("reproducing the Stage-062/063 best-alignment formulas", notes line 55; "Stage-062/063 formula", notes line 178). The Mathematica `.wl` source is already correct (banner `STAGE 064`, wl L26; closing `Stage 064` L204) — only the SymPy source and the stale `.txt` banners carry the old `47`/`45`/`46` labels. This is the known `+17` renumber drift, label-only, `material_change: false`.

**Why this matters:**
Self-labels naming the wrong stage number mislead a reader/auditor about which stage the script verifies and which upstream stages it reproduces. No math is affected.

**Required change:**
Per the in-loop Reading-2 policy (verdict:findings ⇒ fix unambiguous self-labels), update the SymPy source label-only:
- L2 docstring `Stage 47` → `Stage 064`.
- L25 `Stage-45/46 best-alignment formulas` → `Stage-062/063 best-alignment formulas`.
- L122 assertion label `Stage-45 best-alignment formula` → `Stage-062/063 best-alignment formula`.
- L180 `Stage-45/46 best-alignment formulas` → `Stage-062/063 best-alignment formulas`.
- L49 banner `STAGE 64` → `STAGE 064` (3-digit consistency with the `.wl` banner and sibling stages); correspondingly L176 `STAGE 64 AUDIT PASSED` → `STAGE 064 AUDIT PASSED`.
These are comment/string edits only; no expression, symbol, or assertion changes. Then re-run to refresh the `.txt` (folds into F1).

**Verification:**
SymPy source contains no `Stage 47`, `Stage-45`, or `Stage-46` tokens; banner prints `STAGE 064`. Refreshed SymPy `.txt` banner matches.

## Independent-derivation check (Mathematica)

The `.wl` is an **independent** re-derivation, not a transliteration. Three corresponding sections diverge structurally:
1. Closure: SymPy uses `sp.diff(F_loc, sigma_loc)` + `sp.solve` on an algebraic energy (py L63-64); Mathematica uses `VariationalD[energyDensity, sigmaFun[yLoc], yLoc]` + `Solve` on a field functional (wl L37-39). Different machinery (algebraic minimiser vs. variational derivative).
2. Overlap integrals: SymPy uses a **Gaussian** profile `exp(-y^2/(2L^2))` (py L72,98); Mathematica uses a **Lorentzian** `1/(1+z^2/L^2)` (wl L48). Different test profiles is the strongest possible independence signal — both must satisfy the same identities.
3. Cauchy bound: SymPy proves `≥0` via a **discrete two-point** algebraic identity `w1 w2 (H1−H2)^2/(H1^2 H2^2)` (py L129-139); Mathematica proves it via a **continuous variable-H Lorentzian** integral with closed forms `Pi^2 L^2/(64 hw^2)` and `L^2(15Pi^2/128 − 256/225)/hw^3` (wl L106-153). Entirely different routes to `C^2 ≤ 1`.

I independently confirmed the Mathematica hardcoded Cauchy closed forms against the integrals they are checked against: with `chiPhiL^2 = (1+y^2/L^2)^{-2}`, `nppInt = πL/2`, `i1Var = 3πL/(8 hw)`, `i2Var = 5πL/(16 hw^2)` ⇒ `nppInt·i2Var − i1Var^2 = π^2L^2/(64 hw^2)` ✓; and `ffIntegral=3πL/(8hw)`, `ggIntegral=5πL/(16hw^2)`, `fgIntegral=16L/(15 hw^{3/2})` ⇒ `pairGap = L^2/hw^3 (15π^2/128 − 256/225)` ✓. The hardcoded constants are correct closed forms of integrals the script itself computes symbolically — non-tautological, not magic numbers.

## Engine cross-check

Both engines verify the same paper claims and (per the committed transcripts and the residual structure) both report all residuals = 0 / PASS. They use different profiles and different Cauchy routes, so a numeric side-by-side is not meaningful; the shared symbolic deliverables (`O=g_phi I1`, `N_ss=g_phi^2 I2`, `C^2=I1^2/(N_pp I2)`, `C^2=1` matched, `G_eq=g_phi^2 Npp/(KX Hw)`, `Lambda^2/Theta=g_phi^2 I1`, `F_eff=(1/2)(KX−Lambda^2/Theta)phi^2`) agree between A1-A8 (SymPy) and A9-A16 (Mathematica). No engine disagreement. (Note: the *committed* transcripts are stale per F1; this cross-check is read off the current script sources, which agree.)

## Verdict justification

`findings`. The mathematics holds up under attack: the closure law, overlap identities, coherence form, matched-layer reductions, both Cauchy-gap routes (discrete and continuous — I re-derived the Mathematica hardcoded closed forms and they are correct), the general-H softening, and the eliminated-source softening are all non-tautological, well-anchored, exercised on concrete profiles, and faithful to the paper's five deliverables in both engines. Every emitted deliverable value reconciles with the card and notes (0 misaligned). The two findings are non-mathematical: (F1) both committed `.txt` transcripts are stale in mtime AND content (they belong to a prior revision — different sections, obsolete prints, old `STAGE 47/047` banners), and (F2) the SymPy source carries stale `+17` pre-renumber self/cross labels (`Stage 47`, `Stage-45/46`) that should read `064` / `062/063`. F2 is label-only with `material_change: false`; despite the `paper_misalignment` category tag it is a self-label correction routed under the in-loop Reading-2 policy, not a math discrepancy requiring user resolution. Attacks tried that failed: (a) checked A15/the general-H softening for tautology — load is carried by the non-tautological integrand checks A14; (b) checked the Cauchy hardcoded constants for being unanchored magic numbers — they are exact closed forms of the script's own integrals; (c) checked symbol domains (`H,w,Theta,Lambda,KX>0`; `sigma` unrestricted) — all physically justified; (d) checked engine independence (Gaussian vs Lorentzian, discrete vs continuous) — genuinely independent. I confirmed I read the paper card, the notes, and the appendix row, and the scripts' verified claims match them.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 7 deliverable values checked, 0 misaligned.

All values read from the current script sources (`.py`/`.wl`); the committed `.txt` transcripts are stale (F1), so reconciliation is based on the script source + the residual identities they assert, as the augmentation guard permits.

| value | source (py/wl) | .tex/.md location | status |
|---|---|---|---|
| `chi_sigma = g_phi chi_phi/H` | py L68; wl L43-44 | tex L16-19 (eq:source-profile, boxed); md L85 | MATCH |
| `O_(sigma phi) = g_phi I_1` | py L81,85; wl L81,84 | tex L31; md L103 | MATCH |
| `N_(sigma sigma) = g_phi^2 I_2` | py L82,86; wl L82,85 | tex L32; md L105 | MATCH |
| `C_(sigma phi)^2 = I_1^2/(N_pp I_2)` (`≤1`, `=1` matched) | py L87,121; wl L86,93,103 | tex L38-39 (eq:coherence, boxed), L51; md L116,146 | MATCH |
| `Delta K_X^eq = g_phi^2 I_1` | py L88,156-159; wl L87,179-180 | tex L44; md L172 | MATCH |
| `G_eq = g_phi^2 I_1/K_X` | py L88; wl L87 | tex L46 (eq:softening); md L176 | MATCH |
| matched-layer `G_eq = g_phi^2 N_pp/(K_X H_w)` | py L122; wl L104 | tex L53 (eq:matched-layer, Output); md L53,190 | MATCH |

INTERNAL (scaffolding for the `C^2≤1` Cauchy verification and integral cross-checks; not stage deliverables, no finding): Gaussian `Npp_int=√π L`, `I1_int=√π L/H_w`, `I2_int=√π L/H_w^2` (py); Lorentzian `nppInt=πL/2`, `i1Var=3πL/(8hw)`, `i2Var=5πL/(16hw^2)` (wl); discrete two-point gap `w1 w2 (H1−H2)^2/(H1^2 H2^2)` (py); continuous gap closed forms `π^2 L^2/(64 hw^2)` and `L^2(15π^2/128 − 256/225)/hw^3` (wl); `f_eff`/`sigma_stat` intermediate forms.

## Self-test notes

Checked: (1) variable independence — the SymPy `solve`/`diff` minimisers and the Mathematica `VariationalD`/`Solve` act on energies that genuinely depend on `sigma_loc`/`sigmaFun`, so the closure derivative is not identically zero; the general-H `Lambda^2/Theta` collapse is backed by non-tautological integrand checks (A14). (2) Symmetry/parity — all infinite-domain integrals have even integrands (`chiPhi^2`, `chiPhi^2/H^n`, `f g`, `f^2`, `g^2`), so nonzero results are consistent; I re-derived the Lorentzian closed forms and they match the hardcoded constants. (3) Trivial-case — substituting `H1=H2` into the discrete gap gives 0 (equality case, correct); constant-`H` matched layer gives `C^2=1` as asserted. F2 is label-only (`material_change: false`); F1 requires only an output refresh — no assertion edits, so no new paper_misalignment is introduced.
