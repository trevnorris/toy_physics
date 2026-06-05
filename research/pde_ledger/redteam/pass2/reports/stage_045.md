---
unit_id: 045
batch: III.1
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
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage045_coherent_local_tracking.md]
  paper_appendix: present
---

# Audit unit 045 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_045.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 68; cross-refs at 315, 347)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage045_coherent_local_tracking_mathematica_audit.txt`

## What the paper claims

Stage 045 imposes the first concrete "coherent local D/N kernel" hypothesis: the mixed lane `W` and support lane `phi` couple to the *same* local wall/U support density, so the four modal amplitudes take the form `g_W = lam_W/sqrt(mu_eta mu_W)`, `g_R = gamma lam_W/sqrt(mu_U mu_W)`, `g_B = lam_phi/sqrt(mu_eta mu_phi)`, `g_S = gamma lam_phi/sqrt(mu_U mu_phi)`. The card's `\stagefield{Output}` states verbatim: "The exact tracking condition \eqref{eq:app-stage045-tracking} and coherent tracking factor \eqref{eq:app-stage045-Rtr}." The boxed deliverables are (a) the coherence condition `g_B g_R = g_W g_S` (eq app-stage045-coherent-condition), (b) the resulting exact tracking surface `R_phi = R_U = R_tr` (eq -tracking), and (c) the tracking factor `R_tr = [1+chi_0/(1+delta_U)]/(1+chi_0)` with range `1/(1+delta_U) < R_tr < 1` (eq -Rtr). The notes carry four further deliverables beyond the terse card: `rho_0 = sigma_0 =: chi_0`; the total baseline `M_tr = M_mix + M_supp = 8(1+chi_0)^2/[pi^2(1-eps_eta)]·[Z_W/(1-eps_W^split)+Z_phi/(1-eps_phi^split)]`; the collapse of the Stage-044 quadratic branch equation to the one-parameter tracking law `M_tr = G_tr` with the D/N form `G_tr = 9 xi(xi+delta)/[9 delta+(9+2R^2)xi]`; and the normalization collapse `R_target = F_tr` with the D/N form `F_tr = [9δ+(9+2R²)ξ]²[9δ+(9+2R)ξ]²/[81(1-ξ)(9δ²+18δξ+(9+2R²)ξ²)²]`. The card body also boxes the tracking laws `M_tr = G_tr(...)`, `R_target = F_tr(...)` (eq -tracking-laws).

## What the script claims to verify

The docstring lists four checks: (1) the coherent kernel implies `g_B g_R = g_W g_S` exactly; (2) the interference ratios coincide `rho_0 = sigma_0`; (3) `R_tr` satisfies the two exact range identities; (4) the (mislabeled "Stage-27", actually Stage-044) quadratic branch equation collapses to the one-parameter tracking law. The assertions go further than the docstring summary: section 1 extracts `g_*` two independent ways and cross-checks them before testing the coherence identity; section 3 builds `M_mix`/`M_supp`/`M_tr` symbolically (printed, not asserted — carry-forward from upstream); section 4 verifies the quadratic collapse, solves for `M_tr_req`, specializes `G_tr` to the D/N value `lam0=2/9`, and imports the Stage-044 `D_cont`/`F_cont` residual to verify both the generic tracking-F collapse and the D/N `F_tr` form against the notes' closed forms. All `expect_zero`/`expectZero` checks pass with residual 0.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Coherence condition `g_B g_R = g_W g_S` (eq -coherent-condition) | sympy:68 / wl:56 `expect_zero("g_B g_R - g_W g_S", ...)` on extracted g's | match |
| Tracking surface `R_phi = R_U` (eq -tracking) | sympy:74 / wl:59 `rho_0 - sigma_0 == 0` (interference ratios coincide ⇒ directions coincide) | match |
| Tracking factor `R_tr` + range `1/(1+δ_U)<R_tr<1` (eq -Rtr) | sympy:94-101 / wl:62-69 two range identities pinning `1-R_tr` and `R_tr-1/(1+δ_U)` | match |
| `chi_0 := rho_0 = sigma_0` (notes §3) | sympy:74 / wl:59 (the equality is verified; `chi_0` is then a defined symbol) | match |
| `M_tr = M_mix + M_supp` baseline (notes §4) | sympy:114-117 / wl:77-80 symbolic build, printed (carry-forward, not asserted) | match (printed) |
| `M_tr = G_tr`, D/N `G_tr` form (notes §5, card eq -tracking-laws) | sympy:153-169 / wl:107-117 quadratic collapse + `G_tr` D/N specialization | match |
| `R_target = F_tr`, D/N `F_tr` form (notes §6, card eq -tracking-laws) | sympy:189-199 / wl:138-149 Stage-044 F collapse + `F_tr` D/N form | match |

`paper_alignment: aligned` — every boxed/stated deliverable in card and notes has a faithful, non-tautological script-side check in both engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 63-66 | `expect_zero(g_* extracted - reference)` | g_* amplitude pattern (notes §1) | yes (two independent extraction routes) |
| A2 | sympy | 68 | `expect_zero("g_B g_R - g_W g_S")` | coherence condition (eq -coherent-condition) | yes |
| A3 | sympy | 74 | `expect_zero("rho_0 - sigma_0")` | tracking surface R_phi=R_U (eq -tracking) | yes |
| A4 | sympy | 94-101 | two range identities | R_tr factor + range (eq -Rtr) | yes |
| A5 | sympy | 153-156 | `expect_zero("tracking quadratic collapse")` | M_tr=G_tr collapse (notes §5) | yes |
| A6 | sympy | 160-163 | `expect_zero("G_tr formula")` | generic G_tr (notes §5) | yes |
| A7 | sympy | 169 | `expect_zero("G_tr D/N specialization")` | D/N G_tr form (notes §5) | yes |
| A8 | sympy | 189 | `expect_zero("Stage-044 tracking F collapse")` | generic F_tr (notes §6) | yes |
| A9 | sympy | 199 | `expect_zero("F_tr collapse from Stage-044 residual")` | D/N F_tr form (notes §6, eq -tracking-laws) | yes |
| B1-B9 | mathematica | 47-149 | `expectZero[...]` mirroring A1-A9 | same deliverables | yes (independent routes: `D[D[...]]`, `Series`) |

No assertion is orphaned; every load-bearing check traces to a card/notes deliverable. `M_mix`/`M_supp`/`M_tr` (sympy:118-120, wl:81-83) are printed only — legitimate carry-forward from Stages 022/026, flagged as such by inline comment; not a finding (the prefactor structure is verified upstream, and the values reconcile to the notes — see Value Reconciliation).

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.py:9` (docstring), `:125`, `:128` (banner), `:134` (comment)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage045_coherent_local_tracking_sympy_audit.txt:35`

**What's wrong:**
The SymPy script labels the section-4 source equation "the Stage-27 quadratic branch equation" in four places (docstring line 9, section comment line 125, `banner("4. Exact collapse of the Stage-27 branch equation")` line 128, comment line 134), and the committed transcript echoes it at output line 35 ("4. Exact collapse of the Stage-27 branch equation"). The quadratic branch equation this section actually imports and collapses is **Stage 044's** continuum-selected branch equation — the same script's later import comment (line 173) and the `expect_zero` labels (lines 189, 199) correctly call it "Stage-044", and the imported `D_cont_stage044`/`F_cont_stage044` forms (lines 174-182) match Stage 044's `D_cont`/`F_cont` verbatim. The notes (§5) and card (eq -tracking-laws) attribute the quadratic to Stage 044. "Stage-27" is a stale pre-renumber label (the +17 numbering-drift class: 27→044). Note the Mathematica script does **not** carry this label — it has no per-section banner and uses only "Stage-044" in its comment (wl:119,121) and labels.

**Why this matters:**
Pure label drift. The mathematics is correct and the load-bearing labels are right; this is a self-reference/documentation inconsistency, not a math defect. It is informational and matches the known SCRIPT/OUTPUT-band numbering-drift class that the user has explicitly deferred to a dedicated separate pass (`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`). It does not block the audit and no value is affected.

**Required change:**
Defer to the dedicated content-keyed numbering pass. The mechanically-correct content-keyed fix is `Stage-27` → `Stage-044` at sympy lines 9, 125, 128, 134 (and the transcript line 35 will refresh on the next run). This is NOT to be offset-swept; the content (the import is Stage 044's branch equation) keys the fix. No directive is issued for this finding (label-only, deferred class; see Verdict justification).

**Verification:**
After the deferred numbering pass, sympy lines 9/125/128/134 read "Stage-044" and a fresh transcript line 35 reads "4. Exact collapse of the Stage-044 branch equation". No assertion residual changes.

## Independent-derivation check (Mathematica)

The `.wl` is an independent re-derivation, not a transliteration. Two corresponding sections show genuinely different routes:

1. **g_* extraction.** SymPy (lines 49-52) reads bilinear coefficients via the chained `.coeff(W_sym).coeff(eta_sym)` accessor. Mathematica (lines 39-42) instead takes second cross-derivatives `D[D[couplingDensity, Wsym], etasym]`, with an explicit comment (wl:36-38) that this is an "Independent route ... not via the SymPy-style `.coeff(...).coeff(...)` chain." Different mechanism, same coefficients.

2. **Tracking collapse.** SymPy (line 143) collapses by direct substitution `branch_eq.subs(R_phi, R_U)`. Mathematica (lines 97-100) instead expands the numerator as a polynomial in `(rPhi - rU)` via `Series[branchNumRaw, {rPhi, rU, 0}]` and reads off the zeroth-order term, with a comment stating this is the independent route. Both yield the identical `numTrack`.

The Stage-044 `F_cont` import (sympy 174-182 / wl 122-131) is necessarily the same closed form in both — it is a carried-forward residual definition, not an algebra step, so identical form there does not indicate transliteration. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at the level claimed. Every `expect_zero`/`expectZero` residual is 0 in both transcripts. The printed intermediate forms match modulo cosmetic factoring:

- `rho_0 = sigma_0`: sympy `g_U*gamma*sqrt(mu_eta)/(K_U*sqrt(mu_U))` ≡ wl `(gamma*gU*Sqrt[muEta/muU])/KU`.
- `R_tr`: sympy `(chi_0+delta_U+1)/((chi_0+1)*(delta_U+1))` ≡ wl `(1+chi0+deltaU)/(1+chi0+deltaU+chi0*deltaU)` (same after expanding the denominator).
- `M_tr_req`: sympy `xi*(delta+xi)/(R_U**2*lambda_0*xi+delta+xi)` ≡ wl `(xi*(delta+xi))/(delta+xi+lambda0*rU^2*xi)`.
- tracking numerator/denominator, `M_tr`, and the final "coherent normalization residual" `R_target - F_tr` agree (the residual is printed, not asserted to zero, because `R_target` is the free physical target — correct; it is the test quantity, not an identity).

`engines_agree: true`.

## Verdict justification

`findings`, but only the single low-severity `stale_output`/stale-label F1, which is a label-only artifact of the known SCRIPT/OUTPUT-band numbering drift (SymPy "Stage-27" for what is Stage 044's branch equation) — user-deferred to a dedicated numbering pass, so no directive is written and the pipeline is not blocked. The mathematics is sound: I tried to break it on several fronts and failed. The coherence identity `g_B g_R = g_W g_S` is genuinely non-tautological — it holds because of the specific `gamma`-shared, same-source amplitude pattern, and the script extracts the g's two independent ways before testing it, so a sign/coefficient error in the kernel would surface. The `rho_0 = sigma_0` check has different `mu` structure on each side that must cancel (it does). The quadratic collapse and both `G_tr`/`F_tr` D/N specializations match the notes' boxed closed forms exactly, and the imported Stage-044 residual forms match Stage 044's script verbatim. Outputs are fresh (02:04 > 01:43 for both engines). Both engines derive the result via genuinely different routes. Paper card, notes, and appendix row all align with the verified claim.

## Self-test notes

I checked variable-independence (the `Series[branchNumRaw,{rPhi,rU,0}]` zeroth-order extraction is well-posed — `branchNumRaw` genuinely depends on `rPhi`, so the series is non-trivial; the `D[D[coupling,Wsym],etasym]` cross-derivatives are non-zero because `couplingDensity` is bilinear in those symbols). No unbounded-domain integrals appear, so parity is moot. Trivial-case spot check: the quadratic-collapse residual `num_track + collapsed_num` cancels term-by-term by hand (Mmix/Msupp/delta·xi/xi² all pair off), confirming the `assert_zero` is a real identity, not a vacuous pass. The single finding is informational/deferred, so no directive self-test for an applied edit was needed.

## Value Reconciliation (pass-2 augmentation)

Enumeration of every deliverable value/closed-form the scripts emit, located in `.tex` (`/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_045.tex`) and `.md` notes (`/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage045_coherent_local_tracking.md`):

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `g_B g_R - g_W g_S = 0` (coherence identity) | py:68 / wl:56; sympy out:13, wl out:13-14 | .tex:16-18 (eq -coherent-condition `g_B g_R = g_W g_S`); .md:78-91 | MATCH |
| `rho_0 = sigma_0 = g_U gamma sqrt(mu_eta)/(K_U sqrt(mu_U))` ⇒ `R_phi=R_U` | py:70-74 / wl:51-59; sympy out:14-16 | .tex:20-24 (eq -tracking `R_phi=R_U=R_tr`); .md:101,116-121 (`chi_0:=rho_0=sigma_0`) | MATCH |
| `R_tr = (1+chi_0/(1+delta_U))/(1+chi_0)` | py:86 / wl:53; sympy out:21 | .tex:26-33 (eq -Rtr); .md:124-126 | MATCH |
| range id `1-R_tr = chi_0 delta_U/[(1+chi_0)(1+delta_U)]` | py:89,94-97 / wl:62-65; sympy out:22,24 | .md:130-131 | MATCH |
| range id `R_tr-1/(1+delta_U) = delta_U/[(1+chi_0)(1+delta_U)]`; range `1/(1+δ_U)<R_tr<1` | py:90,98-101 / wl:66-69; sympy out:23,25 | .tex:31 (`1/(1+δ_U)<R_tr<1`); .md:133-139 | MATCH |
| `M_mix = 8 Z_W(1+chi_0)^2/[pi^2(1-eps_eta)(1-eps_W^split)]` | py:115 / wl:78; sympy out:30 | .md:149-151 | MATCH |
| `M_supp = 8 Z_phi(1+chi_0)^2/[pi^2(1-eps_eta)(1-eps_phi^split)]` | py:116 / wl:79; sympy out:31 | .md:152-154 | MATCH |
| `M_tr = 8(1+chi_0)^2/[pi^2(1-eps_eta)]·[Z_W/(1-eps_W^split)+Z_phi/(1-eps_phi^split)]` | py:117 / wl:80; sympy out:32 | .tex:40 (`M_tr=G_tr(...)`); .md:158-164 | MATCH |
| `M_tr_req = G_tr = xi(delta+xi)/[delta+(1+lam0 R_U²)xi]` (generic) | py:158-163 / wl:110-112; sympy out:40-41 | .md:199-204 (G_tr; D/N form) | MATCH |
| D/N `G_tr = 9 xi(delta+xi)/[9 delta+(9+2R_U²)xi]` | py:167-169 / wl:115-117; sympy out:42 | .md:202-204 | MATCH |
| D/N `F_tr = [9δ+(9+2R²)ξ]²[9δ+(9+2R)ξ]²/[81(1-ξ)(9δ²+18δξ+(9+2R²)ξ²)²]` | py:191-199 / wl:140-149; sympy out:44 | .tex:42 (`R_target=F_tr(...)`); .md:216-220 | MATCH |
| `lam0 (D/N) = 2/9` | py:166 / wl:114 | .md:202-204 (the `9`/`9+2R²` D/N coefficients encode lam0=2/9; consistent) | MATCH (encoded) |

**Internal scaffolding (no finding):** extracted-vs-reference g_* residuals (verification cross-checks, not deliverables); the printed final "coherent normalization residual = R_target - F_tr" (test quantity built from the matched `F_tr`, with `R_target` the free target — not an emitted constant); pass/PASS flags; the symbolic intermediate `couplingDensity`/`branch_eq`/`numTrack`/`denTrack` (drivers for assertions).

`reconciliation: complete; 12 values checked, 0 misaligned`
