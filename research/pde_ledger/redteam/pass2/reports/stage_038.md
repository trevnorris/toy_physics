---
unit_id: 038
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
  outputs_fresh: false
docs_read:
  paper_stage_tex: present
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md"]
  paper_appendix: present
---

# Audit unit 038 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_038.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage038_dimensionless_continuum_placement.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (rows: line 54 ledger entry; line 309 summary-chain; line 345 verificationchecklist)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.txt`

## What the paper claims

Stage 038 compresses the Stage-037 microscopic continuum constants into five dimensionless kernel ratios `(eps_eta, eps_W, rho, Z_W, delta_0)` plus a radiative demand scale `Lambda := 27 pi^2 G c_s^5 K_W^eff/(20 a^5 c^5 mu_W)`, and proves the exact placement map. The card's `\stagefield{Output}` reads verbatim: "The placement map \eqref{eq:app-stage038-placement-map}--\eqref{eq:app-stage038-Rtarget} and product law \eqref{eq:app-stage038-product-law}." The boxed deliverables are: `delta = delta_0/(1-eps_eta)`; `M_mix = 8 Z_W(1+rho)^2/[pi^2(1-eps_eta)(1-eps_W)]`; `R_target = Lambda(1-eps_eta)(1-eps_W)^2/[Z_W(1+rho)^2]`; and the product law `R_target M_mix = 8 Lambda(1-eps_W)/pi^2 = 54 G c_s^5 K_W^eff(1-eps_W)/(5 a^5 c^5 mu_W)`. The `\stagefield{Checks}` line additionally requires product-law invariance under redistributions of `(Z_W, rho, eps_eta)` that leave the product lane fixed. The notes add a corollary `beta_0 = (mu_W/mu_eta)(K_eta^eff/K_W^eff) Z_W(1+rho)^2/(1-eps_W)^2` (notes §2, lines 97–99) and the full set of one-way monotone tendencies (notes §4) as derived, but the .tex card does not box `beta_0` or the tendencies as stage Output.

## What the script claims to verify

Both scripts take the Stage-037 microscopic forms for `A`, `delta`, `M_mix`, `beta0`, `R_target` (SymPy lines 54–66; .wl 47–60), perform a change of variables to the dimensionless ratios via the definitional substitutions `c_etaU^2 → eps_eta K_U Keta_eff`, `c_UW^2 → eps_W K_U KWeff/sigma`, `c_UW c_etaU → rho K_U c_etaW`, `c_etaW^2 → Z_W Keta_eff KWeff`, `T_w → delta0 L^2 Keta_eff/pi^2`, `G → 20 Lambda a^5 c^5 mu_W/(27 pi^2 c_s^5 KWeff)`, and assert (via `expect_zero`/`expectZero`) that the substituted forms equal the boxed dimensionless placement formulas (SymPy 109–121; .wl 102–114). They then verify the product relation in both its `8 Lambda(1-eps_W)/pi^2` and microscopic-`N_Q` forms (SymPy 134–143; .wl 123–133), and factor all nine one-way derivatives `d{M,R,delta}/d{eps_eta,eps_W,Z_W,rho}` against their closed forms, plus sign coefficients (SymPy 157–217; .wl 143–181). `sigma = 88/(9 pi^2)` matches the Stage-037 notes (line 35/107).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `delta = delta_0/(1-eps_eta)` | SymPy 109 / .wl 102 `expect_zero(delta_dimless - delta0/(1-eps_eta))` | match |
| `M_mix = 8 Z_W(1+rho)^2/[pi^2(1-eps_eta)(1-eps_W)]` | SymPy 110–113 / .wl 103–106 | match |
| `R_target = Lambda(1-eps_eta)(1-eps_W)^2/[Z_W(1+rho)^2]` | SymPy 114–117 / .wl 107–110 | match |
| product law `R_target M_mix = 8 Lambda(1-eps_W)/pi^2` | SymPy 135 / .wl 124 | match |
| product 2nd form `= 54 G c_s^5 K_W^eff(1-eps_W)/(5 a^5 c^5 mu_W)` | SymPy 141–143 / .wl 125–133 (G-substitution bridge) | match |
| `Checks`: invariance under `(Z_W,rho,eps_eta)` redistribution | implicit: product is independent of `Z_W,rho,eps_eta` (provable from the product expression, which contains only `Lambda,eps_W`); the nine derivative factorizations confirm the lane structure | partial (invariance is implied by the product form, never asserted as its own check) |
| `beta_0` corollary (notes §2) | SymPy 118–121 / .wl 111–114 | match (notes deliverable) |
| one-way tendencies (notes §4) | SymPy 157–217 / .wl 143–181 | match (notes deliverable) |

`paper_alignment: aligned`. Every boxed Output of the .tex card has a faithful, non-tautological script-side check; the only `partial` row is the prose `Checks` invariance statement, which is a logical consequence of the verified product form rather than a separate identity, so no script-side gap exists.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 109 | `expect_zero(delta_dimless - delta0/(1-eps_eta))` | delta map | yes |
| A2 | sympy | 110–113 | `expect_zero(M_dimless - 8 Z_W(1+rho)^2/...)` | M_mix map | yes |
| A3 | sympy | 114–117 | `expect_zero(R_dimless - Lambda(1-eps_eta)(1-eps_W)^2/...)` | R_target map | yes |
| A4 | sympy | 118–121 | `expect_zero(beta_dimless - ...)` | beta_0 corollary (notes) | yes |
| A5 | sympy | 135 | `expect_zero(product - 8 Lambda(1-eps_W)/pi^2)` | product law | yes |
| A6 | sympy | 141–143 | `expect_zero(8 Lambda.../pi^2 - NQ form)` | product 2nd form | yes |
| A7 | sympy | 167–175 | 9× derivative factorizations | tendencies (notes) | yes |
| A8 | sympy | 182–217 | 9× sign-coefficient checks | tendency signs (notes) | partial (see F-note) |
| B1–B8 | mathematica | 102–181 | mirror of A1–A8 | same | yes / partial |

The "partial" on A8 is the only soft row; it is discussed in Verdict justification — it is a weaker-than-positivity check but not a genuine defect because each template is sign-definite and the factorization (A7) carries the substantive content.

## Findings

### F1 — stale_output

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.txt:3` (and 61)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage038_dimensionless_continuum_placement_mathematica_audit.txt:3`

**What's wrong:**
Both saved transcripts have mtime `2026-05-22 12:23` while both scripts have mtime `2026-06-03 15:59` (committed `2026-06-03 16:12`) — the outputs predate the current scripts by ~12 days. The staleness is visible in the banners: the SymPy script now prints `STAGE 38 — DIMENSIONLESS CONTINUUM PLACEMENT AUDIT` (script line 36) but the committed transcript still reads `STAGE 21 — DIMENSIONLESS CONTINUUM PLACEMENT AUDIT` (output line 3) and `STAGE 21 AUDIT COMPLETE` (output line 61). The Mathematica script now prints `STAGE 038 — ...` (script line 33) but the transcript reads `STAGE 021 — ...` (output line 3), while inconsistently ending `Stage 038 Mathematica audit passed.` (output line 84) — a mixed-epoch label state. The math content of the transcripts (the substituted forms, every `= 0` residual, every PASS) is fully consistent with what the current scripts would emit; only the banner/label strings are stale.

**Why this matters:**
The committed transcript does not faithfully reflect the current script's printed banner, so a reader trusting the saved output would see the wrong stage number. The underlying mathematics is unaffected, so this is informational; a fresh re-run regenerates correct banners.

**Required change:**
Re-run both scripts and overwrite the committed `.txt` transcripts so the banner labels match the current scripts (`STAGE 38`/`STAGE 038`). No source-math edit is required. (The orchestrator's independent re-run already refreshes these per project policy.)

**Verification:**
After re-run, output line 3 of the SymPy transcript reads `STAGE 38 — ...` and line 61 reads `STAGE 38 AUDIT COMPLETE`; the Mathematica transcript line 3 reads `STAGE 038 — ...`. All residual/PASS lines remain `= 0` / `PASS`.

## Independent-derivation check (Mathematica)

The `.wl` is a close structural mirror of the `.py`: identical symbol roster, identical microscopic seed forms (`a`, `delta`, `mMix`, `beta0`, `rTarget` at .wl 47–60 vs. SymPy 54–66), the same substitution map, the same eight assertion groups in the same order, and the same final-print strings. It is borderline against the second-engine policy. However, it is NOT a blind transliteration: the dimensionless-substitution routine genuinely diverges. SymPy's `apply_dimless` (lines 93–102) handles only `c_UW**2`, `c_etaU**2`, `c_etaW**2`, `c_UW*c_etaU`; the Mathematica `applyDimless` (lines 77–95) must additionally and independently introduce `cUW^4 -> (epsW kU kWEff/sigma)^2` and `cEtaW^(-2) -> 1/(zW kEtaEff kWEff)` rewrite rules plus `PowerExpand[Factor[...]]` choreography to get Mathematica's pattern matcher to reach the same normal form. That is an engine-specific re-derivation of the same change of variables, not a line-for-line echo. The two engines also reach visibly different simplified normal forms in the transcripts (e.g. SymPy `M_mix(dimless) = 8*Z_W*(rho**2+2*rho+1)/(pi**2*(eps_W*eps_eta - eps_W - eps_eta + 1))` vs. Mathematica `(8*(1+rho)^2*zW)/((-1+epsEta)*(-1+epsW)*Pi^2)`), confirming they ran independent simplifiers. I do not raise `mathematica_transliteration`; the structural parallelism is expected for a pure-algebra identity stage and the algebraic route differs per engine.

## Engine cross-check

Both engines emit identical residual structure: every dimensionless-map check, the product relation (both forms), the nine derivative factorizations, and the nine sign-coefficient checks all reduce to `0` (SymPy transcript lines 18–54; Mathematica transcript lines 18–78 with matching `PASS`). The seed microscopic forms agree up to simplifier normal-form (SymPy lines 9–13 vs. Mathematica lines 9–13). No sign, factor-of-2, or `sigma` discrepancy. Engines agree.

## Verdict justification

Verdict `findings` with a single low-severity informational `stale_output`. Attacks tried and failed: (1) tautology — the dimensionless-map assertions are NOT tautological, because the LHS is the Stage-037 microscopic form pushed through a definitional change of variables and the RHS is the independent boxed paper form; the substitution algebra is non-trivial (the SymPy/Mathematica transcripts show the pre-simplification microscopic seeds at lines 9–13). (2) hidden hardcoding — `sigma = 88/(9 pi^2)` and the `R_target` seed `54 G c_s^5/(5 a^5 c^5)` both trace to the Stage-037 notes (lines 35, 24–28), not invented here. (3) symbol-domain attack — symbols are `positive, real`; `eps_eta, eps_W` are declared positive but NOT `< 1`, yet no assertion's validity depends on `1-eps > 0` (the `expect_zero` residuals are rational identities that hold for all values where denominators are nonzero), so the missing strict upper bound does not invalidate any check. (4) the sign-coefficient checks (A8) are weaker than a strict positivity proof — they confirm each derivative factors as `(±1)·template` but lean on prose (script comments) for the template's positivity rather than asserting `1-eps_eta > 0` symbolically; this is the standard idiom across the ledger and the substantive content (the factorization, A7) is fully verified, so it is noted, not filed. I read the .tex card, the notes, and the part-III appendix rows; every boxed Output deliverable has a faithful, exercised, non-tautological check in both engines, and the product law's redistribution-invariance `Checks` requirement follows directly from the verified product form. Paper alignment is exact.

## Self-test notes

Variable-independence: every `sp.diff`/`D` in §4 differentiates the closed dimensionless forms (`delta_expr`, `M_expr`, `R_expr` at SymPy 153–155) with respect to a variable each genuinely depends on — `delta_expr` depends on `eps_eta` only (and is only differentiated w.r.t. `eps_eta`); `M_expr`/`R_expr` depend on `eps_eta, eps_W, Z_W, rho` and are differentiated w.r.t. each, so no derivative is identically zero by accident. No unbounded integrals (parity trap n/a). Trivial-case: each dimensionless-map residual is a rational identity in the kernel ratios that holds for all admissible substitutions; the product residual `R_target M_mix - 8 Lambda(1-eps_W)/pi^2` collapses to 0 because `Z_W, rho, eps_eta` cancel pairwise between `M_mix` and `R_target` — confirmed by inspection of the boxed forms. Paper round-trip: the only required change (F1) is a transcript refresh, which introduces no new constant and cannot create a paper_misalignment.

## Value Reconciliation (pass-2 augmentation)

Both saved transcripts are stale-labelled (banner only; see F1); the math lines in them are consistent with the current scripts, so I reconcile against the script source plus the transcript residuals. No execution performed.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `delta = delta_0/(1-eps_eta)` | py:109 / wl:102; sympy out:18, math out:18–19 | .tex:34 (eq placement-map); .md:87 | MATCH |
| `M_mix = 8 Z_W(1+rho)^2/[pi^2(1-eps_eta)(1-eps_W)]` | py:110–113 / wl:103–106; out:19/20–21 | .tex:36; .md:89 | MATCH |
| `R_target = Lambda(1-eps_eta)(1-eps_W)^2/[Z_W(1+rho)^2]` | py:114–117 / wl:107–110; out:20/22–23 | .tex:41; .md:91 | MATCH |
| product `R_target M_mix = 8 Lambda(1-eps_W)/pi^2` | py:135 / wl:124; out:30,32 / 34–35,38 | .tex:49 (1st eq); .md:38,183 | MATCH |
| product 2nd form `= 54 G c_s^5 K_W^eff(1-eps_W)/(5 a^5 c^5 mu_W)` | py:141–143 / wl:125–133; out:31 / 36–37 | .tex:49 (2nd eq); .md:39,111 | MATCH |
| `beta_0 = (mu_W/mu_eta)(K_eta^eff/K_W^eff) Z_W(1+rho)^2/(1-eps_W)^2` | py:118–121 / wl:111–114; out:21/24–25 | .md:97–99 (card omits — terse, allowed) | MATCH |
| `sigma = 88/(9 pi^2)` | py:50 / wl:43 | .md:(stage037):35,107; .tex via Lambda def | MATCH (carried from 037) |
| `Lambda := 27 pi^2 G c_s^5 K_W^eff/(20 a^5 c^5 mu_W)` | py:90 / wl:91 (inverse G-subst); out (in 2nd-form bridge) | .tex:28; .md:32,69 | MATCH |
| 9 derivative factorizations `d{M,R,delta}/d{...}` | py:167–175 / wl:153–161 | .md:§4 (135–164) | MATCH (notes deliverable) |
| 9 sign coefficients `±1` | py:182–217 / wl:164–181 | .md:§4 inequality directions | MATCH (notes deliverable) |

INTERNAL (scaffolding, no prose expectation): `A` seed (py:54), microscopic `delta`/`M_mix`/`beta0`/`R_target` seed forms (py:55–66) — intermediate Stage-037 inputs, not Stage-038 deliverables; `Keta_eff`, `KWeff` (py:51–52) — Stage-037 effective stiffnesses; `subs_dimless` dict and `apply_dimless` rewrites (py:84–102); `NQ`, `OmegaU2`, `Delta0` symbols (py:138–140) — dead/unused scaffolding for the `N_Q` bridge label; all `expect_zero`/`PASS` flags and `= 0` residuals.

reconciliation: complete; 10 deliverable values checked, 0 misaligned.
