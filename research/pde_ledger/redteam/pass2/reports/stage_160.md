---
unit_id: 160
batch: IV.6
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-08T00:00:00Z
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
  notes_stage_files: ["/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage160_bare_mixed_port_slippage.md"]
  paper_appendix: present
---

# Audit unit 160 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_160.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage160_bare_mixed_port_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (subsec `app-part04-bare-slippage`, lines 1064–1119; result anchor MTDC-T8.5)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.txt`

## What the paper claims

The card's `Derivation ledger` boxed result is the bottom-line claim:
`δγ_W = (δγ_0 − δκ_0/3)/(1+𝔯_{F1}²)`. The notes derive this in three logical steps: (1) the exact compensated-branch identity `δγ_W − (1/3)δκ_W = (δγ_0 − (1/3)δκ_0)/(1+r_{c,*})`, in which the `δr_c` dependence cancels identically; (2) under the Stage-159 canonical-even gate `δκ_W = 0` this collapses to `δγ_W = (δγ_0 − δκ_0/3)/(1+r_{c,*})`; (3) pure-scale harmlessness — if `δγ_0 = δκ_0/3` then `δγ_W = 0` (hence `Δ_Q = 0`, `N_Q−1 = 0`). The notes then push to the final defect law `Δ_Q = −9σ_*Υ_Π/[(1−σ_*)(1+𝔯_{F1}²)]·δΠ_tan` (and `N_Q−1` with `+`), with the tangential transport `δΠ_tan = 0.832409471081635 δΣ_0 − 1.16275838754222 δS`. The card's three `Checks` (deviations about the renormalized point, even-preservation before reading the odd defect, tangent motion gives `δ_⊥=0`) are satisfied conceptually by these steps. The appendix carries the same `δγ_W` formula (eq `app-part04-delta-gammaW-bare-slippage`) and the same `δΠ_tan` coefficients (eq `app-part04-deltaPi-tan`).

## What the script claims to verify

Both engines verify, symbolically in the free parameter `r_c` (= `r_{c,*}` = `𝔯_{F1}²`): (i) the exact compensated-branch slippage identity `δγ_W − (1/3)δκ_W − (δγ_0 − δκ_0/3)/(1+r_c) = 0`, with `δκ_W`/`δγ_W` obtained from `κ_0/(1+r_c)` and `γ_0/(1+r_c)` linearized about the canonical point `κ_{0,*}=(1+r_c)/3`, `γ_{0,*}=(1+r_c)/9`; (ii) the gate-reduced form `δγ_W = (δγ_0 − δκ_0/3)/(1+r_c)`; (iii) pure-scale harmlessness `δγ_W|_{δγ_0=δκ_0/3} = 0`; (iv) assembly of the final defect laws `Δ_Q`/`N_Q−1` after inserting `δΠ_tan` and the `Υ_Π` susceptibility. Steps (i)–(iii) are real `expect_zero`/`expectZero` assertions; step (iv) is print-only assembly.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Card box: `δγ_W = (δγ_0 − δκ_0/3)/(1+𝔯_{F1}²)` | sympy L54 / wl L47 print of `(dg0−dk0/3)/(1+r_c)`, anchored by the L50/L44 identity assertion | match (symbolic, r_c = 𝔯_{F1}²) |
| Notes step 1: exact identity `δγ_W − (1/3)δκ_W = …`, `δr_c` cancels | sympy L50–51 / wl L44–45 `expect_zero` of the full identity (residual 0) | match |
| Notes step 3 / card Check 3: pure-scale harmlessness `δγ_W=0` | sympy L58–59 / wl L49 `expect_zero` (residual 0) | match |
| Notes §5 final laws `Δ_Q`, `N_Q−1` in mouth vars | sympy L66–74 / wl L56–64 print-only assembly | partial (printed, not asserted) |
| Notes §4 `δΠ_tan` coefficients 0.832409471081635, 1.16275838754222 | sympy L66 / wl L56 literals | match |
| Notes numeric collapses (0.240311770175051; 1+𝔯²≈4.161261; 9/(1+𝔯²)≈2.162806; 11.5758539789133; etc.) | not computed by either engine (kept symbolic in r_c) | not script-verified (notes-only narrative) |

`paper_alignment: aligned`. Every load-bearing card/appendix identity has a non-tautological script assertion; the symbolic-vs-numeric gap is a notes-level narrative elaboration, not a card deliverable (the card and appendix box the symbolic form, which the scripts verify exactly).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50–51 | `expect_zero(dgW − dkW/3 − (dg0−dk0/3)/(1+rc))` | notes step-1 exact identity / card box | yes |
| A2 | sympy | 58–59 | `expect_zero(dgW_gate.subs(dg0, dk0/3))` | notes step-3 / card Check 3 | yes |
| A3 | mathematica | 44–45 | `expectZero[same identity]` | notes step-1 / card box | yes |
| A4 | mathematica | 49 | `expectZero[deltaGammaWGate /. dGamma0->dKappa0/3]` | notes step-3 / card Check 3 | yes |
| A5 | both | sympy 66–74 / wl 56–64 | print-only (`Δ_Q`, `N_Q−1`, `δΠ_tan`) | notes §4/§5 | n/a (print) |

A1/A3 are the load-bearing checks; A2/A4 the harmlessness corollary. All four are non-tautological: `δγ_W`/`δκ_W` are built independently (sympy: series-coefficient extraction; wl: hand-coded total differential) and the residual is checked against the separately-written RHS, so a sign/factor error in the linearization would surface a nonzero residual. A5 carries no assertion but mirrors the notes verbatim.

## Findings

None.

## Independent-derivation check (Mathematica)

The `.wl` is a genuinely independent route, not a transliteration:

- SymPy builds `δκ_W`/`δγ_W` by **automated series extraction**: `kW = ((k0_star + eps*dk0)/(1+rc+eps*drc)).series(eps,0,2).removeO()` then `dkW = kW.coeff(eps,1)` (py L41–45). It introduces a perturbation parameter `eps` and a `δr_c` drift `drc`, and lets SymPy compute the linear coefficient.
- Mathematica instead **hand-writes the total differential** in closed form: `deltaKappaW = Together[dKappa0/(1+rStar) − (kappa0Canon/(1+rStar)^2)*deltaR]` (wl L34–39). There is no series machinery, no `eps`, no `removeO`; the derivative is supplied analytically.

These are two different ways to obtain the linearization (machine series vs analytic differential), which is exactly the kind of cross-engine independence the policy wants. The downstream `expectZero` operands also differ in form (`Together[Expand[...]]` + `FullSimplify` vs SymPy `simplify(expand(...))`). Not a transliteration.

## Engine cross-check

Outputs agree at the level claimed. SymPy: `dκ_W = (−dr_c/3 + dκ0)/(r_c+1)`, `dγ_W = (−dr_c/9 + dγ0)/(r_c+1)`; both identity and harmlessness residuals print `0`. Mathematica: `dκ_W = −((deltaR − 3 dKappa0)/(3+3 rStar))` (= identical after factoring out 3), `dγ_W = −((deltaR − 9 dGamma0)/(9+9 rStar))` (identical), both `expectZero` print `0` + `PASS`. The final-law prints match symbolically; the only surface difference is that Mathematica pre-multiplied the `δΠ_tan` coefficients by 9 into the `Δ_Q`/`N_Q−1` numerators (`7.491685…`, `10.464825…` = 9×0.832409…, 9×1.162758…), which is the same expression SymPy leaves factored as `9·(…)`. No disagreement.

## Verdict justification

`clean`. I read the card, the notes, and the appendix subsection before the scripts. The card's boxed result and the appendix's `δγ_W` / `δΠ_tan` equations are exactly the symbolic identities the two engines assert as `expect_zero`/`expectZero` checks, and those checks are non-tautological (the perturbed `δκ_W`/`δγ_W` are constructed independently of the RHS they are subtracted against). Attacks tried and failed: (a) tautology — the LHS linearization and the RHS closed form are built from different expressions, so the residual could be nonzero if the factor-of-3/-9 or the `δr_c`-cancellation were wrong; the output confirms `δr_c` cancels exactly in `δγ_W − δκ_W/3`. (b) transliteration — the two engines use materially different linearization machinery (see above). (c) symbol-domain — all symbols are plain reals; no positivity/branch assumption is leaned on, and the only operations are rational-function simplifications that are valid for real symbols (the `1/(1+r_c)` pole at `r_c=−1` is outside the physical range `r_c=𝔯²>0` and is never evaluated numerically). (d) stale output — both `.txt` are newer than their scripts. (e) the `168π²`/`100π²` and `√(4107−100π²)/(10π)` warning constants do not appear in this stage; the only literals are the `δΠ_tan` coefficients, which match the notes and appendix verbatim. Minor non-finding: the `.wl` second block redeclares the radius symbol as `rc` (L54/L57) while the first block uses `rStar` (L31–49); both are independent free reals, so this is a cosmetic naming inconsistency with no math impact (it does not couple the two blocks, and each block's assertions stand alone).

## Self-test notes

Checked: (1) variable independence — the linearizations differentiate w.r.t. the perturbation directions `dk0`/`dg0`/`dr_c` that actually appear in `κ_0/(1+r_c)`, `γ_0/(1+r_c)`; no identically-zero derivative. (2) trivial-case — substituting `δγ_0=δκ_0/3` collapses the gate form to `0` as asserted (A2/A4), and the identity residual is `0` symbolically (A1/A3). (3) no integrals/parity here. (4) no missing-script directive needed. Conclusion: assertions exercise the card/notes claims faithfully; no directive written.

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 6 deliverable values checked, 0 misaligned.

Note: the scripts deliberately keep everything symbolic in `r_c` (= `𝔯_{F1}²`); they emit no numeric collapses. The notes' numeric values (0.240311770175051, 4.161261012190819, 2.16280593157546, 11.5758539789133, 2.51482073756543, 1.80034014155495, 6.42981496203006) are notes-only narrative elaborations of the symbolic result and are not produced by either engine, so they are not "values the scripts emit" and are excluded from the reconciliation per the procedure (they are not MISSING-DELIVERABLE — the boxed *symbolic* form is the stated card/appendix deliverable and it reconciles).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δγ_W − δκ_W/3 = (δγ_0 − δκ_0/3)/(1+r_c)` (symbolic identity, residual 0) | py L50–51 / wl L44–45; sympy.txt L7, math.txt L7–8 | card L16 box; notes L118–127; appendix L1080–1086 | MATCH |
| `δγ_W = (δγ_0 − δκ_0/3)/(1+r_c)` (gate-reduced) | py L54 / wl L47; sympy.txt L8, math.txt L9 | card L16; notes L149–157; appendix L1082–1086 | MATCH |
| pure-scale harmlessness `δγ_W = 0` | py L58–59 / wl L49; sympy.txt L9, math.txt L10–11 | notes L216–219; appendix L1087 | MATCH |
| `δΠ_tan = 0.832409471081635 δΣ_0 − 1.16275838754222 δS` | py L66 / wl L56; sympy.txt L14, math.txt L16 | notes L252–259; appendix L1113–1118 | MATCH |
| `Δ_Q = −9σ_*Υ_Π δΠ_tan/[(1−σ_*)(1+r_c)]` (symbolic) | py L68/L81 / wl L58; sympy.txt L16/L26, math.txt L18/L28 | notes L292–299; appendix L1104–1110 (in `Ξ_slip` form) | MATCH |
| `N_Q−1 = +9σ_*Υ_Π δΠ_tan/[(1−σ_*)(1+r_c)]` (symbolic) | py L69/L82 / wl L59; sympy.txt L17/L27, math.txt L19/L29 | notes L302–309; appendix L1104–1110 | MATCH |

INTERNAL scaffolding (no finding): `eps` perturbation parameter, `drc`/`deltaR` drift symbol, intermediate `kW`/`gW` series, `dgW_gate`, `dgW_tan`/`dgWTan`, the `banner`/`pass`/`fail`/`expectZero` helpers, the 9×-prefactored numeric coefficients `7.491685239734715`/`10.464825487879981` (= 9× the `δΠ_tan` coefficients; Mathematica `FullSimplify` distributed the 9, same expression SymPy keeps factored).
