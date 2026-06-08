---
unit_id: 163
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
  notes_stage_files: [moving_throat_pde_stage163_off_family_normal_coordinate.md]
  paper_appendix: present
---

# Audit unit 163 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_163.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage163_off_family_normal_coordinate.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows: §subsec:app-part04-off-family-normal-coordinate, lines 1121–1168; stage-range/listing lines 32, 80–86, 1179)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.txt`

## What the paper claims

Stage 163 is a Part IV checkpoint (anchor MTDC-T8.5) that reduces the entire first-order off-family defect of the co-evolving core–mouth branch to a single normal coordinate. The card's `\stagefield{Purpose}` is terse ("a co-evolving core–mouth transport ledger step. Its audit target is the verification output quoted below"); the body quote states "The first true off-family drift is the scalar `\delta_\perp`, with explicit microscopic parent-action formula." The notes are authoritative and enumerate the deliverables: (1) the parent compensation function `F = 1+r² − 4(g−r)²` and its zero set being the Stage 119 compensation family; (2) the lower branch `g_-(r) = r − ½√(1+r²)` with slope `g_-'(r) = 1 − r/(2√(1+r²))`; (3) the normal coordinate `δ_⊥ := δg − g_-'(r*)δr` and the linearization identity `δF = 4√(1+r*²) δ_⊥`; (4) the compensation-ratio transport `δR_q = −δ_⊥/√(1+r*²)`; (5) the exact microscopic formula `δ_⊥ = g* δln(g_qK_s/g_sλ) + [4√(1+r*²)]⁻¹ δln(K_sK_q/λ²)` (the boxed "final output of Part IV"); (6) outlet transport `δC = 16σ*δ_⊥/√(1+r*²)`, `δE2`, `δE4`, `Δ_Q`; (7) the tangent/normal split of `δΠ`; and (8) numeric Family-1 coefficients. The card `\stagefield{Checks}` requires: deviations taken about the renormalized canonical point; even-preservation imposed before reading the odd defect; tangent motion on the parent family gives `δ_⊥ = 0`.

## What the script claims to verify

The SymPy script (docstring, 5 numbered checks) verifies: (1) `δF − 4s·δ_⊥ = 0` and `δR + δ_⊥/s = 0` plus that tangent motion `dg = g_-'·dr` keeps `δF = δR = 0`; (2) the microscopic `δ_⊥` identity (the linearized `δg − g_-'·δr` equals the boxed two-channel log form); (3) the outlet transport identity `δC − 4σ*·δF/(1+r²) = 0` with `δC, δE2, δE4, Δ_Q` printed as closed forms; (4) the `δΠ` tangent/normal split `δΠ − δΠ_expected = 0`; (5) Family-1 numeric coefficient readbacks. The Mathematica script verifies the same identities AND adds two independent routes the SymPy lacks: the implicit-function-theorem derivation of the slope `g_-'(r) = −(∂F/∂r)/(∂F/∂g)|_{g=g_-}` cross-checked against `D[gMinus,r]`, and a full `Series[..., {eps,0,1}]` perturbation derivation of `δr`/`δg` from the parent log-variations, re-deriving `δ_⊥` independently.

## Paper ↔ script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| `δF = 4√(1+r*²) δ_⊥` (notes 76–81; appendix eq deltaF-deltaperp) | py:60 / wl:45 `δF − 4s δ_⊥ = 0` | match |
| `δR_q = −δ_⊥/√(1+r*²)` (notes 99–105) | py:61 / wl:46 `δR + δ_⊥/s = 0` | match |
| Tangent motion ⇒ `δ_⊥=0` ⇒ `δF=δR=0` (card check 3) | py:62–63 (subs dg→gp·dr) | match |
| slope `g_-'(r) = 1 − r/(2√(1+r²))` (notes 55–60) | py:64 print; wl:38–40,47 implicit-fn cross-check | match |
| microscopic `δ_⊥` two-channel form (notes 165–173; appendix boxed) | py:79–83 / wl:55–58 + series route wl:61–81 | match |
| `δC = 16σ*δ_⊥/√(1+r*²) = 4σ*δF/(1+r*²)` (notes 227–235) | py:91,103–106 / wl:88,97 | match |
| `δE2, δE4, Δ_Q` outlet forms (notes 248–281) | py:94–96 / wl:89–91 (printed closed forms) | match |
| `δΠ` tangent/normal split (notes 314–327) | py:115–117 / wl:104–106 | match |
| Family-1 numeric coefficients (notes 112,119,126,193,196,242,337,339,352) | py:129–141 / wl:116–121 | match |

paper_alignment: aligned.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 60 | `expect_zero(δF − 4s δ_⊥)` | deliverable 3 (δF identity) | yes |
| A2 | sympy | 61 | `expect_zero(δR + δ_⊥/s)` | deliverable 4 (δR_q) | yes |
| A3 | sympy | 62 | `expect_zero(δF|dg=gp·dr)` | card check 3 (tangent ⇒ 0) | yes |
| A4 | sympy | 63 | `expect_zero(δR|dg=gp·dr)` | card check 3 | yes |
| A5 | sympy | 83 | `expect_zero(δ_⊥micro − δ_⊥expected)` | deliverable 5 (microscopic form) | yes |
| A6 | sympy | 103 | `expect_zero(δC − 4σ*δF/(1+r²))` | deliverable 6 (outlet) | yes |
| A7 | sympy | 117 | `expect_zero(δΠ − δΠ_expected)` | deliverable 7 (δΠ split) | yes |
| B1 | math | 40 | `expectZero(gPrime − gPrimeImplicit)` | slope (independent implicit-fn route) | yes |
| B2 | math | 45 | `expectZero(δF − 4s δ_⊥)` | deliverable 3 | yes |
| B3 | math | 46 | `expectZero(δR + δ_⊥/s)` | deliverable 4 | yes |
| B4 | math | 58 | `expectZero(δ_⊥micro − δ_⊥expected)` | deliverable 5 | yes |
| B5 | math | 76 | `expectZero(δr series − hand form)` | deliverable 5 (independent series route) | yes |
| B6 | math | 77 | `expectZero(δg series − hand form)` | deliverable 5 (independent series route) | yes |
| B7 | math | 79 | `expectZero(δ_⊥ series − two-channel form)` | deliverable 5 (independent series route) | yes |
| B8 | math | 97 | `expectZero(δC − 4σ*δF/(1+r²))` | deliverable 6 | yes |
| B9 | math | 106 | `expectZero(δΠ − δΠ_expected)` | deliverable 7 | yes |

No print-only deliverable is left unguarded: the `δE2/δE4/Δ_Q` closed forms (py:99–101, wl:93–96) are derived from `δC`, which is itself pinned by the asserted `δC − 4σ*δF/(1+r²) = 0`; their `5·deltaC` / `16` / `27`-factor structure matches the notes coefficients (80/√, 16/√, 16/√ respectively). Numeric readbacks (Family-1 block) are summaries, not deliverables proven in-script.

## Findings

None. No directive written.

### Higher-bar checkpoint analysis (load-bearing constant re-derivation)

The load-bearing constant is the linearization coefficient `4√(1+r*²)` in `δF = 4√(1+r*²) δ_⊥`, together with the slope `g_-'(r) = 1 − r/(2√(1+r²))`. I re-derived both from the script logic without trusting any literal:

- On the branch `g = g_- = r − s/2` (`s = √(1+r²)`): `g − r = −s/2`, so `F = 1+r² − 4·(s²/4) = 1+r² − (1+r²) = 0` — the compensation family holds (non-vacuous: it is exactly zero, not an asserted literal).
- `∂F/∂g = −8(g−r)`; at `g=g_-`: `−8·(−s/2) = 4s`. `∂F/∂r = 2r + 8(g−r)`; at `g=g_-`: `2r − 4s`. Hence `δF = 4s·dg + (2r−4s)·dr`.
- `4s·δ_⊥ = 4s(dg − g_-'·dr) = 4s·dg − 4s(1 − r/(2s))dr = 4s·dg − (4s − 2r)dr = 4s·dg + (2r−4s)dr = δF`.

So `δF − 4s·δ_⊥ ≡ 0` is a genuine algebraic identity, not `X − X`. The coefficient `4√(1+r*²)` is genuinely derived from `∂F/∂g` on the branch, and the printed Family-1 value `8.15966765224253` equals `4·√(1+1.77799353547498²)` (confirmed against notes:126). The slope is independently re-derived in Mathematica by the implicit-function theorem (wl:38–40), an algebraically distinct route. The checkpoint clears the higher bar.

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. It re-uses the same final identities (as both must, since they target the same physics) but reaches `δ_⊥` by two routes absent from the SymPy script:

1. Slope via implicit-function theorem (wl:38–40):
   `gPrimeImplicit = -((D[fComp, r])/(D[fComp, g])) /. g -> gMinus;`
   `expectZero["gPrime matches implicit-function route", gPrime - gPrimeImplicit];`
   SymPy only computes `gp = sp.diff(gminus, r)` (py:51) — direct differentiation, never the `−F_r/F_g` route.

2. Microscopic `δ_⊥` via series perturbation (wl:61–81):
   `pertRule = {Ks -> Ks (1 + eps dlnKs), ...}; deltaRSeries = Coefficient[Series[rExpr /. pertRule, {eps, 0, 1}] // Normal, eps];`
   This derives `δr`, `δg` by expanding `r = λ/√(K_sK_q)` and `g = g_q√K_s/(g_s√K_q)` to first order in a perturbation parameter, then re-derives `δ_⊥` (wl:79). SymPy hardcodes the linearized `delta_r`, `delta_g` directly (py:75–76) and never builds them from the defining ratios. The Mathematica series route therefore independently confirms the linearization the SymPy script assumes — a strictly stronger second engine. This is the opposite of the transliteration defect.

## Engine cross-check

Both engines agree on every asserted identity (all `= 0` / all PASS in transcripts). The printed symbolic forms agree up to cosmetic ordering:
- `g_-'(r)`: SymPy `−r/(2√(r²+1)) + 1`; Mathematica `1 − r/(2√(1+r²))` — identical.
- `δC`: SymPy `8σ*(2dg√(r²+1) + dr(r − 2√(r²+1)))/(r²+1)`; Mathematica `16(dg + (dr(−2 + r/√(1+r²))/2))σ*/√(1+r²)` — algebraically identical (factor `16/√` = `8·2/√` after rationalizing).
- Family-1 coefficients agree to 20 sig figs: `4√(1+r*²) = 8.1596676522425…`, `16/√(1+r*²) = 7.8434567102020…`, `Σ0_can·S_can/√ = 1.5284331782324…`. No `engine_disagreement`.

## Verdict justification

Clean. I attacked: (a) tautology — the seven `expect_zero`/nine `expectZero` checks each subtract two independently-constructed expressions (`δF` built from `∂F/∂g,∂F/∂r` vs `4s·δ_⊥` built from the slope; microscopic `δ_⊥` from substituting linearized ratios vs the boxed two-channel form), so none is `X−X`; the Mathematica series and implicit-function routes make the microscopic and slope checks doubly non-tautological. (b) hardcoded result — the load-bearing `4√(1+r*²)` is derived from `∂F/∂g` on the branch, not pinned as a literal; the numeric `8.159…` is a printed readback consistent with `4·√(1+r*²)`. (c) symbol assumptions — `r, g, dg, dr` real with `r>0, g>0` is correct for the parent ratios; positivity does not over-strengthen any simplify (the identities are polynomial/rational, not branch-sensitive). (d) tangent-motion check (card check 3) — `dF.subs(dg, gp*dr)` genuinely substitutes the tangency condition and yields 0, exercising the `δ_⊥=0` claim. (e) numbering — notes line 45 correctly cites "Stage 119" (the Step-0 fix held); the card `\stagefield{Purpose}` correctly reads "Stage~163" with no +17 (=180) self-label garble. (f) stale constants — no `168π²`/`100π²`; the Family-1 radius `1.77799353547498` equals the canonical `√(4107−100π²)/(10π)`. Outputs are fresher than their scripts. Paper card + notes + appendix and the script's verified claim match exactly. The checkpoint clears the higher bar on all four counts: both engines present and independent, assertions substantive and non-tautological, paper-aligned, and the load-bearing constant re-derived from script logic.

## Self-test notes

Checked: (1) variable independence — every `sp.diff`/`D` is taken w.r.t. a variable the expression genuinely depends on (`F`, `R`, `gMinus`, `rExpr`, `gExpr` all depend on `r`/`g`; the `Series` is in `eps` on which all perturbed symbols depend), so no derivative is identically zero and no `expect_zero` passes vacuously. (2) No unbounded integrals in this stage — parity trap N/A. (3) Trivial-case: on the branch `g−r=−s/2` makes `F=0` and `δF−4s·δ_⊥` collapses to `0` as re-derived above; tangent substitution `dg=gp·dr` zeroes `δ_⊥` and hence `δF,δR`. (4) Paths confirmed: `.py` in `scripts/`, `.wl` in `mathematica/`. No fix prescribed, so no paper round-trip needed.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value the scripts emit was located in the docs. All reconcile.

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `δF = 4√(1+r*²) δ_⊥` (symbolic) | py:60, wl:45 / out py:5, wl:7 | notes:76–81; appendix:1154 | MATCH |
| `δR_q = −δ_⊥/√(1+r*²)` (symbolic) | py:61, wl:46 / out py:6, wl:9 | notes:99–105 | MATCH |
| `g_-'(r) = 1 − r/(2√(1+r²))` | py:64, wl:47 / out py:9, wl:11 | notes:55–60 | MATCH |
| microscopic `δ_⊥` two-channel form | py:84, wl:59 / out py:11, wl:14 | notes:165–173, 380–388; appendix:1159–1166 | MATCH |
| `δC = 16σ*δ_⊥/√(1+r*²)` (symbolic) | py:98, wl:93 / out py:12, wl:21 | notes:227–235 | MATCH |
| `δE2, δE4, Δ_Q` outlet closed forms | py:99–101, wl:93–96 / out py:13–15, wl:22–24 | notes:248–281 | MATCH |
| `δΠ` tangent/normal split (symbolic) | py:118, wl:107 / out py:18, wl:29 | notes:314–327 | MATCH |
| `4√(1+r*²) = 8.15966765224253` | py:136, wl:116 / out py:23, wl:34 | notes:126 | MATCH |
| `−1/√(1+r*²) = −0.490216044387626` | py:137, wl:117 / out py:24, wl:35 | notes:119 | MATCH |
| `g_* = 0.758035078944663` | py:138, wl:118 / out py:25, wl:36 | notes:193, 203 | MATCH |
| `1/(4√(1+r*²)) = 0.122554011096907` | py:139, wl:119 / out py:26, wl:37 | notes:196, 206 | MATCH |
| `Σ0_can·S_can/√(1+r*²) = 1.52843317823248` | py:140, wl:120 / out py:27, wl:38 | notes:352 | MATCH |
| `16/√(1+r*²) = 7.84345671020202` | py:141, wl:121 / out py:28, wl:39 | notes:242 | MATCH |
| `r_* = 1.77799353547498` (Family-1) | py:123, wl:109 | notes:112, 341 | MATCH |
| `Σ0_can = 4.651033550168876` | py:125, wl:111 | notes:337 | MATCH |
| `S_can = 0.6703621156734617` | py:126, wl:112 | notes:339 | MATCH |

INTERNAL (verification scaffolding, no prose expected): `s = √(1+r²)`, `gminus`, `gp`, `F`, `R`, intermediate `dF`/`dR`, `delta_r`/`delta_g`, `pertRule`/`eps` series intermediates, `Rstar = 1/4`, the four "even-preserving tangent motion keeps δF/δR = 0" residuals, and all PASS flags.

reconciliation: complete; 16 deliverable values checked, 0 misaligned
