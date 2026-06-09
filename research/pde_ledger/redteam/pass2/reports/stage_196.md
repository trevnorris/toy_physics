---
unit_id: 196
batch: V.3
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-09T00:00:00Z
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
  notes_stage_files: [moving_throat_pde_stage196_higher_odd_irrelevance_theorem.md]
  paper_appendix: present
---

# Audit unit 196 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_196.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage196_higher_odd_irrelevance_theorem.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part05.tex` (rows at lines 123, 1318, 1467, 1487, 1569)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_mathematica_audit.txt`

## What the paper claims

The card's `\stagefield{Output}` states: *"Proves all extra isotropic odd tails beginning at \(O(\omega^7)\) are invisible to the point-particle \(2.5\)PN theorem."* The notes (§7) enumerate five concrete deliverables: (1) an exact response-side higher-odd difference identity showing any extra retarded tail `O(ω^7)` changes the grouped response only at `O(ω^7)`; (2) the same at the DtN operator level (`Yhat_2^def,≥7 − Yhat_2^def,5 = −L0·L_{≥7}/(D5·D_{≥7})`); (3) that the Stage-194 outgoing-normalization compiler `χ_Q = 3(Sβ^5+9Σ_5)/(3S−Σ_0)` is unchanged by all higher-odd DtN data beginning at `O(z^7)`; (4) the Packet-A consequence that the Stage-195 source-map-reduced residual `Δ_norm = P0^target(1/χ_Q − 1)` is likewise unchanged; (5) the final reduced statement that the only live retarded obstruction at 2.5PN is `Δ_Q := χ_Q − 1`. The card is intentionally terse ("original-stage audit card … records the claim boundary, not a second independent proof"); the notes are authoritative on the symbolic content.

## What the script claims to verify

The SymPy script (banner "STAGE 196 — EXACT HIGHER-ODD IRRELEVANCE THEOREM") verifies all five deliverables symbolically: §I builds the grouped one-pole module with an explicit `i·τ_Q·ω^7` tail and asserts the closed-form difference identity plus that the difference series starts exactly at `ω^7` (response unchanged through `ω^5`); §II does the same at the DtN level with an `i·L7·z^7` tail; §III imposes the Stage-194 canonical-even matching and asserts `χ_Q` equals `3(Sβ^5+9Σ_5)/(3S−Σ_0)` and `∂χ_Q/∂L7 = 0`; §IV asserts `∂N_Q/∂L7 = 0` and `∂Δ_norm/∂L7 = 0`; §V asserts the L7 coefficient first enters the normalized response at `z^7` with value `−i/L0`. The Mathematica script asserts the identical deliverables, but §III re-derives the canonical-even slots from the native `SphericalHankelH1[2,x]` outgoing operator and `Solve`s for `Σ_2,Σ_4` instead of hardcoding them.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| (1) response-side diff identity, change first at ω^7 | py:59,72,80 / wl:90,91,98 | match |
| (2) DtN-side diff identity, change first at z^7 | py:103,104,109 / wl:123,124,125 | match |
| (3) χ_Q unchanged by O(z^7) data | py:158,159 / wl:181,182,184 | match |
| (4) Δ_norm, N_Q unchanged | py:175,176,177 / wl:200,201,203 | match |
| (5) Δ_Q := χ_Q − 1 as the sole obstruction | ledger prints py:207-208 / wl:212 (stated, not separately asserted; follows from (3)) | match (the operative content is χ_Q stability, fully asserted) |

`paper_alignment: aligned`. The card text "Mathematica audit: none yet" (stage_196.tex:11) is contradicted by a present, passing `.wl` — a paper-side card-text lag (paper-cleanup class), noted in F2.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 59 | `simplify(Y7−Y5) − Ydiff_expected == 0` | claim 1 | yes |
| A2 | sympy | 72 | `series(Y7−Y5, ω, 0, 7) == 0` | claim 1 (unchanged through ω^5) | yes |
| A3 | sympy | 78-81 | diff series `− i·τ_Q·ω^7/4 == 0` | claim 1 (starts at ω^7) | yes |
| A4 | sympy | 103 | DtN diff `− Ydef_diff_expected == 0` | claim 2 | yes |
| A5 | sympy | 104 | `series(Ydef diff, z, 0, 7) == 0` | claim 2 | yes |
| A6 | sympy | 158 | `chi_from_series − chi_from_def == 0` | claim 3 | yes |
| A7 | sympy | 159 | `diff(chi_from_series, L7) == 0` | claim 3 (L7-independence) | yes (L7 present in matchedDtn z^7) |
| A8 | sympy | 175,176 | `diff(N_Q,L7)==0`, `diff(Δ_norm,L7)==0` | claim 4 | yes |
| A9 | sympy | 187 | `L7_coeff_in_Y + i/L0 == 0` | claim 1/2 (first entry at z^7) | yes |
| B1 | mathematica | 90 | `retDifference − retDifferenceTarget == 0` | claim 1 | yes |
| B2 | mathematica | 91 | `retDifferenceLowCoeffs == {0..0}` | claim 1 | yes |
| B3 | mathematica | 123 | `dtnDifference − dtnDifferenceTarget == 0` | claim 2 | yes |
| B4 | mathematica | 146-154 | `Solve` even rules unique; B5/B6 below verify | claim 3 (derives, not asserts, Σ_2,Σ_4) | yes |
| B5 | mathematica | 177,178 | `(σ2/.rule) − σ2Target == 0`, σ4 likewise | claim 3 | yes |
| B6 | mathematica | 181 | `chiFromSeries − chiSympyTarget == 0` | claim 3 | yes |
| B7 | mathematica | 182,184-186 | `D[chiFromSeries,l7]==0`; `D[Coefficient[matchedDtn,z,k],l7]==0` for k=0..5 | claim 3 (L7-independence) | yes (l7 present in matchedDtn z^7) |
| B8 | mathematica | 200,201 | `D[nNatural,l7]==0`, `D[deltaNatural,l7]==0` | claim 4 | yes |
| B9 | mathematica | 128 | `Coefficient[...,l7] + I/l0 == 0` | claim 1/2 | yes |

No tautological or unanchored rows. Every `assert_zero`/`assert_nonzero`-shaped derivative check has its differentiation variable (`L7`/`l7`) genuinely present in the parent expression at `z^7` (verified below), so none is trivially-true.

## Findings

### F1 — stale_output

**Severity:** low (informational)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.py` (mtime 2026-06-03 15:59:11)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.txt` (mtime 2026-06-01 12:02:48)

**What's wrong:**
The saved SymPy output is older than the `.py` and its content disagrees with the current script: the committed `.txt` banner reads `STAGE 179 — EXACT HIGHER-ODD IRRELEVANCE THEOREM` (lines 3, 165) while the current `.py` banner is `STAGE 196 …` (py:35, py:189). This is the known +17 renumber drift (179+17=196): the output was captured before the banner relabel. The Mathematica output is current (banner "STAGE 196", math output line 3). Only the banner string differs; every numeric/symbolic result in the stale `.txt` still matches what the current script computes (the algebra was untouched by the relabel — confirmed by comparing the §I–§V result lines against the current `.py` deliverables).

**Why this matters:**
Cosmetic/provenance only — the captured transcript carries a stale stage label. The verifier re-runs the script independently, which will regenerate the correct "STAGE 196" banner and refresh the output. Non-blocking.

**Required change:**
None for Codex. Orchestrator's independent re-run of `python3 scripts/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_sympy_audit.py` refreshes the committed output and clears the stale 179 banner.

**Verification:**
After re-run, `scripts/output/…_sympy_audit.txt` line 3 reads `STAGE 196 — EXACT HIGHER-ODD IRRELEVANCE THEOREM` and its mtime is newer than the `.py`.

### F2 — paper_misalignment (subtype: paper_missing_script_claim — card-text lag)

**Severity:** low
**Files:**
- paper side: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_196.tex:11` quote: *"Mathematica audit: none yet."*
- script side: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_mathematica_audit.wl` (present, all checks PASS — math output lines 13-69)

**What's wrong:**
The card's `\stagefield{Verification}` line still says the Mathematica audit does not exist, but a complete, passing dual-engine `.wl` was created in batch V.3 (commit 27e393b). The card text lags the actual verification state.

**Why this matters:**
The published card understates verification coverage; a reader citing this card would believe only SymPy verifies the stage. Purely paper-prose; no math impact.

**Required change:**
None for Codex (paper.tex is off-limits to the red-team). Routes to the user / paper-cleanup tracker: update stage_196.tex:11 to cite `mathematica/moving_throat_pde_stage196_higher_odd_irrelevance_theorem_mathematica_audit.wl` once the V.3 retrofit batch is signed off. This is a paper-cleanup-class item, not a script defect.

**Verification:**
stage_196.tex:11 names the `.wl` path; no script change.

## Independent-derivation check (Mathematica)

The `.wl` is **GENUINELY INDEPENDENT on the load-bearing claim (χ_Q stability, §III)** and borderline-but-acceptable on the bookkeeping identities (§I/§II). Evidence:

- **§III (the heart — genuinely independent).** The `.py` *hardcodes* the DtN parameterization and the matching laws:
  `L0_stage194 = -3*S + Sigma0` … `L5_stage194 = S*beta**5/9 + Sigma5` (py:123-126), `Sigma2_match = -(3*S*beta**2 - 3*S + Sigma0)/9` (py:133), then *verifies* them. The `.wl` instead *regenerates* the canonical-even window from the native outgoing operator:
  `lambdaOut = FunctionExpand[x*D[SphericalHankelH1[2, x], x]/SphericalHankelH1[2, x]]` (wl:134) → `lambdaWindow5` (wl:135), maps `x → stretch*z`, and **`Solve`s** for the slots:
  `evenRules = Solve[{Coefficient[normalDtn,z,2]==1/9, Coefficient[normalDtn,z,4]==4/81}, {sig2,sig4}, Reals]` (wl:146-153), checking uniqueness (wl:154). The `.py` never imports a Hankel function (`grep` for `SphericalHankel`/`FunctionExpand` in the `.py` returns nothing). Both routes land on the same `χ_Q = 3(scale·stretch^5+9·sig5)/(3·scale−sig0)` (py output line 130; wl output line 41). This is the strongest possible independence on the one claim that actually matters here.

- **§I/§II (borderline — forced-form identities).** Both engines build the same module `3/4 + 1/(4(1−frontCore−tail))` (py:49-50 / wl:70-73) and the same DtN module `l0/(l0+l2 z²+l4 z⁴+i l5 z⁵+tail)` (py:91-96 / wl:104-108), and both compare the difference to the *same* closed form `H/(4(1−X)(1−X−H))` (py:51 `Ydiff_expected` / wl:75 `retDifferenceTarget`) and `−L0·L_{≥7}/(D5·D7)` (py:97 / wl:110). The `.wl` even labels these "- SymPy target". This is the same algebraic construction in both engines. It is **not** flagged as `mathematica_transliteration` because the difference of two rational functions is *algebraically forced* to that unique closed form — there is no alternate "route" to a difference-of-fractions identity; both engines necessarily produce it. The coefficient-extraction mechanism does differ (`sp.series().removeO()` py:61-63 vs `Series`/`Coefficient` tables wl:76-78,113). Per the consistent-threshold rule, a forced-unique identity reached by both engines is not a port of "echoed algebra"; the load-bearing independence (§III) carries the dual-engine requirement.

Conclusion: independence is real where it counts. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree on every deliverable:
- `σ_Q^can = 9/(8 Ω^5)` (py out L9-13 / wl out L9).
- response through ω^5: `1 + ω²/4Ω² + ω⁴/4Ω⁴ + 9 i χ ω⁵/32Ω⁵` (py out L33-36 / wl out L10).
- difference starts at `i τ ω⁷/4` (py out L44-47 / wl out L18... wl out L18 line "first response-side higher-odd correction" PASS).
- `χ_Q = 3(Sβ⁵+9Σ₅)/(3S−Σ₀)` (py out L121-130 / wl out L41).
- `N_Q = (3S−Σ₀)/(Sβ⁵+9Σ₅)`·(scaled) and `Δ_norm` (py out L138-150 / wl out L60-61).
- `∂χ_Q/∂L7 = ∂N_Q/∂L7 = ∂Δ_norm/∂L7 = 0` (py out L133,151,152 / wl out, all PASS).
- L7 enters the normalized response first at z⁷ with coefficient `−i/L0` (py out L159-162 / wl out L31-32).
No sign, factor, or residual disagreement. `engines_agree: true`.

## Verdict justification

`findings` — but only two low-severity, non-math items: a stale SymPy output (`.txt` carries the pre-renumber "STAGE 179" banner; refreshed by the verifier's independent re-run) and a paper-side card-text lag ("Mathematica audit: none yet" despite a present passing `.wl`). The mathematics holds up under attack. Attacks tried that failed: (a) the `l7`/`L7` independence self-test trap — I confirmed `l7` genuinely appears in `matchedDtn` at the `z^7` slot (sympy output L96: `729·i·L7·z^7/(2187S−729Σ₀)`; `.wl` builds it explicitly at wl:143), so `D[chiFromSeries,l7]=0` and the z^0..z^5 independence checks are *non-vacuous* can-fail statements, not trivially-true; (b) transliteration on the load-bearing claim — the `.wl` §III re-derives the canonical-even slots from the native `SphericalHankelH1[2,x]` operator and `Solve`s for `Σ_2,Σ_4`, a route the `.py` does not contain, so it is genuinely independent; (c) value reconciliation — every emitted deliverable (`χ_Q`, `Σ_2,Σ_4`, `N_Q`, `Δ_norm`, `P0_target`, `Ω_Q`, `σ_Q^can`, `Δ_Q`) reconciles against the notes. I read the card, notes (full), and the appendix rows; the script's verified claim matches the paper's stated claim.

## Self-test notes

(1) Variable independence: `L7`/`l7` is genuinely present in `Ydef7`/`matchedDtn`/`dtnWithTail` at the `z^7` slot (py:129, wl:143; sympy out L96 shows the `729·i·L7·z^7` term), so every `diff(·,L7)`/`D[·,l7]` check is non-vacuous — the derivative is a real `0` because `χ_Q,N_Q,Δ_norm` depend only on `z^0..z^5` data, not because `L7` is absent. The `assert_nonzero`-flavored counterpart is implicitly covered: `L7_coeff_in_Y + i/L0 == 0` (py:187, wl:128) confirms `L7` enters with nonzero coefficient `−i/L0` at `z^7`. (2) Parity/symmetry: no unbounded-domain integrals here; all checks are power-series coefficient identities in `ω`/`z` about 0 — even/odd power structure (`ω^5` odd reaction slot, `ω^7` next odd) is respected and matches the notes. (3) Trivial-case: substituting `τ_Q=0`/`L7=0` collapses the tailful module to the tailless one and every difference identity reduces to `0` consistently. No directive is written: F1 is an orchestrator re-run (no Codex edit) and F2 is a paper-side prose lag (Codex cannot touch paper.tex); there is no script-side fix for Codex to apply.

## Value Reconciliation (pass-2 augmentation)

Deliverable-level reconciliation. The committed SymPy `.txt` is stale only in its banner label (STAGE 179); all numeric/symbolic result lines still reflect the current script, so I base reconciliation on script source + saved outputs as the augmentation permits, noting the `stale_output` signal (F1).

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `σ_Q^can = 9/(8 Ω_Q^5)` | py:45 / wl:68 (py out L9, wl out L9) | notes:90 `σ_Q^can = 9/(8Ω_Q^5)` | MATCH |
| response through ω^5 `1+ω²/4Ω²+ω⁴/4Ω⁴+9iχω⁵/32Ω⁵` | py:73 / wl:80 (py out L33-36, wl out L10) | notes:141-153 | MATCH |
| response-side diff `= H_Q/(4(1−X)(1−X−H))` | py:51,59 / wl:75,90 | notes:120-125 | MATCH |
| DtN-side diff `= −L0 L_{≥7}/(D5 D_{≥7})` | py:97,103 / wl:110,123 | notes:206-211 | MATCH |
| DtN response through z^5 `1−L2z²/L0+(L2²/L0²−L4/L0)z⁴−iL5z⁵/L0` | py:112 / wl:114-119 (py out L67-78) | notes:221-233 | MATCH |
| `Σ_2 = −(3Sβ²−3S+Σ_0)/9` | py:133 / wl:157,177 (py out L84-88, wl out L43) | notes:263 | MATCH |
| `Σ_4 = −(3Sβ⁴−3S+Σ_0)/27` | py:134 / wl:158,178 (py out L89-93, wl out L45) | notes:265 | MATCH |
| `χ_Q = 3(Sβ⁵+9Σ_5)/(3S−Σ_0)` | py:141 / wl:163,181 (py out L121-130, wl out L41) | notes:50-51, 274 | MATCH |
| `N_Q = 1/χ_Q = (3S−Σ_0)/(3(Sβ⁵+9Σ_5))` | py:168 / wl:192,194 (py out L138-144, wl out L60) | notes:300 (`N_Q=χ_Q^{-1}`) | MATCH |
| `P0_target = 54 G c_s^5/(5 a^5 c^5)` | py:167 / wl:191 (in py out L147, wl out L61) | notes:297 | MATCH |
| `Δ_norm = P0_target(1/χ_Q − 1)` | py:169,177 / wl:193,203 (py out L145-150, wl out L61) | notes:290-291 | MATCH |
| L7 first-entry coeff `−i/L0` at z^7 | py:184,187 / wl:128 (py out L158-162) | notes:67-78 (tail first at z^7), :234 | MATCH |
| `Δ_Q := χ_Q − 1` (sole obstruction) | ledger py:208 / wl:212 | notes:365 (boxed) / card:13,15 | MATCH |

Internal scaffolding (accounted for, no finding): `X_Q`, `H_Q`, `frontCore`, `Ydiff_expected`/`retDifferenceTarget`, `Ydef_diff_expected`/`dtnDifferenceTarget`, `D5`/`D7`, `evenRules`/`evenRule`, intermediate `Y5`/`Y7`/`Ydef5`/`Ydef7`, series-truncation residuals, all PASS/FAIL flags.

reconciliation: complete; 13 deliverable values checked, 0 misaligned.
