---
unit_id: 134
batch: IV.4
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-06T00:00:00Z
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
  notes_stage_files: [/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md]
  paper_appendix: present
---

# Audit unit 134 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_134.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage134_family1_mouth_fixedpoint.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (subsection "Coupled mouth-layer fixed point", lines 674-725; `\Pi_*` at line 663; `\input{stages/stage_134}` at line 1302)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.txt`

## What the paper claims

Stage 134 specializes the general Stage 133 coupled mouth-layer fixed-point law to the first explicit Family-1 branch: one static shell-compliance lane (`\kappa_s=0`, contributing the exact static limit `\mathcal S(\Pi,0)=1`) and one first mixed D/N half-wave lane (`\kappa_q=\pi/2`). The card's quoted `\stagefield{Output}`-equivalent (the boxed result block in the card body) is: "First explicit reduction gives `\(\Pi=M_s+M_q\mathcal S_q(\Pi)\)`." The notes add the closed-form kernel `\mathcal S_q(\Pi)=\Pi[\tfrac{\pi}{2}\tanh\tfrac{\pi}{2}+\Pi(e^{-\Pi}\operatorname{sech}\tfrac{\pi}{2}-1)]/[(1-e^{-\Pi})(\tfrac{\pi^2}{4}-\Pi^2)]`, the canonical compensation value `\mathcal S_q(\Pi_*)\approx0.658075937605428` at the imported `\Pi_*\approx1.50882951349316`, and the canonical gain line `M_s\approx1.50882951349316-0.658075937605428\,M_q`. Distinct deliverables: (D1) the fixed-point law `\Pi=M_s+M_qS_q(\Pi)`; (D2) the static shell limit `\mathcal S(\Pi,0)=1`; (D3) the closed-form `S_q(\Pi)`; (D4) `S_q(\Pi_*)\approx0.658075937605428`; (D5) the canonical gain line. The card's `\stagefield{Checks}` explicitly defers outlet consistency of `(M_s,M_q)` to Stage 135 and susceptibility closure to Stage 137 ("carried forward here"), and states the fixed points are "numerically located, not closed-form constants."

## What the script claims to verify

The SymPy script builds the D/N kernel `S(Pi,kappa)` from the parent form, takes the `kappa->0` limit (static shell channel), substitutes `kappa=pi/2` to form `S_q`, and forms the fixed-point RHS `M_s + M_q*S_q`. It then makes three substantive assertions: (1) `S_shell - 1 == 0` exactly (static shell limit, D2); (2) `S_q(p)` at `p=1/2,1,2` matches three external 30-digit literals (claimed by comment to be from an independent mpmath run) within 1e-12, exercising the closed form D3; (3) `S_q(Pi_star)` matches the notes literal `0.658075937605428` within 1e-12 (D4). The gain line (D5) is printed symbolically but intentionally NOT asserted, with a comment explaining that re-asserting `intercept == Pi_*` / `slope == -S_q(Pi_*)` would be an X-X tautology and that outlet consistency is verified downstream at Stage 135. The Mathematica script mirrors this exact choreography step-for-step.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| D1 — fixed-point law `Π = M_s + M_q S_q(Π)` | Both engines form `M_s + M_q*S_q` and PRINT it (`.py:31`, `.wl:41`); structure of `S_q` is exercised by D3/D4 checks | match (the law's shape is the carrier; the load-bearing content is `S_q`, checked) |
| D2 — static shell limit `S(Π,0)=1` | `.py:48-50` `assert simplify(S_shell-1)==0`; `.wl:44` `expectZero["static shell channel", sShell-1]` | match |
| D3 — closed-form `S_q(Π)` | `.py:64-72` three external-literal spot checks; `.wl:57-62` same | match |
| D4 — `S_q(Π_*) ≈ 0.658075937605428` | `.py:76-78` assert vs notes literal; `.wl` computes `sStar` but only PRINTS (no assert) | match (sympy) / partial (mathematica prints, does not assert) |
| D5 — canonical gain line | printed only, not asserted (both engines), by design (X-X tautology avoidance; deferred to St.135) | match (intentional non-assertion, paper Checks downgrade aligns) |

`paper_alignment: aligned`. Every paper deliverable maps to a script-side check or to a documented, paper-card-aligned deferral. The Checks-downgrade (Cluster C of `redteam/resolutions/batch_IV4_paper_alignment.md`) is reflected faithfully in both the card and the script comments.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 50 | `assert simplify(S_shell - 1) == 0` | D2 (static shell limit) | yes (real `kappa->0` limit; gives 1 only if kernel correct) |
| A2 | sympy | 64-71 | `assert abs(S_q(p)-target) < 1e-12` at p=1/2,1,2 | D3 (closed-form S_q) | yes (external literals; wrong kernel fails) |
| A3 | sympy | 77-78 | `assert abs(S_q(Pi_star)-0.658075937605428) < 1e-12` | D4 | yes (notes literal anchor) |
| A4 | mathematica | 44 | `expectZero["static shell channel", sShell-1]` | D2 | yes |
| A5 | mathematica | 57-62 | `expectClose["S_q at p=...", N[sQ/.p->...], literal, 1e-12]` at p=1/2,1,2 | D3 | yes |
| A6 | mathematica | (none) | `sStar` PRINTED only at :70, not asserted | D4 | no (print-only; SymPy carries the D4 assert) |

A1-A5 are "yes". A6 is a print, not an assertion — but D4 is asserted on the SymPy side (A3), so the deliverable is covered. The gain line (D5) is print-only on both sides by design (documented X-X avoidance), which the paper Checks downgrade endorses.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:31-62`
- (corresponding SymPy) `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:21-72`

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script rather than an independent second-engine derivation. The kernel functions are character-for-character identical algebra:

- SymPy `.py:22-24`: `Pi * (kappa * sp.tanh(kappa) + Pi * (sp.exp(-Pi) / sp.cosh(kappa) - 1)) / ((1 - sp.exp(-Pi)) * (kappa**2 - Pi**2))`
- Mathematica `.wl:32-33`: `p*(k*Tanh[k] + p*(Exp[-p]/Cosh[k] - 1))/((1 - Exp[-p])*(k^2 - p^2))`

Every subsequent step is the same choreography in the same order: shell limit via `Limit[..., k0->0]` (`.wl:36-39`) vs `sp.limit(S(Pi,kk), kk, 0)` (`.py:29`); `S_q` via substituting `kappa=Pi/2` (`.wl:40` vs `.py:30`); `fixedPointLaw = Ms + Mq*sQ` (`.wl:41` vs `.py:31`); and the three numeric spot-checks at `p=1/2,1,2` using the **identical** 30-digit literal targets `0.608336415687717065435990381419`, `0.633127670034487546375729566676`, `0.681366857005321783286541952613` (`.wl:57-62` vs `.py:61-63`). The two engines do not derive the result by structurally distinct routes; the `.wl` simply re-types the `.py` algebra in Mathematica syntax.

This transliteration was identified during pass-1 (batch IV.4): `redteam/resolutions/batch_IV4_paper_alignment.md:91` (group M5) prescribed for stage 134 "use `Series` + `Limit` in Mathematica to extract `S_q` as the `kappa -> pi/2` boundary, rather than retyping the hand-derived expression." That fix was never applied. The batch-5 commit `0a9b203` only removed the X-X gain-line print/assert; the kernel transliteration was retained. The project's own tracker confirms this: `notes/MATHEMATICA_MIRROR_POLICY.md:32` records "**134** introduced NO new mirror — ... the corroborated shell-limit + S_q spot-checks were kept" — i.e., the pre-existing transliterated kernel was carried forward, not de-mirrored. (Contrast IV.4's de-mirrored stages 127/133/137, `MIRROR_POLICY.md:48`.)

**Why this matters:**
The dual-engine policy exists so a transcription or sign error in the hand-derived kernel would surface as a cross-engine disagreement. Because both engines type the same expression and check it against the same external literals, an error introduced into the kernel form would be replicated identically in both — the second engine provides no independent confirmation of the kernel's algebraic structure (only of the numeric literals, which both already share). The only genuine cross-engine corroboration here is numeric agreement at three sample points, which is weaker than an independent symbolic re-derivation.

**Required change:**
Replace the Mathematica `sQ` construction (`.wl:40`) with an independent symbolic route that does not re-type the SymPy closed form. Per the pass-1 M5 plan: derive `S_q` as the `kappa -> pi/2` boundary of the D/N response via `Series`/`Limit` in Mathematica (the `kappa = pi/2` point is a removable singularity of `sKernel` because `(kappa^2 - p^2)` is finite there, so `S_q = Limit[sKernel[p, kap], kap -> Pi/2]` extracted symbolically is a structurally distinct evaluation route from substituting `Pi/2` into the SymPy-form expression). Then `expectZero["S_q independent-route", sQindependent - sQ]` to cross-check the independent route against the substituted form. Keep the existing shell-limit and spot-check assertions. See directive for the exact construction and the self-test.

**Verification:**
After the fix, `.wl` should contain a new `Limit`/`Series`-based `sQ` derivation distinct from `sKernel[p, Pi/2]`, plus an `expectZero` (or `pass`) line cross-checking it; the Mathematica output should show a new PASS line for that cross-check, and the existing three spot-checks and shell limit must still PASS, exit 0.

## Independent-derivation check (Mathematica)

The `.wl` is NOT an independent derivation; it is a transliteration of the `.py` (see F1, with the three corresponding sections quoted). The kernel, the limit-extraction strategy, the `Pi/2` substitution, the `Ms + Mq*sQ` assembly, and the three numeric literal targets are all reproduced step-for-step. The only Mathematica-distinct elements are cosmetic (the `expectZero`/`expectClose` harness and the `PiStar` comment-naming to dodge pitfall #13). This is the basis for the F1 `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree at the level they assert:

| quantity | SymPy output | Mathematica output |
|---|---|---|
| `S_shell - 1` | `0` (line 27) | `0`, PASS (line 6-7) |
| `S_q(1/2)` | `0.608336415687717065435990381419` (line 29) | `0.6083364156877170654359903814193...` PASS (line 9-10) |
| `S_q(1)` | `0.633127670034487546375729566676` (line 30) | `0.6331276700344875463757295666760...` PASS (line 11-12) |
| `S_q(2)` | `0.681366857005321783286541952613` (line 31) | `0.6813668570053217832865419526134...` PASS (line 13-14) |
| `S_q(Pi_star)` | `0.658075937605428494269581645208` (line 24) | `0.65807593760542948674050367268...` (line 16) |

`S_q(Pi_star)` agrees to ~15 digits (`...605428494...` vs `...605429486...`), as expected since `Pi_star` is supplied as a 15-digit literal (`1.50882951349316`), bounding the propagated precision. No engine disagreement at the asserted tolerance (1e-12). The shell limit and three spot-checks are bit-for-bit consistent.

## Verdict justification

The math holds up: the static-shell limit `S(Π,0)=1` is a real (non-tautological) `kappa->0` computation and is correct; the `S_q` closed form matches both the notes (lines 51-60) and the appendix kernel `eq:app-part04-S-kernel` exactly; the three external-literal spot checks and the `S_q(Π_*)` anchor are non-tautological (a wrong kernel would miss the literals) and pass; the gain line is correctly print-only by design and the paper Checks downgrade endorses the downstream deferral to Stages 135/137; the value reconciliation is fully clean. I attacked the kernel sign/factor (the `(κ²-Π²)` denominator and the `e^{-Π}sech κ` numerator term — both match notes and appendix verbatim), the static limit (re-derived by hand to 1, correct), and the spot-check literals (hand-estimated `S_q(1)≈0.6333`, matching the literal). The single defect is structural: the Mathematica script is a line-by-line transliteration of the SymPy script (F1), a fix that pass-1 planned (M5) but never applied. This is not stop-cold — it does not affect the math result, no downstream constant changes, and it is mechanically fixable by adding an independent symbolic route in the `.wl`. Hence `verdict: findings`, `stop_cold: null`, `paper_alignment: aligned`.

## Value Reconciliation (pass-2 augmentation)

Every RESULT/deliverable value the scripts emit, reconciled against the `.tex` card and `.md` notes:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `S_shell = 1` (static shell limit) | py:48-51 / wl:44; out py:5,27, wl:5-7 | notes:25-27 (boxed `S(Π,0)=1`); appendix `S(Π,κ)` kernel l.690-694 | MATCH |
| `S_q(Π)` closed form | py:30 / wl:40; out py:6-14, wl:8 | notes:50-61 (boxed); appendix `eq:app-part04-S-kernel` l.688-695 | MATCH |
| Fixed-point law `Π = M_s + M_q S_q(Π)` | py:31 / wl:41; out py:15-23, wl:15 | notes:41-47 & 132-136 (boxed); stage_134.tex:16; appendix `eq:app-part04-F1-mouth-fixedpoint` l.702-707 | MATCH |
| `S_q(Π_*) ≈ 0.658075937605428` | py:40,76-79 / wl:67,70; out py:24,33, wl:16 | notes:86 (`\mathcal S_q(\Pi_*)\approx0.658075937605428`) | MATCH |
| `Π_* = 1.50882951349316` (imported) | py:39 / wl:66; out py:26, wl:18 | notes:72; appendix:663; (owned by stage 130/131) | MATCH |
| Canonical gain line `M_s ≈ 1.50882951349316 - 0.658075937605428 M_q` | py:43 / wl:68; out py:26, wl:18 | notes:91-94 & 138-142 (boxed); stage Result | MATCH |

INTERNAL (verification scaffolding, no finding, correctly absent from prose): the three spot-check anchor literals `0.608336415687717065435990381419` (S_q at 1/2), `0.633127670034487546375729566676` (S_q at 1), `0.681366857005321783286541952613` (S_q at 2); tolerance `1e-12`; the printed symbolic `S_shell`/`S_q`/`fixedPointLaw` expression dumps; PASS/FAIL flags and residual diffs.

reconciliation: complete; 6 deliverable values checked, 0 misaligned

## Self-test notes

I checked: (1) variable independence — there are no `sp.diff`/`D[]` derivative-vs-wrong-variable traps in this script; the only "derivative-like" step is the `kappa->0` limit, which genuinely depends on `kappa` and is non-vacuous (yields 1). (2) Symmetry/parity — no unbounded integrals here. (3) Trivial-case pre-check — I hand-verified the static-shell `kappa->0` limit reduces to exactly 1 from the kernel, and hand-estimated `S_q(1)≈0.6333` matching the literal target, confirming the asserts are non-trivially satisfiable and would fail on a wrong kernel. (4)/(5) For the F1 directive, the proposed `Limit[sKernel[p, kap], kap -> Pi/2]` route targets `mathematica/` (correct dir), and is a removable-singularity limit (denominator `(kap^2 - p^2)` is finite and nonzero at `kap=pi/2` for generic `p`, with the `p=pi/2` pole excluded by `p>0` generic), so the limit equals direct substitution — the `expectZero[sQindependent - sQ]` cross-check is a genuine zero (not vacuous) and introduces no new paper_misalignment (no new constant). I confirmed I read the paper card, the notes, and the relevant appendix subsection before opening the scripts, and the script's verified claim matches the paper's claim.
