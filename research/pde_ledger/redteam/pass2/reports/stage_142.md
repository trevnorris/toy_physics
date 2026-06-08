---
unit_id: 142
batch: IV.5
auditor_model: claude-opus-4-8[1m]
audit_date: 2026-06-07T00:00:00Z
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
  notes_stage_files: [notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md]
  paper_appendix: present
---

# Audit unit 142 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_142.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage142_selfconsistent_mouth_branch.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_142}` at line 1318; no extra narrative row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.txt`

## What the paper claims

Stage 142 closes the mouth problem by identifying the free core branch label with the explicit mouth-source factor, `g_c = g_Π`, which collapses the gain pair `(M_s, M_q)` to a single scalar branch law. The card's bottom-line claim (quoted verbatim): *"Identifying \(\mathfrak g_c=\mathfrak g_\Pi\) gives \(\Pi=\Sigma_0[1-R_q(\Pi)\mathcal S_q(\Pi)]\)."* (stage_142.tex:16). The notes add the explicit deliverables: the gain ratio `R_q(Π) = (g_Π − r_F1)²/(1+r_F1²)`, the shell-gain function `Σ_0(Π) = Π/(1 − R_q S_q)`, the self-matched traction `T̂_m(Π) = sqrt(9Π/(20(1−R_q S_q)))`, and the canonical Family-1 point: `g_-^F1 ≈ 0.758035078944663`, `Π_* ≈ 1.50882951349316`, `R_q(Π_*) = 1/4`, `S_q(Π_*) ≈ 0.658075937605429`, `Σ_0(Π_*) ≈ 1.80594111095636`, `T̂_m(Π_*) ≈ 0.901484054174205` (md:96–121). The card's Checks item 3 requires fixed points be "numerically located, not closed-form constants."

## What the script claims to verify

Both engines build the same chain — `r_F1 = sqrt(4107−100π²)/(10π)`, the hardcoded mouth-source factor `g_Π`, the carried-forward `S_q(Π)` tanh/sech closed form (Stage 140 closure), then `R_q`, `Σ_0`, `T̂_m`. The load-bearing checks are: (a) `R_q(g_-) − 1/4 = 0` (a definitional/branch-blind solver-consistency identity, explicitly labeled as such in the comments); (b) a NON-tautological cross-stage anchor evaluating `R_q` at Stage 131's independently-derived `Π_ext = 1.50882951349315558300555075595` and confirming it lands on `1/4`; (c) `expect_close` on all six canonical-point deliverable values. The Mathematica script additionally (d) re-derives `g_Π` by direct symbolic integration of the Stage-130 mouth-source law `∫₀¹ σ_Π(z) cos(πz/2) dz` and confirms it equals the hardcoded `g_Π`, and (e) checks `100π²(1+r²) = 4107`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| Branch law `Π = Σ_0[1 − R_q S_q]` (card:16) | Ledger print + `Σ_0 = Π/(1−R_q S_q)` construction (py:41, wl:55) | match |
| `R_q(Π) = (g_Π−r)²/(1+r²)` (md:36–40) | `Rq` def + `R_q(g_-)=1/4` + ext anchor (py:40,54,91; wl:54,94,118) | match |
| `Σ_0(Π) = Π/(1−R_q S_q)` (md:59–63) | `Sigma0` def + canonical-point value (py:41,100; wl:55,123) | match |
| `T̂_m(Π) = sqrt(9Π/(20(1−R_q S_q)))` (md:78–82) | `That = sqrt((9/20)Σ_0)` + value (py:42,101; wl:56,124) | match |
| `g_-^F1 ≈ 0.758035078944663` (md:97) | `expect_close g_-^{F1}` (py:97; wl:120) | match |
| `Π_* ≈ 1.50882951349316` (md:106) | `expect_close Pi_* value` (py:98; wl:121) | match |
| `R_q(Π_*) = 1/4` (md:112) | solver-consistency + independent anchor (py:82,94; wl:114,119) | match |
| `S_q(Π_*) ≈ 0.658075937605429` (md:114) | `expect_close S_q(Pi_*)` (py:99; wl:122) | match |
| `Σ_0(Π_*) ≈ 1.80594111095636` (md:119) | `expect_close Sigma_0(Pi_*)` (py:100; wl:123) | match |
| `T̂_m(Π_*) ≈ 0.901484054174205` (md:121) | `expect_close That(Pi_*)` (py:101; wl:124) | match |
| Checks item 3: "numerically located" | nsolve (py:57) / FindRoot (wl:97) | match |
| `g_Π` source-projection origin (Stage 130) | wl re-derives via Integrate (wl:70–76) | match (extra, in-scope) |

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 54 | `expect_zero(R_q(g_-) − 1/4)` | `R_q(Π_*)=1/4` (definitional) | partial (branch-blind; labeled) |
| A2 | sympy | 82 | `R_q(Π_*) nsolve − 1/4 < 1e-15` | `R_q(Π_*)=1/4` (solver-consistency) | partial (redundant; labeled) |
| A3 | sympy | 94 | `R_q(Π_ext) − 1/4 < 1e-12` at Stage-131 Π_* | `R_q(Π_*)=1/4` (non-tautological) | yes |
| A4 | sympy | 97–101 | `expect_close` ×5 on canonical values | all 5 canonical deliverables | yes |
| A5 | mathematica | 76 | `expectZero(g_Π,integral − g_Π,closed)` | `g_Π` source-projection (Stage 130) | yes (independent re-derivation) |
| A6 | mathematica | 82 | `expectZero(100π²(1+r²) − 4107)` | `r_F1` constant consistency | no (tautological for r-def) |
| A7 | mathematica | 94 | `expectZero(R_q(g_-) − 1/4)` | `R_q(Π_*)=1/4` (definitional) | partial (mirrors A1) |
| A8 | mathematica | 119 | `R_q(Π_ext) = 1/4` at Stage-131 Π_* | `R_q(Π_*)=1/4` (non-tautological) | yes |
| A9 | mathematica | 120–124 | `expectApprox` ×5 on canonical values | all 5 canonical deliverables | yes |

A6 is tautological-by-construction (`1+r² = 4107/(100π²)` for the given `r`), but it is a labeled sanity check on the carried-forward constant, not a stand-in for a physics claim; the load-bearing independent work is A3/A5/A8. No standalone finding.

## Findings

None. All attacks below failed.

## Independent-derivation check (Mathematica)

PARTIAL → the call is **INDEPENDENT for the load-bearing piece**. The `S_q`, `R_q`, `Σ_0`, `T̂_m` closed forms in the `.wl` (lines 53–56) are the same hardcoded expressions as the `.py` (lines 39–42), so that algebra is shared. BUT the most load-bearing object, `g_Π`, is genuinely re-derived independently in the `.wl`:

- `.wl:69–76`: `sigmaPi = piM*Exp[-piM*zVar]/(1 - Exp[-piM])` then `Integrate[sigmaPi*Cos[Pi*zVar/2], {zVar,0,1}]` and `expectZero[g_Π,integral − g_Π,closed]`. This is built only from the Stage-130 source profile and the cos projection — it does NOT reuse the hardcoded `gPi` algebra, so it cannot share a typo with it. Output (wl.txt:13–15) shows the integral evaluated to exactly the hardcoded closed form (residual 0). The `.py` has no analogous derivation.
- This matches the heads-up: the projection integral `∫₀¹ σ_Π(z) cos(πz/2) dz` has replaced an earlier `Normal[Series[...]]` self-check; `grep` confirms zero `Series`/`Normal` remain.
- The independent anchor uses Stage 131's `Π_*` (1.50882951349315558300555075595), NOT 142's own nsolve output (1.5088295134931555274…, which diverges at digit ~16) — confirmed at `.wl:117` and `.py:90`, both equal to Stage-131's `.txt` value (output line 2).
- The `g_Π` denominator sign is `(e^Π − 1)` (correct) in both `g_Π` (py:35, wl:49) and Stage 131's nsolve target (stage131.py:23), not `(1 − e^Π)`. The `(1 − e^Π)` that appears squared inside `R_q` (py.txt:12) is sign-irrelevant because it is squared.

Because one engine independently re-derives the central hardcoded form via a structurally distinct route (symbolic integration), this is NOT a `mathematica_transliteration` finding; the shared `S_q`/`R_q` algebra is acceptable given both engines anchor that algebra to the independent Stage-131 `Π_*`.

## Engine cross-check

Both engines pass and agree to the precision they claim:

| value | sympy (.txt) | mathematica (.txt) | agreement |
|---|---|---|---|
| g_-^F1 | 0.758035078944662826… | 0.758035078944662826919680890414… | full (≥27 digits) |
| Π_* | 1.50882951349315552747… (nsolve) | 1.50882951349315558300… (FindRoot) | digits 1–16; FindRoot matches Stage-131 |
| R_q(Π_*) | 0.250000…0019 | 0.250000…0002 | both → 1/4 |
| S_q(Π_*) | 0.658075937605429271930… | 0.658075937605429274616… | digits 1–17 |
| Σ_0(Π_*) | 1.80594111095635380721… | 1.80594111095635387237… | digits 1–17 |
| T̂_m(Π_*) | 0.901484054174204022702… | 0.901484054174204038964… | digits 1–17 |

The ~16–17-digit divergence on Π_*-dependent values is exactly the known nsolve-vs-FindRoot solver gap (FindRoot at WorkingPrecision 80 is the more precise and reproduces Stage 131; nsolve diverges at digit ~16). Both `expect_close`/`expectApprox` tolerances (1e-12) comfortably absorb it, and both engines independently re-anchor `R_q=1/4` at Stage-131's `Π_*`. No `engine_disagreement`.

## Verdict justification

Clean. I attacked: (1) the hardcoded `g_Π` — the `.wl` independently re-derives it from the Stage-130 source integral (residual 0), and the `.py` independent anchor at Stage-131's `Π_*` would break on any `g_Π` typo; both held. (2) the `R_q(g_-)=1/4` check — it IS definitional/branch-blind, but the script explicitly labels it as solver-consistency and supplies the genuine non-tautological anchor A3/A8, so it is not a disguised tautology. (3) the carried `r_F1 = sqrt(4107−100π²)/(10π)` — verified `4107 = 3·37²` so `g_- = (2sqrt(4107−100π²) − 37sqrt(3))/(20π)` matches Stage 131's exact form; no stale `168π²`. (4) freshness — both outputs (18:02) are newer than both scripts (16:51), so no `stale_output`. (5) all six canonical deliverable values reconcile with the notes. I read the card, notes, and appendix; the script's verified claim matches the paper's stated claim.

## Self-test notes

- **Variable independence**: no `diff`/`D` derivative-vanishing trap in this stage (the one derivative-style object, the `g_Π` source integral, is over `zVar` and the integrand genuinely depends on `zVar`).
- **Symmetry/parity**: the `.wl` `g_Π` integral is over `[0,1]` (bounded, asymmetric), no parity-cancellation trap; integrand `σ_Π·cos(πz/2)` is positive on (0,1), integral legitimately nonzero.
- **Trivial-case / pole pre-check**: confirmed `S_q` denominator `((π/2)²−Π²)` is nonzero at `Π_*≈1.509` (≈0.19; pole at Π=π/2≈1.571 is above Π_*), so the canonical-point evaluations are well-defined. Confirmed `A6` is exactly-zero-by-construction (acceptable labeled sanity check, not load-bearing).

## Value Reconciliation (pass-2 augmentation)

reconciliation: complete; 10 values checked, 0 misaligned

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `Π = Σ_0[1 − R_q S_q]` (branch law) | py:105 / wl:128; py.txt:39, wl.txt:62 | stage_142.tex:16; md:43–48 | MATCH |
| `Σ_0(Π) = Π/(1 − R_q S_q)` | py:41,106 / wl:55,129 | md:59–63 | MATCH |
| `T̂_m(Π) = sqrt(9Π/(20(1−R_q S_q)))` | py:42,107 / wl:56,130 | md:78–82 | MATCH |
| `g_-^F1 = 0.758035078944663` | py.txt:20 / wl.txt:34 (0.7580350789446628…) | md:97 (≈0.758035078944663) | MATCH |
| `Π_* = 1.50882951349316` | py.txt:21 (1.50882951349315552747…) / wl.txt:35 (1.50882951349315558300…) | md:106 (≈1.50882951349316) | MATCH (both round to md value) |
| `R_q(Π_*) = 1/4` | py.txt:23 (0.2500…0019) / wl.txt:37 (0.2500…0002) | md:112 (=1/4) | MATCH |
| `S_q(Π_*) = 0.658075937605429` | py.txt:24 (0.658075937605429271…) / wl.txt:38 (…274616…) | md:114 (≈0.658075937605429) | MATCH |
| `Σ_0(Π_*) = 1.80594111095636` | py.txt:25 (1.80594111095635380…) / wl.txt:39 (…387237…) | md:119 (≈1.80594111095636) | MATCH (rounds to md value) |
| `T̂_m(Π_*) = 0.901484054174205` | py.txt:26 (0.901484054174204022…) / wl.txt:40 (…204038…) | md:121 (≈0.901484054174205) | MATCH (rounds to md value) |
| `g_Π = 2Π(2Πe^Π+π)/((4Π²+π²)(e^Π−1))` | py.txt:10 / wl.txt:13,23 | md:9–13; stage_130 md:34–39 | MATCH |

INTERNAL (scaffolding, no finding expected in prose): `r_F1 = sqrt(4107−100π²)/(10π)` (carried-forward upstream constant, derived/anchored in Stage 131); `Π_ext = 1.50882951349315558300555075595` (a cross-stage anchor borrowed from Stage 131, correctly identical to Stage 131's `.txt` output); per-check residuals/tolerances; the `100π²(1+r²)=4107` identity value; the `g_Π`-integral residual.

All ten emitted deliverable values reconcile with the notes (the `.tex` card legitimately carries only the symbolic branch law and routes the numeric values to the `.md`, consistent with Checks item 3 requiring numerically-located fixed points). No MISMATCH and no MISSING-DELIVERABLE.
