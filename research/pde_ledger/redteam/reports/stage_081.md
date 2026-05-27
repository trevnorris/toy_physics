---
unit_id: 081
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage081_family1_pi_thresholds.md
  paper_appendix: present
---

# Audit unit 081 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_081.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage081_family1_pi_thresholds.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 140, `\input` at line 280)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.txt`

## What the paper claims

Stage 081 expresses the Family-1 demand thresholds in the branch-product variable `Pi_tr`, including the blocking parameter `eps_blk`. The card's `\stagefield{Output}` is the pair of equations: (i) the support-ratio inversion `zeta_req = (Pi_tr - C_mix) / [C_mix - eps_blk (2 C_mix - Pi_tr)]` (eq:app-stage081-Pi-zeta), and (ii) the unblocked limit `Pi_tr / C_mix = 1 + zeta_req` at `eps_blk = 0` (eq:app-stage081-unblocked). The notes go further and name the closed-form map `Q(zeta; eps_blk) := [1 + (1 - 2 eps_blk) zeta] / [1 - eps_blk zeta]`, list the anchor values `Q(0) = 1`, `Q(1) = 2`, the strictly-positive derivative `dQ/dzeta = (1 - eps_blk)/(1 - eps_blk zeta)^2`, the explicit numerical thresholds (`Pi_suff^(chi)/C_mix ~= 3.46622291347846`, etc.) at `lambda_mu = 1`, `eps_blk = 0`, and the hard ceiling `Pi_max^(F1)` with denominator-positivity condition `eps_blk < 1/zeta_max^(F1) ~= 0.405263689711371`. Part-III appendix row 140 summarises: "Thresholds in `Pi_tr` with blocking parameter."

## What the script claims to verify

Both scripts symbolically invert the Stage-35 support-demand law to obtain `Pi_tr = C_mix Q(zeta; eps_blk)`, exhibit the closed-form `Q`, assert the anchor values `Q(0) = 1` and `Q(1) = 2`, substitute the five carry-forward numerical thresholds from Stage 63 (`zeta_suff^(chi)`, `zeta_fail^(chi)`, `zeta_suff^(J)`, `zeta_fail^(J)`, `zeta_max^(F1)`), and confirm that in the unblocked limit `eps_blk = 0` each maps to `1 + zeta` to 10^-14 precision. The Mathematica script additionally cross-checks the symbolic `qq` (built from `piOfZeta/cMix` via `Solve`) against the explicit closed form `(1 + zeta - 2 epsBlk zeta)/(1 - epsBlk zeta)` and checks the blocking-ceiling reciprocal `epsCeiling * zetaMaxF1 - 1 == 0`. The SymPy script also prints the derivative `dQ/dzeta` (not asserted) and the blocking ceiling (not asserted).

## Paper <-> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `zeta_req = (Pi_tr - C_mix)/[C_mix - eps_blk(2 C_mix - Pi_tr)]` (eq:app-stage081-Pi-zeta) | both scripts define `zeta_expr` from this and solve for `Pi`; Mathematica additionally asserts `qq == (1 + zeta - 2 eps_blk zeta)/(1 - eps_blk zeta)` against the Solve-derived form | match |
| Unblocked limit `Pi_tr/C_mix = 1 + zeta_req` (eq:app-stage081-unblocked) | Mathematica asserts `(Pi_*/C_mix - (1 + zeta_*)) /. epsBlk -> 0 == 0` for each of the five thresholds; SymPy prints `sp.N(..., 20)` without asserting | match (Mathematica) / print-only (SymPy) |
| Anchor values `Q(0)=1`, `Q(1)=2` | both `expectZero[Q(0)-1]`, `expectZero[Q(1)-2]` | match |
| Notes-side `Q` closed form | Mathematica `expectZero["Q matches closed form", qq - (1 + zeta - 2 epsBlk zeta)/(1 - epsBlk zeta)]`; SymPy prints but does not assert | match (Mathematica) |
| Notes-side derivative `dQ/dzeta = (1 - eps_blk)/(1 - eps_blk zeta)^2` | SymPy prints `sp.diff(Q, zeta)`; Mathematica omits | partial (print-only; not load-bearing for paper Output, which is only the two equations) |
| Blocking ceiling `eps_blk < 1/zeta_max^(F1) ~= 0.40526...` | both compute `1/zetaMaxF1`; Mathematica asserts the reciprocal relation `epsCeiling * zetaMaxF1 - 1 == 0` | match (Mathematica) / print-only (SymPy) |
| Stage identification "Stage 081" | SymPy docstring at line 3 says "Stage 64"; SymPy banner at line 28 says "STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS"; Mathematica banner correctly says "STAGE 081" | mismatch (sympy script self-identifies as Stage 64; paper card, notes, appendix, and Mathematica script all use 081) |

`paper_alignment: partial` — the math content is faithful to the paper Output; the only finding is a stale stage label in the SymPy docstring/banner. The "partial" rather than "aligned" reflects that one row in the cross-check is `mismatch`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 40 | `expect_zero("Q(0)-1", Q.subs(zeta, 0) - 1)` | anchor `Q(0)=1` (notes-side anchor that pins eq:app-stage081-Pi-zeta at zeta=0) | yes |
| A2 | sympy | 41 | `expect_zero("Q(1)-2", Q.subs(zeta, 1) - 2)` | anchor `Q(1)=2` (notes-side anchor) | yes |
| M1 | mathematica | 54-55 | `expectZero["Q matches closed form", qq - (1 + zeta - 2 epsBlk zeta)/(1 - epsBlk zeta)]` | ties the Solve-derived inversion (eq:app-stage081-Pi-zeta) to the notes-side closed form | yes |
| M2 | mathematica | 59 | `expectZero["Q(0)-1", (qq /. zeta -> 0) - 1]` | anchor `Q(0)=1` | yes |
| M3 | mathematica | 60 | `expectZero["Q(1)-2", (qq /. zeta -> 1) - 2]` | anchor `Q(1)=2` | yes |
| M4 | mathematica | 80 | `expectApprox["Pi_suff^(chi)/C_mix at eps=0 matches 1+zeta", (piSuffChiOverC - (1 + zetaSuffChi1)) /. epsBlk -> 0, 0, 10^-14]` | unblocked limit (eq:app-stage081-unblocked) at zeta = zeta_suff^(chi) | yes |
| M5 | mathematica | 81 | analogous for `Pi_fail^(chi)` | unblocked limit at zeta = zeta_fail^(chi) | yes |
| M6 | mathematica | 82 | analogous for `Pi_suff^(J)` | unblocked limit at zeta = zeta_suff^(J) | yes |
| M7 | mathematica | 83 | analogous for `Pi_fail^(J)` | unblocked limit at zeta = zeta_fail^(J) | yes |
| M8 | mathematica | 84 | analogous for `Pi_max^(F1)` | unblocked limit at zeta = zeta_max^(F1) | yes |
| M9 | mathematica | 88 | `expectApprox["blocking ceiling reciprocal", N[epsCeiling*zetaMaxF1 - 1, 40], 0, 10^-14]` | reciprocal identity for the blocking ceiling (notes section 3 denominator positivity) | yes |

The Mathematica engine does the heavy assertion work. SymPy only asserts the two anchor values; everything else on the SymPy side is print-only. Per the paper card the load-bearing identities are the two `\stagefield{Output}` equations, both of which the Mathematica script exercises with assertions.

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** target_mismatch (script docstring/banner self-identifies as a different stage)
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:3` — docstring line `"""SymPy audit for Stage 64."""`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py:28` — banner `banner("STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS")`
- saved output `scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt:3` echoes the same "STAGE 64" header

**What's wrong:**
The SymPy script's docstring and banner self-identify the stage as "Stage 64", but the paper card (`paper/stages/stage_081.tex:1`, heading "Stage 081: Family--1 Demand Window in the Branch-Product Variable"), the source notes (`notes/stages/moving_throat_pde_stage081_family1_pi_thresholds.md:1`, title "Stage 081"), the part-III appendix table (`paper/appendices/stage_appendix_part03.tex:140`, row "081 & Branch-product demand window"), and the companion Mathematica script (`mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:38`, banner "STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS") all identify it as Stage 081. The mismatch is a leftover from the global stage renumbering (commit `0d09ef6 fully reorder the pde ledger`); the filename and Mathematica banner were updated but the SymPy in-file string literals were not. The math is correct — the script faithfully audits the Stage 081 claims — but its self-identification disagrees with the paper.

**Why this matters:**
A transcript header that points to a non-existent (or different) stage number is a debugging hazard. Anyone reading `scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt` sees the line `STAGE 64 — FAMILY-1 PRODUCT THRESHOLDS` and is left to guess whether the transcript was misfiled, the script audits the wrong stage, or the label is stale. It also weakens the cross-engine sanity check, since the two engines disagree on their own stage number. Per v2 policy, the auditor must not rewrite either side silently — even when the resolution is "obviously" script-side.

**Required change:**
See `## Resolve before fix_loop` in the directive. Codex must not auto-edit either side until the user confirms direction.

**Verification:**
After resolution, the SymPy script's line 3 docstring reads "SymPy audit for Stage 081." and line 28 banner reads `"STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS"`. A re-run of `redteam exec-sympy 081` produces a transcript whose third line reads `STAGE 081 — FAMILY-1 PRODUCT THRESHOLDS`.

## Independent-derivation check (Mathematica)

Both engines start from the same physical premise `zeta_expr = (Pi - Cmix)/(Cmix - eps_blk(2 Cmix - Pi))` and solve symbolically for `Pi`, but they diverge in assertion content:

- SymPy (lines 33-36): defines `zeta_expr`, calls `Pi_of_zeta = sp.solve(sp.Eq(zeta, zeta_expr), Pi)[0]`, simplifies, and forms `Q = sp.simplify(Pi_of_zeta / Cmix)`. The only assertions are the two anchor values `Q(0)=1` and `Q(1)=2`.
- Mathematica (lines 45-55): defines `zetaExpr`, calls `piOfZeta = piTr /. First[Solve[zeta == zetaExpr, piTr]]`, strips the `ConditionalExpression` wrapper that `Solve` introduces under non-trivial assumptions (line 51, retrofit; also line 53 on `qq`), forms `qq = FullSimplify[piOfZeta/cMix]`, and then asserts `expectZero["Q matches closed form", qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)]`. It then runs the eight further assertions M2-M9 listed in the inventory.

This is not a transliteration. The Mathematica script tests an identity that the SymPy script does not (the closed-form match against the notes), exercises the unblocked-limit and blocking-ceiling identities as assertions rather than prints, and uses `Solve` + `ConditionalExpression` strip rather than a plain `sp.solve`. The two engines arrive at the same `Q` form from the same premise via independent symbolic machinery.

`ConditionalExpression` strip retrofit check: lines 26, 51, and 53 of the .wl carry the strip. Line 46's `Solve` is followed by the strip on line 51 before any substitution at lines 59-60 and 68-72. No new `Solve` / `Reduce` adjacent to substitutions is added without the strip. The retrofit is in place.

## Engine cross-check

Both engines agree on the symbolic form of Q:

- SymPy output line 6: `Q(zeta;eps_blk) = (2*eps_blk*zeta - zeta - 1)/(eps_blk*zeta - 1)`
- Mathematica output line 8: `Q(zeta;eps_blk) = (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)`

These differ only by multiplying numerator and denominator by `-1`; algebraically identical.

Both engines agree numerically at `eps_blk = 0`:

| Threshold | SymPy (output lines 16-20) | Mathematica (output lines 18-27) |
|---|---|---|
| `Pi_suff^(chi)/C_mix` | `3.4662229134784601214` | `0` diff vs `1 + zetaSuffChi1` (PASS at 10^-14) |
| `Pi_fail^(chi)/C_mix` | `3.4675291327386998930` | PASS |
| `Pi_suff^(J)/C_mix`   | `3.4425757147717899187` | PASS |
| `Pi_fail^(J)/C_mix`   | `3.4675273685505798582` | PASS |
| `Pi_max^(F1)/C_mix`   | `3.4675292294560100537` | PASS |
| blocking ceiling      | `0.40526368971137149977` | `0.40526368971137148173360806443142537808` (PASS reciprocal) |

The 17th-digit difference in the blocking-ceiling values reflects SymPy's `sp.N(..., 20)` vs Mathematica's 40-digit `N[]`; both round to `0.405263689711371...` and the symbolic reciprocal check holds. `engines_agree: true`.

`outputs_fresh: true` — sympy script mtime Apr 1 12:39 < sympy output Apr ... wait: sympy script Apr 1 12:39; sympy output May 22 23:40 (newer). Mathematica script May 23 10:36; mathematica output May 25 00:21 (newer). Both transcripts post-date their scripts.

## Verdict justification

The math content of Stage 081 is correctly verified: both scripts independently invert the Stage-35 support-demand law to recover the branch-product map `Q(zeta; eps_blk)`, both confirm the anchor values, and the Mathematica engine additionally asserts the closed-form match, the unblocked-limit identity for each carry-forward threshold, and the blocking-ceiling reciprocal. I attempted to break the inversion by checking whether `Solve` could return a `ConditionalExpression` that the script forgets to strip before substituting `zeta -> 0` or `zeta -> 1` — but the retrofit on lines 51 and 53 of the .wl handles this, and no new Solve/Reduce was added without the strip. I attempted to find a hardcoded numerical value the paper does not state — but every numerical threshold (`zetaSuffChi1 = 2.46622291347846`, etc.) is a Stage-63 carry-forward explicitly listed in the notes file at the values used. I attempted to confirm that the v1 findings (hardcoded `qq` divorced from `piOfZeta`, tautological derivative check, eps=0 literal targets, blocking-ceiling-against-its-own-digits) were genuinely repaired — they were. The single remaining finding is a script-side stage-label residue ("Stage 64") in the SymPy docstring/banner from the pre-renumber commit; the paper card, notes, appendix, and Mathematica script all call this Stage 081. Verdict: `findings` (one low-severity `paper_misalignment`); `stop_cold: null`.

## Self-test notes

Traps checked: (1) `Solve` branch ambiguity — none; `zetaExpr` is linear in `piTr` so the inversion is unique. (2) `simplify` under aggressive assumptions — `$Assumptions` restricts to `epsBlk < 1, piTr > 0, cMix > 0, zeta >= 0`, none of which silently kills the branch; the SymPy script's `real=True` declaration is consistent. (3) Symmetry/parity traps — not applicable, no integrals in this unit. (4) Trivial-case pre-check — substituting `zeta = 0` into `(1 + zeta - 2 epsBlk zeta)/(1 - epsBlk zeta)` literally gives `1`; substituting `zeta = 1` gives `(2 - 2 epsBlk)/(1 - epsBlk) = 2`; both confirmed by inspection before reading the script outputs. (5) The `ConditionalExpression` retrofit was preemptively applied on lines 26 / 51 / 53 and no new bare `Solve` was introduced. (6) The "Resolve before fix_loop" question for F1 has only one mathematically-sensible answer (update script to "Stage 081") but per v2 policy I still route through the user gate rather than auto-fix.
