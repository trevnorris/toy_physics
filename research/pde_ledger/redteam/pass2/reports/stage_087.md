---
unit_id: 087
batch: III.5
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
  notes_stage_files: [moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish.md]
  paper_appendix: present
---

# Audit unit 087 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_087.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row 152 + `\input` line 292)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.txt`

## What the paper claims

Stage 087 is an exact-closure / open consolidation statement. Its `\stagefield{Output}` reads: "The one-ratio finish line \eqref{eq:app-stage087-final-condition}", where the boxed condition (eq. `app-stage087-final-condition`) is `\rho_\alpha=\frac{\alpha_{\rm req}}{\alpha_{\rm mix}}` "is the only remaining support-side input from the outgoing branch." The card's body adds that, given `\rho_\alpha`, the Family-1 support/source verdict follows from Stage 086, and that "No separate dependence on `N_Q^{target}`, `\beta_0`, `s_-`, `\lambda_-`, or `\widehat m_-` remains." `Inputs` are Stages 085-086. The notes go further and report the concrete unblocked Family-1 window at `\lambda_\mu=1`: `rho_alpha <= 3.46622291347846` (guaranteed success), `rho_alpha >= 3.46752913273870` (guaranteed failure), `rho_alpha < 3.46752922945601` (absolute constructive ceiling). The appendix row 152 paraphrases: "Remaining Family-1 support/source gap reduced to one outgoing-branch loading ratio." Distinct deliverables: (D1) the symbolic one-ratio finish line `rho_alpha = alpha_req/alpha_mix`; (D2) the three numeric window literals (notes only).

## What the script claims to verify

Both engines define the blocked criterion `zeta_req(rho_alpha, eps_blk) = (rho_alpha - 1)/(1 - eps_blk*(2 - rho_alpha))`, print it and its derivative, then assert (i) the unblocked reduction `zeta_req(rho_alpha, 0) == rho_alpha - 1`, exercising the "one-ratio" reduction at the natural unblocked limit. The Mathematica script additionally asserts the derivative equals the hand-supplied closed form `(1 - eps_blk)/(1 - eps_blk*(2 - rho_alpha))^2`. Both scripts then "cross-check" the three Family-1 window literals (`rho_suff`, `rho_fail`, `rho_max`) which the docstrings/comments describe as anchoring "against the upstream stage-086 quoted values to catch renumber or transcription drift." Mathematica also numerically checks the three `zeta = rho-1` window values. The scripts' docstrings explicitly state Stage 087 is a "checkpoint-consolidation statement, not a fresh derivation," with the actual cancellation chain verified upstream in stages 081-086.

## Paper ↔ script cross-check

| paper deliverable | script-side check | status |
|---|---|---|
| D1 — `rho_alpha = alpha_req/alpha_mix` is the one remaining support input; criterion reduces to a single ratio | `unblocked zeta_req` assert (py:55, wl:44) + `d zeta_req exact formula` (wl:43); FINAL LEDGER print restates `rho_alpha = alpha_req/alpha_mix` | match (reduction `zeta_req(·,0)=rho_alpha-1` is the substantive, non-tautological exercise of the one-ratio criterion) |
| D2 — window literals 3.46622291347846 / 3.46752913273870 / 3.46752922945601 (notes §2) | `rho_*` cross-checks (py:73-75, wl:61-63); `zeta_*` numeric checks (wl:73-75) | partial — values are present and CORRECT vs the notes, but the assertion compares each literal against the SAME in-file literal, so it cannot perform the drift-detection its comment claims (see F1) |

`paper_alignment: aligned` — both deliverables are present and numerically consistent with the card/notes; the defect is a hollow verification mechanism, not a value or target mismatch.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 55 | `expect_zero(zeta_req(eps_blk=0) - (rho_alpha-1))` | D1 (one-ratio reduction) | yes |
| A2 | sympy | 73 | `expect_close(rho_suff, Float("3.46622291347846",30), 1e-13)` | D2 | no — self-comparison vs identical literal |
| A3 | sympy | 74 | `expect_close(rho_fail, Float("3.46752913273870",30), 1e-13)` | D2 | no — self-comparison |
| A4 | sympy | 75 | `expect_close(rho_max, Float("3.46752922945601",30), 1e-13)` | D2 | no — self-comparison |
| A5 | mathematica | 43 | `expectZero[dZeta - dZetaExpected]` | D1 (derivative form) | yes |
| A6 | mathematica | 44 | `expectZero[(zetaReq/.epsBlk->0) - (rhoAlpha-1)]` | D1 (one-ratio reduction) | yes |
| A7 | mathematica | 61 | `expectApprox[rhoSuff, 3.46622291347846, 1e-14]` | D2 | no — self-comparison vs identical literal |
| A8 | mathematica | 62 | `expectApprox[rhoFail, 3.46752913273870, 1e-14]` | D2 | no — self-comparison |
| A9 | mathematica | 63 | `expectApprox[rhoMax, 3.46752922945601, 1e-14]` | D2 | no — self-comparison |
| A10 | mathematica | 73 | `expectApprox[zetaSuff, 2.466222913478460121439184, 1e-14]` | D2 (zeta window) | no — target is just `rhoSuff-1` re-typed |
| A11 | mathematica | 74 | `expectApprox[zetaFail, 2.467529132738699892968270, 1e-14]` | D2 | no — target is `rhoFail-1` re-typed |
| A12 | mathematica | 75 | `expectApprox[zetaMax, 2.467529229456010053667114, 1e-14]` | D2 | no — target is `rhoMax-1` re-typed |

A1/A5/A6 are the genuine, load-bearing assertions. A2-A4 and A7-A12 are self-referential and cannot fail (see F1).

## Findings

### F1 — insufficient_verification

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:58-60,69-75`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:46-75`

**What's wrong:**
Both scripts comment that the `rho_*` checks anchor the window literals "against the upstream stage-086 quoted values to catch renumber/transcription drift." SymPy comment (py:69-72): "Cross-check the Family-1 ratio-window literals against the upstream stage 086 quoted values to catch renumber/transcription drift." Mathematica comment (wl:53-55): "The rho_X literals below are carried from stage 086; the cross-check below anchors them against the upstream stage-086 paper values."

But the assertions compare each literal against *the identical literal re-typed in the same file*:
- py:58 `rho_suff = sp.Float("3.46622291347846")` vs py:73 `expect_close(..., rho_suff, sp.Float("3.46622291347846", 30), ...)` — same digit string on both sides; the ~1e-16 "diff" in the output (`scripts/output/...txt:11-13`) is only the 15-digit-vs-30-digit float-parse rounding, not a real comparison.
- wl:57 `rhoSuff = ToExpression["3.46622291347846\`20"]` vs wl:61 `expectApprox["rho_suff^(chi) vs stage-086", rhoSuff, 3.46622291347846, 10^-14]` — identical digits; the Mathematica output records `diff = 0.` (`mathematica/output/...txt:11-15`), confirming a literal self-comparison.
- The Mathematica `zeta_*` checks (wl:73-75) are likewise self-referential: at `epsBlk->0`, `zetaSuff = rhoSuff - 1`, and the target `2.466222913478460121439184` is exactly `3.46622291347846 - 1`, so the check is `(rhoSuff-1) - (rhoSuff-1) ≈ 0`.

No part of either file actually reads, imports, or reproduces stage-086's output, so the check cannot detect the renumber/transcription drift its own comment names as its purpose. The only substantive verifications in this stage are A1/A5/A6 (the `zeta_req(·,0) = rho_alpha-1` reduction and the derivative form), which do correctly exercise the paper's one-ratio claim.

This is NOT a value mismatch: the three literals are present and CORRECT relative to the notes (`notes/...stage087...md:60-63`), so the deliverable values reconcile. The defect is purely that the drift-guard is hollow — it asserts a number against a copy of itself.

**Why this matters:**
The script advertises a guard against exactly the renumber/transcription drift that the project's numbering history is plagued by, but the guard is a no-op: if someone copy-edited `rho_suff` to a wrong value, both sides of the comparison would change together and the assert would still pass. The stage's protection against the failure mode it claims to protect against is illusory.

**Required change:**
Make the window-literal checks actually exercise an independent relation instead of comparing a literal to a copy of itself. Two acceptable, in-scope options (Codex picks the cleaner one and applies to BOTH engines for symmetry):
- (a) Drop the redundant `expect_close`/`expectApprox` self-comparisons (py:73-75, wl:61-63) and instead assert the *ordering/consistency* relations that the window must satisfy and that would actually break under bad transcription: `rho_suff < rho_fail < rho_max` AND `rho_fail < rho_max` with the tight gap `rho_max - rho_fail < 1e-6` (the constructive-ceiling property), i.e. structural relations among the three literals rather than each-against-itself. These can fail if any literal is mistyped.
- (b) If an upstream anchor is genuinely wanted, leave a single literal block but reword the comments (py:69-72, wl:46-55) to drop the false "cross-check against upstream stage-086" framing and call them what they are ("window literals carried from the notes; ordering sanity check below"). Then add the ordering asserts from (a).
Either way, the load-bearing A1/A5/A6 reduction asserts stay untouched.

NOTE: this is low-severity scaffolding hygiene, not a correctness or paper-alignment defect. If the orchestrator prefers to leave it (consolidation stage, values verified correct against notes, genuine reduction asserts present), that is defensible; the finding documents the hollow guard so it is on record.

**Verification:**
After fix, the sympy/mathematica outputs should show new ordering/consistency checks (e.g. `rho_suff < rho_fail`) that pass, and the self-comparison lines (py:73-75 / wl:61-63 with `diff = 0.` / `diff = 1e-16`) should either be gone or accompanied by a relation that could fail. The A1/A5/A6 reduction lines must still print and pass. The three literals must remain `3.46622291347846 / 3.46752913273870 / 3.46752922945601` (unchanged, they match the notes).

## Independent-derivation check (Mathematica)

The `.wl` is NOT a transliteration of the `.py`. It shares the same `zeta_req` definition and the same unblocked reduction assert (unavoidable — that is the single physics object of the stage), but it independently (a) supplies and verifies a hand-written closed form for the derivative `dZetaExpected = (1 - epsBlk)/(1 - epsBlk*(2 - rhoAlpha))^2` (wl:39,43), which the SymPy script only prints, never asserts; and (b) adds numeric `zeta_*` window checks (wl:73-75) absent from SymPy. The variable choreography and helper structure differ (`FullSimplify[Together[Expand[...]]]` vs `sp.simplify(sp.expand(...))`; `expectApprox`/`Assumptions` vs `expect_close`). The shared steps reflect a shared trivial physics object, not line-by-line porting. No `mathematica_transliteration` finding.

## Engine cross-check

Both engines agree:
- `zeta_req` printed in algebraically equivalent forms: SymPy `(1 - rho_alpha)/(-eps_blk*(rho_alpha-2)-1)` (`scripts/output/...txt:5`) ≡ Mathematica `(-1 + rhoAlpha)/(1 + epsBlk*(-2 + rhoAlpha))` (`mathematica/output/...txt:5`). Same rational function (numerator/denominator both negated).
- `unblocked zeta_req = 0` in both (sympy txt:7, math txt:9-10).
- Window values identical to printed precision: `zeta at success ratio = 2.46622291347846...` in both (sympy txt:8, math txt:17). Same for failure/max.
- Cross-check diffs both ~1e-16/0 (sympy txt:11-13, math txt:11-15). No `engine_disagreement`.

## Verdict justification

The stage's load-bearing claim — that the Family-1 support/source side reduces to the single ratio `rho_alpha = alpha_req/alpha_mix` — is genuinely exercised by the `zeta_req(rho_alpha, 0) = rho_alpha - 1` reduction (A1/A6) and the derivative-form check (A5), both non-tautological and engine-agreed. The three window literals are present and numerically CORRECT against the notes, so there is no value/target paper_misalignment and the value reconciliation is complete. The one real defect (F1, low severity) is that the `rho_*`/`zeta_*` "cross-checks against upstream stage-086" compare each literal against an identical copy of itself in the same file (output diffs of `0.`/1e-16 confirm), so the advertised renumber/transcription-drift guard is hollow. Because the values themselves are correct and the substantive reduction asserts hold, the verdict is `findings` (one low-severity `insufficient_verification`), not `clean` and not `stop_cold`. Read the paper card, notes, and appendix row; the script's verified claim matches the paper's one-ratio Output.

DEFERRED (out of scope, noted not flagged): cross-references to other stages — sympy:12 "former stage 65", sympy:13 "former stage 69 closure", sympy:80 / output line 18 "Stage-69 ratio window" (the window is now stage 086 per appendix row 150). These are stale CROSS-references to other stages, deferred per the numbering policy. The scripts' own canonical "Stage 087" self-labels are correct.

## Self-test notes

Checked: (1) variable-independence — the only derivative `D[zetaReq, rhoAlpha]` is taken w.r.t. a symbol `zetaReq` genuinely depends on, so it is non-trivially nonzero (output shows the rational form, not 0); no identically-zero-derivative trap. (2) Symmetry/parity — no integrals in this stage. (3) Trivial-case — substituting `eps_blk=0` into `(rho_alpha-1)/(1-eps_blk*(2-rho_alpha))` gives `rho_alpha-1`, so A1/A6 residual is genuinely 0 and the assert is meaningful, not vacuous. (4) Paper round-trip — the proposed F1 fix (ordering/consistency asserts) introduces no new constant and keeps the three literals at their notes-matching values, so it cannot create a new paper_misalignment. (5) The self-comparison trap (A2-A4, A7-A12) is the finding itself, confirmed by the `diff = 0.`/1e-16 output lines.

## Value Reconciliation (pass-2 augmentation)

Outputs are FRESH: sympy `.py` mtime 2026-05-27 10:18:19, output 10:24:52; mathematica `.wl` 10:18:25, output 10:26:20. Both outputs are newer than their scripts. Reconciliation uses script source + committed outputs; no scripts were run.

Deliverable-level reconciliation:

| value | source (py/wl + output line) | .tex/.md location | status |
|---|---|---|---|
| `rho_alpha = alpha_req/alpha_mix` (symbolic one-ratio finish line) | py:79 FINAL LEDGER, wl FINAL print; sympy txt:17 | `.tex:17` (boxed eq) + notes:15,38 | MATCH |
| `rho_suff = 3.46622291347846` | py:58, wl:57; sympy txt:11, math txt:11 | notes:60 | MATCH |
| `rho_fail = 3.46752913273870` | py:59, wl:58; sympy txt:12, math txt:13 | notes:62 | MATCH |
| `rho_max = 3.46752922945601` | py:60, wl:59; sympy txt:13, math txt:15 | notes:63 | MATCH |
| `zeta_req(rho_alpha,eps_blk)` closed form `(rho_alpha-1)/(1-eps_blk*(2-rho_alpha))` | py:50, wl:37; sympy txt:5, math txt:5 | (probe form, restated implicitly in notes §2 criterion) | MATCH (consistent with notes one-ratio criterion) |

INTERNAL (genuine scaffolding / derived intermediates, no prose finding expected):
- `d zeta_req / d rho_alpha` printed form and its expected closed form (derivative-verification intermediate).
- `zeta at success/failure/max ratio` = 2.466222913478460.../2.467529132738699.../2.467529229456010... (these are `rho_X - 1`; the zeta-variable window, derived in-script; the notes report the rho-window, not the zeta-window).
- cross-check `diff` values (~1e-16, 0.), tolerances `1e-13`/`1e-14`, PASS flags.

reconciliation: complete; 5 deliverable values checked, 0 misaligned.

All stated deliverable values (the symbolic one-ratio finish line and the three numeric window literals) reconcile exactly with the `.tex` card and/or `.md` notes. The F1 finding is an `insufficient_verification` (hollow self-comparison guard), NOT a value mismatch — no paper_misalignment, no user resolution required.
