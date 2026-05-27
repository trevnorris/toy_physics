---
unit_id: 087
batch: III.5
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
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
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part03.tex` (row line 152)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.txt`

## What the paper claims

Stage 087's `\stagefield{Output}` is the boxed equation `rho_alpha = alpha_req / alpha_mix` is the only remaining support-side input from the outgoing branch. The card and notes state that the explicit Family-1 support theorem has been reduced from a multivariable closure problem (involving `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, `Pe_req`) to a single normalized loading ratio. The notes additionally tabulate three numeric thresholds at `lambda_mu = 1` and unblocked limit: `rho_alpha <= 3.46622291347846` (guaranteed success), `rho_alpha >= 3.46752913273870` (guaranteed failure), and `rho_alpha < 3.46752922945601` (absolute constructive ceiling). The stage is explicitly a checkpoint/finish-line consolidation of prior reduction work; the appendix row (line 152) labels it "Outgoing branch loading-ratio finish line".

## What the script claims to verify

Both scripts define an auxiliary function `zeta_req(rho_alpha, eps_blk) = (rho_alpha - 1)/(1 - eps_blk*(2 - rho_alpha))`, take a derivative, then substitute eps_blk -> 0 and rho_alpha -> {rho_suff, rho_fail, rho_max} to check numeric values. The Mathematica script additionally checks that the symbolic derivative matches a hand-written closed form `(1 - eps_blk)/(1 - eps_blk*(2 - rho_alpha))^2`. The script banner is labeled "STAGE 70" rather than 087. No assertion in either script exercises the paper's actual `Output` (that the Family-1 support theorem's other variables have cancelled).

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side coverage |
|---|---|
| (D1) `rho_alpha = alpha_req/alpha_mix` is the *only* remaining support-side input; `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, `Pe_req` have all dropped out | **missing** — neither script exhibits or checks any of these cancellations |
| (D2) Numeric thresholds `rho_suff = 3.46622291347846`, `rho_fail = 3.46752913273870`, `rho_max = 3.46752922945601` are correctly carried | **partial** — the literals are present in the script as inputs to `zeta_req`, but the script only verifies `zeta = rho - 1` (a tautology) rather than deriving these threshold values from any upstream criterion |
| (D3) Functional form / monotonicity of the supporting `zeta_req` map | **extra** — script tests `zeta_req` properties (Mathematica gives a non-trivial derivative match); the paper card does not mention `zeta_req` for stage 087 |

`paper_alignment` = partial (most deliverables missing or tautologically covered; one extra not in paper).

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored? |
|---|---|---|---|---|---|
| A1 | sympy | 27 | `expect_zero("unblocked zeta_req", zeta_req.subs(eps_blk, 0) - (rho_alpha - 1))` | none — pure simplification identity | no (tautological) |
| A2 | sympy | 41 | `expect_zero("zeta_suff - 2.46622291347846", zeta_suff - 2.46622291347846)` | D2 numerically | no (tautological: zeta_suff = rho_suff - 1 = 2.46622291347846 by construction) |
| A3 | sympy | 42 | `expect_zero("zeta_fail - 2.46752913273870", ...)` | D2 numerically | no (same tautology) |
| A4 | sympy | 43 | `expect_zero("zeta_max - 2.46752922945601", ...)` | D2 numerically | no (same tautology) |
| A5 | math | 43 | `expectZero["d zeta_req exact formula", dZeta - dZetaExpected]` | D3 (not a paper claim) | yes (non-trivial derivative residual) |
| A6 | math | 44 | `expectZero["unblocked zeta_req", (zetaReq /. epsBlk -> 0) - (rhoAlpha - 1)]` | none | no (tautological) |
| A7 | math | 58 | `expectApprox["zeta_suff numeric check", zetaSuff, 2.466222913478460121439184, 10^-14]` | D2 | no (tautological: zetaSuff = rhoSuff - 1) |
| A8 | math | 59 | `expectApprox["zeta_fail numeric check", ...]` | D2 | no (same tautology) |
| A9 | math | 60 | `expectApprox["zeta_max numeric check", ...]` | D2 | no (same tautology) |

## Findings

### F1 — paper_misalignment

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_087.tex:14-20` (Output / boxed equation)
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish.md:22-38` (list of variables that drop out)
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:19-48` (entire body)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:32-65` (entire body)

**Subtype:** `script_missing_paper_claim` (combined with target_mismatch on the banner label).

**What's wrong:**
Paper card (lines 14-20) frames the stage's `Output` as the boxed identity that `rho_alpha = alpha_req/alpha_mix` is the *only* remaining support-side input — i.e. the support theorem has lost dependence on `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, and `Pe_req` (notes lines 22-38 enumerate these). The script bodies do not introduce any of those symbols and do not perform or check any cancellation. Instead they reuse a stage-70-era auxiliary `zeta_req(rho_alpha, eps_blk) = (rho_alpha - 1)/(1 - eps_blk*(2 - rho_alpha))` and inspect it at eps_blk = 0. Both script banners read "STAGE 70" / "STAGE 070", confirming the scripts were authored for a different (earlier) stage and were not rewritten for stage 087's `Output`.

**Why this matters:**
A `clean` verdict on these scripts is meaningless: their `PASS` says nothing about whether stage 087's actual reduction claim is true. The reader cannot tell whether the support theorem actually became independent of `s_-`, `lambda_-`, etc., or whether the script was simply written against the wrong target.

**Required change:**
ROUTED TO USER (see directive's `## Resolve before fix_loop`). Either: (a) accept that stage 087 is a status/checkpoint consolidation whose cancellations were already verified in stages 068-069/085-086 — in that case update the script banners to "STAGE 087" and either keep zeta_req as a downstream consistency check or replace these scripts with explicit pointers to the upstream cancellation scripts; (b) keep this stage as load-bearing and add explicit cancellation checks to the scripts.

**Verification:**
After user picks a direction, the new check (or the new banner/docstring) appears in both scripts and the saved outputs reflect them.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:27,41-43`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:44,58-60`

**What's wrong:**
The check at sympy line 27 is `zeta_req.subs(eps_blk, 0) - (rho_alpha - 1)`. Since `zeta_req = (rho_alpha - 1)/(1 - eps_blk*(2 - rho_alpha))`, substituting eps_blk = 0 gives literally `(rho_alpha - 1)/1 = (rho_alpha - 1)`; the residual is 0 by direct algebraic simplification, not by any physical input. Identical for mathematica line 44.

Worse: the numeric checks at sympy lines 41-43 are
- `zeta_suff = zeta_req.subs({rho_alpha: rho_suff, eps_blk: 0}) = rho_suff - 1 = 3.46622291347846 - 1 = 2.46622291347846`
- the assertion then checks `zeta_suff - 2.46622291347846 == 0`.

That is `(rho_suff - 1) - (rho_suff - 1) == 0` and has no chance of failing. The hand-written literal on the RHS is just `rho_suff - 1` typed out. The Mathematica `expectApprox` calls at lines 58-60 do the same comparison with a 1e-14 tolerance — equally tautological at the level of physics; only floating-point arithmetic could ever make them fail.

**Why this matters:**
Five out of six assertions in the sympy script and four of five "expect" calls in Mathematica are guaranteed by construction. They consume verifier compute and produce green PASS lines that imply no real check was done.

**Required change:**
ROUTED TO USER (see directive). If F1's resolution adopts direction (b) — keep stage as load-bearing — these tautologies must be replaced with substantive checks. If direction (a) — stage is a status/checkpoint — these can be removed or replaced with explicit `print` statements that don't pretend to assert.

**Verification:**
After fix, the assertions either (i) are removed, or (ii) reference values derived from an upstream source rather than re-deriving and comparing against themselves.

### F3 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_sympy_audit.py:19-48`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish_mathematica_audit.wl:32-63`

**What's wrong:**
Both scripts use the same variable choreography:
- SymPy lines 21-24: declare `rho_alpha, eps_blk`, define `zeta_req`, compute `sp.diff(zeta_req, rho_alpha)`, print.
- Mathematica lines 34-42: declare `rhoAlpha, epsBlk`, define `zetaReq`, compute `D[zetaReq, rhoAlpha]`, print.

- SymPy line 27: `expect_zero("unblocked zeta_req", zeta_req.subs(eps_blk, 0) - (rho_alpha - 1))`.
- Mathematica line 44: `expectZero["unblocked zeta_req", (zetaReq /. epsBlk -> 0) - (rhoAlpha - 1)]`.

- SymPy lines 30-35 and Mathematica lines 46-52: declare the same three numeric ratios and substitute them.

The Mathematica adds one independent check (`d zeta_req exact formula` against a hand-written derivative), which is genuinely non-tautological. But the rest of the file is a direct port. Given the unit's tiny scope, this is borderline; flagging as `low` so the reviewer can choose to address it together with F1's resolution.

**Why this matters:**
The two-engine policy is meant to catch algebra errors via independent re-derivation; a port catches only typing errors. For stage 087 specifically, this is mostly cosmetic because the bulk of the algebra is tautological anyway (see F2). After F1/F2 are resolved, the engines should be re-derived independently.

**Required change:**
Once F1 resolution is known, re-author the Mathematica script (if kept) to derive the kept claims using Mathematica-native primitives (`Solve`, `Reduce`, `Series`) rather than mirroring sympy step-for-step.

**Verification:**
The two scripts no longer share a one-to-one assertion correspondence.

## Independent-derivation check (Mathematica)

The Mathematica script is largely a transliteration of the sympy script (see F3). The single exception is the `expectZero["d zeta_req exact formula", dZeta - dZetaExpected]` check at line 43, which is a non-trivial residual: it compares `D[zetaReq, rhoAlpha]` against the hand-written closed form `(1 - epsBlk)/(1 - epsBlk*(2 - rhoAlpha))^2`. This closed form is correct (verified by mental algebra) and provides genuine value-add over the sympy script. The remaining checks mirror sympy line-for-line.

## Engine cross-check

Both engines pass and report consistent intermediate values:
- `zeta_req` form: sympy gives `(1 - rho_alpha)/(-eps_blk*(rho_alpha - 2) - 1)`; mathematica gives `(-1 + rhoAlpha)/(1 + epsBlk*(-2 + rhoAlpha))`. Multiply numerator and denominator of each by -1 and they are identical.
- Numeric `zeta_suff` / `zeta_fail` / `zeta_max`: agree to better than 1e-14.

Engines agree, but as noted they're agreeing on tautologies.

## Verdict justification

The scripts run clean and the engines agree, so on a purely operational reading, stage 087 has green output. But the audit's job is to attack: the load-bearing assertions are tautological (F2), and none of them exercise the paper's stated `Output` (F1). The banner labels "STAGE 70" / "STAGE 070" are smoking-gun evidence that these scripts were written for an earlier stage's purpose and copy-promoted. The mathematica script is a partial transliteration (F3). Because the resolution of F1 depends on whether stage 087 is intended to be a status/checkpoint consolidation (in which case the missing cancellation checks are upstream) or a load-bearing claim (in which case they belong here), F1 routes to the user, and the directive holds the orchestrator.

## Self-test notes

- Variable independence: confirmed by mental algebra — `(rho_alpha - 1)/(1 - eps_blk*(2 - rho_alpha))` at eps_blk = 0 collapses to `(rho_alpha - 1)`, so `expect_zero(... - (rho_alpha - 1))` is identically zero (F2).
- Symmetry/parity: not applicable; no integrals in this unit.
- Trivial-case substitution: confirmed `zeta_suff = rho_suff - 1 = 2.46622291347846` exactly when eps_blk = 0; the literal on the RHS of the expect_zero is bit-identical (F2 tautology).
- Path specifications: both scripts exist in the standard locations (scripts/, mathematica/); the directive does not need to create new files, only flag for user resolution.
- Paper round-trip: confirmed the boxed `\stagefield{Output}` at stage_087.tex:14-20 maps to a *cancellation* claim, not to a `zeta_req` identity; the scripts therefore do not address it.
