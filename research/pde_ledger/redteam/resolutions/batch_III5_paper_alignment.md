---
batch: III.5
created: 2026-05-27
status: applied_and_scripts_pass
items:
  - Q1: 087 F1 — is stage 087 a status/checkpoint consolidation or load-bearing? → APPROVED (a) status/checkpoint
  - Q2: 087 F2 — what to do with the tautological window-value sanity checks → APPROVED (i) expect_close vs upstream
  - Q3: 089 F1 — Pe_req=0 link verification → APPROVED (a) strengthen scripts
notes:
  - All 6 sympy + 6 mathematica scripts re-run sequentially (single-seat) and exit 0.
  - One orchestrator hot-fix on 088 sympy: `Y_rho.subs(omega**2/Omega_Q**2, u)` failed after sp.simplify reshaped the denominator into `(Omega_Q**2 - omega**2)` form. Fixed by substituting `omega**2 -> u * Omega_Q**2` then sp.simplify.
---

# Batch III.5 paper-alignment resolutions

Two stages flagged `paper_misalignment` requiring user direction (087 F1 and 089 F1). All other findings (086 F1/F2, 087 F3, 088 F1/F2/F3, 089 F2/F3/F4, 090 F1/F2/F3) are non-paper_misalignment and route directly to the standard fix loop.

Twelve banner relabels (all 6 stages × both engines have stale `STAGE 68`/`069`/`70`/`070`/`71`/`071`/`72`/`072`/`73`/`073` banners from the pre-renumber numbering — global-renumber leftovers, same pattern as III.4) are orchestrator-direct script-side label fixes with unambiguous direction; they are applied directly by the orchestrator and do not appear here.

The orchestrator recommendation pass writes `## Recommendation` blocks below each question with `direction:` and rationale (Codex stalled on III.4 Q1 so we kept the orchestrator-direct fallback proven there). The user reviews and approves.

---

## Q1 — Stage 087 F1: scripts verify auxiliary `zeta_req(rho_alpha, eps_blk)` but the paper claims a *cancellation* theorem

**Auditor finding (full):** `redteam/reports/stage_087.md` F1 (subtype `script_missing_paper_claim` with banner-level `target_mismatch`).

### Three pieces of evidence

- `paper/stages/stage_087.tex:6-8` — `\stagefield{Purpose}{Stage~087 records that the explicit Family-1 support/source side has been reduced to a single outgoing-branch loading ratio.}` and `\stagefield{Inputs}{Stages~085--086.}`
- `paper/stages/stage_087.tex:13-21` boxes `\rho_\alpha = \alpha_{\rm req}/\alpha_{\rm mix}` as the only remaining support-side input and states `No separate dependence on N_Q^target, beta_0, s_-, lambda_-, or mhat_- remains.`
- `notes/stages/moving_throat_pde_stage087_outgoing_branch_loading_ratio_finish.md:5-17` — "Stages 65-69 compressed the full reduced moving-throat PDE to one surviving quadrupole residual, then to one explicit product window, and finally to one pure loading-ratio criterion. **This stage states the cleanest final verdict of that reduction chain.**"

The notes explicitly say the cancellations were performed in stages 65-69 (now post-renumber 081-086). Stage 087 itself introduces no new symbols — its derivation is a one-paragraph restatement that the prior chain has compressed to one number.

### Options

- **(a) Status/checkpoint consolidation.** Stage 087 is what its `\stagefield{Purpose}` says it is: a *records that* paragraph carrying forward the cancellation theorem proved in stages 081-086. Direction: relabel banners "STAGE 70/070" → "STAGE 087" (this happens orthogonally for all 6 stages in the orchestrator-direct sweep), add a docstring/comment line in both scripts naming the upstream verification sources for the cancellation, and keep the `zeta_req(rho_alpha, eps_blk)` evaluation as downstream consistency. F2 below handles the tautological window-value checks.
- **(b) Load-bearing.** Rewrite both scripts to introduce `s_-`, `lambda_-`, `beta_0`, `mhat_-`, `Pi_tr`, `C_mix`, `Pe_req` as free symbols and exhibit the explicit support-theorem expression in those symbols, then assert that under the closure `mhat_-^2 Gamma_{5,-} = 2 G / (5 c^5)` the dependence collapses to `rho_alpha` only. This duplicates work already done in stages 068-069/081-086.
- **(c) Mixed.** Keep as carry-forward (a) but add a minimal symbolic check that `zeta_req` has the right shape for the upstream window. Adds shape coverage but introduces redundancy.

### Destination-verification grep (run locally by orchestrator)

- `paper/stages/stage_087.tex` `\stagefield{Inputs}` line explicitly names stages 085-086, not 081-086 (i.e., paper itself is treating 087 as carry-forward from the immediately upstream consolidation pair).
- Notes file paragraph at line 17 reads "Everything else on the explicit support/source side is now fixed." and section 4 "Best current expert verdict" frames stage 087 as a status statement.
- No downstream stage (088-096) references stage 087 for any new identity; downstream `\stagefield{Inputs}` lines name stages 086 (for 088) and stage 087 only for the loading ratio statement carried forward.

## Recommendation

- direction: **a** (status/checkpoint consolidation)
- rationale: The paper's `\stagefield{Purpose}` language ("Stage 087 records that...") and `\stagefield{Inputs}{Stages 085--086.}` both designate 087 as a consolidation. The notes file is unambiguous: cancellations live in stages 65-69 (= post-renumber 081-086). The stage introduces no new symbols. Treating 087 as load-bearing would double the verification cost without checking anything not already verified upstream.
- downstream_impact: None mathematically. Cosmetically: scripts get a comment block naming `scripts/output/moving_throat_pde_stage068*` (= post-renumber 081), `scripts/output/moving_throat_pde_stage069*` (= 082), and the just-verified 085/086 outputs as the upstream sources for the cancellation chain. (Note: the parenthetical "(now post-renumber 081-086)" labels in this resolution are explanatory — the codex apply just adds the comment naming the post-renumber filenames it can find.)
- notes: Per `feedback_script_review_depth` and `feedback_redteam_priorities`, the script's math correctness is what matters; checkpoint-consolidation stages that *only restate* upstream results don't need a second derivation. F2's tautological residuals get a separate (sub)question because they are independent of F1's direction once F1 = (a) is fixed.

---

## Q2 — Stage 087 F2: the `zeta_suff`/`zeta_fail`/`zeta_max` window-value "anchors" are tautological

**Auditor finding (full):** `redteam/reports/stage_087.md` F2 (`tautological_check`).

### What the script does

- `scripts/.../stage087..._sympy_audit.py:27` — `expect_zero("unblocked zeta_req", zeta_req.subs(eps_blk, 0) - (rho_alpha - 1))`. Reduces algebraically to `(rho_alpha - 1) - (rho_alpha - 1) == 0`.
- `scripts/.../stage087..._sympy_audit.py:41-43` — three `expect_zero("zeta_X - <literal>", sp.N(zeta_X - sp.Float("<literal>"), 18))` calls, each comparing `rho_X - 1` to a literal that equals `rho_X - 1`.
- `mathematica/.../stage087..._mathematica_audit.wl:44, 58-60` — identical pattern.

The three numeric literals (`2.46622291347846`, `2.46752913273870`, `2.46752922945601`) are quoted in the notes file (`notes/stages/moving_throat_pde_stage087_*.md:60-63`) but their derivation lives in upstream stage 069 (= post-renumber 082).

### Options

- **(i) Replace with upstream-anchored cross-checks.** Replace the four tautological `expect_zero` lines with `expect_close(zeta_X, sp.Float("<value from upstream stage 069 output>"), tol=1e-12)` style comparisons. This catches the failure mode where a renumber miscopies the literal value, or where an upstream change shifts the window.
- **(ii) Demote to `print` statements.** The values are useful to display but the "anchor" check is illusory. Convert to `print(f"zeta_suff = {sp.N(zeta_suff, 18)}")` etc. — pure logging.

### Destination-verification grep

- The same three literals also appear in `notes/.../moving_throat_pde_stage087_*.md:59-63` as the prose statement of the window. These are the canonical narrative form; the script should compare against them, not against `rho_X - 1`.
- Stage 069/082's `.txt` output is the upstream source. The verifier wave will need to confirm those values still match if we go with (i).

## Recommendation

- direction: **i** (replace with `expect_close` cross-checks against upstream stage 069/082 quoted values)
- rationale: Direction (i) is strictly stronger than (ii) without adding meaningful complexity — it converts illusory "anchors" into real numeric anchors. This is the same approach 089 F3 prescribes for its analogous tautologies and the same approach we used for 086 F1's `expect_close` helper.
- downstream_impact: None.
- notes: Codex's apply block writes a small `expect_close(name, val, target, tol)` helper (similar to 086's planned helper from its F1 directive) at module top, then four calls comparing `zeta_suff/zeta_fail/zeta_max` and the unblocked `zeta_req(eps_blk=0)` against the notes-quoted upstream targets. Verifier confirms residuals are now `~1e-15` (real numerical equality), not `0` (definitional).

---

## Q3 — Stage 089 F1: scripts verify `zeta_min < A_F1` but paper boxes `Pe_req = 0` (CHECKPOINT)

**Auditor finding (full):** `redteam/reports/stage_089.md` F1 (subtype `script_missing_paper_claim`). **Stage 089 is a CHECKPOINT** (`is_checkpoint: true`) so the higher bar applies.

### The link chain on the paper side

`paper/stages/stage_089.tex:13-29`:
- Eq. `app-stage089-zero-bias-supply` (line 13-17): `At zero transport bias, Omega_{Pe=0}=1, so Family-1 supplies zeta_F1(0) = A_F1 ≈ 1.00005192880220.`
- Eq. `app-stage089-success` (line 18-23) — boxed: `zeta_req^min = 1/3 < A_F1`.
- Eq. `app-stage089-Pe-zero` (line 25-29) — boxed Output: `Pe_req = 0`.

The chain is: `Omega(Pe=0) = 1` ⟹ `zeta_F1(0) = A_F1` ⟹ (`zeta_req^min = 1/3 < A_F1`) ⟹ `Pe_req = 0`. The script only verifies the third link (`zeta_min < A_F1`).

### What the script omits

The Omega expression in `scripts/.../stage089..._sympy_audit.py:41-44`:

```
Omega = simplify( pi * Pe * (2*Pe*exp(Pe) + pi) / ((4*Pe^2 + pi^2) * (exp(Pe) - 1)) )
```

At `Pe = 0` this is `0/0`. The `Omega(Pe=0) = 1` limit is a non-trivial l'Hospital — not verifiable by substitution. The Mathematica script has the same structure. Neither engine computes the limit; neither asserts `Pe_req = 0` directly.

Note: `Omega(Pe=0) = 1` is a Stage-062 result (= post-renumber stage 075? — needs upstream check), but per the v2 audit doctrine, a checkpoint stage that states the link should verify the link, not hand-wave it.

### Options

- **(a) Strengthen scripts.** Add explicit symbolic checks in both engines:
  1. `sp.limit(Omega, Pe, 0) == 1` (and `Limit[omegaPe[pe], pe -> 0] == 1` in Mathematica)
  2. `sp.limit(zeta_F1, Pe, 0) - A_F1 == 0` (and Mathematica analogue)
  3. Define `Pe_req = 0 if zeta_req^min < A_F1 else sp.nsolve(...)` and assert `Pe_req == 0` when the precondition holds.
  
  This closes the chain end-to-end and meets the checkpoint bar.

- **(b) Carry-forward comment.** Declare `Omega(Pe=0) = 1` as upstream (Stage 062/075) and add an in-script comment naming the exact upstream verification file/line. Anchors the link rather than verifies it. Weaker for a checkpoint.

- **(c) Update paper.** Weaken the boxed `Pe_req = 0` to a corollary in prose, with `zeta_req^min < A_F1` as the actual `Output`. Removes the link from the script's responsibility. Likely a paper-side regression — `Pe_req = 0` is a load-bearing downstream input.

### Destination-verification grep

- `paper/stages/stage_090.tex` references stage 089 directly. (Need to verify what 090 carries forward — if 090 uses `Pe_req = 0` as input, then option (c) cascades downstream.)
- `notes/stages/moving_throat_pde_stage089_*.md:99-109` quotes `Pe_req = 0` as the explicit minimal-branch result. Notes already treat it as established.
- The Omega expression matches the standard Family-1 transport map shape from stage 062/075. Confirming `Omega(Pe→0) = 1` is a one-line `sp.limit` call.

## Recommendation

- direction: **a** (strengthen scripts to close the chain end-to-end)
- rationale: Stage 089 is a CHECKPOINT — the higher bar requires substantive assertions on the boxed Output. The missing links (`Omega(Pe=0) = 1`, `zeta_F1(0) = A_F1`, `Pe_req = 0`) are all one-liner `sp.limit` checks in SymPy and `Limit[..., pe -> 0]` in Mathematica. The Omega expression is `0/0` at `Pe = 0` and Mathematica's hybrid solver / SymPy's series-based `Limit` both handle it. Adding these costs ~6 lines per engine and converts an unsupported boxed Output into a verified one.
- downstream_impact: None paper-side. Script-side: directive prescribes the new `sp.limit` calls + a `Pe_req = 0` assertion gated on `zeta_req^min < A_F1`. Per pitfall #10, do NOT use `sp.nsolve` for the conditional branch; the precondition holds at construction so the assignment `Pe_req = sp.Integer(0)` is unconditional in this stage.
- notes: F2 (mathematica_transliteration), F3 (tautological_check), F4 (hardcoded_result) are independent of F1's direction and route to Codex regardless. F4's auditor-suggested `sp.nsolve` rederivation of `Pe_suff_chi`, `Pe_fail_chi` is at risk for pitfall #10 — flag in the Codex apply prompt to use `mpmath.findroot` with bracketing, or just keep the literals with the provenance-comment path from the directive's option (b).

---

## Banner relabel sweep (orchestrator-direct, no user gate)

All 6 III.5 stages have stale "STAGE 68/069/70/070/71/071/72/072/73/073" banners — global-renumber leftovers. Same pattern as III.4. Applied directly without user gate:

| Stage | Engine | Current banner | Target |
|---|---|---|---|
| 085 | .py | `STAGE 68 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS` | `STAGE 085 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS` |
| 085 | .wl | `STAGE 068 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS` | `STAGE 085 — EXACT CANCELLATION OF OUTGOING-NORMALIZATION FACTORS` |
| 086 | .py | `STAGE 69 — FAMILY-1 PURE LOADING-RATIO WINDOW` | `STAGE 086 — FAMILY-1 PURE LOADING-RATIO WINDOW` |
| 086 | .wl | `STAGE 069 — FAMILY-1 PURE LOADING-RATIO WINDOW` | `STAGE 086 — FAMILY-1 PURE LOADING-RATIO WINDOW` |
| 087 | .py | `STAGE 70 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE` | `STAGE 087 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE` |
| 087 | .wl | `STAGE 070 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE` | `STAGE 087 — FINAL REDUCED FINISH-LINE IN THE LOADING-RATIO VARIABLE` |
| 088 | .py | `STAGE 71 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE` | `STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE` |
| 088 | .wl | `STAGE 071 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE` | `STAGE 088 — LOADING-RATIO EXTRACTION FROM THE MINIMAL ISOTROPIC MODULE` |
| 089 | .py | `STAGE 72 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH` | `STAGE 089 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH` |
| 089 | .wl | `STAGE 072 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH` | `STAGE 089 — EXPLICIT FAMILY-1 VERDICT FOR THE MINIMAL ISOTROPIC BRANCH` |
| 090 | .py | `STAGE 073 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION` | `STAGE 090 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION` |
| 090 | .wl | `STAGE 073 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION` | `STAGE 090 — UPDATED REDUCED STATUS AFTER THE LOADING-RATIO EXTRACTION` |

Also: stage 089 sympy docstring (line 3) and stage 090 wl/py docstrings reference stale stage numbers — those get fixed orchestrator-direct alongside the banner.
