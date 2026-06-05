# Batch III.1 (stages 037-048) — pde_ledger SECOND PASS

Date: 2026-06-05

Part III.1 — Continuum kernel, generalized branch, rank-2.

## Method

- v2 paper-grounded auditor PLUS the **exhaustive script→doc value-reconciliation
  augmentation** (`redteam/pass2/RECONCILIATION_AUGMENTATION.md`) — every audit
  agent read both its rendered prompt and the augmentation doc.
- 12 clean per-stage audit agents in parallel (audit agents execute nothing → 0
  Mathematica seats).
- Independent exec reliability gate — orchestrator re-ran BOTH engines for every
  touched stage (sequential, ≤1 seat) and refreshed the committed transcripts via
  `sed '1,/^---$/d;/^# exit_code:/d' <log> > <out>`. The exec re-run + an **arbiter
  grep** on the source are ground truth — the audit agents applied an inconsistent
  threshold to the numbering labels (some wrote label-fix directives, some deferred),
  and the arbiter grep resolved the scope consistently (and caught 3 inline-comment
  self-labels the audits missed; see below).
- Codex as the sole fix-applier. 5 Codex waves of 2 (≤1 seat per wave; the two
  `.wl`-editing math stages 039/043 paired last; exec deferred until builds finished).
- Clean per-stage verify agent for each stage (incl. 045).

## Result

All 12 stages reached `verified` at v2 depth + value-reconciliation augmentation.
**12/12 verified, `material_change=false` on all 12. No stop-cold, no blocked,
0 Codex deviations, all iter-1 exit 0.**

| Outcome | Stages |
| --- | --- |
| Real script-side math findings → fixed → verified (2) | 039, 043 |
| Label-only self-label fix + output refresh → verified (9) | 037, 038, 040, 041, 042, 044, 046, 047, 048 |
| Deferred-clean (source already canonical + outputs fresh) → verified (1) | 045 |

No checkpoints in III.1. **No genuine `paper_misalignment` anywhere** (037's audit
mislabeled a script self-label as `paper_misalignment` — it is label-only, no value
dispute; no user value-gate was triggered). The full value reconciliation came back
0-misaligned (below), confirming no script value disagrees with the card/notes.

## Resolutions

### Real script-side findings (2)

- **039** — TWO genuine findings + a self-label fix:
  - **F1 `tautological_check`** — the only check intended to exercise the boxed
    direction factor `R_U` was identically true by construction: `z0=κ0·g_W·(1+ρ0)`,
    `z1=κ1·g_W·(1+ρ0/(1+δU))`, so `z1·(1+ρ0) − (κ1/κ0)·z0·(1+ρ0/(1+δU)) ≡ 0` (the κ0
    cancels; `R_U` never referenced). Replaced (both engines) with the falsifiable
    `z1/z0 − (κ1/κ0)·R_U == 0`, where `R_U=(1+ρ0/(1+δU))/(1+ρ0)` is independently
    defined from `ρ0,δU` (NOT via z0/z1) → fails if R_U's closed form is wrong.
  - **F2 `insufficient_verification`** — the surviving exact product law
    `R_target·M_mix = 8Λ(1−ε_W,split)/π²` (the headline "factorization survives"
    deliverable, notes §5) was only `print`ed. Added `expect_zero("product law
    survives", …)` in both engines (LHS from independently-simplified
    M_mix_split·R_target_split; RHS the notes' closed form via ε_W,split).
  - **F3 label-only** — docstring + subbanners `Stage 22`/`22.k`→`Stage 39`/`39.k`.
  - Both new checks PASS in both engines; engines agree; no regression.
    `material_change=false` (assertions added over already-correct values + a
    de-tautologization; no verified RESULT value moved).
- **043** — ONE genuine finding + a self-label fix:
  - **F2 `insufficient_verification`** — the `M_supp` baseline-value check
    re-substituted the SAME literal `B=8/π²` into both sides of an identity already
    proven in the free symbol `B` (line 141), so it could not fail and never derived
    the number. Replaced (both engines) with a check that DERIVES the baseline from
    the stage's frozen overlap constant: `B=κ0²=(9/11)·σ`, `σ=88/(9π²)` ⇒ `8/π²`,
    asserted via `expect_zero("baseline B = 8/pi^2 from frozen sigma", …)`. Now the
    value `8/π²` is exercised (fails if κ0² were mis-stated); structural-form check
    retained.
  - **F1 label-only** — docstring + subbanners `Stage 26`/`26.k`→`Stage 43`/`43.k`.
  - New check PASS in both engines; matches notes:149+151 exactly; `material_change=false`.

### Label-only self-label fixes (9) — numbering-drift interim policy

Per the user-confirmed interim policy (Reading 2, decided 2026-06-05) and
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` §Interim policy: for a `verdict:findings`
stage, fix the audit-flagged UNAMBIGUOUS self-labels (matching each file's canonical
banner) + refresh its outputs; defer clean-stage stale outputs, cross-refs, and
ambiguous refs to the dedicated pass. Every III.1 stage was `verdict:findings` (all
flagged stale_output), so each got its unambiguous SymPy self-labels canonicalized:

- module docstring filename + header (`stage{20,21}`/`Stage {20,21}`→canonical on
  037/038; `Stage {23,24,25,27,29,30,31}`→canonical on 040/041/042/044/046/047/048),
- `.py` sub-stage section indices `NN.k`→canonical `NN.k` (039/040/041/042/043/044),
- closing `All Stage-NN … passed` pass-lines (046/047),
- **+ 3 inline-comment self-refs the audits missed but the arbiter grep caught**:
  040:136 `(section 23.2)`→`40.2`, 041:107 `section 24.1`→`41.1` (applied via a Codex
  iter-2 on those two units).

Every `.py` diff is **strip-the-number identical to HEAD** (label-only; verified
mechanically). The Mathematica `.wl` source was already canonical on all 12 (banners
3-digit, plain `1.`-`k.` subbanners, canonical closings) → **no `.wl` source edit on
any label-only stage**; the stale `.wl` `.txt` outputs were cured by the orchestrator
re-run alone. Variable names (`F_stage18`, `G_stage19`, `F_stage23`) and all
cross-references to OTHER stages were left untouched.

### 045 — deferred-clean

Source self-labels already canonical (docstring `Stage 045`, closing `All Stage-045`,
banner `STAGE 045`) and outputs already fresh; the only finding was a low-severity
`stale_output` on 4 `Stage-27` CROSS-refs (to upstream stage 044), deferred. No fix;
verified clean (math re-confirmed: coherence identity `g_B g_R=g_W g_S`,
quadratic→tracking collapse, D/N `G_tr`/`F_tr` forms all non-tautological + aligned).

## Value reconciliation (pass-2 augmentation)

Applied on all 12 stages; **143 deliverable values checked batch-wide, 0 misaligned.**
Per stage: 037=12, 038=10, 039=12, 040=10, 041=8, 042=14, 043=16, 044=9, 045=12,
046=13, 047=14, 048=13. No MISMATCH, no MISSING-DELIVERABLE anywhere.

## Dual-engine / mirror status

All 12 already carried both independent engines from pass 1; the re-audit re-confirmed
genuine independence on every stage (e.g. 039 `.wl` derives δ_split via
`a1Direct/a0Expected−1` vs SymPy's postulate-and-confirm; 043 `.wl` uses `Det[…]`,
`Limit[…,δU→∞]`, `Series`, and three numeric sign test-points with no SymPy
counterpart). **No new `.wl`; 039/043 `.wl` were EDITED to add the F-fix assertions
(still independent routes — verify agents confirmed not transliterations); 0 mirror
reclassification, 0 sanctioned mirrors.**

## Deferred to the dedicated SCRIPT/OUTPUT-band numbering pass

(`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, PENDING — content-keyed, never
offset-sweep.) Left untouched by this batch, as policy dictates:
- **CROSS-refs** to other stages in the `.py`/`.wl` source (e.g. 037 `Stage-17/19`;
  038 `Stage-20`; 039/041 `Stage-21`/`Stage-23`; 040 `Stage-18`/`Stage 22`; 042/044
  `Stage-24`/`Stage-25`/`Stage-23`; 045 `Stage-27`; 047 `Stage 28`). These map by
  content (the +17 epoch), not by sweep.
- **One AMBIGUOUS self-vs-cross ref**: 047:121 comment `the exact Stage-30
  support-loading coefficient` — `M_supp` is established upstream (043/048), so this
  is genuinely ambiguous (stale self-label vs cross-ref); it is a comment (never
  reaches the transcript), so deferral leaves the committed output canonical.

## Trackers

6 prose trackers synced (PAPER_CLEANUP **P5-05** = no new paper/notes items — III.1
made ZERO paper/.tex/notes edits; all numbering deferrals route to the dedicated
SCRIPT/OUTPUT-band plan; no checkpoints, no EM stages, no new constants).
