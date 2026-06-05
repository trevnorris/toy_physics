# Batch III.2 (stages 049-060) — pde_ledger SECOND PASS

Date: 2026-06-05

Part III.2 — Tracking, zeta thresholds, asymmetry, boost (continuum kernel,
transport, reachability, microclosure tail).

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
  and the arbiter grep resolved the scope consistently (and caught the stale
  committed Mathematica outputs the audits missed; see the INFRA NOTE below).
- Codex as the sole fix-applier, with the `.py`-only-vs-`.wl` seat distinction: the
  11 label-only stages touch only `.py` source (no seat), while the one `.wl`-editing
  math stage (060) consumes a Mathematica seat for its build; exec deferred until the
  build finished.
- Clean per-stage verify agent for each stage.

## Result

All 12 stages reached `verified` at v2 depth + value-reconciliation augmentation.
**12/12 verified, `material_change=false` on all 12. No stop-cold, no blocked,
0 Codex deviations, all iter-1 exit 0.** Every one of the 12 was `verdict:findings`
(each flagged at least one stale self-label, so each got an unambiguous self-label
canonicalized + outputs refreshed).

| Outcome | Stages |
| --- | --- |
| Real script-side math finding → fixed → verified (1) | 060 |
| Label-only self-label fix + output refresh → verified (11) | 049, 050, 051, 052, 053, 054, 055, 056, 057, 058, 059 |

One checkpoint in range — **051** — re-verified clean (no certified constant moved).
**No genuine `paper_misalignment` anywhere.** The full value reconciliation came
back 0-misaligned (below), confirming no script value disagrees with the card/notes.

## Resolutions

### Real script-side finding (1)

- **060** — ONE genuine finding + a self-label fix:
  - **F2 `tautological_check`** — the Mathematica `wl:140` assertion
    `expectZero["Xi_micro - Lambda^2 L^2/(Theta T_X)", xiMicro - lambdaPhi^2*ell^2/(theta*tX)]`
    subtracted `xiMicro` from a verbatim copy of its own `wl:132` definition
    (identically 0, could not fail). Replaced with
    `expectZero["Xi_micro chi-route equals D/M-route", xiMicroFromChi - xiMicroFromDM]`
    — a genuine cross-check between two INDEPENDENTLY constructed routes (the
    susceptibility route `chiSigma→1/theta` at `wl:133` vs the Einstein/diffusion
    route `dSigma→mSigma*theta` at `wl:136`, where `mSigma` must cancel). The check
    is now non-tautological; the value `Xi_micro = Lambda_phi² L²/(Theta T_X)` is
    UNCHANGED, and the deliverable stays covered (also by the retained
    `wl:134-138`/`wl:141-142`). The 060 `.wl` remains an independent engine (verify
    agent confirmed it is not a transliteration of the SymPy path).
  - **F1 label-only** — self-banner number canonicalized to stage 060.
  - `material_change=false` (de-tautologization over an already-correct value; no
    verified RESULT value moved).

### Label-only self-label fixes (11 + 060's F1) — numbering-drift interim policy

Per the user-confirmed interim policy (Reading 2) and
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` §Interim policy: for a
`verdict:findings` stage, fix the audit-flagged UNAMBIGUOUS self-labels (matching
each file's canonical banner) + refresh its outputs; defer cross-refs and ambiguous
refs to the dedicated pass. Every III.2 stage was `verdict:findings`, so each got at
least one numbering self-label corrected. The rule (matching the III.1 precedent)
corrects the stale NUMBER only and PRESERVES each label's existing format:

- 2-digit docstrings / closing-prints (e.g. `Stage 32`→`Stage 49`) kept 2-digit,
- 3-digit filename-docstrings on 051/052 kept 3-digit (the on-disk filename is
  3-digit, `stage34`→`stage051`),
- already-correct 2-digit `STAGE NN` banners LEFT UNPADDED (no cosmetic padding —
  that belongs to the dedicated plan),
- 049 also fixed its `.wl:93` closing `Stage 32`→`Stage 49`.

All `.py`/`.wl` source diffs are **strip-the-number identical to HEAD** except
060.wl's one math line (the F2 de-tautologization above). Same class as the I.2
stage-021 / II.1 032-036 / III.1 9-stage fixes.

## Value reconciliation (pass-2 augmentation)

Applied on all 12 stages; **126 deliverable values checked batch-wide, 0 misaligned.**
Per stage: 049=10, 050=11, 051=9, 052=12, 053=11, 054=14, 055=9, 056=12, 057=7,
058=11, 059=10, 060=10. No MISMATCH, no MISSING-DELIVERABLE anywhere.

## Dual-engine / mirror status

All 12 already carried both independent engines from pass 1; independence
re-confirmed on every stage. **060.wl was EDITED** (tautology → independent
cross-route check) and **remains an independent engine** (verify agent confirmed,
not a transliteration). **0 new `.wl`, 0 mirror reclassification, 0 sanctioned
mirrors.**

## INFRA NOTE — exec_logs vs committed outputs

`$RT exec-sympy` / `exec-mathematica` write the fresh transcript to
`redteam/pass2/exec_logs/` but do NOT auto-overwrite the committed `output/*.txt`;
the orchestrator must sed-refresh the committed transcripts
(`sed '1,/^---$/d;/^# exit_code:/d' <log> > <out>`). In III.2 the Mathematica
committed outputs for 050–059 were stale (Codex only re-ran the `.wl` for 049/060),
caught by the **arbiter grep** and fixed by the refresh. After the refresh the
arbiter grep confirmed every committed output is free of stale old-epoch
self-banners; all 24 reliability-gate runs (12 SymPy + 12 Mathematica) exit 0.

## Deferred to the dedicated SCRIPT/OUTPUT-band numbering pass

(`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, PENDING — content-keyed, never
offset-sweep.) Left untouched by this batch, as policy dictates — the CROSS-refs and
compound dual-epoch refs in the `.py`/`.wl` source:
- 050 `py:61` `Stage 32`→stage049;
- 051 `py:20-21` `Stage 050/034` + `py:126`/`wl:87` `Stage 047/030` (compound
  dual-epoch);
- 055 `py:73` `Stage-35`→stage052;
- 056 `py:7` `Stage-36`→stage053;
- 059 `py:6` `Stage-41`→stage058 + `py:9`/`75` `Stage-39`→stage056;
- 060 `py:159` `Stage-39`→stage056.

These map by content (the +17 epoch), not by sweep.

## Trackers

6 prose trackers synced (PAPER_CLEANUP **P5-06** = no new paper/notes items — III.2
made ZERO paper/.tex/notes edits; the only source change is a script (060.wl); all
numbering deferrals route to the dedicated SCRIPT/OUTPUT-band plan; checkpoint 051
re-verified clean with no constant change, no EM stages, no new constants).
