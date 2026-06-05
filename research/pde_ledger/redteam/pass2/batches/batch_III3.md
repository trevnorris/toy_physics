# Batch III.3 (stages 061-072) — pde_ledger SECOND PASS

Date: 2026-06-05

Part III.3 — Microclosure, gain thresholds, equilibrium, walls (microscopic gain
thresholds, parent-action gain, parent thresholds, equilibrium alignment, thin-wall
confinement, wall figure of merit, sech-Gaussian resonance, resonance thresholds,
final reduced verdict [checkpoint], GNLS wall-shell, tanh-wall branch, explicit
branch thresholds).

## Method

- v2 paper-grounded auditor PLUS the **exhaustive script→doc value-reconciliation
  augmentation** (`redteam/pass2/RECONCILIATION_AUGMENTATION.md`) — every audit
  agent read both its rendered prompt and the augmentation doc.
- 12 clean per-stage audit agents in parallel (audit agents execute nothing → 0
  Mathematica seats).
- Independent exec reliability gate — orchestrator re-ran BOTH engines for every
  stage (sequential, ≤1 seat) and refreshed the committed transcripts via
  `sed '1,/^---$/d;/^# exit_code:/d' <log> > <out>`. The exec re-run + an **arbiter
  grep** on the source AND the refreshed committed outputs are ground truth — the
  audit agents again applied an inconsistent threshold to numbering labels (some
  wrote full label-fix directives, some under-called; one missed a filename-style
  `stageNN_` self-label that the `\b`-less re-grep caught on 068), so the orchestrator
  authored every directive and made the self-vs-cross call.
- Codex as the sole fix-applier, with the `.py`-only-vs-`.wl` seat distinction: the
  8 label-only stages touch only `.py` (0 seat → run concurrently in two waves of 4),
  the one `.wl`-editing math stage (070) consumes a seat and ran solo; the 3
  refresh-only stages needed no Codex (orchestrator re-run + refresh only).
- Clean per-stage verify agent for each stage.

## Result

All 12 stages reached `verified` at v2 depth + value-reconciliation augmentation.
**12/12 verified, `material_change=false` on all 12. No stop-cold, no blocked,
0 Codex deviations, all iter-1 exit 0.** Every one of the 12 was `verdict:findings`
(stale committed outputs everywhere; 9 also carried at least one unambiguous stale
self-label).

| Outcome | Stages |
| --- | --- |
| Real script-side math finding (de-taut) + self-label fix + refresh → verified (1) | 070 |
| Label-only self-label fix + output refresh → verified (8) | 061, 062, 063, 064, 065, 066, 068, 071 |
| Refresh-only (stale_output; source already canonical) → verified (3) | 067, 069, 072 |

One checkpoint in range — **069** (`final_reduced_verdict`) — cleared the higher bar
(below). **No genuine `paper_misalignment` anywhere.** Full value reconciliation came
back 0-misaligned (below).

## Resolutions

### Real script-side finding (1)

- **070** — ONE genuine `tautological_check` (F1) + two self-labels (F2):
  - **F1 `tautological_check`** — the "anchoring cross-check" `I_1/J_1 = 4πa²ℓ` was a
    self-cancelling tautology in BOTH engines: `(4πa²ℓ·If/Hw)/(If/Hw) ≡ 4πa²ℓ`, the
    `If/Hw` factor cancels so the result is independent of the sech moment value
    (SymPy `py:61-63`; Mathematica `wl:92-94`, whose own comment says "independent of
    I_f's value"). The genuinely informative quantity (the sech profile yields
    `I_f = 2/3`, `I_g = 14/15`) was only PRINTED, never asserted. **Fix:** made the
    sech moments load-bearing — SymPy now asserts `expect_zero("sech-profile moment
    I_f = 2/3", If_sym - 2/3)`; Mathematica now asserts `IfNum ≈ 2/3` and
    `IgNum ≈ 14/15` against the already-computed NIntegrate values (`tol=10^-10`).
    The structural `I_1/J_1` ratio line is RETAINED as documentation alongside.
  - **⚠️ Orchestrator false-positive guard caught a WRONG constant in the audit's
    proposed fix:** the audit (and the script's own stale print annotation at
    `wl:86`) claimed the sech second moment is `I_g = 8/15`. The orchestrator computed
    it independently — `I_g = ∫(f'')²dξ = 2 − 16/3 + 64/15 = 14/15` (≈0.9333), NOT
    8/15 (≈0.5333) — and directed Codex to assert `14/15`. The arbiter exec confirmed:
    the 30-digit NIntegrate yields `0.93333…`, `PASS: sech-profile moment I_g = 14/15`.
    Had `8/15` been asserted, the script would have FAILED. The directive also
    corrected the wrong `wl:86` print annotation `8/15`→`14/15`.
  - The symbolic `kappa`/`W_wall`/`Xi` deliverable assertions are UNCHANGED and still
    pass; the sech moments are the stage's OWN internal anchor (not a deliverable, not
    consumed downstream — Stage 071 separately uses the tanh profile with `I_f=1/3`).
    The 070 `.wl` remains an independent engine (verify agent confirmed; the numeric
    profile cross-check is genuinely independent of the SymPy path).
  - **F2 self-labels** — `py:3` filename `stage53`→`stage070` (3-digit) + `py:5` prose
    `Stage 53`→`Stage 70` (2-digit).
  - `material_change=false` (de-tautologization over an already-correct value + a
    print-annotation correction; no verified RESULT value moved).

### Checkpoint 069 — higher-bar scrutiny

`final_reduced_verdict` is a consolidation checkpoint. Higher-bar findings: it pins
**no** numeric constant — `C_res² = 0.994418…`/`P_res = 1.005612…` are upstream
(Stages 067/068) deliverables carried as the FREE symbols `Cres2Prim` (Mathematica) /
`Pres_gap` (SymPy), never re-asserted against themselves. The substantive deliverables
(three-zone verdict, matched-window width, exact side-band widths,
`P_res-1=(1-C²)/C²`, strict interior-point ordering) are genuine symbolic identities
that would fail if a threshold form were wrong, and the two engines build them by
independent routes (SymPy `Pres`-first/multiplication; Mathematica `Cres2`-first with
`PresGap` recovered via `Solve`). One secondary check (`Pres_from_ratio`, SymPy A6) is
tautological-by-construction but NON-load-bearing, and the same identity gets genuine
cross-route content from the Mathematica `Solve` inversion — so no certification rests
on it; left as honest documentation (its overstated comment noted as an observation,
NOT a finding — prose-only, out of red-team script-math scope). The only 069 finding
was `stale_output` (committed banner `STAGE 052`); resolved by refresh. Checkpoint
cleared.

### Label-only self-label fixes (8 + 070's F2) — numbering-drift interim policy

Per the user-confirmed interim policy (Reading 2) and
`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md` §Interim policy: for a
`verdict:findings` stage, fix the UNAMBIGUOUS self-labels (matching each file's
canonical banner) + refresh outputs; defer cross-refs / variable names / ambiguous
refs to the dedicated pass. The rule corrects the stale NUMBER only and PRESERVES each
label's existing format (III.1/III.2 precedent):

- 2-digit prose docstrings / closing-prints kept 2-digit (`Stage 44`→`Stage 61`,
  `All Stage 45 … passed`→`Stage 62`, …),
- 3-digit filename-style docstrings kept 3-digit (`stage51`→`stage068`,
  `stage53`→`stage070`, `stage54`→`stage071`),
- already-correct 2-digit `STAGE NN` banners LEFT UNPADDED (no cosmetic padding —
  deferred to the dedicated plan; the audit agents' directives proposed padding and
  cross-ref fixes, which the orchestrator stripped out before applying).

Self-labels fixed: 061 `py:3`; 062 `py:3`+`py:112`; 063 `py:3`+`py:124`; 064 `py:3`;
065 `py:3`; 066 `py:3`; 068 `py:4`; 070 `py:3`+`py:5`; 071 `py:3`+`py:5`. 067, 069,
072 had no stale self-labels in source (canonical already) → refresh-only. All
`.py`/`.wl` source diffs are **strip-the-number identical to HEAD** except 070's three
added math lines (the F1 de-tautologization + annotation fix). Same class as I.2
stage-021 / II.1 032-036 / III.1 / III.2 fixes.

## Value reconciliation (pass-2 augmentation)

Applied on all 12 stages; **116 deliverable values checked batch-wide, 0 misaligned.**
Per stage: 061=12, 062=9, 063=8, 064=7, 065=11, 066=9, 067=8, 068=7, 069=9, 070=10,
071=13, 072=13. No MISMATCH, no MISSING-DELIVERABLE anywhere. (Note: 066 = wall
figure-of-merit — no stale "168"/"100" figure-of-merit prefactor exists at this stage;
its only constant is the `4π` shell prefactor, which matches card/notes/both engines.)

## Dual-engine / mirror status

All 12 already carried both independent engines from pass 1; independence
re-confirmed on every stage. **070.wl was EDITED** (added the two numeric sech-moment
asserts + corrected the I_g annotation) and **remains an independent engine** (verify
agent confirmed; its concrete numeric profile cross-check is a genuinely independent
route the SymPy script does not perform). 072.wl is transliteration-leaning but
accepted per the established mirror policy (a pure substitution+limit-identity stage
whose algebra/limits are independently re-executed). **0 new `.wl`, 0 mirror
reclassification, 0 sanctioned mirrors.**

## INFRA NOTE — exec_logs vs committed outputs (III.2 lesson re-confirmed)

`$RT exec-sympy` / `exec-mathematica` write the fresh transcript to
`redteam/pass2/exec_logs/` but do NOT auto-overwrite the committed `output/*.txt`; the
orchestrator must sed-refresh. In III.3, all 12 committed transcripts (mtime
2026-05-11/05-22) predated the scripts (mtime 2026-06-03) and most predated June-3
content additions (e.g. 064's GENERAL-H / CONTINUOUS-CAUCHY sections, 066's
`dW/dJ1`/`dW/dT_X` monotonicity lines, 068's `(C-form)/(P-form)` bands + `P_res
numeric residual`, 070's entire numeric-profile section). The orchestrator re-ran all
24 engines (12 SymPy + 12 Mathematica, all exit 0) and sed-refreshed every committed
output; the refreshed transcripts adopt the canonical bare pass-2 format (matching the
III.2-verified outputs — the old `# … Audit Output / # Status: PASS / EXIT_CODE: 0`
wrapper was the pre-pass-2 runner format). An **arbiter grep** on the refreshed
committed outputs confirmed no stale old-epoch (NNN−17) appears as a self-banner or
closing — the only stale-band remnants are CROSS-references (064 `Stage-45/46`
assertion labels, 069 `Stage-49/51` family prints, 070 `Stage-47/48` normalization
prints), all deferred.

## Deferred to the dedicated SCRIPT/OUTPUT-band numbering pass

(`redteam/NUMBERING_SCRIPT_OUTPUT_BAND_PLAN.md`, PENDING — content-keyed, never
offset-sweep.) Left untouched by this batch, as policy dictates — the CROSS-refs,
variable names, and ambiguous refs in `.py`/`.wl` source:
- 063 `py:76` `Stage-44`→stage061;
- 064 `py:25`/`122`/`180` + `wl:104` `Stage-45/46`→stages062/063 (one inside an
  assertion label);
- 065 `py:22` `Stage-44`→stage061;
- 066 `py:14`/`59` `Stage-48`→stage065;
- 069 `py:9`/`11`/`26`/`94`/`112`/`120`/`176`/`179` + `wl:99`/`116`/`164`/`167`
  `Stage-49/51`→stages066/068 (`py:26` "Stage 049" is the same ref in 3-digit form —
  ambiguous self-vs-cross, kept for the dedicated pass);
- 070 `py:57`/`59` + `wl:87`/`88` `Stage-47/48`→stages064/065, plus the variable names
  `J1_stage48` (py:62-63) / `J1Stage48` (wl:78).

These map by content (the +17 epoch), not by sweep.

## Trackers

6 prose trackers synced (PAPER_CLEANUP **P5-07** = no new paper/notes items — III.3
made ZERO paper/.tex/notes edits; the only source changes are scripts (8 label-only
`.py`, plus 070's `.py`+`.wl` math/annotation lines); all numbering deferrals route to
the dedicated SCRIPT/OUTPUT-band plan; checkpoint 069 cleared the higher bar with no
constant change; no EM stages in range; no new pinned constants — 070's sech moments
`I_f=2/3`, `I_g=14/15` are an internal anchor, now asserted).
