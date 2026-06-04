# Numbering Broad-Sweep — execution plan (the deferred ~1,200-mention pool)

**Date:** 2026-06-03
**Status:** PREP COMPLETE — awaiting user gate to start band 1. Read with
`redteam/NUMBERING_BROAD_SWEEP_WORKLIST.md` (machine inventory) and
`redteam/NUMBERING_RECON_CLASSD_MAP.md` (the completed Class-D pass this extends).

## What this is
The third and final stage-number reconciliation pass. Phase 1 (deterministic A/B/C,
`e2a4780`) and Phase 2 (Class-D high-signal body cross-refs in the 57 citing notes +
deferred tail, `28784b7`) are done. This pass covers the **remaining body cross-refs
in ALL ~236 stage notes** — predominantly the notes that were NEVER in the 57-file
Class-D citing set (e.g. stage198/199/201/204/211, the 188–234 endgame band), plus
special-form refs (`\mathcal S_{NNN}`) and low-band stragglers (stage119/121/123…).

## Scope (mechanical inventory: `broad_sweep_scan.py`)
- **PRIMARY = 300 candidates / 58 notes** — forward bare `Stage NNN` refs (cited > self)
  + all special forms. This is the sweep's target.
- **Band 001–090 = 0 candidates.** The drift is entirely in the post-EM-extension
  content (091→253); the early EM/linearized band is clean. Big bands: 164–200 (131),
  201–234 (111).
- **SECONDARY = 2354 backward bare refs** — mostly GENUINE upstream citations. NOT in
  this sweep. A handful may be stale (scan saw −31/−67 offsets) but yield is low and
  noise is high. Held for a separate explicit decision after PRIMARY lands.
- Scope is NOTES ONLY (`notes/stages/*.md`). Script (.py/.wl) body cross-refs are out
  (their self-banners were Class C; the one special `.wl` cross-ref set was done in P2).

## Ground rules (identical to Class D — do not relax)
1. **Canonical = filename / `\label{stage:NNN}`.** NEVER an offset sweep. Confirmed
   offsets in this neighborhood run +17/+34/+51/+68/+85 and the SAME stale number maps
   to DIFFERENT canonicals by content — content decides every row.
2. **Map by deliverable, not by number.** Read the full citing sentence, extract the
   described deliverable, OPEN the candidate target note, confirm its H1/boxed result
   matches. Evidence = the matching line in the target.
3. **Self-references are common here** (e.g. stage199:263 boxes "Theorem (Stage 250 …
   pairwise orbit-transport law)" = stage199's OWN theorem → 199). When a note states
   its own result with a stale number, map to the note's own canonical.
4. **Genuine forward refs may exist** (a note legitimately deferring to a later stage).
   If content shows the cited number really is the later stage's deliverable, LEAVE it.
   Class D found none among its flagged rows, but the broad pool is looser — verify.
5. **False-positive guard:** a longer stem that extends a shorter canonical stem is a
   real distinct stage (e.g. stage132_*_status ≠ 129). Confirm the cited file exists.
6. **Label-only application:** swap ONLY the number token; every other byte on the line
   stays identical (preserve hyphen/en-dash/slash/`\mathcal S_{}` bracket exactly). No
   prose/section/value/`\label` changes.
7. **WHO APPLIES:** doc-only label edits → orchestrator / fresh-context agent applies
   DIRECTLY (NOT Codex; Codex is for script CODE/MATH). Per [[numbering-drift-root-cause]].

## Per-band protocol (one band = one gated batch)
1. **Classify (read-only agents, fan-out by band).** Each agent takes its band's
   worklist rows, and for EACH row returns: cited token, described deliverable, proposed
   canonical (or LEAVE/genuine-forward, or CANNOT-CONFIRM), and the target-file
   `file:line` evidence. Agents also re-scan their notes for any candidate the inventory
   missed (adjacent residuals: `## k.` self-titles, "what Stage N changes" transitions).
2. **Orchestrator review.** Spot-check every self-ref, every unusual offset, every
   CANNOT-CONFIRM and LEAVE verdict against the files before applying.
3. **Apply (label-only).** A fresh applier agent (or orchestrator) applies the confirmed
   rows; CANNOT-CONFIRM rows stay unchanged and roll to a documented tail.
4. **Verify mechanically.** `git diff` must be byte-identical except number tokens
   (reuse the Class-D `verify_labelonly.py` normalizer). New targets all in 1–253.
5. **Record + gate.** Append results to a per-band section; HALT for explicit user go
   before the next band ([[feedback-sequential-audit-chunks]] — never roll forward).

## Recommended batch order (small→large to calibrate, then bulk)
| batch | band | candidates | notes | rationale |
| --- | --- | ---: | ---: | --- |
| 1 (pilot) | 136–163 | 10 | 4 | smallest; validates format + methodology |
| 2 | 091–135 | 17 | 9 | low-band stragglers (119/121/123/124/126/127/131/133/135/137/138) |
| 3 | 235–253 | 31 | 13 | endgame sympy-audit notes |
| 4 | 164–200 | 131 | 13 | bulk A (split 2 agents if needed) |
| 5 | 201–234 | 111 | 19 | bulk B (split 2 agents if needed) |

Bands 4/5 are large; each can be 2 parallel classifier agents (≤7 notes each) — but
agents only READ here, so concurrency is unconstrained (no Mathematica seat limit).

## Batch results

### Batch 1 (pilot, band 136–163) — ✅ DONE 2026-06-03
Classifier agent + orchestrator spot-verify, then applier agent (label-only), then
mechanical label-only verify. **8 of 10 PRIMARY mapped+applied; 2 were false positives.**
- stage137: `Stage 236`→**134** (Family-1 fixed-point law, box `Π=M_s+M_q𝒮_q(Π)`),
  `Stage 237`→**135** (Outlet-Consistent Mouth Closure). [+102 epoch]
- stage138: `Stage 239`→**137** ×4 (Explicit Core-to-Mouth Gain Map), `Stage 237`→**135**
  ×2, and the RANGE `Stages 170–174`→**119–123** (parent-normalized 𝔯/𝔤/Σ₀ — a DIFFERENT
  −51 epoch in the same file, content-confirmed; 170–174 = outlet map, wrong content).
- **False positives (LEAVE):** stage152:114 / stage153:57 `\mathcal S_1` = the Family-1
  source-correction SCALAR (≈0.6163, boxed beside 𝔤_1≈0.684), NOT "Stage 1". → New
  false-positive class: `\mathcal S_n` math symbols. The `mathcalS` token kind needs a
  content check every time (cf. the genuinely-stale stage193 `\mathcal S_{244}`).
- Diff: 2 files / 8 lines, byte-identical except number tokens; en-dash preserved. Uncommitted.

### Batch 2 (band 091–135) — ✅ DONE 2026-06-03
Two parallel read-only classifiers (119–126 / 127–135) + orchestrator spot-verify of the
load-bearing/repeated maps, then applier + robust multiset label-only verify. **24 edits / 9
notes** (17 forward = all the PRIMARY worklist candidates, MAP with 0 LEAVE/0 false-positive
among forward; + 7 backward stale refs found in-note). Multi-epoch, all content-confirmed:
- **+102 forward:** `Stage 220`→118 (parent-action core params, ×4 across 119/123),
  `Stage 221`→119 (𝔯/𝔤 defs, ×2), `Stage 231`→129 (mouth boundary layer: σ_Π + Π_m — NOT 130,
  which is the bias-map; content-distinguished, ×3 across 131/133), `Stage 233`→131 (parent
  micro-threshold, ×2 stage135), `Stage 223`→121 (self/geometric-r, ×2 stage121/124).
- **+51 ranges:** `Stages 172–173`→121–122 (123/126), `Stages 172–174`→121–123 (124),
  `Stages 176–177`→125–126 (127, positive-source theorem/families).
- **+17 backward:** `Stage-98`→115 (×3 stage119, core-balance), `Stage 99`→116 (×2 stage121,
  D/N tube), `Stages 97–99`→114–116 + `Stage 98`→115 (stage135, two-channel core/core-balance).
- **False positives LEFT:** `\mathcal S(Π,κ)` / `\mathcal S_q(Π)` (D/N mouth-response KERNEL,
  not "Stage n") in stage133/135. Confirmed canonical 220/221/231/233 stems are unrelated
  (dynamic-port / survival-window / continuum-pullback / orbit-lock) → genuinely stale, not forward.
- Verify: removed-multiset ≡ added-multiset (strictly label-only); 24/24 balanced. Committed.
- FYI (out of these 9 notes, for a later batch): canonical stage115 itself carries internal
  stale "Stage 97"/"Stage 95" self-labels.

### Batch 3 (band 235–253) — ✅ DONE 2026-06-03
Two parallel read-only classifiers (235–243 / 244–252) + orchestrator verify of every MAP.
**Only 6 of 31 candidates were stale; 25 were GENUINE and correctly LEFT** — this endgame band
is a sequential pipeline where most "Stage N" mentions are real next-stage/upstream continuation
refs (future-tense "Stage 245 should…", "Stage 253 must map…"). The genuine-vs-stale test
(future-looking + content-matches-cited-stage = genuine; "already proved" + cited content
mismatches = stale) was decisive.
- **6 stale (MAP), all +17 old-epoch:** stage237 L5/23/31/323/391 `Stage 253`→**236** (the
  equal-drift dependent-plane ray `-q_η(0,1,1)` / rigid-mouth quotient-failure carrier = stage236's
  deliverable; canonical 253 = physical calibration, unrelated); stage239 L499 range
  `Stages 252–237`→**235–237** (static-blind LINE-owner 235 → CURVE-owner 237; 252→235 at +17, 237 canonical).
- **25 LEFT:** 4 genuine next-stage continuation (235→236, 236→237, 243→244) + all 21 of band
  244–252 (genuine pipeline refs, 0 stale) + 3 false-positive `\mathcal S_{242}`/`\mathcal P_{242}`
  math-object subscripts in stage243 (242 correctly owns the packet/source pair).
- Verify: strictly label-only (6/6). Committed. **Lesson: the endgame band's default is GENUINE;
  bands 4/5 (164–234) are mid-derivation and expected to be mostly STALE like batches 1/2.**

### Calibration finding — SECONDARY is NOT "mostly genuine"
While in the 4 notes the classifier found content-confirmed STALE **backward-number** refs
at the dominant **+17** epoch (canonical = cited+17), which the forward-only PRIMARY filter
had placed in SECONDARY: stage137:10 `Stages 97–99`→**114–116**, stage137:58 `Stage 97`→**114**
(canonical 097 = single_normalization_defect, NOT the two-channel core), stage138:60 `Stage 98`
→**115**. Conversely stage152/153 `Stage 147/148` are genuine (current-epoch). **Implication:**
the 2354 backward refs are a MIX of stale (+17) and genuine — they need the same per-ref
content classification, not a blanket defer. RECOMMEND revising the plan so classifier agents
ALSO classify backward refs in each note they already read (near-zero marginal cost, avoids a
second pass over the same files). HELD pending user decision (the original gate deferred SECONDARY).

## Definition of done
All 300 PRIMARY candidates classified + applied-or-documented; mechanical label-only
verification clean per band; a final tail doc lists any CANNOT-CONFIRM rows with leads.
SECONDARY (2354 backward) explicitly deferred to a separate decision.
