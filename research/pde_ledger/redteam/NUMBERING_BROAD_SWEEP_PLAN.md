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

### Batch 4 (band 164–200) — ✅ DONE 2026-06-03
4 parallel read-only classifiers (198 / 199 / 190-191-193-194 / 170-177-178-179-187-189-192) +
orchestrator anchor-verify (188/189/191/192/193/197/198 H1s + 159/168/169) + 5 parallel appliers
(distinct files). **119 label-only edits / 11 notes.** Strip-all-digits proof: only numbers changed;
every per-file stale source-number fully gone (0 residual).
- **Uniform +51 ladder** (the 188–194 cluster, all content-anchored): 238→187, 239→188, 240→189,
  241→190, 242→191, 243→192, 244→193, 245→194. stage198 (242/243/248/249→191/192/197/198, 27 lines)
  & stage199 (243/248/249/250→192/197/198/199-SELF, 27 lines) were the dense ones; self-refs include
  boxed theorems ("Theorem (Stage 250 …)"→199, "(Stage 242 home-stretch)"→191).
- **Off-epoch (content-anchored):** stage190:588 `Stage 253`→**168** (−85, off-bundle slippage decomp);
  stage170 `Stage 237`→**169** (−68, no-go feed-down, ×4), `Stage 244`→**159** (−85, hybrid-outlet δE₂/δE₄, ×2),
  `Stage 238`→**170** (SELF).
- **Range refs (−51):** stage177 225–227→174–176, stage178 226–228→175–177, stage179 228–229→177–178,
  stage187 235–237→184–186 (corroborates the deferred-tail A2 fix: same monomials, same +51), plus
  239–241→188–190 / 238–239→187–188 in 190/191.
- **User's `\mathcal S_{244}` example RESOLVED:** stage193:332 `\mathcal S_{244}`→`\mathcal S_{193}` (SELF object-subscript);
  also 2 inventory-missed `post-243`→`post-192`. Inventory regex also missed multi-token lines (line 5 with
  3 refs, etc.) — the classifiers' full in-note grep caught them (lesson: the worklist undercounts dense lines).
- **GENUINE-left:** stage189 (2× "Stage 190") and stage192 (2× "Stage 193") are real next-stage refs → 0 edits.

### Batch 5 (band 201–234) — ✅ DONE 2026-06-03
5 parallel read-only classifiers (211 / 201 / 204+202 / 234+223 / 216+12-singletons) + orchestrator
anchor-verify (incl. the mixed/unusual offsets) + 5 parallel appliers. **80 label-only edits / 7 notes.**
Strip-all-digits proof: only numbers changed; per-file residuals 0; titles intact.
- **stage211** (21 lines, uniform −34): 243→209, 244→210, 245→211(SELF), 246→212 (cross-confirmed by stage212's canonical convention).
- **stage201** (16, +51): 238→187, 243→192, 249→198, 250→199, 251→200, 252→201(SELF); compound "249/250"→"198/199".
- **stage204** (18, MIXED, content-mapped not swept): 237→203 (×10), 253→202 (×4), 238→204 (×4 SELF), 243→192.
- **stage223** (13, −17): 239→222 (×9), 240→223 (×4 SELF).
- **stage234** (10, MIXED): 250→233 (×8), 251→234 (×2 SELF), 241→224, 249→232, 239→188 — three distinct offsets in one note, content-mapped.
- **stage232** 249→232 (SELF, +17); **stage233** "Stages 241 and 242"→"226 and 227" (+15, unusual — corroborated by stage233's own status line citing canonical "Stage-226").
- **GENUINE-left (batch-5E):** stage202 (0 edits) + the simplex-ladder notes 207/208/210/212/213/214/215/216/217/222/228 — all forward refs are real next-stage "should build…" pipeline refs; backward ranges (210→206–209, 228→223–227, 216→215/213/212) content-match; stage216 `\mathcal S_5^{can}` = math object (support cardinality), not a stage. **Lesson reaffirmed: the simplex/ray ladder (207–218) cites neighbors genuinely; only 2 of those 13 notes had any stale ref.**

## ✅ BROAD SWEEP (PRIMARY) COMPLETE — 2026-06-03
All 5 gated batches done. **~240 label-only line-edits across ~31 notes**, every one content-mapped
(never offset-swept) and mechanically verified label-only. Commits: batch1 `7f14506`, batch2 `172b657`,
batch3 `9f6bb7b`, batch4 `506fcbb`, batch5 `<this>`. Multi-epoch throughout (−17/−34/−51/−68/−85 forward
+ +15/+17 backward + self-refs), often several epochs in ONE note. Confirmed GENUINE-left where the cited
stage is the real next/neighbor (the 235–253 endgame pipeline + the 207–218 simplex ladder) — never swept.
Inventory lesson: `broad_sweep_scan.py`'s regex undercounts multi-token lines / `post-NNN` / compound
"N/M" tokens; the classifiers' full per-note grep caught the remainder. **REMAINING = SECONDARY only**
(2354 backward bare refs, a MIX of stale-+17 and genuine) — held for a separate decision.

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

## SECONDARY pass — execution plan (next session)
Worklist = `redteam/NUMBERING_SECONDARY_WORKLIST.md` (row-level, regenerate with `broad_sweep_scan.py`
against the then-current tree first). **2600 backward bare refs, MOSTLY GENUINE — default verdict LEAVE.**

**Honest expectation (low yield):** there is NO clean mechanical staleness signal — the `dist=self−cited`
column is citation distance, not staleness (a genuine stage022→"Stage 3" has dist +19; a stale
stage119→"Stage 98" also has large dist). Staleness is `cited→canonical` and needs content. From the
PRIMARY pass, stale backward refs cluster in the **+17 renumber zone, self∈[100,180]**; band 001–090 was
not renumbered (≈all genuine), and the 219–253 endgame cites neighbors genuinely. The worklist count is
inflated by self-refs and line-2 H1 titles (all → LEAVE). Realistic yield: a few dozen stale refs at most.

**Known stale classes already confirmed (fast-path when content matches):** `97→114` (two-channel core),
`98→115` (core-balance), `99→116` (D/N tube). Expect a few more to surface; log each new (cited→canonical).

**Recommended fan-out (agents):**
1. **Triage low band + endgame cheaply.** One agent confirms a SAMPLE of band 001–090 (and 219–253) is
   genuine, then LEAVE the whole band (don't per-ref-classify 1000+ genuine refs). Document the sample.
2. **Focus the real work on self∈[091,218].** Partition those notes into ~6–8 clusters; one read-only
   classifier per cluster content-classifies every backward ref (default LEAVE; MAP only on content
   mismatch + matching upstream). Same per-row evidence format as PRIMARY.
3. **Apply + verify per cluster** exactly like PRIMARY: orchestrator reviews every MAP, applier agent
   applies label-only, strip-all-digits proof (removed==added after `sed 's/[0-9]//g'`) + per-file
   residual-0, titles intact, commit per cluster, gate between clusters ([[feedback-sequential-audit-chunks]]).

The verified methodology, agent-prompt shapes, and the label-only proof are all in the PRIMARY batch
sections above — reuse them verbatim.

## Definition of done
All 300 PRIMARY candidates classified + applied-or-documented; mechanical label-only
verification clean per band; a final tail doc lists any CANNOT-CONFIRM rows with leads.
SECONDARY (2354 backward) explicitly deferred to a separate decision.

---
## SECONDARY pass — Cluster 1 (focus 091–218 + 2 targeted 219–253) — ✅ DONE 2026-06-03
Worklist regenerated against current tree (2600 backward refs). 8 read-only classifiers over focus self∈[091,218] + 2 sampling agents (001–090, 219–253). **25 label-only edits / 12 notes**, every one content-mapped (never offset-swept), strip-all-digits verified.
- **091–135 (13):** stage114:59/103 `Stage-95`→112, stage115:6 `Stage 97`→114 & `Stage 95`→112, :10/:74/:103 `Stage-95`→112, :31 `Stage-97`→114, stage116:40 `Stage 98`→115, :70 `Stage-97`→114, stage117:26 `Stages 97–99`→114–116, stage118:201/240 `Stage-97`→114.
- **136–163 (6):** stage158:344 `Stage 95`→112, stage160:64 `Stage 97`→114, :73 `Stage 98`→115, stage161:16/110/405 `Stage 99`→116.
- **179–189 (4):** stage180:141/252 `Stage-21`→38, stage182:83/138 `Stage-30`→47.
- **219–253 targeted (2):** stage220:396 `Stage-004`→021 (transfer factor; old-004=projected_maxwell), stage251:15 `Stage-12/13`→29/30 (selected quadrupole coeff; old-012/013=projected_maxwell).
- **NEW stale classes confirmed (verified vs target H1 + old-canonical H1):** `95→112` (Robin–mixed compensation law; old-095=geometry contamination), `21→038` (continuum placement map R_target/β₀; old-021=Maxwell one-port), `30→047` (coherent D/N kernel couplings; old-030=selected-mode response), `004→021`, `012/013→029/030`. (All +17.)
- Padding rule: preserve source token width (Stage-004→021 stays 3-digit; Stage-21→38, Stage-30→47 stay 2-digit; expand only when forced, e.g. 95→112).

### ⚠️ SCOPE FINDING — band 001–090 is NOT leave-able (deferred to Cluster 2)
The "001–090 was never renumbered → leave wholesale" premise is **FALSE for backward refs**. Note files jump 003→021; canonical 004–020 are the inserted projected-Maxwell / parent-throat-action scripts; the original wall/BdG chain was pushed +17 but bare `Stage N` prose refs were never renumbered. Sample: ~12/16 stale at +17 (e.g. stage024:5 `Stage 6`→023 grouped bundle; stage085:31 `Stage 13`→030; stage089:78 `Stage 35`→052). GENUINE survivors = refs to Stage 1/2/3 (never moved) and zero-padded `Stage 021+` forms. ~662 backward refs / 72 notes need full per-ref content classification — same methodology, NEVER offset-sweep (offsets split by content).
Also: the **164–178 classifier (C3a) was UNRELIABLE** — it waved through `Stage 6`/sub-21 refs as "genuine no-note-file canonical refs"; those are the same +17 stale class → re-run 164–178 in Cluster 2.

### Deferred tail (Cluster 2+)
- 136–163: 6 CANNOT-CONFIRM compounds — `Stages 99 and 170` (×3: 162:19/221, 163:45; "99"-half→119, "170" garbled), `Stages 95–98` (159:16; lead→112), `Stages 98–99` (161:61, 162:58; →115/116/119) — need range-aware handling.
- 219–253: in-band (cited ≥219) refs genuine → LEFT wholesale; uncertain `stage229:5 "Stage 143/093"` flagged.
- Focus zero-MAP bands 190–200, 201–218 verified genuine.

---
## SECONDARY pass — Cluster 2 (band 001–090 + 164–178) — ✅ DONE 2026-06-04
The "001–090 never renumbered" premise was FALSE: per notes/LINEAR_STAGE_RENUMBERING_MANIFEST.json the original chain was shifted +17 (001–003 unchanged; 004–020 = inserted projected-Maxwell scripts, no notes; 021 = former Stage 4; 022–253 = former 005–236 at +17). Band 001–090's bare prose `Stage N` refs were never updated → overwhelmingly stale.
- **8 read-only classifiers** (001-029/030-039/040-049/050-059/060-069/070-079/080-090 + a 164-178 re-run of the C3a `Stage 6` blind spot) → ~549 maps to redteam/secondary_maps/c2_*.md.
- **7 adversarial verifiers** (one per 001-090 cluster) re-checked every map vs the manifest +17 rule AND the old-form-vs-new-form trap (a genuine current-canonical ref mistaken for an old number). **0 wrong maps**; surfaced coverage gaps + 1 CANNOT.
- **Final cleanup classifier** over the remaining bare refs: +33 maps (21 stale SELF-refs the cluster agents missed in body prose; 12 genuine missed cross-refs incl. stage080 title `Stage-61`→078, stage084 `Stage 65`→082, stage045 `Stage-27`→044, stage081 `Stage-35`→052, stage170 `Stage-6`→023) and **49 genuine canonical 2-digit neighbor-refs LEFT** (the 070-089 explicit-branch band cites neighbors by their real canonical numbers).
- **Applied via redteam/apply_secondary_maps.py** (deterministic; per-edit uniqueness + label-only guard): **591 changed lines / 68 notes**. Strip-all-digits proof PASS (removed≡added after marker+digit strip → strictly label-only); numstat balanced 591/591.
- **NEW stale classes (all +17):** 5→022, 6→023, 7→024, 8→025 … 22→039; cross classes 16→033, 18→035, 27→044, 35→052, 61→078, 65→082; SELF-refs map to own canonical.
- **LEFT (documented):** stage062:233 `Stage 43` = author mis-citation (content home 061, NOT +17 060 → separate content fix, not label-only); stage039:256 `Stage 23` = forward ref (→040, forward pass); ~49 genuine canonical 2-digit neighbor refs in 070-089 (correct as-is; 3-digit padding is a separate cosmetic concern).
- Cluster 1 padding normalized to 3-digit and amended into its commit (e177262).
### Still deferred (Cluster 3+)
- 136–163: 6 CANNOT-CONFIRM compounds (`Stages 99 and 170`, `95–98`, `98–99`) — range-aware handling.
- 219–253 in-band refs LEFT genuine; uncertain stage229:5 `Stage 143/093` flagged.
- stage062:233 author mis-citation.
